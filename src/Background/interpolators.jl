#=
Part of this work is derived from astropy/photutils and astropy/astropy. The relevant derivations
are considered under a BSD 3-clause license. =#

using Interpolations: InterpolationType, CubicSplineInterpolation, AbstractInterpolation
using ImageTransformations: imresize!
using NearestNeighbors
using Distances

const MinkowskiMetric = Union{Euclidean,Chebyshev,Cityblock,Minkowski}

"""
    ZoomInterpolator(factors)

Use a cubic-spline interpolation scheme to increase resolution of a mesh.

`factors` represents the level of "zoom", so an input mesh of size `(10, 10)` with factors `(2, 2)` will have an output size of `(20, 20)`. If only an integer is provided, it will be used as the factor for every axis.

# Examples
```jldoctest
julia> ZoomInterpolator(2)([1 0; 0 1])
4×4 Array{Float64,2}:
  1.0          0.75   0.25   -2.77556e-17
  0.75         0.625  0.375   0.25
  0.25         0.375  0.625   0.75
 -5.55112e-17  0.25   0.75    1.0

julia> ZoomInterpolator(3, 1)([1 0; 0 1])
6×2 Array{Float64,2}:
  1.0          -2.77556e-17
  1.0          -2.77556e-17
  0.666667      0.333333
  0.333333      0.666667
 -5.55112e-17   1.0
 -5.55112e-17   1.0

```
"""
struct ZoomInterpolator <: BackgroundInterpolator
    factors::NTuple{2,<:Integer}
end

ZoomInterpolator(factor::Integer) = ZoomInterpolator((factor, factor))
ZoomInterpolator(f1::Integer, args...) = ZoomInterpolator((f1, args...))

function (z::ZoomInterpolator)(mesh::AbstractArray{T}) where T
    itp = CubicSplineInterpolation(axes(mesh), mesh)
    out = similar(mesh, float(T), size(mesh) .* z.factors)
    return imresize!(out, itp)
end
​

"""
    IDWInterpolator(coordinates, values, weights = nothing, leafsize = 8.0)

This performs Inverse Distance Weighted Interpolation

`coordinates` represent the sample points filled column wise and `values` is the value at those sample points, `weight` are an additional multiplicative factor in the standard
inverse distance weight and `leafsize` determines at what number of points to stop splitting the KDTree further.

The interpolator can be called with some additional parameter being, `n_neighbors` which is the number of nearest neighbors to be considered,
`power` is the exponent for distance in the weighing factor, `reg` is the offset for the weighing factor in denominator,
`conf_dist` is the distance below which two points would be considered as the same point.

!!! warning
    0 < `n_neighbor` <= number of elements in tree

# Examples
```jldoctest


julia> IDWInterpolator(coordinates, value)(positions)
1×2 Array{Float64,2}:
 9.57581  2.75894

```
"""
struct IDWInterpolator <: BackgroundInterpolator
    factors::NTuple{2,<:Integer}
    leafsize::Integer
    n_neighbors::Integer
    power::Real
    reg::Real
    conf_dist::Real
end

IDWInterpolator(factors; leafsize = 8, n_neighbors = 8, power = 1.0, reg = 0.0, conf_dist = 1e-12) = IDWInterpolator(factors, leafsize, n_neighbors, power, reg, conf_dist)
# convenience constructors
IDWInterpolator(factor::Integer; kwargs...) = IDWInterpolator((factor, factor); kwargs...)
IDWInterpolator(factor::Integer, args...; kwargs...) = IDWInterpolator((factor, args...); kwargs...)

function (IDW::IDWInterpolator)(mesh::AbstractArray{T}) where T
    knots = Array{Float64}(undef, 2, length(mesh))
    idxs = CartesianIndices(mesh)
    for (i, idx) in enumerate(CartesianIndices(mesh))
        @inbounds knots[:, i] .= Tuple(idx)
    end

    itp = ShepardIDWInterpolator(knots, mesh, IDW.leafsize, IDW.n_neighbors, IDW.power, IDW.reg, IDW.conf_dist)
    out = similar(mesh, float(T), size(mesh) .* IDW.factors)
    return imresize!(out, itp)
end


###############################################################################
# Interface with Interpolations.jl

struct IDW <: InterpolationType end

struct ShepardIDWInterpolator{T,N} <: AbstractInterpolation{T,N,IDW}
    tree::KDTree{<:AbstractVector,<:MinkowskiMetric,T}
    values::Array{T,N}
    n_neighbors::Integer
    power::Real
    reg::Real
    conf_dist::Real
end

Base.axes(itp::ShepardIDWInterpolator) = axes(itp.values)

function ShepardIDWInterpolator(knots,
    values::AbstractArray{T},
    leafsize = 8,
    n_neighbors = 8,
    power = 1,
    reg = 0,
    conf_dist = 1e-12) where T

    tree = KDTree(knots, leafsize = leafsize)
    return ShepardIDWInterpolator(tree, values, n_neighbors, power, reg, conf_dist)
end

function (itp::ShepardIDWInterpolator{T,N})(points::Vararg{T,N}) where {T,N}

    _points = Array{T}(undef, length(points), 1)
    _points .= points

    # find the n-closest indices and distances
    idxs, dist = knn(itp.tree, _points, itp.n_neighbors, true)

    dist[1][1] <= itp.conf_dist && return itp.values[idxs[1][1]]

    # no-allocation loop calculating using Shepard's scheme
    num = den = zero(T)
    for (i, d) in zip(idxs[1], dist[1])
        w = 1 / (itp.reg + d^itp.power)
        num += w * itp.values[i]
        den += w
    end

    return num / den
end
