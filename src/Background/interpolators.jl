#=
Part of this work is derived from astropy/photutils and astropy/astropy. The relevant derivations
are considered under a BSD 3-clause license. =#

using Interpolations: CubicSplineInterpolation
using ImageTransformations: imresize!
using NearestNeighbors

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

"""
    IDWInterpolator(coordinates, values, weights = nothing, leafsize = 8.0)

This performs Inverse Distance Weighted Interpolation

`coordinates` represent the sample points filled column wise and `values` is the value at those sample points, `weight` are an additional multiplicative factor in the standard
inverse distance weight and `leafsize` determines at what number of points to stop splitting the KDTree further.

!!! warning
    0 < `leafsize` < number of elements in tree

# Examples
```jldoctest
julia> coordinates = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0; 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];

julia> weights = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];

julia> value = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];

julia> positions = [9.9 1.0; 10.0 1.0];

julia> alg = IDWInterpolator(coordinates,value, weights, 8);

julia> alg(positions)
1×2 Array{Float64,2}:
 9.57581  2.75894

julia> IDWInterpolator(coordinates, value)(positions)
1×2 Array{Float64,2}:
 9.57581  2.75894

```
"""
struct IDWInterpolator <: BackgroundInterpolator
    coordinates::AbstractArray
    values::AbstractArray
    weights::Union{AbstractArray, Nothing}
    leafsize::Real

    function IDWInterpolator(coordinates::AbstractArray, values::AbstractArray, weights::Union{AbstractArray, Nothing}, leafsize::Real)
        leafsize <= 0 && error("Invalid leafsize=$leafsize. leafsize must be greater than 0")
        new(coordinates, values, weights, leafsize)
    end
end

IDWInterpolator(coordinates, value; weight = nothing, leafsize = 8) = IDWInterpolator(coordinates, value, weight, leafsize)

function idw_util(alg::IDWInterpolator, idxs::AbstractArray, dist::AbstractArray, power::Real, reg::Real, conf_dist::Real)

    dist[1] <= conf_dist && return alg.values[idxs[1]]

    num = den = zero(eltype(dist))

    for (i,j) = zip(1:length(idxs), 1:length(dist))
        w_i = 1 / (reg + dist[i])^power
        num += w_i * alg.values[idxs[i]] * (alg.weights === nothing ? 1 : alg.weights[idxs[i]])
        den += w_i
    end
    return num / den
end


function idw_interpolator(alg::IDWInterpolator, positions::AbstractArray, n_neighbors::Real, power::Real, reg::Real, conf_dist::Real)
    tree = KDTree(alg.coordinates, leafsize = alg.leafsize)

    idxs, dist = knn(tree, positions, n_neighbors, true)

    _out = zeros(eltype(dist[1]), 1, size(positions, 2))
    for i = 1:length(_out)
        _out[i] = idw_util(alg, idxs[i], dist[i], power, reg, conf_dist)
    end

    return _out
end

(alg::IDWInterpolator)(positions; n_neighbors = 8,
                                  power = 1.0,
                                  reg = 0.0,
                                  conf_dist = 1e-12) = idw_interpolator(alg, positions, n_neighbors, power, reg, conf_dist)
