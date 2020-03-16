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

The interpolator can be called with some additional parameter being, `n_neighbors` which is the number of nearest neighbors to be considered,
`power` is the exponent for distance in the weighing factor, `reg` is the offset for the weighing factor in denominator,
`conf_dist` is the distance below which two points would be considered as the same point.

!!! warning
    0 < `n_neighbor` <= number of elements in tree

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
function shepherd_util(idxs, dist, values, weights, power, reg, conf_dist)

    dist[1] <= conf_dist && return values[idxs[1]]

    num = den = zero(eltype(dist))

    for i in eachindex(idxs)
        w_i = 1 / (reg + (dist[i])^power)
        num += w_i * values[idxs[i]] * (weights === nothing ? 1 : weights[idxs[i]])
        den += w_i
    end
    return num / den
end

function ShepherdInterpolator(coordinates, values, points, weights, leafsize, n_neighbors, power, reg, conf_dist)
    # generating the KDTree
    tree = KDTree(coordinates, leafsize = leafsize)

    # querying the n-closest indexes and distances
    idxs, dist = knn(tree, points, n_neighbors, true)

    # calculate values at every desired point, can be done by some broadcasting trick
    out = Array{Float64}(undef, size(points,2))
    for i in eachindex(out)
        out[i] = shepherd_util(idxs[i], dist[i], values, weights, power, reg, conf_dist)
    end

    return out
end

ShepherdInterpolator(coordinates, values, points; weights = nothing, leafsize = 8, n_neighbors = 8, power = 1.0, reg = 0.0, conf_dist = 1e-12) =
        ShepherdInterpolator(coordinates, values, points, weights, leafsize, n_neighbors, power, reg, conf_dist)


struct IDWInterpolator <: BackgroundInterpolator
    factors::NTuple{2,<:Integer}
    leafsize
    n_neighbors
    power
    reg
    conf_dist
end

IDWInterpolator(factors; leafsize = 8, n_neighbors = 8, power = 1.0, reg = 0.0, conf_dist = 1e-12) = IDWInterpolator(factors, leafsize, n_neighbors, power, reg, conf_dist)

function (IDW::IDWInterpolator)(mesh::AbstractArray{T}, weights::AbstractArray{T}) where T
    # factors have to be integer tupules only
    # initialising the output array
    out = similar(mesh, float(T), size(mesh) .* IDW.factors)

    # Resizing the mesh and getting the know points

    # Creating an array of all points in modified mesh for query

    # Think something for the weights, as in how to access them

    # Call the ShepherdInterpolator on the corresponding data generated above

    # Fill out the final output array

    # return the output array
end

# finally make a constructor which takes weights as an optional argument



#####################################################################################
#
# # supports IDWInterpolation on 2-D images only
# # leafsize = 8, n_neighbors = 8, power = 1.0, reg = 0.0, conf_dist = 1e-12
# struct IDWInterpolator <: BackgroundInterpolator
#     factors::NTuple{2,<:Integer}
#     leafsize::Real
#     n_neighbors::Real
#     power::Real
#     reg::Real
#     conf_dist::Real
# end
#
# IDWInterpolator(factors; leafsize = 8, n_neighbors = 8, power = 1.0, reg = 0.0, conf_dist = 1e-12) = IDWInterpolator(factors, leafsize, n_neighbors, power, reg, conf_dist)


# function (IDW::IDWInterpolator)(mesh::AbstractArray{T,2}) where T
#     # this is for the final output
#     out = similar(mesh, float(T), size(mesh) .* IDW.factors)
#     display(out)
#     # generate the coefficients for the known meshpoints points where the values are known to form the KDTree
#     # current implementation only considers a 2-D mesh
#     coordinates = Array{Union{Missing, Int}}(missing, 2, length(mesh))
#     for i=1:size(mesh,1), j=1:size(mesh,2)
#         coordinates[1,(i - 1)*size(mesh,2) + j] = i * IDW.factor[1]
#         coordinates[2,(i - 1)*size(mesh,2) + j] = j *
#
#         if i == 1 ||  j == 1
#
#     end
#
#     # scale up coordinates before passing into tree creation
#     coordinates = coordinates .* IDW.factors
#
#     # declaring _coordinates because KDTree requires floating point input
#     _coordinates = 1.0 * coordinates
#
#     # filling the known points of array out
#     for i=1:size(coordinates,2)
#         out[coordinates[1,i], coordinates[2,i]] = mesh[convert(Int,coordinates[1,i] / IDW.factors[1]), convert(Int,coordinates[2,i] / IDW.factors[2])]
#     end
#
#     display(out)
#     # generate the KDTree
#     tree = KDTree(_coordinates, leafsize = IDW.leafsize)
#
#     # generate points to be queried on the tree for interpolation
#     positions = Array{Number}(undef, 2, length(out))
#     for i=1:size(out,1), j=1:size(out,2)
#         positions[1,(i - 1)*size(out,2) + j] = i * 1.0
#         positions[2,(i - 1)*size(out,2) + j] = j * 1.0
#     end
#
#
#     # query on tree and get distances returned in sorted manner
#     idxs, dist = knn(tree, positions, IDW.n_neighbors, true)
#
#     # filling the out array with interpolated values
#     for i in eachindex(out)
#         out[i] = idw_util(IDW, idxs[i], dist[i], out)
#     end
#
#     # returning the desired array
#     return out
# end


# alg = IDWInterpolator((2,2), n_neighbors = 3)([1 0; 0 1])
#
# display([1 0; 0 1])
#
#
# # I need to change the spacing in fabric
#
# x = range(0, 10, length=11)
