#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

module Aperture

using TypedTables
using Transducers

export photometry,
       Subpixel,
       CircularAperture,
       CircularAnnulus,
       EllipticalAperture,
       EllipticalAnnulus,
       RectangularAperture,
       RectangularAnnulus

"""
    AbstractAperture{T} <: AbstractMatrix{T}

The abstract super-type for Apertures.

Apertures can be thought of as a cutout or stamp of a geometric shape with
shading applied. For example, a circular aperture with a diameter of 3 pixels
will require a 5x5 pixel grid (when perfectly on-grid) to represent.

```jldoctest ap1
julia> ap = CircularAperture(3, 3, 2.5)
5×5 CircularAperture{Float64} with indices 1:5×1:5:
 0.136857  0.769325  0.983232  0.769325  0.136857
 0.769325  1.0       1.0       1.0       0.769325
 0.983232  1.0       1.0       1.0       0.983232
 0.769325  1.0       1.0       1.0       0.769325
 0.136857  0.769325  0.983232  0.769325  0.136857
```

This is a useful way of thinking about apertures: if we have some data,
we can weight the data with the aperture.

```jldoctest ap1
julia> data = fill(2, 5, 5);

julia> idxs = map(intersect, axes(ap), axes(data)) |> CartesianIndices;

julia> weighted_cutout = data[idxs] .* ap[idxs]
5×5 Matrix{Float64}:
 0.273713  1.53865  1.96646  1.53865  0.273713
 1.53865   2.0      2.0      2.0      1.53865
 1.96646   2.0      2.0      2.0      1.96646
 1.53865   2.0      2.0      2.0      1.53865
 0.273713  1.53865  1.96646  1.53865  0.273713
```

Performing aperture photometry is merely summing the weighted cutout shown above.

```jldoctest ap1
julia> flux = sum(weighted_cutout)
39.269908169872416

julia> flux ≈ (π * 2.5^2) * 2 # area of circle times intensity of 2
true
```

What's interesting about the implementation of apertures, though,
is they are lazy. This means there is no stored matrix of aperture values;
rather, they are calculated on the fly as needed.

```jldoctest ap1
julia> axes(ap)
(1:5, 1:5)

julia> ap[-10, -10] # out-of-bounds, but calculated on the fly
0.0

julia> ap .* ones(5, 7) # broadcasts to eachindex(data), regardless of ap bound
5×7 Matrix{Float64}:
 0.136857  0.769325  0.983232  0.769325  0.136857  0.0  0.0
 0.769325  1.0       1.0       1.0       0.769325  0.0  0.0
 0.983232  1.0       1.0       1.0       0.983232  0.0  0.0
 0.769325  1.0       1.0       1.0       0.769325  0.0  0.0
 0.136857  0.769325  0.983232  0.769325  0.136857  0.0  0.0
```

This allows extremely efficient computation of aperture photometry from small to
medium sized apertures.

```julia-repl
julia> using BenchmarkTools

julia> @btime sum(idx -> \$ap[idx] * \$data[idx], \$idxs)
  1.097 μs (0 allocations: 0 bytes)
39.26990816987243
```

This is essentially the full implementation of [`photometry`](@ref), save for
the packing of additional information into a tabular form.
"""
abstract type AbstractAperture{T}  <: AbstractMatrix{T} end

"""
    bounds(::AbstractAperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture.
"""
bounds(::AbstractAperture)

center(ap::AbstractAperture) = ap.x, ap.y # greedy

"""
    size(::AbstractAperture)

Return (`ny`, `nx`) of the aperture.
"""
function Base.size(ap::AbstractAperture)
    xmin, xmax, ymin, ymax = bounds(ap)
    return xmax - xmin + 1, ymax - ymin + 1
end

function Base.axes(ap::AbstractAperture)
    xmin, xmax, ymin, ymax = bounds(ap)
    return xmin:xmax, ymin:ymax
end

"""
    Subpixel(ap, N=1) <: AbstractAperture

Use a subpixel quadrature approximation for pixel shading instead of exact
geometric methods.

For any pixel laying on the border of `ap`, this alters the shading algorithm by
breaking the border pixel up into `(N, N)` subpixels. The shading value is the
fraction of these subpixels within the geometric border of `ap`.

Using a subpixel shading method is sometimes faster than exact methods at the
cost of accuracy. For [`CircularAperture`](@ref) the subpixel method is only
faster than the exact method for `N` ~ 7. for [`EllipticalAperture`](@ref) the
cutoff is `N` ~ 12, and for [`RectangularAperture`](@ref) the cutoff is `N` ~ 20.


# Examples

```jldoctest
julia> ap = CircularAperture(3, 3, 2.5)
5×5 CircularAperture{Float64} with indices 1:5×1:5:
 0.136857  0.769325  0.983232  0.769325  0.136857
 0.769325  1.0       1.0       1.0       0.769325
 0.983232  1.0       1.0       1.0       0.983232
 0.769325  1.0       1.0       1.0       0.769325
 0.136857  0.769325  0.983232  0.769325  0.136857

julia> sub_ap = Subpixel(ap, 5)
5×5 Subpixel{Float64, CircularAperture{Float64}} with indices 1:5×1:5:
 0.12  0.76  1.0  0.76  0.12
 0.76  1.0   1.0  1.0   0.76
 1.0   1.0   1.0  1.0   1.0
 0.76  1.0   1.0  1.0   0.76
 0.12  0.76  1.0  0.76  0.12
```

!!! note
    `photutils` offers a `center` shading method which is equivalent to using
    the `Subpixel` method with 1 subpixel. To avoid unneccessary namespace
    cluttering, we simply instruct users to use `Subpixel(ap)` instead.
"""
struct Subpixel{T,AP<:AbstractAperture{T}} <: AbstractAperture{T}
    ap::AP
    N::Int
end

function Base.show(io::IO, c::Subpixel)
    print(io, "Subpixel(")
    print(io, c.ap)
    print(io, ", $(c.N))")
end

center(sp_ap::Subpixel) = center(sp_ap.ap)
bounds(sp_ap::Subpixel) = bounds(sp_ap.ap)
overlap(sp_ap::Subpixel, i, j) = overlap(sp_ap.ap, i, j)

Subpixel(ap::AbstractAperture) = Subpixel(ap, 1)

@enum OverlapFlag Inside Outside Partial

function Base.getindex(ap::AbstractAperture{T}, idx::Vararg{Int,2}) where T
    i, j = idx
    flag = overlap(ap, i, j)
    # TODO: revisit a better way to handle subpixel apertures
    #flag === Outside && return zero(T)
    flag === Inside && return one(T) / one(T)
    cx, cy = center(ap)
    return partial(ap, i - cx, j - cy)
end

# This bypasses checking aperture axes for broadcasting
Broadcast.combine_axes(ap::AbstractAperture, arrs...) = Broadcast.combine_axes(arrs...)
Broadcast.combine_axes(arr, ap::AbstractAperture) = axes(arr)

###########

"""
    photometry(::AbstractAperture, data::AbstractMatrix, [error]; [f = sum])
    photometry(::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, [error]; [f = sum])

Perform aperture photometry on `data` given aperture(s). If `error` (the
pixel-wise standard deviation) is provided, will calculate sum error. If a list
of apertures is provided the output will be a `TypedTables.Table`, otherwise a
`NamedTuple`. An optional function `f` can be passed to return additional statistics
within each aperture. This can be useful for, e.g., computing the PSF of each source. By default, just the sum within each aperture is returned.

!!! tip
    This code is automatically multi-threaded. To take advantage of this please
    make sure `JULIA_NUM_THREADS` is set before starting your runtime.
"""
function photometry(ap::AbstractAperture, data::AbstractMatrix, error; f = sum)
    cx, cy = center(ap)
    meta = (xcenter = cx, ycenter = cy)
    idxs = map(intersect, axes(ap), axes(data), axes(error))
    if any(isempty, idxs)
        if f == sum
            return (meta..., aperture_sum = 0.0, aperture_sum_err = NaN)
        else
            return (meta..., aperture_sum = 0.0, aperture_sum_err = NaN, aperture_f = 0.0)
        end
    end
    img_ap = CartesianIndices(idxs) |> Map(idx -> ap[idx] * data[idx])
    img_ap_var = CartesianIndices(idxs) |> Map(idx -> ap[idx] * error[idx]^2)

    aperture_sum = sum(img_ap)
    aperture_sum_var = sum(img_ap_var)
    aperture_sum_err = sqrt(aperture_sum_var)

    if f == sum
        return (; meta..., aperture_sum, aperture_sum_err)
    else
        aperture_f = f(img_ap)
        # Note: not supported if f returns a Tuple
        #aperture_f_var = f(img_ap_var)
        #aperture_f_err = sqrt(aperture_f_var)
        return (; meta..., aperture_sum, aperture_sum_err, aperture_f)
    end
end

function photometry(ap::AbstractAperture, data::AbstractMatrix; f = sum)
    cx, cy = center(ap)
    meta = (xcenter = cx, ycenter = cy)
    idxs = map(intersect, axes(ap), axes(data))
    if any(isempty, idxs)
        if f == sum
            return (meta..., aperture_sum = 0.0)
        else
            return (meta..., aperture_sum = 0.0, aperture_f = 0.0)
        end
    end

    img_ap = CartesianIndices(idxs) |> Map(idx -> ap[idx] * data[idx])
    aperture_sum = sum(img_ap)

    if f == sum
        return (; meta..., aperture_sum)
    else
        aperture_f = f(img_ap)
        return (; meta..., aperture_sum, aperture_f)
    end
end

function photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, error; f = sum)
    rows = tcollect(aps |> Map(ap -> photometry(ap, data, error; f)))
    return Table(rows)
end

function photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix; f = sum)
    rows = tcollect(aps |> Map(ap -> photometry(ap, data; f)))
    return Table(rows)
end

include("circular.jl")
include("elliptical.jl")
include("rectangle.jl")
include("overlap.jl")
include("plotting.jl")

end
