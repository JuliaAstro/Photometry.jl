#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

module Aperture

using TypedTables

export photometry,
       Subpixel,
       CircularAperture,
       CircularAnnulus,
       EllipticalAperture,
       EllipticalAnnulus,
       RectangularAperture,
       RectangularAnnulus

"""
The abstract super-type for Apertures
"""
abstract type AbstractAperture end

"""
    bounds(::AbstractAperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture.
"""
bounds(::AbstractAperture)

"""
    size(::AbstractAperture)

Return (`ny`, `nx`) of the aperture.
"""
Base.size(ap::AbstractAperture) = map(length, axes(ap))
Base.size(ap::AbstractAperture, dim) = length(axes(ap, dim))

function cutout_indices(ap::AbstractAperture, data::AbstractArray)
    # slices for indexing the larger array
    xmin, xmax, ymin, ymax = bounds(ap)

    x_slice = max(xmin, firstindex(data, 2)):min(xmax, lastindex(data, 2))
    y_slice = max(ymin, firstindex(data, 1)):min(ymax, lastindex(data, 1))
    return CartesianIndices((y_slice, x_slice))
end

Base.getindex(ap::AbstractAperture, idx::CartesianIndex) = getindex(ap, idx.I...)
Base.getindex(ap::AbstractAperture, idxs) = map(idx -> ap[idx], idxs)
Base.eachindex(ap::AbstractAperture) = CartesianIndices(axes(ap))
Base.collect(ap::AbstractAperture) = ap[eachindex(ap)]

function Base.axes(ap::AbstractAperture)
    xmin, xmax, ymin, ymax = bounds(ap)
    return ymin:ymax, xmin:xmax
end
Base.axes(ap::AbstractAperture, dim) = axes(ap)[dim]


struct Subpixel{AP<:AbstractAperture} <: AbstractAperture
    ap::AP
    N::Int
end

function Base.show(io::IO, c::Subpixel)
    print(io, "Subpixel(")
    print(io, c.ap)
    print(io, ", $(c.N))")
end

bounds(sp_ap::Subpixel) = bounds(sp_ap.ap)
overlap(sp_ap::Subpixel, args...) = overlap(sp_ap.ap, args...)

Subpixel(ap::AbstractAperture) = Subpixel(ap, 1)


function Base.getproperty(ap::Subpixel, key::Symbol)
    key in (:N, :ap) && return getfield(ap, key)
    return getfield(ap.ap, key)
end

@enum OverlapFlag Inside Outside Partial

function Base.getindex(ap::AbstractAperture, i, j)
    flag = overlap(ap, i, j)
    flag === Outside && return 0.0
    flag === Inside && return 1.0

    return partial(ap, j - ap.x, i - ap.y)
end


###########

"""
    photometry(::AbstractAperture, data::AbstractMatrix, [error])
    photometry(::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, [error])

Perform aperture photometry on `data` given aperture(s). If `error` (the pixel-wise standard deviation) is provided, will calculate sum error. If a list of apertures is provided the output will be a `TypedTables.Table`, otherwise a `NamedTuple`.

# Methods
* `:exact` - Will calculate the exact geometric overlap
* `:center` - Will only consider full-pixel overlap (equivalent to subpixel method with 1 subpixel)
* `(:subpixel, n)` - Use `n^2` subpixels to calculate overlap

!!! note
    The `:exact` method is slower than the subpixel methods by at least an order of magnitude, so if you are dealing with large images and many apertures, we recommend using `:subpixel` with some reasonable `n`, like 10.

!!! tip
    This code is automatically multi-threaded. To take advantage of this please make sure `JULIA_NUM_THREADS` is set before starting your runtime.
"""
function photometry(ap::AbstractAperture, data::AbstractMatrix, error)
    idxs = cutout_indices(ap, data)
    isempty(idxs) && return (xcenter = ap.x, ycenter = ap.y, aperture_sum = 0.0, aperture_sum_err = NaN)
    aperture_sum = sum(idx -> ap[idx] * data[idx], idxs)
    aperture_sum_err = sqrt(sum(idx -> ap[idx] * error[idx]^2, idxs))

    return (xcenter = ap.x, ycenter = ap.y, aperture_sum = aperture_sum, aperture_sum_err = aperture_sum_err)
end

function photometry(ap::AbstractAperture, data::AbstractMatrix)
    idxs = cutout_indices(ap, data)
    isempty(idxs) && return (xcenter = ap.x, ycenter = ap.y, aperture_sum = 0.0)
    aperture_sum = sum(idx -> ap[idx] * data[idx], idxs)
    return (xcenter = ap.x, ycenter = ap.y, aperture_sum = aperture_sum)
end

function photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, error)
    rows = similar(aps, NamedTuple{(:xcenter, :ycenter, :aperture_sum, :aperture_sum_err)})
    Threads.@threads  for idx in eachindex(rows)
        rows[idx] = photometry(aps[idx], data, error)
    end
    return Table(rows)
end

function photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix)
    rows = similar(aps, NamedTuple{(:xcenter, :ycenter, :aperture_sum)})
    Threads.@threads for idx in eachindex(rows)
        rows[idx] = photometry(aps[idx], data)
    end
    return Table(rows)
end

include("circular.jl")
include("elliptical.jl")
include("rectangle.jl")
include("overlap.jl")
include("plotting.jl")

end

