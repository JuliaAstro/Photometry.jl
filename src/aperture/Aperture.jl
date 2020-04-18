#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

module Aperture

using DataFrames: DataFrame

export mask,
       cutout,
       aperture_photometry,
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
    mask(::AbstractAperture; method=:exact)

Return an array of the weighting of the aperture in the minimum bounding box. For an explanation of the different methods, see [`aperture_photometry`](@ref).
"""
mask(::AbstractAperture)

"""
    bbox(::AbstractAperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture.
"""
bbox(::AbstractAperture)

"""
    size(::AbstractAperture)

Return (`ny`, `nx`) of the aperture.
"""
function Base.size(a::AbstractAperture)
    box = bbox(a)
    return (box[4] - box[3] + 1, box[2] - box[1] + 1)
end

function edges(ap::AbstractAperture)
    ibox = bbox(ap)
    xmin = ibox[1] - ap.x - 0.5
    xmax = ibox[2] - ap.x + 0.5
    ymin = ibox[3] - ap.y - 0.5
    ymax = ibox[4] - ap.y + 0.5
    return (xmin, xmax, ymin, ymax)
end

function overlap_slices(c::AbstractAperture, shape::Tuple)
    xmin, xmax, ymin, ymax = bbox(c)

    # No overlap
    if xmin ≥ shape[2] || ymin ≥ shape[1] || xmax ≤ 1 || ymax ≤ 1
        return nothing, nothing
    end

    # slices for indexing the larger array
    slices_large = (max(ymin, 1):min(ymax, shape[1]),
                    max(xmin, 1):min(xmax, shape[2]))

    # slices for indexing the smaller array
    slices_small = (max(2 - ymin, 1):min(ymax - ymin, shape[1] - ymin) + 1,
                    max(2 - xmin, 1):min(xmax - xmin, shape[2] - xmin) + 1)

    return slices_large, slices_small
end

"""
    cutout(::AbstractAperture, data)

Get the cutout of the aperture from the `data`. This will handle partial overlap by padding the data with zeros.
"""
function cutout(c::AbstractAperture, data::AbstractMatrix{T}) where T
    box = bbox(c)
    maxy, maxx = size(data)
    ny, nx = size(c)

    # Check if x or y are less than our minimum index
    partial_overlap = box[1] < 1 || box[3] < 1 || box[2] > maxx || box[4] > maxy

    if !partial_overlap
        cutout = data[box[3]:box[4], box[1]:box[2]]
    end

    if partial_overlap || size(cutout) != size(c)
        slices_large, slices_small = overlap_slices(c, size(data))
        # no overlap
        slices_small === nothing && return nothing

        cutout = zeros(T, ny, nx)
        cutout[slices_small...] .= data[slices_large...]
    end
    return cutout
end

function apply(a::AbstractAperture, data::AbstractMatrix; method = :exact)
    cut = cutout(a, data)
    cut === nothing && return similar(data, 0, 0)
    return cut .* mask(a, method = method)
end

"""
    aperture_photometry(::AbstractAperture, data::AbstractMatrix, [error]; method=:exact)
    aperture_photometry(::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, [error]; method=:exact)

Perform aperture photometry on `data` given aperture(s). If `error` (the pixel-wise standard deviation) is provided, will calculate sum error. If a list of apertures is provided the output will be a `DataFrame`, otherwise a `NamedTuple`.

# Methods
* `:exact` - Will calculate the exact geometric overlap
* `:center` - Will only consider full-pixel overlap (equivalent to subpixel method with 1 subpixel)
* `(:subpixel, n)` - Use `n^2` subpixels to calculate overlap

!!! note
    The `:exact` method is slower than the subpixel methods by at least an order of magnitude, so if you are dealing with large images and many apertures, we recommend using `:subpixel` with some reasonable `n`, like 10.

!!! tip
    This code is automatically multi-threaded. To take advantage of this please make sure `JULIA_NUM_THREADS` is set before starting your runtime.
"""
function aperture_photometry(a::AbstractAperture, data::AbstractMatrix, error; method = :exact)
    data_weighted = apply(a, data, method = method)
    aperture_sum = sum(data_weighted)
    variance_weighted = apply(a, error.^2, method = method)
    aperture_sum_err = sqrt(sum(variance_weighted))

    return (xcenter = a.x, ycenter = a.y, aperture_sum = aperture_sum, aperture_sum_err = aperture_sum_err)
end

function aperture_photometry(a::AbstractAperture, data::AbstractMatrix; method = :exact)
    data_weighted = apply(a, data, method = method)
    aperture_sum = sum(data_weighted)

    return (xcenter = a.x, ycenter = a.y, aperture_sum = aperture_sum)
end

function aperture_photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, error; method = :exact)
    rows = similar(aps, NamedTuple{(:xcenter, :ycenter, :aperture_sum, :aperture_sum_err)})
    Threads.@threads  for idx in eachindex(rows)
        rows[idx] = aperture_photometry(aps[idx], data, error; method = method)
    end
    return DataFrame(rows)
end

function aperture_photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix; method = :exact)
    rows = similar(aps, NamedTuple{(:xcenter, :ycenter, :aperture_sum)})
    Threads.@threads for idx in eachindex(rows)
        rows[idx] = aperture_photometry(aps[idx], data; method = method)
    end
    return DataFrame(rows)
end



include("circular.jl")
include("elliptical.jl")
include("rectangle.jl")
include("overlap.jl")
include("plotting.jl")
end
