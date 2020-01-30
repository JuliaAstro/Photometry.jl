module Aperture

using DataFrames: DataFrame

export mask,
       cutout,
       aperture_photometry

abstract type AbstractAperture end

"""
    mask(::AbstractAperture)

Return an array of the weighting of the aperture in the minimum bounding box.
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


function overlap_slices(c::AbstractAperture, shape::Tuple)
    xmin, xmax, ymin, ymax = bbox(c)

    # No overlap
    if xmin ≥ shape[2] || ymin ≥ shape[1] || xmax ≤ 1 || ymax ≤ 1
        return nothing, nothing
    end

    # TODO something is wrong here, indexing is screwy
    slices_large = (max(ymin, 1):min(ymax, shape[1]), 
                    max(xmin, 1):min(xmax, shape[2]))
    slices_small = (max(-ymin, 1):min(ymax - ymin, shape[1] - ymin), 
                    max(-xmin, 1):min(xmax - xmin, shape[2] - xmin))

    return slices_large, slices_small
end

function cutout(c::AbstractAperture, data::AbstractMatrix{T}) where {T}
    box = bbox(c)
    maxy, maxx = size(data)
    ny, nx = size(c)

    # Check if x or y are less than our minimum index
    partial_overlap = box[1] < 1 || box[3] < 1 || box[2] > maxx || box[4] > maxy
    
    if !partial_overlap
        cut = data[box[3]:box[4], box[1]:box[2]]
    end
    
    if partial_overlap || size(cut) != size(c)
        slices_large, slices_small = overlap_slices(c, size(data))
        # no overlap
        slices_small === nothing && return nothing
        
        cut = zeros(T, ny, nx)
        @show slices_small
        @show slices_large
        cut[slices_small...] .= data[slices_large...]
        return cut
    end
    return cut
end

"""
    *(::AbstractAperture, data::AbstractMatrix)
    *(data::AbstractMatrix, ::AbstractAperture)

Return the data weighted by the aperture. If there is no overlap, will return an empty Array
"""
function Base.:*(c::AbstractAperture, data::AbstractMatrix{T}) where {T}
    cut = cutout(c, data)
    cut === nothing && return typeof(data)(undef, 0, 0)
    return cut .* mask(c)
end

Base.:*(data::AbstractMatrix, c::AbstractAperture) = c * data

"""
    aperture_photometry(::AbstractAperture, data::AbstractMatrix)
    aperture_photometry(::AbstractVector{<:AbstractAperture}, data::AbstractMatrix)

Perform aperture photometry on `data` given aperture(s). If a list of apertures is provided the output will be a `DataFrame`, otherwise a `NamedTuple`.

"""
function aperture_photometry(a::AbstractAperture, data::AbstractMatrix)
    data_weighted = data * a
    aperture_sum = sum(data_weighted)
    return (xcenter = a.x, ycenter = a.y, aperture_sum = aperture_sum)
end

aperture_photometry(a::AbstractVector{<:AbstractAperture}, data::AbstractMatrix) = DataFrame(aperture_photometry.(a, Ref(data)))

include("circular.jl")
include("overlap.jl")

end