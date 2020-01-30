module Photometry

using DataFrames: DataFrame

export mask,
       cutout,
       aperture_photometry

abstract type Aperture end

"""
    mask(::Aperture)

Return an array of the weighting of the aperture in the minimum bounding box.
"""
mask(::Aperture)

"""
    bbox(::Aperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture.
"""
bbox(::Aperture)

"""
    size(::Aperture)

Return (`ny`, `nx`) of the aperture.
"""
function Base.size(a::Aperture)
    box = bbox(a)
    return (box[4] - box[3] + 1, box[2] - box[1] + 1)
end


function overlap_slices(c::Aperture, shape::Tuple)
    xmin, xmax, ymin, ymax = bbox(c)

    # No overlap
    if xmin ≥ shape[2] || ymin ≥ shape[1] || xmax ≤ 1 || ymax ≤ 1
        return nothing, nothing
    end

    slices_large = (vcat(max(ymin, 1), 2:min(ymax, shape[1])), 
                    vcat(max(xmin, 1), 2:min(xmax, shape[2])))
    slices_small = (vcat(max(-ymin, 1), 2:min(ymax - ymin, shape[1] - ymin)), 
                    vcat(max(-xmin, 1), 2:min(xmax - xmin, shape[2] - xmin)))

    return slices_large, slices_small
end

function cutout(c::Aperture, data::AbstractMatrix{T}) where {T}
    box = bbox(c)

    # Check if x or y are less than our minimum index
    partial_overlap = box[1] < 1 || box[3] < 1

    if !partial_overlap
        cutout = data[box[3]:box[4], box[1]:box[2]]
    end
    ny, nx = size(c)
    
    if partial_overlap || size(cutout) != size(c)
        slices_large, slices_small = overlap_slices(c, size(data))

        # no overlap
        slices_small === nothing && return nothing
        
        cutout = zeros(T, ny, nx)
        cutout[slices_small...] .= data[slices_large...]
    end
    return cutout
end

"""
    *(::Aperture, data::AbstractMatrix)

Return the data weighted by the aperture. If there is no overlap, will return an empty Array
"""
function Base.:*(c::Aperture, data::AbstractMatrix{T}) where {T}
    cut = cutout(c, data)
    cut === nothing && return typeof(data)(undef, 0, 0)
    return cut .* mask(c)
end

Base.:*(data::AbstractMatrix, c::Aperture) = c * data

"""
    aperture_photometry(::Aperture, data::AbstractMatrix)
    aperture_photometry(::AbstractVector{<:Aperture}, data::AbstractMatrix)

Perform aperture photometry on `data` given aperture(s). If a list of apertures is provided the output will be a `DataFrame`, otherwise a `NamedTuple`.

"""
function aperture_photometry(a::Aperture, data::AbstractMatrix)
    data_weighted = data * a
    aperture_sum = sum(data_weighted)
    return (xcenter = a.x, ycenter = a.y, aperture_sum = aperture_sum)
end

aperture_photometry(a::AbstractVector{<:Aperture}, data::AbstractMatrix) = DataFrame(aperture_photometry.(a, Ref(data)))

include("overlap.jl")
include("circular.jl")

end
