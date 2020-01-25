module AperturePhotometry

export area,
       mask,
       cutout,
       aperture_photometry

abstract type Aperture end

"""
    area(::Aperture)

Returns the geometric area of the aperture
"""
area(::Aperture)

"""
    mask(::Aperture)

Return a `BitArray` of the masked aperture values in the minimum bounding box
"""
mask(::Aperture)

"""
    bbox(::Aperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture
"""
bbox(::Aperture)


function overlap_slices(c::Aperture, shape::Tuple)
    xmin, ymin, xmax, ymax = round.(Int, bbox(c))

    # No overlap
    if xmin ≥ shape[2] || ymin ≥ shape[1] || xmax ≤ 0 || ymax ≤ 0
        return nothing, nothing
    end

    slices_large = (max(ymin, 0) + 1:min(ymax, shape[1]) - 1, max(xmin, 0) + 1:min(xmax, shape[2]) - 1)
    slices_small = (max(-ymin, 0) + 1:min(ymax - ymin, shape[1] - ymin) - 1, max(-xmin, 0) + 1:min(xmax - xmin, shape[2] - xmin) - 1)

    return slices_large, slices_small
end

function cutout(c::Aperture, data::AbstractMatrix{T}) where {T}
    box = bbox(c)
    partial_overlap = box[1] < 0 || box[2] < 0

    if !partial_overlap
        cutout = data[box[2]:box[4], box[1]:box[3]]
    end

    if partial_overlap || size(cutout) != (box[4] - box[2], box[3] - box[1])
        slices_large, slices_small = overlap_slices(c, size(data))
        # no overlap
        slices_small === nothing && return T[]
        
        cutout = zeros(T, Int(box[4] - box[2]), Int(box[3] - box[1]))
        cutout[slices_small...] .= data[slices_large...]

    end
    return cutout
end

function Base.:*(c::Aperture, data::AbstractMatrix{T}) where {T}
    cut = cutout(c, data)
    cut === nothing && return T[]
    return cut .* mask(c)
end

Base.:*(data::AbstractMatrix, c::Aperture) = c * data

function aperture_photometry(a::Aperture, data::AbstractMatrix)
    data_weighted = data * a
    aperture_sum = sum(data_weighted)
    return (xcenter = a.x, ycenter = a.y, aperture_sum = aperture_sum)
end

include("overlap.jl")
include("circular.jl")

end
