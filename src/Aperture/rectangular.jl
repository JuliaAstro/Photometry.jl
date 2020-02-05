export RectangularAperture,
       RectangularAnnulus

 
"""
    RectangularAperture(x, y, w, h)
    RectangularAperture([x, y], w, h)

A rectangular aperture.

# Examples
```jldoctest
julia> ap = RectangularAperture(0, 0, 10, 5)
RectanguRelarAperture(0, 0, w=10, h=5)

```
"""

struct RectangularAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    w::T
    h::T
end

# RectangularAperture(center::AbstractVector, w, h) = RectangularAperture(center..., w, h)
function RectangularAperture(x, y, w, h) 
    if (w < 0 || h < 0)
        error("Invalid parameters" * "Require w >= 0.0 , h >= 0.0")
    end
    
    half_width = w / 2
    half_height = h / 2


    return RectangularAperture(promote(x, y, w, h)...)
end

function Base.show(io::IO, rect::RectangularAperture)
    print(io, "RectangularAperture($(rect.x), $(rect.y), w=$(rect.w), h=$(rect.h))")
end

# half_width = w / 2
# half_width = h / 2

bbox(rect::RectangularAperture{<:Integer}) = (rect.x - rect.half_width, rect.x + rect.half_width, rect.y - rect.half_height, rect.y + rect.half_height)

function bbox(rect::RectangularAperture{<:AbstractFloat})
    xmin = floor(Int, rect.x - rect.half_width)
    xmax = ceil(Int, rect.x + rect.half_width)
    ymin = floor(Int, rect.y - rect.half_height)
    ymax = ceil(Int, rect.y + rect.half_height)
    return (xmin, xmax, ymin, ymax)
end

function mask(rect::RectangularAperture, theta = 0, subpixels = 5; method = :exact)
    bounds = edges(rect)
    box = bbox(rect)
    ny, nx = size(rect)
    return rectangular_overlap(bounds..., nx, ny, rect.w, rect.h, theta, subpixels, method = method) 
end

######################################################

"""
    RectangularAnnulus(x, y, w_in, w_out, h_in, h_out)
    RectangularAnnulus([x, y], w_in, w_out, h_in, h_out)

A rectangular aperture.

# Examples
```jldoctest
julia> ap = RectangularAnnulus(0, 0, 5, 10, 6, 12)
RectangularAnnulus(0, 0, w_in=5, w_out=10, h_in=6, h_out=12)

```
"""
struct RectangularAnnulus{T <: Number} <: AbstractAperture
    x::T
    y::T
    w_in::T
    w_out::T
    h_in::T
    h_out::T

end

# RectangularAnnulus(center::AbstractVector, w_in, w_out, h_in, h_out) = RectangularAnnulus(center..., w_in, w_out, h_in, h_out)
function RectangularAnnulus(x, y, w_in, w_out, h_in, h_out)
    if(w_out <= w_in || h_out <= h_in)
        error("Invalid parameters" * "Require w_out > w_in, h_out > h_in")

        half_width_out = w_out / 2
        half_height_out = h_out / 2
    
    end
    return RectangularAnnulus(promote(x, y, w_in, w_out, h_in, h_out)...)
end


function Base.show(io::IO, rect::RectangularAnnulus)
    print(io, "RectangularAnnulus($(rect.x), $(rect.y), w_in=$(rect.w_in), w_out=$(rect.w_out), h_out=$(rect.h_out))")
end


bbox(rect::RectangularAnnulus{<:Integer}) = (rect.x - rect.half_width_out, rect.x + rect.half_width_out, rect.y - rect.half_height, rect.y + rect.half_height)


function bbox(rect::RectangularAnnulus{<:AbstractFloat})
    xmin = floor(Int, rect.x - rect.half_width_out)
    xmax = ceil(Int, rect.x + rect.half_width_out)
    ymin = floor(Int, rect.y - rect.half_height_out)
    ymax = ceil(Int, rect.y + rect.half_height_out)
    return (xmin, xmax, ymin, ymax)
end

function mask(rect::RectangularAnnulus, theta = 0, subpixels = 5; method = :exact)
    bounds = edges(rect)
    box = bbox(rect)
    ny, nx = size(rect)
    out = rectangular_overlap(bounds..., nx, ny, rect.w_out, rect.h_out, theta, subpixels, method = method)
    out .-= rectangular_overlap(bounds..., nx, ny,  rect.w_out, rect.h_out, theta, subpixels, method = method)
end

function edges(rect::Union{RectangularAperture,RectangularAnnulus})
    ibox = bbox(rect)
    xmin = ibox[1] - rect.w - 0.5
    xmax = ibox[2] - rect.w + 0.5
    ymin = ibox[3] - rect.h - 0.5
    ymax = ibox[4] - rect.h + 0.5
    return (xmin, xmax, ymin, ymax)
end
