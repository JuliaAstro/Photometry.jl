export RectangularAperture,
       RectangularAnnulus

 
"""
    RectangularAperture(x, y, w, h, theta)
    RectangularAperture([x, y], w, h, theta)

A rectangular aperture.

# Examples
```jldoctest
julia> ap = RectangularAperture(0, 0, 10, 5, 50)
RectanguRelarAperture(0, 0, w=10, h=5)

```
"""

struct RectangularAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    w::T
    h::T
    theta::T
end

function RectangularAperture(x, y, w, h, theta = 0.0) 
    if (w < 0 || h < 0)
        error("Invalid parameters" * "Require w >= 0.0 , h >= 0.0")
    end
    
    return RectangularAperture(promote(x, y, w, h, theta = 0.0)...)
end

function Base.show(io::IO, rect::RectangularAperture)
    print(io, "RectangularAperture($(rect.x), $(rect.y), w=$(rect.w), h=$(rect.h), theta=$(rect.theta))")
end


function bbox(rect::RectangularAperture{<:Number})
    
    w = rect.w
    h = rect.h

    half_width = w / 2
    half_height = h / 2

    xmin = floor(Int, rect.x - half_width)
    xmax = ceil(Int, rect.x + half_width)
    ymin = floor(Int, rect.y - half_height)
    ymax = ceil(Int, rect.y + half_height)
    return (xmin, xmax, ymin, ymax)
end

# function mask(rect::RectangularAperture, theta = 0, subpixels = 5; method = :exact)
#     bounds = edges(rect)
#     box = bbox(rect)
#     ny, nx = size(rect)
#     return rectangular_overlap(bounds..., nx, ny, rect.w, rect.h, theta, subpixels, method = method) 
# end

######################################################

"""
    RectangularAnnulus(x, y, w_in, w_out, h_in, h_out, theta)
    RectangularAnnulus([x, y], w_in, w_out, h_in, h_out, theta)

A rectangular annulus.

# Examples
```jldoctest
julia> ap = RectangularAnnulus(0, 0, 5, 10, 6, 12, 50)
RectangularAnnulus(0, 0, w_in=5, w_out=10, h_in=6, h_out=12, theta = 50)

```
"""
struct RectangularAnnulus{T <: Number} <: AbstractAperture
    x::T
    y::T
    w_in::T
    w_out::T
    h_in::T
    h_out::T
    theta::T
end

function RectangularAnnulus(x, y, w_in, w_out, h_in, h_out, theta = 0.0)
    if(w_out <= w_in || h_out <= h_in)
        error("Invalid parameters" * "Require w_out > w_in, h_out > h_in")

        
    end
    return RectangularAnnulus(promote(x, y, w_in, w_out, h_in, h_out)...)
end


function Base.show(io::IO, rect::RectangularAnnulus)
    print(io, "RectangularAnnulus($(rect.x), $(rect.y), w_in=$(rect.w_in), w_out=$(rect.w_out), h_in=$(rect.h_in), h_out=$(rect.h_out), theta=$(rect.theta))")
end




function bbox(rect::RectangularAnnulus{<:Number})

    
    half_width_out = rect.w_out / 2
    half_height_out = rect.h_out / 2
    
    xmin = floor(Int, rect.x - half_width_out)
    xmax = ceil(Int, rect.x + half_width_out)
    ymin = floor(Int, rect.y - half_height_out)
    ymax = ceil(Int, rect.y + half_height_out)
    return (xmin, xmax, ymin, ymax)
end

