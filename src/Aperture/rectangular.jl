export RectangularAperture,
       RectangularAnnulus

 
"""
    RectangularAperture(x, y, w, h, theta)
    RectangularAperture([x, y], w, h, theta)

A rectangular aperture.

w : float
    The full width of the rectangle in pixels. For theta=0 the width side is along the x axis.

h : float
    The full height of the rectangle in pixels. For theta=0 the height side is along the y axis.
    
theta : float (optional)
    The rotation angle in degrees of the rectangle “width” side from the positive x axis. 
    The rotation angle increases counterclockwise. The default is 0.

Raises Error
    If either width (w) or height (h) is negative.

# Examples
```jldoctest
julia> ap = RectangularAperture(0, 0, 10, 5, 50)
RectanguRelarAperture(0, 0, w=10, h=5, theta=50)

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
        error("Invalid parameters. Require w >= 0.0 , h >= 0.0")
    end
    
    return RectangularAperture(promote(x, y, w, h, theta)...)
end

function Base.show(io::IO, rect::RectangularAperture)
    print(io, "RectangularAperture($(rect.x), $(rect.y), w=$(rect.w), h=$(rect.h), theta=$(rect.theta))")
end


function bbox(rect::RectangularAperture{<:Number})
    
    half_width = rect.w / 2
    half_height = rect.h / 2

    xmin = floor(Int, rect.x - half_width)
    xmax = ceil(Int, rect.x + half_width)
    ymin = floor(Int, rect.y - half_height)
    ymax = ceil(Int, rect.y + half_height)
    return (xmin, xmax, ymin, ymax)
end

function mask(rect::RectangularAperture; method = :exact)
    bounds = edges(rect)
    box = bbox(rect)
    ny, nx = size(rect)
    return rectangular_overlap(bounds..., nx, ny, rect.w, rect.h, rect.theta, method = method) 
end

######################################################

"""
    RectangularAnnulus(x, y, w_in, w_out, h_in, h_out, theta)
    RectangularAnnulus([x, y], w_in, w_out, h_in, h_out, theta)

A rectangular annulus.

w_in : float
    The inner full width of the rectangular annulus in pixels. For theta=0 the width side is along the x axis.

w_out : float
    The outer full width of the rectangular annulus in pixels. For theta=0 the width side is along the x axis.

h_out : float
    The outer full height of the rectangular annulus in pixels.
    
    The "inner full height" (h_in) is calculated as:

        h_in = h_out * (w_in / w_out)

    For theta=0 the height side is along the y axis.

theta : float (optional)
    The rotation angle in degrees of the rectangle “width” side from the positive x axis.
    The rotation angle increases counterclockwise. The default is 0.

Raises Error
        If inner width (w_in) is greater than outer width (w_out).
        If either the inner width (w_in) or the outer height (h_out) is negative.


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
    if(w_out <= w_in || w_in < 0 || h_out < 0)
        error("Invalid parameters. Require w_out > w_in, h_in > 0, h_out > 0")

        
    end
    return RectangularAnnulus(promote(x, y, w_in, w_out, h_in, h_out, theta)...)
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

function mask(rect::RectangularAnnulus; method = :exact)
    bounds = edges(rect)
    box = bbox(rect)
    ny, nx = size(rect)
    out = rectangular_overlap(bounds..., nx, ny, rect.w_out, rect.h_out, rect.theta, method = method)
    out .-= rectangular_overlap(bounds..., nx, ny, rect.w_in, rect.h_in, rect.theta, method = method)
end

function edges(rect::Union{RectangularAperture,RectangularAnnulus})
    ibox = bbox(rect)
    xmin = ibox[1] - rect.x
    xmax = ibox[2] - rect.x
    ymin = ibox[3] - rect.y
    ymax = ibox[4] - rect.y
    return (xmin, xmax, ymin, ymax)
end

