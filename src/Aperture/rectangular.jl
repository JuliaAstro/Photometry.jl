export RectangularAperture
# RectangularAnnulus


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

RectangularAperture(center::AbstractVector, w, h) = RectangularAperture(center..., w, h)
RectangularAperture(x, y, w, h) = RectangularAperture(promote(x, y, w, h)...)

function Base.show(io::IO, c::RectangularAperture)
    print(io, "RectangularAperture($(c.x), $(c.y), w=$(c.w), h=$(c.h))")
end

bbox(c::RectangularAperture{<:Integer}) = (c.x - c.(w/2), c.x + c.(w/2), c.y - c.(h/2), c.y + c.(h/2))

function bbox(c::RectangularAperture{<:AbstractFloat})
    xmin = floor(Int, c.x - c.(w/2))
    xmax = ceil(Int, c.x + c.(w/2))
    ymin = floor(Int, c.y - c.(h/2))
    ymax = ceil(Int, c.y + c.(h/2))
    return (xmin, xmax, ymin, ymax)
end

function mask(c::RectangularAperture; method = :exact)
    bounds = edges(c)
    box = bbox(c)
    ny, nx = size(c)
    return Rectangular_overlap(bounds..., nx, ny, c.r, method = method)
end

# #######################################################

# """
#     CircularAnnulus(x, y, r_in, r_out)
#     CircularAnnulus([x, y], r_in, r_out)

# A circular aperture.

# # Examples
# ```jldoctest
# julia> ap = CircularAnnulus(0, 0, 5, 10)
# CircularAnnulus(0, 0, r_in=5, r_out=10)

# ```
# """
# struct CircularAnnulus{T <: Number} <: AbstractAperture
#     x::T
#     y::T
#     r_in::T
#     r_out::T
# end

# CircularAnnulus(center::AbstractVector, r_in, r_out) = CircularAnnulus(center..., r_in, r_out)
# CircularAnnulus(x, y, r_in, r_out) = CircularAnnulus(promote(x, y, r_in, r_out)...)

# function Base.show(io::IO, c::CircularAnnulus)
#     print(io, "CircularAnnulus($(c.x), $(c.y), r_in=$(c.r_in), r_out=$(c.r_out))")
# end

# bbox(c::CircularAnnulus{<:Integer}) = (c.x - c.r_out, c.x + c.r_out, c.y - c.r_out, c.y + c.r_out)

# function bbox(c::CircularAnnulus{<:AbstractFloat})
#     xmin = floor(Int, c.x - c.r_out + 0.5)
#     xmax = ceil(Int, c.x + c.r_out + 0.5)
#     ymin = floor(Int, c.y - c.r_out + 0.5)
#     ymax = ceil(Int, c.y + c.r_out + 0.5)
#     return (xmin, xmax, ymin, ymax)
# end

# function mask(c::CircularAnnulus; method = :exact)
#     bounds = edges(c)
#     box = bbox(c)
#     ny, nx = size(c)
#     out = circular_overlap(bounds..., nx, ny, c.r_out, method = method)
#     out .-= circular_overlap(bounds..., nx, ny, c.r_in, method = method)
# end

# function edges(c::Union{CircularAperture,CircularAnnulus})
#     ibox = bbox(c)
#     xmin = ibox[1] - c.x - 0.5
#     xmax = ibox[2] - c.x + 0.5
#     ymin = ibox[3] - c.y - 0.5
#     ymax = ibox[4] - c.y + 0.5
#     return (xmin, xmax, ymin, ymax)
# end
