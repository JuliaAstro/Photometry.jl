
export CircularAperture,
       CircularAnnulus


"""
    CircularAperture(x, y, r)
    CircularAperture([x, y], r)

A circular aperture.

# Examples
```jldoctest
julia> ap = CircularAperture(0, 0, 10)
CircularAperture(0, 0, r=10)

```
"""
struct CircularAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    r::T
end

CircularAperture(center::AbstractVector, r) = CircularAperture(center..., r)
CircularAperture(x, y, r) = CircularAperture(promote(x, y, r)...)

function Base.show(io::IO, c::CircularAperture)
    print(io, "CircularAperture($(c.x), $(c.y), r=$(c.r))")
end

bbox(c::CircularAperture{<:Integer}) = (c.x - c.r, c.x + c.r, c.y - c.r, c.y + c.r)

function bbox(c::CircularAperture{<:AbstractFloat})
    xmin = floor(Int, c.x - c.r + 0.5)
    xmax = ceil(Int, c.x + c.r + 0.5)
    ymin = floor(Int, c.y - c.r + 0.5)
    ymax = ceil(Int, c.y + c.r + 0.5)
    return (xmin, xmax, ymin, ymax)
end

function mask(c::CircularAperture; method = :exact)
    bounds = edges(c)
    ny, nx = size(c)
    return circular_overlap(bounds..., nx, ny, c.r, method = method)
end

#######################################################

"""
    CircularAnnulus(x, y, r_in, r_out)
    CircularAnnulus([x, y], r_in, r_out)

A circular aperture.

# Examples
```jldoctest
julia> ap = CircularAnnulus(0, 0, 5, 10)
CircularAnnulus(0, 0, r_in=5, r_out=10)

```
"""
struct CircularAnnulus{T <: Number} <: AbstractAperture
    x::T
    y::T
    r_in::T
    r_out::T
end

CircularAnnulus(center::AbstractVector, r_in, r_out) = CircularAnnulus(center..., r_in, r_out)
CircularAnnulus(x, y, r_in, r_out) = CircularAnnulus(promote(x, y, r_in, r_out)...)

function Base.show(io::IO, c::CircularAnnulus)
    print(io, "CircularAnnulus($(c.x), $(c.y), r_in=$(c.r_in), r_out=$(c.r_out))")
end

bbox(c::CircularAnnulus{<:Integer}) = (c.x - c.r_out, c.x + c.r_out, c.y - c.r_out, c.y + c.r_out)

function bbox(c::CircularAnnulus{<:AbstractFloat})
    xmin = floor(Int, c.x - c.r_out + 0.5)
    xmax = ceil(Int, c.x + c.r_out + 0.5)
    ymin = floor(Int, c.y - c.r_out + 0.5)
    ymax = ceil(Int, c.y + c.r_out + 0.5)
    return (xmin, xmax, ymin, ymax)
end

function mask(c::CircularAnnulus; method = :exact)
    bounds = edges(c)
    ny, nx = size(c)
    out = circular_overlap(bounds..., nx, ny, c.r_out, method = method)
    out .-= circular_overlap(bounds..., nx, ny, c.r_in, method = method)
end

function edges(c::Union{CircularAperture,CircularAnnulus})
    ibox = bbox(c)
    xmin = ibox[1] - c.x - 0.5
    xmax = ibox[2] - c.x + 0.5
    ymin = ibox[3] - c.y - 0.5
    ymax = ibox[4] - c.y + 0.5
    return (xmin, xmax, ymin, ymax)
end
