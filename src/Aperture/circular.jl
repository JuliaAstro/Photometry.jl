
export CircularAperture,
       CircularAnnulus


"""
    CircularAperture(x, y, r)
    CircularAperture([x, y], r)

A circular aperture.

A circular aperture with radius `r`. `r` must be greater than or equal to 0. 

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

    function CircularAperture(x::T, y::T, r::T) where T <: Number
        r ≥ 0 || error("Invalid radius r=$r. r must be greater than or equal to 0.")
        new{T}(x, y, r)
    end
end

CircularAperture(center::AbstractVector, r) = CircularAperture(center..., r)
CircularAperture(x, y, r) = CircularAperture(promote(x, y, r)...)

function Base.show(io::IO, c::CircularAperture)
    print(io, "CircularAperture($(c.x), $(c.y), r=$(c.r))")
end

function bbox(c::CircularAperture)
    xmin = floor(Int, c.x - c.r)
    xmax = floor(Int, c.x + c.r)
    ymin = floor(Int, c.y - c.r)
    ymax = floor(Int, c.y + c.r)
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

A circular annulus with inner radius `r_in` and outer radius `r_out`. 0 ≤ `r_in` ≤ `r_out`.

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

    function CircularAnnulus(x::T, y::T, r_in::T, r_out::T) where T <: Number
        0 ≤ r_in ≤ r_out || error("Invalid radii ($r_in, $r_out). r_out must be greater than r_in which must be greater than or equal to 0.")
        new{T}(x, y, r_in, r_out)
    end
end

CircularAnnulus(center::AbstractVector, r_in, r_out) = CircularAnnulus(center..., r_in, r_out)
CircularAnnulus(x, y, r_in, r_out) = CircularAnnulus(promote(x, y, r_in, r_out)...)

function Base.show(io::IO, c::CircularAnnulus)
    print(io, "CircularAnnulus($(c.x), $(c.y), r_in=$(c.r_in), r_out=$(c.r_out))")
end

function bbox(c::CircularAnnulus)
    xmin = floor(Int, c.x - c.r_out)
    xmax = floor(Int, c.x + c.r_out)
    ymin = floor(Int, c.y - c.r_out)
    ymax = floor(Int, c.y + c.r_out)
    return (xmin, xmax, ymin, ymax)
end

function mask(c::CircularAnnulus; method = :exact)
    bounds = edges(c)
    ny, nx = size(c)
    out = circular_overlap(bounds..., nx, ny, c.r_out, method = method)
    out .-= circular_overlap(bounds..., nx, ny, c.r_in, method = method)
end
