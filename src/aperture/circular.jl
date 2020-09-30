#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

###############

"""
    CircularAperture(x, y, r)
    CircularAperture(position, r)

A circular aperture.

A circular aperture with radius `r`. `r` must be greater than or equal to 0.

# Examples
```jldoctest
julia> ap = CircularAperture(0, 0, 10)
CircularAperture(0, 0, r=10)
```
"""
struct CircularAperture{T<:Number} <: AbstractAperture
    x::T
    y::T
    r::T
end

function overlap(ap::CircularAperture, i, j)
    dist = sqrt((i - ap.y)^2 + (j - ap.x)^2)
    dr = sqrt(2) / 2 # corner-center distance of pixel
    dist > ap.r + dr && return Outside
    dist < ap.r - dr && return Inside
    return Partial
end

partial(ap::CircularAperture, x, y) = circular_overlap_single_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r)
partial(ap::Subpixel{<:CircularAperture}, x, y) = circular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r, ap.N)


CircularAperture(center, r) = CircularAperture(center..., r)
CircularAperture(x, y, r) = CircularAperture(promote(x, y, r)...)

function Base.show(io::IO, c::CircularAperture)
    print(io, "CircularAperture($(c.x), $(c.y), r=$(c.r))")
end


function bounds(c::CircularAperture)
    xmin = ceil(Int, c.x - c.r - 0.5)
    ymin = ceil(Int, c.y - c.r - 0.5)
    xmax = ceil(Int, c.x + c.r - 0.5)
    ymax = ceil(Int, c.y + c.r - 0.5)
    return (xmin, xmax, ymin, ymax)
end

#######################################################


"""
    CircularAnnulus(x, y, r_in, r_out)
    CircularAnnulus(position, r_in, r_out)

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
end


CircularAnnulus(center, r_in, r_out) = CircularAnnulus(center..., r_in, r_out)
CircularAnnulus(x, y, r_in, r_out) = CircularAnnulus(promote(x, y, r_in, r_out)...)

function Base.show(io::IO, c::CircularAnnulus)
    print(io, "CircularAnnulus($(c.x), $(c.y), r_in=$(c.r_in), r_out=$(c.r_out))")
end


function overlap(ap::CircularAnnulus, i, j)
    dist = sqrt((i - ap.y)^2 + (j - ap.x)^2)
    dr = sqrt(2) / 2 # corner-center distance of pixel
    ap.r_in - dr < dist < ap.r_out + dr || return Outside
    ap.r_in + dr < dist < ap.r_out - dr && return Inside
    return Partial
end

function partial(ap::CircularAnnulus, x, y)
    f1 = circular_overlap_single_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_out)
    f2 = circular_overlap_single_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_in)
    return f1 - f2
end

function partial(ap::Subpixel{<:CircularAnnulus}, x, y)
    f1 = circular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_out, ap.N)
    f2 = circular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_in, ap.N)
    return f1 - f2
end

function bounds(c::CircularAnnulus)
    xmin = ceil(Int, c.x - c.r_out - 0.5)
    ymin = ceil(Int, c.y - c.r_out - 0.5)
    xmax = ceil(Int, c.x + c.r_out - 0.5)
    ymax = ceil(Int, c.y + c.r_out - 0.5)
    return (xmin, xmax, ymin, ymax)
end
