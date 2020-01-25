
export CircularAperture,
       CircularAnnulus

struct CircularAperture{T <: Number} <: Aperture
    center::Tuple{T,T}
    r::T
end

function Base.show(io::IO, c::CircularAperture)
    x, y = c.center
    print(io, "CircularAperture($x, $y, r=$(c.r))")
end

area(c::CircularAperture) = π * c.r^2

struct CircularAnnulus{T <: Number} <: Aperture
    center::Tuple{T,T}
    r_in::T
    r_out::T
end

function Base.show(io::IO, c::CircularAnnulus)
    x, y = c.center
    print(io, "CircularAnnulus($x, $y, r_in=$(c.r_in), r_out=$(c.r_out))")
end

area(c::CircularAnnulus) = π * (c.r_out^2 - c.r_in^2)
