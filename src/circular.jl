
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

bbox(c::CircularAperture) = (-c.r, c.r, -c.r, c.r)

function mask(c::CircularAperture)
    box = bbox(c)
    nx = ny = Int(2c.r)
    out = BitArray(false, nx, ny)
    xs = box[1]:box[2]
    ys = hcat(box[3]:box[4])
    out[sqrt.(xs.^2 .+ ys.^2) .< c.r] .= true
    return out
end

#######################################################

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
bbox(c::CircularAnnulus) = (-c.r, c.r, -c.r, c.r)

function mask(c::CircularAnnulus)
    box = bbox(c)
    nx = ny = Int(2c.r)
    out = BitArray(false, nx, ny)
    xs = box[1]:box[2]
    ys = hcat(box[3]:box[4])
    out[c.r_in .< sqrt.(xs.^2 .+ ys.^2) .< c.r_out] .= true
    return out
end
