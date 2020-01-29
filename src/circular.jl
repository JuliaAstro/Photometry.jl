
export CircularAperture,
       CircularAnnulus

struct CircularAperture{T <: Number} <: Aperture
    x::T
    y::T
    r::T
end

CircularAperture(center::AbstractVector, r) = CircularAperture(center..., r)

function Base.show(io::IO, c::CircularAperture)
    print(io, "CircularAperture($(c.x), $(c.y), r=$(c.r))")
end

bbox(c::CircularAperture) = (c.x - c.r, c.y - c.r, c.x + c.r, c.y + c.r)

function mask(c::CircularAperture; method = :exact)
    bounds = edges(c)
    box = bbox(c)
    nx = Int(box[4] - box[2])
    ny = Int(box[3] - box[1])
    return circular_overlap(bounds..., nx, ny, c.r, method = method)
end

#######################################################

struct CircularAnnulus{T <: Number} <: Aperture
    x::T
    y::T
    r_in::T
    r_out::T
end

CircularAnnulus(center::AbstractVector, r_in, r_out) = CircularAnnulus(center..., r_in, r_out)

function Base.show(io::IO, c::CircularAnnulus)
    print(io, "CircularAnnulus($(c.x), $(c.y), r_in=$(c.r_in), r_out=$(c.r_out))")
end

bbox(c::CircularAnnulus) = (c.x - c.r_out, c.y - c.r_out, c.x + c.r_out, c.y + c.r_out)

function mask(c::CircularAnnulus; method = :exact)
    bounds = edges(c)
    box = bbox(c)
    nx = box[3] - box[1]
    ny = box[4] - box[2]
    out = circular_overlap(bounds..., nx, ny, c.r_out, method = method)
    out .-= circular_overlap(bounds..., nx, ny, c.r_in, method = method)
end

function edges(c::Union{CircularAperture,CircularAnnulus})
    ibox = bbox(c) .- 0.5
    xmin = ibox[1] - c.x
    ymin = ibox[2] - c.y
    xmax = ibox[3] - c.x
    ymax = ibox[4] - c.y
    return (xmin, ymin, xmax, ymax)
end
