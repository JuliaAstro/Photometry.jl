#= 
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#


"""
    RectangularAperture(x, y, w, h, θ)
    RectangularAperture(position, w, h, θ)

A rectangular aperture.

A rectangular aperture with width `w`, height `h`, and position angle `θ` in degrees.

# Examples
```jldoctest
julia> ap = RectangularAperture(0, 0, 10, 4, 0)
RectangularAperture(0, 0, w=10, h=4, θ=0°)
```
"""
struct RectangularAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    w::T
    h::T
    theta::T
    function RectangularAperture(x::T, y::T, w::T, h::T, θ::T) where T <: Number
        (w < 0  || h < 0) && error("Invalid side lengths w=$w, h=$h. Both must be positive.")
        new{T}(x, y, w, h, mod(θ, 360))
    end
end

RectangularAperture(center::AbstractVector, w, h, θ) = RectangularAperture(center..., w, h, θ)
RectangularAperture(x, y,  w, h, θ) = RectangularAperture(promote(x, y, w, h, θ)...)

function Base.show(io::IO, ap::RectangularAperture)
    print(io, "RectangularAperture($(ap.x), $(ap.y), w=$(ap.w), h=$(ap.h), θ=$(ap.theta)°)")
end

function bbox(ap::RectangularAperture{T}) where T
    w2 = ap.w / 2
    h2 = ap.h / 2
    sint, cost = sincos(deg2rad(ap.theta))
    
    dx1 = abs(w2 * cost - h2 * sint)
    dy1 = abs(w2 * sint + h2 * cost)
    dx2 = abs(w2 * cost + h2 * sint)
    dy2 = abs(w2 * sint - h2 * cost)
    
    dx = max(dx1, dx2)
    dy = max(dy1, dy2)
    
    xmin = ceil(Int, ap.x - dx - 0.5)
    ymin = ceil(Int, ap.y - dy - 0.5)
    xmax = ceil(Int, ap.x + dx - 0.5)
    ymax = ceil(Int, ap.y + dy - 0.5)
    return xmin, xmax, ymin, ymax
end

function mask(ap::RectangularAperture; method = :exact)
    bounds = edges(ap)
    ny, nx = size(ap)
    return rectangular_overlap(bounds..., nx, ny, ap.w, ap.h, ap.theta, method = method)
end

#######################################################

"""
    RectangularAnnulus(x, y, w_in, w_out, h_out, θ)
    RectangularAnnulus(position, w_in, w_out, h_out, θ)

A rectangular annulus with inner width `w_in`, outer width `w_out`, outer height `h_out`, and position angle `θ` in degrees. `h_in` is automatically calculated from `w_in / w_out * h_out`. Note that `w_out ≥ w_in > 0`.

# Examples
```jldoctest
julia> ap = RectangularAnnulus(0, 0, 5, 10, 8, 45)
RectangularAnnulus(0.0, 0.0, w_in=5.0, w_out=10.0, h_in=4.0, h_out=8.0, θ=45.0°)
```

!!! warning
    The `:exact` method is not implemented for `RectangularAnnulus`
"""
struct RectangularAnnulus{T <: Number} <: AbstractAperture
    x::T
    y::T
    w_in::T
    w_out::T
    h_in::T
    h_out::T
    theta::T

    function RectangularAnnulus(x::T, y::T, w_in::T, w_out::T, h_out::T, θ::T) where T <: Number
        0 < w_in ≤ w_out || error("Invalid sides ($w_in, $w_out). `w_out` must be greater than or equal to `w_in` which must be greater than 0.")
        h_out ≤ 0 && error("Invalid side($h_out). `h_out` must be greater than 0.")
        P = typeof(one(T) / 1)
        new{P}(x, y, w_in, w_out, w_in / w_out * h_out, h_out, mod(θ, 360))
    end
end

RectangularAnnulus(center::AbstractVector, w_in, w_out, h_out, θ) = RectangularAnnulus(center..., w_in, w_out, h_out, θ)
RectangularAnnulus(x, y, w_in, w_out, h_out, θ) = RectangularAnnulus(promote(x, y, w_in, w_out, h_out, θ)...)

function Base.show(io::IO, ap::RectangularAnnulus)
    print(io, "RectangularAnnulus($(ap.x), $(ap.y), w_in=$(ap.w_in), w_out=$(ap.w_out), h_in=$(ap.h_in), h_out=$(ap.h_out), θ=$(ap.theta)°)")
end

function bbox(ap::RectangularAnnulus)
    w2 = ap.w_out / 2
    h2 = ap.h_out / 2
    sint, cost = sincos(deg2rad(ap.theta))

    dx1 = abs(w2 * cost - h2 * sint)
    dy1 = abs(w2 * sint + h2 * cost)
    dx2 = abs(w2 * cost + h2 * sint)
    dy2 = abs(w2 * sint - h2 * cost)

    dx = max(dx1, dx2)
    dy = max(dy1, dy2)

    xmin = ceil(Int, ap.x - dx - 0.5)
    ymin = ceil(Int, ap.y - dy - 0.5)
    xmax = ceil(Int, ap.x + dx - 0.5)
    ymax = ceil(Int, ap.y + dy - 0.5)
    return (xmin, xmax, ymin, ymax)
end

function mask(ap::RectangularAnnulus; method = :exact)
    bounds = edges(ap)
    ny, nx = size(ap)
    out = rectangular_overlap(bounds..., nx, ny, ap.w_out, ap.h_out, ap.theta, method = method)
    out .-= rectangular_overlap(bounds..., nx, ny,  ap.w_in, ap.h_in, ap.theta, method = method)
end
