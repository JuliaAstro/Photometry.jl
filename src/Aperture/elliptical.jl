#= 
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

export EllipticalAperture,
       EllipticalAnnulus

"""
    EllipticalAperture(x, y, a, b, θ)
    EllipticalAperture([x, y], a, b, θ)

An elliptical aperture with semi-major axis `a`, semi-minor axis `b`, and angle `θ`. `a` and `b` must be ≥ 0, `θ` is measured in degrees counter-clockwise the standard x-axis.

# Examples
```jldoctest
julia> ap = EllipticalAperture(0, 0, 4, 2, 35)
EllipticalAperture(0, 0, a=4, b=2, θ=35°)
```
"""
struct EllipticalAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    a::T
    b::T
    theta::T

    function EllipticalAperture(x::T, y::T, a::T, b::T, theta::T) where T <: Number
        a < 0 && error("Invalid axis a=$a. a must be greater than or equal to 0")
        b < 0 && error("Invalid axis b=$b. a must be greater than or equal to 0")
        new{T}(x, y, a, b, mod(theta, 360))
    end
end


EllipticalAperture(x, y, a, b, theta) = EllipticalAperture(promote(x, y, a, b, theta)...)
EllipticalAperture(center::AbstractVector, a, b, theta) = EllipticalAperture(center..., a, b, theta)

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), θ=$(e.theta)°)")
end

function oblique_coefficients(a, b, theta)
    sintheta, costheta = sincos(deg2rad(theta))
    a2 = a^2
    b2 = b^2
    cxx = costheta^2 / a2 + sintheta^2 / b2
    cyy = sintheta^2 / a2 + costheta^2 / b2
    cxy = 2 * costheta * sintheta * (1 / a2 - 1 / b2)
    return cxx, cyy, cxy
end

function bbox(e::EllipticalAperture)

    sintheta, costheta = sincos(deg2rad(e.theta))

    t = atan((-e.b * tand(e.theta)) / e.a)

    sint, cost = sincos(t)
    xmin = e.x + e.a * cost * costheta - e.b * sint * sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmin = min(xmin, e.x + e.a * cost * costheta - e.b * sint * sintheta)
    end

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmax = max(xmax, e.x + e.a * cost * costheta - e.b * sint * sintheta)
    end

    t = atan((e.b * cotd(e.theta)) / e.a)

    sint, cost = sincos(t)
    ymin = e.y + e.b * sint * costheta + e.a * cost * sintheta
    ymax = ymin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        ymin = min(ymin, e.y + e.b * sint * costheta + e.a * cost * sintheta)
    end

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        ymax = max(ymax, e.y + e.b * sint * costheta + e.a * cost * sintheta)
    end

    xmin = floor(Int, xmin)
    xmax = ceil(Int, xmax)
    ymin = floor(Int, ymin)
    ymax = ceil(Int, ymax)

    return xmin, xmax, ymin, ymax
end

function mask(e::EllipticalAperture; method = :center)
    bounds = edges(e)
    ny, nx = size(e)
    return elliptical_overlap(bounds..., nx, ny, e.a, e.b, e.theta, method = method)
end


#######################################################

"""
    EllipticalAnnulus(x, y, a_in, b_in, theta_in, a_out, b_out, theta_out)
    EllipticalAnnulus([x, y], a_in, b_in, theta_in, a_out, b_out, theta_out)

An elliptical annulus.

# Examples
```jldoctest
julia> ap = EllipticalAnnulus(0, 0, 10, 5, 25, 16, 10, 45)
EllipticalAnnulus(0, 0, a_in=10, b_in=5, theta_in=25°, a_out=16, b_out=10, theta_out=45°)

```
"""
struct EllipticalAnnulus{T <: Number} <: AbstractAperture
    x::T
    y::T
    a_in::T
    b_in::T
    theta_in::T
    a_out::T
    b_out::T
    theta_out::T
end

function EllipticalAnnulus(x, y, a_in, b_in, theta_in, a_out, b_out, theta_out)
    if (b_in < 0.0 || a_in < b_in || b_out < 0.0 || a_out < b_out)
            error("illegal ellipse parameters. Require a_in >= b_in > 0.0, a_out >= b_out > 0.0")
    end
    theta_in = mod(theta_in, 360)
    theta_out = mod(theta_out,360);
    return EllipticalAperture(promote(x, y, a_in, b_in, theta_in, a_out, b_out, theta_out)...)
end

function Base.show(io::IO, e::EllipticalAnnulus)
    print(io, "EllipticalAnnulus($(e.x), $(e.y), a_in=$(e.a_in), b_in=$(e.b_in), theta_in=$(e.theta_in)°, a_out=$(e.a_out), b_out=$(e.b_out), theta_out=$(e.theta_out)°)")
end
