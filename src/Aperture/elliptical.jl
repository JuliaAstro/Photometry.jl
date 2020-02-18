export EllipticalAperture,
       EllipticalAnnulus

"""
    EllipticalAperture(x, y, a, b, theta)

This is an elliptical aperture. The structure stores the basic elements required to define an elliptical aperture.
The elements being x, y, a, b, theta where they correspond to x and y-cordinates
of center of ellipse, length of semi-major axis of ellipse, length of semi-minor axis of ellipse and angle of rotation from
positive x-axis in counter-clockwise sense respectively, here a >= b > 0 and theta is in degrees.

# Examples
```jldoctest
julia> EllipticalAperture(2,2,2,2,1.2)
EllipticalAperture(2.0, 2.0, a=2.0, b=2.0, theta=1.2°)

```
"""
struct EllipticalAperture{T <: Number} <: AbstractAperture
    x::T
    y::T
    a::T
    b::T
    theta::T
end

function EllipticalAperture(x, y, a, b, theta)
    if (b < 0.0 || a < b )
            error("Invalid ellipse parameters. Require a >= b > 0.0")
    end
    theta = mod(theta, 360)
    return EllipticalAperture(promote(x, y, a, b, theta)...)
end

EllipticalAperture(center::AbstractVector, a, b, theta) = EllipticalAperture(center..., a, b, theta)

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), theta=$(e.theta)°)")
end

function oblique_coefficients(a, b, theta)
    sintheta, costheta = sincos(deg2rad(theta))
    a2 = a^2
    b2 = b^2
    cxx = costheta^2 / a2 + sintheta^2 / b2
    cyy = sintheta^2 / a2 + costheta^2 / b2
    cxy = 2*costheta*sintheta * (1/a2 - 1/b2)
    return cxx, cyy, cxy
end

function bbox(e::EllipticalAperture)

    sintheta, costheta = sincos(deg2rad(e.theta))

    t = atan((-e.b*tand(e.theta))/e.a)

    sint, cost = sincos(t)
    xmin = e.x + e.a*cost*costheta - e.b*sint*sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n*pi)
        xmin = min(xmin, e.x + e.a*cost*costheta - e.b*sint*sintheta)
    end

    for n in -2:2
        sint, cost = sincos(t + n*pi)
        xmax = max(xmax, e.x + e.a*cost*costheta - e.b*sint*sintheta)
    end

    t = atan((e.b*cotd(e.theta))/e.a)

    sint, cost = sincos(t)
    ymin = e.y + e.b*sint*costheta + e.a*cost*sintheta
    ymax = ymin

    for n in -2:2
        sint, cost = sincos(t + n*pi)
        ymin = min(ymin, e.y + e.b*sint*costheta + e.a*cost*sintheta)
    end

    for n in -2:2
        sint, cost = sincos(t + n*pi)
        ymax = max(ymax, e.y + e.b*sint*costheta + e.a*cost*sintheta)
    end

    xmin = floor(Int, xmin)
    xmax = ceil(Int, xmax)
    ymin = floor(Int, ymin)
    ymax = ceil(Int, ymax)

    return xmin, xmax, ymin, ymax
end


function edges(e::EllipticalAperture)
    ibox = bbox(e)
    xmin = ibox[1] - e.x - 0.5
    xmax = ibox[2] - e.x + 0.5
    ymin = ibox[3] - e.y - 0.5
    ymax = ibox[4] - e.y + 0.5
    return (xmin, xmax, ymin, ymax)
end

function mask(e::EllipticalAperture; method = :center)
    bounds = edges(e)
    ny, nx = size(e)
    return elliptical_overlap(bounds..., nx, ny, e.a, e.b, e.theta, method = method)
end
