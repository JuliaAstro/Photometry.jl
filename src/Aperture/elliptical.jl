export EllipticalAperture,
       EllipticalAnnulus

"""
   EllipticalAperture(x, y, a, b, theta)

A elliptical aperture.

Where a >= b > 0 and theta in degrees

# Examples
```jldoctest
julia> EllipticalAperture(2,2,2,2,1.2)
EllipticalAperture(2.0, 2.0, a=2.0, b=2.0, theta=1.2)

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
            error("illegal ellipse parameters. Require a >= b > 0.0")
    end
    theta = mod(theta, 360)
    return EllipticalAperture(promote(x, y, a, b, theta)...)
end

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), theta=$(e.theta)°)")
end

function oblique_coefficients(a, b, theta)
    sintheta, costheta = sincos(deg2rad(theta))
    asq = a*a
    bsq = b*b
    cxx = costheta*costheta*inv(asq) + sintheta*sintheta*inv(bsq)
    cyy = sintheta*sintheta*inv(asq) + costheta*costheta*inv(bsq)
    cxy = 2*costheta*sintheta * (inv(asq) - inv(bsq))
    return cxx, cyy, cxy
end

function bbox(e::EllipticalAperture)

    sintheta, costheta = sincos(deg2rad(e.theta))

    t = atan((-e.b*tan(deg2rad(e.theta)))/e.a)

    xmin = e.x + e.a*cos(t)*costheta - e.b*sin(t)*sintheta
    xmax = xmin

    for n in -2:2
        xmin = min(xmin, e.x + e.a*cos(t + n*pi)*costheta - e.b*sin(t + n*pi)*sintheta)
    end

    for n in -2:2
        xmax = max(xmax, e.x + e.a*cos(t + n*pi)*costheta - e.b*sin(t + n*pi)*sintheta)
    end

    t = atan((e.b*cot(deg2rad(e.theta)))/e.a)

    ymin = e.y + e.b*sin(t)*costheta + e.a*cos(t)*sintheta
    ymax = ymin

    for n in -2:2
        ymin = min(ymin, e.y + e.b*sin(t + n*pi)*costheta + e.a*cos(t + n*pi)*sintheta)
    end

    for n in -2:2
        ymax = max(ymax, e.y + e.b*sin(t + n*pi)*costheta + e.a*cos(t + n*pi)*sintheta)
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

function mask(e::EllipticalAperture; method = :exact)
    bounds = edges(e)
    box = bbox(e)
    ny, nx = size(e)
#     return elliptical_overlap(bounds..., nx, ny, rect.w, rect.h, rect.theta, method = method)
end
