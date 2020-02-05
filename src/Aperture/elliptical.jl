export EllipticalAperture,
       EllipticalAnnulus

"""
   EllipticalAperture(x, y, a, b, theta)

A elliptical aperture.

Where a >= b > 0 and -pi/2 <= theta <= pi/2

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
    cxx::T
    cxy::T
    cyy::T
end

function EllipticalAperture(x, y, a, b, theta)
    if (b < 0.0 || a < b || theta < -pi/2 || theta > pi/2)
            error("illegal ellipse parameters. " * "Require a >= b > 0.0, -pi/2 <= theta <= pi/2")
    end
    cxx, cyy, cxy = oblique_coefficients(a, b, theta)
    return EllipticalAperture(promote(x, y, a, b, theta, cxx, cxy, cyy)...)

end

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), theta=$(e.theta))")
end

function oblique_coefficients(a, b, theta)
    sintheta, costheta = sincos(theta)
    asq = a*a
    bsq = b*b
    cxx = costheta*costheta*inv(asq) + sintheta*sintheta*inv(bsq)
    cyy = sintheta*sintheta*inv(asq) + costheta*costheta*inv(bsq)
    cxy = 2*costheta*sintheta * (inv(asq) - inv(bsq))
    return cxx, cyy, cxy
end

function bbox(e::EllipticalAperture{<:AbstractFloat})

    x = e.x
    y = e.y
    a = e.a
    b = e.b
    theta = e.theta

    sintheta, costheta = sincos(theta)

    t = atan((-b*tan(theta))/a)

    xmin = x + a*cos(t)*costheta - b*sin(t)*sintheta
    xmax = xmin

    for n in -2:2
        xmin = min(xmin, x + a*cos(t + n*pi)*costheta - b*sin(t + n*pi)*sintheta)
    end

    for n in -2:2
        xmax = max(xmax, x + a*cos(t + n*pi)*costheta - b*sin(t + n*pi)*sintheta)
    end

    t = atan((b*cot(theta))/a)

    ymin = y + b*sin(t)*costheta + a*cos(t)*sintheta
    ymax = ymin

    for n in -2:2
        ymin = min(ymin, y + b*sin(t + n*pi)*costheta + a*cos(t + n*pi)*sintheta)
    end

    for n in -2:2
        ymax = max(ymax, y + b*sin(t + n*pi)*costheta + a*cos(t + n*pi)*sintheta)
    end

    return xmin, xmax, ymin, ymax
end
