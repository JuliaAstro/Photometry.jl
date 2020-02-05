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
    if (b < 0.0 || a < b || theta < -pi/2.0 || theta > pi/2.0)
            error("illegal ellipse parameters. " * "Require a >= b > 0.0, -pi/2 <= theta <= pi/2")
    end
    cxx, cyy, cxy = oblique_coefficients(a, b, theta)
    return EllipticalAperture(promote(x, y, a, b, theta, cxx, cxy, cyy)...)

end

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), theta=$(e.theta))")
end

function oblique_coefficients(a, b, theta)
    costheta = cos(theta)
    sintheta = sin(theta)
    cxx = costheta*costheta/(a*a) + sintheta*sintheta/(b*b)
    cyy = sintheta*sintheta/(a*a) + costheta*costheta/(b*b)
    cxy = (2.0)*costheta*sintheta * ((1.0)/(a*a) - (1.0)/(b*b))
    return cxx, cyy, cxy
end

function bbox(e::EllipticalAperture{<:AbstractFloat})

    x = e.x
    y = e.y
    a = e.a
    b = e.b
    theta = e.theta

    t = atan((1.0)*(-b*tan(theta))/a)

    xmin = x + a*cos(t)*cos(theta) - b*sin(t)*sin(theta)
    xmax = xmin

    for n in -2:2
        xmin = min(xmin, x + a*cos(t + n*pi)*cos(theta) - b*sin(t + n*pi)*sin(theta))
    end

    for n in -2:2
        xmax = max(xmax, x + a*cos(t + n*pi)*cos(theta) - b*sin(t + n*pi)*sin(theta))
    end

    t = atan((1.0)*(b*cot(theta))/a)

    ymin = y + b*sin(t)*cos(theta) + a*cos(t)*sin(theta)
    ymax = ymin

    for n in -2:2
        ymin = min(ymin, y + b*sin(t + n*pi)*cos(theta) + a*cos(t + n*pi)*sin(theta))
    end

    for n in -2:2
        ymax = max(ymax, y + b*sin(t + n*pi)*cos(theta) + a*cos(t + n*pi)*sin(theta))
    end

    return xmin, xmax, ymin, ymax
end
