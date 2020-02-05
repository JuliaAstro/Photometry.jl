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
