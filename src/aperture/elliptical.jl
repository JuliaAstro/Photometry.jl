#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

"""
    EllipticalAperture(x, y, a, b, θ=0)
    EllipticalAperture(position, a, b, θ=0)

An elliptical aperture with semi-major axis `a`, semi-minor axis `b`, and position angle `θ`. `a` and `b` must be ≥ 0, `θ` is measured in degrees counter-clockwise the standard x-axis.

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
end


EllipticalAperture(x::Number, y::Number, a, b, theta=0) = EllipticalAperture(promote(x, y, a, b, theta)...)
EllipticalAperture(center, a, b, theta=0) = EllipticalAperture(center..., a, b, theta)

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), θ=$(e.theta)°)")
end

oblique_coefficients(ap::EllipticalAperture) = oblique_coefficients(ap.a, ap.b, ap.theta)
oblique_coefficients(ap::Subpixel{<:EllipticalAperture}) = oblique_coefficients(ap.ap)

function oblique_coefficients(a, b, theta)
    sintheta, costheta = sincos(deg2rad(theta))
    a2 = a^2
    b2 = b^2
    cxx = costheta^2 / a2 + sintheta^2 / b2
    cyy = sintheta^2 / a2 + costheta^2 / b2
    cxy = 2 * costheta * sintheta * (1 / a2 - 1 / b2)
    return cxx, cyy, cxy
end

function overlap(ap::EllipticalAperture, i, j)
    cxx, cyy, cxy = oblique_coefficients(ap)
    flags = (
        inside_ellipse(j - 0.5, i - 0.5, ap.x, ap.y, cxx, cyy, cxy),
        inside_ellipse(j - 0.5, i + 0.5, ap.x, ap.y, cxx, cyy, cxy),
        inside_ellipse(j + 0.5, i - 0.5, ap.x, ap.y, cxx, cyy, cxy),
        inside_ellipse(j + 0.5, i + 0.5, ap.x, ap.y, cxx, cyy, cxy)
    )
    all(flags) && return Inside
    all(!, flags) && return Outside

    return Partial
end

function bounds(e::EllipticalAperture)
    iszero(e.a) && return e.x, e.x, e.y, e.y
    sintheta, costheta = sincos(deg2rad(e.theta))

    t = atan(-e.b * tand(e.theta), e.a)
    
    sint, cost = sincos(t)
    xmin = e.x + e.a * cost * costheta - e.b * sint * sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmin = min(xmin, e.x + e.a * cost * costheta - e.b * sint * sintheta)
        xmax = max(xmax, e.x + e.a * cost * costheta - e.b * sint * sintheta)
    end

    t2 = atan(e.b * cotd(e.theta),  e.a)

    sint, cost = sincos(t2)
    ymin = e.y + e.b * sint * costheta + e.a * cost * sintheta
    ymax = ymin

    for n in -2:2
        sint, cost = sincos(t2 + n * pi)
        ymin = min(ymin, e.y + e.b * sint * costheta + e.a * cost * sintheta)
        ymax = max(ymax, e.y + e.b * sint * costheta + e.a * cost * sintheta)
    end

    xmin = floor(Int, xmin)
    xmax = ceil(Int, xmax)
    ymin = floor(Int, ymin)
    ymax = ceil(Int, ymax)

    return xmin, xmax, ymin, ymax
end

partial(ap::EllipticalAperture, x, y) = elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a, ap.b, ap.theta)
partial(ap::Subpixel{<:EllipticalAperture}, x, y) = elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, oblique_coefficients(ap)..., ap.N)


#######################################################

"""
    EllipticalAnnulus(x, y, a_in, a_out, b_out, θ=0)
    EllipticalAnnulus(position, a_in, a_out, b_out, θ=0)
An elliptical annulus with inner semi-major axis `a_in`, outer semi-major axis `a_out`, outer semi-minor axis `b_out`, and position angle `θ`.
`a_out` ≥ `a_in` ≥ 0 and `b_out` must be ≥ 0, `θ` is measured in degrees counter-clockwise the standard x-axis.

`b_in` will automatically be calculated from `(a_in / a_out) * b_out`. Note this may cause a type instability.
# Examples
```jldoctest
julia> ap = EllipticalAnnulus(0, 0, 4, 10, 5, 45)
EllipticalAnnulus(0.0, 0.0, a_in=4.0, a_out=10.0, b_in=2.0, b_out=5.0, θ=45.0°)
```
"""
struct EllipticalAnnulus{T <: Number} <: AbstractAperture
    x::T
    y::T
    a_in::T
    b_in::T
    a_out::T
    b_out::T
    theta::T
end

EllipticalAnnulus(x::Number, y::Number, a_in, a_out, b_out, theta=0) = EllipticalAnnulus(promote(x, y, a_in, a_in / a_out * b_out, a_out, b_out, theta)...)
EllipticalAnnulus(center, a_in, a_out, b_out, theta=0) = EllipticalAnnulus(center..., a_in, a_out, b_out, theta)

function Base.show(io::IO, e::EllipticalAnnulus)
    print(io, "EllipticalAnnulus($(e.x), $(e.y), a_in=$(e.a_in), a_out=$(e.a_out), b_in=$(e.b_in), b_out=$(e.b_out), θ=$(e.theta)°)")
end


function overlap(ap::EllipticalAnnulus, i, j)
    coeffs_out = oblique_coefficients(ap.a_out, ap.b_out, ap.theta)
    flags_out = (
        inside_ellipse(j - 0.5, i - 0.5, ap.x, ap.y, coeffs_out...),
        inside_ellipse(j - 0.5, i + 0.5, ap.x, ap.y, coeffs_out...),
        inside_ellipse(j + 0.5, i - 0.5, ap.x, ap.y, coeffs_out...),
        inside_ellipse(j + 0.5, i + 0.5, ap.x, ap.y, coeffs_out...)
    )

    coeffs_in = oblique_coefficients(ap.a_in, ap.b_in, ap.theta)
    flags_in = (
        inside_ellipse(j - 0.5, i - 0.5, ap.x, ap.y, coeffs_in...),
        inside_ellipse(j - 0.5, i + 0.5, ap.x, ap.y, coeffs_in...),
        inside_ellipse(j + 0.5, i - 0.5, ap.x, ap.y, coeffs_in...),
        inside_ellipse(j + 0.5, i + 0.5, ap.x, ap.y, coeffs_in...)
    )

   all(flags_out) && all(!, flags_in) && return Inside
   all(flags_in) || all(!, flags_out) && return Outside

    return Partial
end

function bounds(e::EllipticalAnnulus)

    sintheta, costheta = sincos(deg2rad(e.theta))

    t = atan(-e.b_out * tand(e.theta), e.a_out)

    sint, cost = sincos(t)
    xmin = e.x + e.a_out * cost * costheta - e.b_out * sint * sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmin = min(xmin, e.x + e.a_out * cost * costheta - e.b_out * sint * sintheta)
        xmax = max(xmax, e.x + e.a_out * cost * costheta - e.b_out * sint * sintheta)
    end

    t2 = atan(e.b_out * cotd(e.theta), e.a_out)

    sint, cost = sincos(t2)
    ymin = e.y + e.b_out * sint * costheta + e.a_out * cost * sintheta
    ymax = ymin

    for n in -2:2
        sint, cost = sincos(t2 + n * pi)
        ymin = min(ymin, e.y + e.b_out * sint * costheta + e.a_out * cost * sintheta)
        ymax = max(ymax, e.y + e.b_out * sint * costheta + e.a_out * cost * sintheta)
    end

    xmin = floor(Int, xmin)
    xmax = ceil(Int, xmax)
    ymin = floor(Int, ymin)
    ymax = ceil(Int, ymax)

    return xmin, xmax, ymin, ymax
end


function partial(ap::EllipticalAnnulus, x, y)
    f1 = elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_out, ap.b_out, ap.theta)
    f2 = elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_in, ap.b_in, ap.theta)
    return f1 - f2
end

function partial(ap::Subpixel{<:EllipticalAnnulus}, x, y)
    coeffs_out = oblique_coefficients(ap.a_out, ap.b_out, ap.theta)
    f1 = elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, coeffs_out..., ap.N)
    coeffs_in = oblique_coefficients(ap.a_in, ap.b_in, ap.theta)
    f2 = elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, coeffs_in..., ap.N)
    return f1 - f2
end
