#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

"""
    EllipticalAperture(x, y, a, b, θ)
    EllipticalAperture(position, a, b, θ)

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

    function EllipticalAperture(x::T, y::T, a::T, b::T, theta::T) where T <: Number
        a < 0 && error("Invalid axis a=$a. a must be greater than or equal to 0")
        b < 0 && error("Invalid axis b=$b. b must be greater than or equal to 0")
        new{T}(x, y, a, b, mod(theta, 360))
    end
end


EllipticalAperture(x, y, a, b, theta) = EllipticalAperture(promote(x, y, a, b, theta)...)
EllipticalAperture(center, a, b, theta) = EllipticalAperture(center..., a, b, theta)

function Base.show(io::IO, e::EllipticalAperture)
    print(io, "EllipticalAperture($(e.x), $(e.y), a=$(e.a), b=$(e.b), θ=$(e.theta)°)")
end

oblique_coefficients(ap::EllipticalAperture) = oblique_coefficients(ap.a, ap.b, ap.theta)

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
    n_in = count(flags)
    n_in === 4 && return Inside
    n_in === 0 && return Outside

    return Partial
end

function bounds(e::EllipticalAperture)
    sintheta, costheta = sincos(deg2rad(e.theta))

    t = atan((-e.b * tand(e.theta)) / e.a)

    sint, cost = sincos(t)
    xmin = e.x + e.a * cost * costheta - e.b * sint * sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmin = min(xmin, e.x + e.a * cost * costheta - e.b * sint * sintheta)
        xmax = max(xmax, e.x + e.a * cost * costheta - e.b * sint * sintheta)
    end

    t2 = atan((e.b * cotd(e.theta)) / e.a)

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

partial(ap::EllipticalAperture) = (x, y) -> elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a, ap.b, ap.theta)
partial(ap::Subpixel{<:EllipticalAperture}) = (x, y) -> elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a, ap.b, ap.theta, ap.N)


#######################################################

"""
    EllipticalAnnulus(x, y, a_in, a_out, b_out, θ)
    EllipticalAnnulus(position, a_in, a_out, b_out, θ)
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

    function EllipticalAnnulus(x::T, y::T, a_in::T, b_in::T, a_out::T, b_out::T, theta::T) where T <: Number
        0 ≤ a_in ≤ a_out || error("Invalid axis a_in=$a_in or a_out=$a_out. a_out must be greater than a_in which must be greater than or equal to 0")
        0 ≤ b_in ≤ b_out || error("Invalid axis b_in=$b_in or b_out=$b_out. b_out must be greater than b_in which must be greater than or equal to 0")
        new{T}(x, y, a_in, b_in, a_out, b_out, mod(theta, 360))
    end
end

EllipticalAnnulus(x, y, a_in, a_out, b_out, theta) = EllipticalAnnulus(promote(x, y, a_in, a_in / a_out * b_out, a_out, b_out, theta)...)
EllipticalAnnulus(center::AbstractVector, a_in, a_out, b_out, theta) = EllipticalAnnulus(center..., a_in, a_out, b_out, theta)

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

    t = atan((-e.b_out * tand(e.theta)) / e.a_out)

    sint, cost = sincos(t)
    xmin = e.x + e.a_out * cost * costheta - e.b_out * sint * sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmin = min(xmin, e.x + e.a_out * cost * costheta - e.b_out * sint * sintheta)
        xmax = max(xmax, e.x + e.a_out * cost * costheta - e.b_out * sint * sintheta)
    end

    t2 = atan((e.b_out * cotd(e.theta)) / e.a_out)

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


function partial(ap::EllipticalAnnulus)
    (x, y) -> elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_out, ap.b_out, ap.theta) -
              elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_in, ap.b_in, ap.theta)
end

function partial(ap::Subpixel{<:EllipticalAnnulus})
    (x, y) -> elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_out, ap.b_out, ap.theta, ap.N) -
              elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_in, ap.b_in, ap.theta, ap.N)
end
