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
5×7 EllipticalAperture{Int64} with indices -2:2×-3:3:
 0.873382  1.0       1.0       0.796137  0.23968   0.0       0.0
 0.844185  1.0       1.0       1.0       0.990119  0.435284  0.0
 0.324917  0.997821  1.0       1.0       1.0       0.997821  0.324917
 0.0       0.435284  0.990119  1.0       1.0       1.0       0.844185
 0.0       0.0       0.23968   0.796137  1.0       1.0       0.873382
```
"""
struct EllipticalAperture{T <: Number} <: AbstractAperture{T}
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
oblique_coefficients(ap::Subpixel{T,<:EllipticalAperture}) where {T} = oblique_coefficients(ap.ap)

function oblique_coefficients(a, b, theta)
    sintheta, costheta = sincosd(theta)
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

function elliptical_bounds(cx, cy, a, b, theta)
    iszero(a) && return cx, cx, cy, cy
    sintheta, costheta = sincosd(theta)

    t = atan(-b * tand(theta), a)
    
    sint, cost = sincos(t)
    xmin = cx + a * cost * costheta - b * sint * sintheta
    xmax = xmin

    for n in -2:2
        sint, cost = sincos(t + n * pi)
        xmin = min(xmin, cx + a * cost * costheta - b * sint * sintheta)
        xmax = max(xmax, cx + a * cost * costheta - b * sint * sintheta)
    end

    t2 = atan(b * cotd(theta),  a)

    sint, cost = sincos(t2)
    ymin = cy + b * sint * costheta + a * cost * sintheta
    ymax = ymin

    for n in -2:2
        sint, cost = sincos(t2 + n * pi)
        ymin = min(ymin, cy + b * sint * costheta + a * cost * sintheta)
        ymax = max(ymax, cy + b * sint * costheta + a * cost * sintheta)
    end

    xmin = ceil(Int, xmin)
    xmax = floor(Int, xmax)
    ymin = ceil(Int, ymin)
    ymax = floor(Int, ymax)

    return xmin, xmax, ymin, ymax
end


bounds(e::EllipticalAperture) = elliptical_bounds(e.x, e.y, e.a, e.b, e.theta)

partial(ap::EllipticalAperture, x, y) = elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a, ap.b, ap.theta)
partial(sub_ap::Subpixel{T,<:EllipticalAperture}, x, y) where {T} = elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, oblique_coefficients(sub_ap.ap)..., sub_ap.N)


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
15×15 EllipticalAnnulus{Float64} with indices -7:7×-7:7:
 0.594853   1.0       1.0       1.0         …  0.0       0.0       0.0
 1.0        1.0       1.0       1.0            0.0       0.0       0.0
 1.0        1.0       1.0       1.0            0.0       0.0       0.0
 1.0        1.0       1.0       1.0            0.0       0.0       0.0
 1.0        1.0       1.0       1.0            0.0       0.0       0.0
 0.814163   1.0       1.0       1.0         …  0.414163  0.0       0.0
 0.369432   1.0       1.0       1.0            0.975704  0.193728  0.0
 0.0112571  0.809079  1.0       1.0            1.0       0.809079  0.0112571
 0.0        0.193728  0.975704  1.0            1.0       1.0       0.369432
 0.0        0.0       0.414163  1.0            1.0       1.0       0.814163
 0.0        0.0       0.0       0.546165    …  1.0       1.0       1.0
 0.0        0.0       0.0       0.00252321     1.0       1.0       1.0
 0.0        0.0       0.0       0.0            1.0       1.0       1.0
 0.0        0.0       0.0       0.0            1.0       1.0       1.0
 0.0        0.0       0.0       0.0            1.0       1.0       0.594853
```
"""
struct EllipticalAnnulus{T <: Number} <: AbstractAperture{T}
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

bounds(e::EllipticalAnnulus) =  elliptical_bounds(e.x, e.y, e.a_out, e.b_out, e.theta)


function partial(ap::EllipticalAnnulus, x, y)
    f1 = elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_out, ap.b_out, ap.theta)
    f2 = elliptical_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.a_in, ap.b_in, ap.theta)
    return f1 - f2
end

function partial(sub_ap::Subpixel{T,<:EllipticalAnnulus}, x, y) where T
    ap = sub_ap.ap
    coeffs_out = oblique_coefficients(ap.a_out, ap.b_out, ap.theta)
    f1 = elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, coeffs_out..., sub_ap.N)
    coeffs_in = oblique_coefficients(ap.a_in, ap.b_in, ap.theta)
    f2 = elliptical_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, coeffs_in..., sub_ap.N)
    return f1 - f2
end
