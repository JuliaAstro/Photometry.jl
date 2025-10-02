#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

"""
    EllipticalAperture(x, y, a, b, θ=0)
    EllipticalAperture(position, a, b, θ=0)

An elliptical aperture with semi-major axis `a`, semi-minor axis `b`, and
position angle `θ`. `a` and `b` must be ≥ 0, `θ` is measured in degrees
counter-clockwise the standard x-axis.

# Examples
```jldoctest
julia> ap = EllipticalAperture(0, 0, 4, 2, 35)
7×5 EllipticalAperture{Int64} with indices -3:3×-2:2:
 0.873382  0.844185  0.324917  0         0
 1         1         0.997821  0.435284  0
 1         1         1         0.990119  0.23968
 0.796137  1         1         1         0.796137
 0.23968   0.990119  1         1         1
 0         0.435284  0.997821  1         1
 0         0         0.324917  0.844185  0.873382
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
        inside_ellipse(i - 0.5, j - 0.5, ap.x, ap.y, cxx, cyy, cxy),
        inside_ellipse(i - 0.5, j + 0.5, ap.x, ap.y, cxx, cyy, cxy),
        inside_ellipse(i + 0.5, j - 0.5, ap.x, ap.y, cxx, cyy, cxy),
        inside_ellipse(i + 0.5, j + 0.5, ap.x, ap.y, cxx, cyy, cxy)
    )
    all(flags) && return Inside
    !any(flags) && return Outside

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

function partial(ap::EllipticalAperture, x, y)
    return elliptical_overlap_exact(
        x - 0.5, y - 0.5, x + 0.5, y + 0.5,
        ap.a, ap.b, ap.theta,
    )
end
function partial(sub_ap::Subpixel{T,<:EllipticalAperture}, x, y) where {T}
    return elliptical_overlap_single_subpixel(
        x - 0.5, y - 0.5, x + 0.5, y + 0.5,
        oblique_coefficients(sub_ap.ap)...,
        sub_ap.N,
    )
end


#######################################################

"""
    EllipticalAnnulus(x, y, a_in, a_out, b_out, θ=0)
    EllipticalAnnulus(position, a_in, a_out, b_out, θ=0)
An elliptical annulus with inner semi-major axis `a_in`, outer semi-major axis
`a_out`, outer semi-minor axis `b_out`, and position angle `θ`.
`a_out` ≥ `a_in` ≥ 0 and `b_out` must be ≥ 0, `θ` is measured in degrees
counter-clockwise the standard x-axis.

`b_in` will automatically be calculated from `(a_in / a_out) * b_out`. Note
this may cause a type instability.
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

function EllipticalAnnulus(x::Number, y::Number, a_in, a_out, b_out, theta=0)
    return EllipticalAnnulus(
        promote(x, y, a_in, a_in / a_out * b_out, a_out, b_out, theta)...
    )
end
function EllipticalAnnulus(center, a_in, a_out, b_out, theta=0)
    return EllipticalAnnulus(
        center..., a_in, a_out, b_out, theta
    )
end

function Base.show(io::IO, e::EllipticalAnnulus)
    print(io, "EllipticalAnnulus($(e.x), $(e.y), a_in=$(e.a_in), a_out=$(e.a_out), b_in=$(e.b_in), b_out=$(e.b_out), θ=$(e.theta)°)")
end


function overlap(ap::EllipticalAnnulus, i, j)
    coeffs_out = oblique_coefficients(ap.a_out, ap.b_out, ap.theta)
    flags_out = (
        inside_ellipse(i - 0.5, j - 0.5, ap.x, ap.y, coeffs_out...),
        inside_ellipse(i - 0.5, j + 0.5, ap.x, ap.y, coeffs_out...),
        inside_ellipse(i + 0.5, j - 0.5, ap.x, ap.y, coeffs_out...),
        inside_ellipse(i + 0.5, j + 0.5, ap.x, ap.y, coeffs_out...)
    )

    coeffs_in = oblique_coefficients(ap.a_in, ap.b_in, ap.theta)
    flags_in = (
        inside_ellipse(i - 0.5, j - 0.5, ap.x, ap.y, coeffs_in...),
        inside_ellipse(i - 0.5, j + 0.5, ap.x, ap.y, coeffs_in...),
        inside_ellipse(i + 0.5, j - 0.5, ap.x, ap.y, coeffs_in...),
        inside_ellipse(i + 0.5, j + 0.5, ap.x, ap.y, coeffs_in...)
    )

   all(flags_out) && !any(flags_in) && return Inside
   all(flags_in) || !any(flags_out) && return Outside

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

@testsnippet elliptical begin
    using Photometry.Aperture: EllipticalAperture, EllipticalAnnulus,
                               bounds, center, oblique_coefficients, photometry
end

@testitem "aperture/elliptical: Apertures" setup=[elliptical] begin
    ap_ellipse = EllipticalAperture(0, 0, 20, 10, 0)
    @test center(ap_ellipse) == (0, 0)
    @test bounds(ap_ellipse) == (-20, 20, -10, 10)
    @test size(ap_ellipse) == (41, 21)
    @test size(ap_ellipse, 1) == 41
    @test all(axes(ap_ellipse) .== (-20:20, -10:10))
    @test eachindex(ap_ellipse) == CartesianIndex(-20, -10):CartesianIndex(20, 10)

    @test EllipticalAperture([0, 0], 20, 10, 0) == ap_ellipse

    ap_ellipse = EllipticalAperture(0, 0, 2, 1, 45)
    @test bounds(ap_ellipse) == (-1, 1, -1, 1)
end

@testitem "aperture/elliptical: Elliptical Aperture" setup=[elliptical] begin
    e = EllipticalAperture(0, 0, 20, 10, 0)
    @test sprint(show, e) == "EllipticalAperture(0, 0, a=20, b=10, θ=0°)"
end

@testitem "aperture/elliptical: oblique_coefficients" setup=[elliptical] begin
    @test all(oblique_coefficients(2, 2, 0) .≈ (0.25, 0.25, 0.0))
    @test all(oblique_coefficients(2, 2, 90) .≈ (0.25, 0.25, 0.0))
    @test all(oblique_coefficients(2, 1, 30) .≈ (7 / 16, 13 / 16, -6sqrt(3) / 16))
end

@testitem "aperture/elliptical: Elliptical Annulus" setup=[elliptical] begin
    e0 = EllipticalAnnulus(0, 0, 8, 16, 4, 45)
    @test sprint(show, e0) == "EllipticalAnnulus(0.0, 0.0, a_in=8.0, a_out=16.0, b_in=2.0, b_out=4.0, θ=45.0°)"
    @test EllipticalAnnulus([0, 0], 8, 16, 4, 45) == EllipticalAnnulus(0, 0, 8, 16, 4, 45)
end

@testitem "aperture/elliptical: Elliptical Annulus bounding box" setup=[elliptical] begin
    e = EllipticalAnnulus(0, 0, 8, 16, 4, 0)
    @test center(e) == (0, 0)
    @test bounds(e) == (-16, 16, -4, 4)
    @test size(e) == (33, 9)
    @test all(axes(e) .== (-16:16, -4:4))
    @test eachindex(e) == CartesianIndex(-16, -4):CartesianIndex(16, 4)
end

@testitem "aperture/elliptical: regression circular ellipse" setup=[elliptical] begin
    # some weird bug where centered on-grid past a certain size (and rotated) would fail
    e = EllipticalAperture(5.5, 5.5, 4, 4, 20)
    @test sum(e) ≈ photometry(e, ones(9, 9)).aperture_sum # just a test that no errors occur
end
