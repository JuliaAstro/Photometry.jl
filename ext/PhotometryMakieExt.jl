"""
Makie plotting support for apertures, mirroring the Plots recipes in
`src/aperture/plotting.jl`: `plot`/`lines` draw an aperture's outline, annuli
draw both rings, and a vector of apertures draws every outline in a single
plot call. An optional trailing argument controls the number of outline
samples, e.g. `lines(ap, 25)`.

One deliberate difference from the Plots recipes: outlines are centered on the
aperture's own `(x, y)` rather than shifted by +0.5, matching Makie's `heatmap`
convention of centering cell `(i, j)` at `(i, j)`.
"""
module PhotometryMakieExt

using Makie: Makie, Point2d
using Photometry.Aperture: AbstractAperture, Subpixel,
    CircularAperture, CircularAnnulus,
    EllipticalAperture, EllipticalAnnulus,
    RectangularAperture, RectangularAnnulus

Makie.plottype(::AbstractAperture) = Makie.Lines
Makie.plottype(::AbstractVector{<:AbstractAperture}) = Makie.Lines

_point(ap::CircularAperture, θ) = Point2d(ap.x + ap.r * cos(θ), ap.y + ap.r * sin(θ))

function _ellipsepoint(x, y, a, b, θ, t)
    u = x + a * cos(θ) * cosd(t) - b * sin(θ) * sind(t)
    v = y + a * cos(θ) * sind(t) + b * sin(θ) * cosd(t)
    return Point2d(u, v)
end
_point(ap::EllipticalAperture, θ) = _ellipsepoint(ap.x, ap.y, ap.a, ap.b, θ, ap.theta)

# Exact rectangle perimeter: for θ ∈ [0, π/2], |c|c + |s|s ≡ 1 pins u to the
# right edge while |c|c - |s|s sweeps v across it, and so on around the sides.
function _rectanglepoint(x, y, w, h, θ, t)
    sinth, costh = sincos(θ)
    u = w / 2 * (abs(costh) * costh + abs(sinth) * sinth)
    v = h / 2 * (abs(costh) * costh - abs(sinth) * sinth)
    sint, cost = sincos(deg2rad(t))
    return Point2d(x + u * cost - v * sint, y + u * sint + v * cost)
end
_point(ap::RectangularAperture, θ) = _rectanglepoint(ap.x, ap.y, ap.w, ap.h, θ, ap.theta)

# Odd default sample count closes the outline exactly at θ = 2π.
_outline(ap, n = 101) = [_point(ap, θ) for θ in range(0, 2π; length = n)]

# NB: the NaN line break must be spliced in as a 1-element vector.
# `vcat` treats a bare Point2d as a static vector and concatenates its scalars,
# degrading the result to Vector{Any}.
const _BREAK = [Point2d(NaN, NaN)]

_rings(ap::CircularAnnulus) = (
    CircularAperture(ap.x, ap.y, ap.r_out),
    CircularAperture(ap.x, ap.y, ap.r_in),
)
_rings(ap::EllipticalAnnulus) = (
    EllipticalAperture(ap.x, ap.y, ap.a_out, ap.b_out, ap.theta),
    EllipticalAperture(ap.x, ap.y, ap.a_in, ap.b_in, ap.theta),
)
_rings(ap::RectangularAnnulus) = (
    RectangularAperture(ap.x, ap.y, ap.w_out, ap.h_out, ap.theta),
    RectangularAperture(ap.x, ap.y, ap.w_in, ap.h_in, ap.theta),
)

function _outline(ap::Union{CircularAnnulus, EllipticalAnnulus, RectangularAnnulus}, n = 101)
    outer, inner = _rings(ap)
    return vcat(_outline(outer, n), _BREAK, _outline(inner, n))
end

# Subpixel is a sampling strategy; its footprint is the wrapped aperture's.
_outline(ap::Subpixel, n = 101) = _outline(ap.ap, n)

function Makie.convert_arguments(::Makie.PointBased, ap::AbstractAperture, npoints::Integer = 101)
    return (_outline(ap, npoints),)
end

function Makie.convert_arguments(::Makie.PointBased, aps::AbstractVector{<:AbstractAperture})
    isempty(aps) && return (Point2d[],)
    return (reduce((a, b) -> vcat(a, _BREAK, b), map(_outline, aps)),)
end

end
