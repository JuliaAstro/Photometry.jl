#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#


"""
    RectangularAperture(x, y, w, h, θ=0)
    RectangularAperture(position, w, h, θ=0)

A rectangular aperture.

A rectangular aperture with width `w`, height `h`,
and position angle `θ` in degrees.

# Examples
```jldoctest
julia> ap = RectangularAperture(0, 0, 10, 4, 0)
11×5 RectangularAperture{Int64} with indices -5:5×-2:2:
 0.25  0.5  0.5  0.5  0.25
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.5   1    1    1    0.5
 0.25  0.5  0.5  0.5  0.25
```
"""
struct RectangularAperture{T <: Number} <: AbstractAperture{T}
    x::T
    y::T
    w::T
    h::T
    theta::T
    function RectangularAperture(x::T, y::T, w::T, h::T, θ::T) where T <: Number
        (w < 0  || h < 0) && error("Invalid side lengths w=$w, h=$h. Both must be positive.")
        new{T}(x, y, w, h, mod(θ, 360))
    end
end

RectangularAperture(center, w, h, θ) = RectangularAperture(center..., w, h, θ)
RectangularAperture(x, y, w, h, θ) = RectangularAperture(promote(x, y, w, h, θ)...)

function Base.show(io::IO, ap::RectangularAperture)
    print(io, "RectangularAperture($(ap.x), $(ap.y), w=$(ap.w), h=$(ap.h), θ=$(ap.theta)°)")
end

function bounds(ap::RectangularAperture)
    w2 = ap.w / 2
    h2 = ap.h / 2
    sint, cost = sincos(deg2rad(ap.theta))

    dx1 = abs(w2 * cost - h2 * sint)
    dy1 = abs(w2 * sint + h2 * cost)
    dx2 = abs(w2 * cost + h2 * sint)
    dy2 = abs(w2 * sint - h2 * cost)

    dx = max(dx1, dx2)
    dy = max(dy1, dy2)

    xmin = ceil(Int, ap.x - dx - 0.5)
    ymin = ceil(Int, ap.y - dy - 0.5)
    xmax = ceil(Int, ap.x + dx - 0.5)
    ymax = ceil(Int, ap.y + dy - 0.5)
    return xmin, xmax, ymin, ymax
end

function overlap(ap::RectangularAperture, i, j)
    x = i - ap.x
    y = j - ap.y
    flags = (
        inside_rectangle(x - 0.5, y - 0.5, ap.w, ap.h, ap.theta),
        inside_rectangle(x - 0.5, y + 0.5, ap.w, ap.h, ap.theta),
        inside_rectangle(x + 0.5, y - 0.5, ap.w, ap.h, ap.theta),
        inside_rectangle(x + 0.5, y + 0.5, ap.w, ap.h, ap.theta)
    )
    all(flags) && return Inside
    !any(flags) && return Outside

    return Partial
end

partial(ap::RectangularAperture, x, y) = rectangular_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.w, ap.h, ap.theta)
partial(sub_ap::Subpixel{T,<:RectangularAperture}, x, y) where {T} = rectangular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, sub_ap.ap.w, sub_ap.ap.h, sub_ap.ap.theta, sub_ap.N)

#######################################################

"""
    RectangularAnnulus(x, y, w_in, w_out, h_out, θ=0)
    RectangularAnnulus(position, w_in, w_out, h_out, θ=0)

A rectangular annulus with inner width `w_in`, outer width `w_out`,
outer height `h_out`, and position angle `θ` in degrees. `h_in` is automatically
calculated from `w_in / w_out * h_out`. Note that `w_out ≥ w_in > 0`.

# Examples
```jldoctest
julia> ap = RectangularAnnulus(0, 0, 5, 10, 8, 45)
13×13 RectangularAnnulus{Float64} with indices -6:6×-6:6:
 0.0       0.0       0.0         …  0.0         0.0       0.0
 0.0       0.0       0.0            0.0         0.0       0.0
 0.0       0.0       0.00252532     0.0         0.0       0.0
 0.0       0.0       0.568542       0.0         0.0       0.0
 0.0       0.568542  1.0            0.215729    0.0       0.0
 0.528175  1.0       1.0         …  1.0         0.215729  0.0
 0.215729  1.0       1.0            1.0         1.0       0.215729
 0.0       0.215729  1.0            1.0         1.0       0.528175
 0.0       0.0       0.215729       1.0         0.568542  0.0
 0.0       0.0       0.0            0.568542    0.0       0.0
 0.0       0.0       0.0         …  0.00252532  0.0       0.0
 0.0       0.0       0.0            0.0         0.0       0.0
 0.0       0.0       0.0            0.0         0.0       0.0
```
"""
struct RectangularAnnulus{T <: Number} <: AbstractAperture{T}
    x::T
    y::T
    w_in::T
    w_out::T
    h_in::T
    h_out::T
    theta::T
end

RectangularAnnulus(center, w_in, w_out, h_out, θ=0) = RectangularAnnulus(center..., w_in, w_out, h_out, θ)
RectangularAnnulus(x::Number, y::Number, w_in, w_out, h_out, θ=0) = RectangularAnnulus(promote(x, y, w_in, w_out, h_out * w_in / w_out, h_out, θ)...)
# RectangularAnnulus(x, y, w_in, w_out, h_out) = RectangularAnnulus(promote(x, y, w_in, w_out, h_out, 0)...)

function Base.show(io::IO, ap::RectangularAnnulus)
    print(io, "RectangularAnnulus($(ap.x), $(ap.y), w_in=$(ap.w_in), w_out=$(ap.w_out), h_in=$(ap.h_in), h_out=$(ap.h_out), θ=$(ap.theta)°)")
end


function overlap(ap::RectangularAnnulus, i, j)
    x = i - ap.x
    y = j - ap.y
    flags_out = (
        inside_rectangle(x - 0.5, y - 0.5, ap.w_out, ap.h_out, ap.theta),
        inside_rectangle(x - 0.5, y + 0.5, ap.w_out, ap.h_out, ap.theta),
        inside_rectangle(x + 0.5, y - 0.5, ap.w_out, ap.h_out, ap.theta),
        inside_rectangle(x + 0.5, y + 0.5, ap.w_out, ap.h_out, ap.theta)
    )

    flags_in = (
        inside_rectangle(x - 0.5, y - 0.5, ap.w_in, ap.h_in, ap.theta),
        inside_rectangle(x - 0.5, y + 0.5, ap.w_in, ap.h_in, ap.theta),
        inside_rectangle(x + 0.5, y - 0.5, ap.w_in, ap.h_in, ap.theta),
        inside_rectangle(x + 0.5, y + 0.5, ap.w_in, ap.h_in, ap.theta)
    )

    all(flags_out) && !any(flags_in) && return Inside
    all(flags_in) || !any(flags_out) && return Outside

    return Partial
end

function bounds(ap::RectangularAnnulus)
    w2 = ap.w_out / 2
    h2 = ap.h_out / 2
    sint, cost = sincos(deg2rad(ap.theta))

    dx1 = abs(w2 * cost - h2 * sint)
    dy1 = abs(w2 * sint + h2 * cost)
    dx2 = abs(w2 * cost + h2 * sint)
    dy2 = abs(w2 * sint - h2 * cost)

    dx = max(dx1, dx2)
    dy = max(dy1, dy2)

    xmin = ceil(Int, ap.x - dx - 0.5)
    ymin = ceil(Int, ap.y - dy - 0.5)
    xmax = ceil(Int, ap.x + dx - 0.5)
    ymax = ceil(Int, ap.y + dy - 0.5)
    return (xmin, xmax, ymin, ymax)
end

function partial(ap::RectangularAnnulus, x, y)
    f1 = rectangular_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5,
                                   ap.w_out, ap.h_out, ap.theta)
    f2 = rectangular_overlap_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5,
                                   ap.w_in, ap.h_in, ap.theta)
    return f1 - f2
end

function partial(sub_ap::Subpixel{T,<:RectangularAnnulus}, x, y) where T
    ap = sub_ap.ap
    f1 = rectangular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5,
                                             ap.w_out, ap.h_out, ap.theta, sub_ap.N)
    f2 = rectangular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5,
                                             ap.w_in, ap.h_in, ap.theta, sub_ap.N)
    return f1 - f2
end
