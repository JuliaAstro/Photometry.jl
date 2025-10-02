#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

###############

"""
    CircularAperture(x, y, r)
    CircularAperture(position, r)

A circular aperture.

A circular aperture with radius `r`. `r` must be greater than or equal to 0.

# Examples
```jldoctest
julia> ap = CircularAperture(0, 0, 10)
 21×21 CircularAperture{Int64} with indices -10:10×-10:10:
 0.0        0.0       0.0         …  0.0         0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.00571026     0.00571026  0.0       0.0
 0.0        0.0       0.491844       0.491844    0.0       0.0
 0.0        0.170878  0.982952    …  0.982952    0.170878  0.0
 0.0        0.659735  1.0            1.0         0.659735  0.0
 0.0590655  0.975524  1.0            1.0         0.975524  0.0590655
 0.293527   1.0       1.0            1.0         1.0       0.293527
 0.445643   1.0       1.0            1.0         1.0       0.445643
 0.495832   1.0       1.0         …  1.0         1.0       0.495832
 0.445643   1.0       1.0            1.0         1.0       0.445643
 0.293527   1.0       1.0            1.0         1.0       0.293527
 0.0590655  0.975524  1.0            1.0         0.975524  0.0590655
 0.0        0.659735  1.0            1.0         0.659735  0.0
 0.0        0.170878  0.982952    …  0.982952    0.170878  0.0
 0.0        0.0       0.491844       0.491844    0.0       0.0
 0.0        0.0       0.00571026     0.00571026  0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.0         …  0.0         0.0       0.0
```
"""
struct CircularAperture{T<:Number} <: AbstractAperture{T}
    x::T
    y::T
    r::T
end

CircularAperture(center, r) = CircularAperture(center..., r)
CircularAperture(x, y, r) = CircularAperture(promote(x, y, r)...)

function overlap(ap::CircularAperture, i, j)
    dist = sqrt((i - ap.x)^2 + (j - ap.y)^2)
    dr = sqrt(2) / 2 # corner-center distance of pixel
    dist > ap.r + dr && return Outside
    dist < ap.r - dr && return Inside
    return Partial
end

function partial(ap::CircularAperture, x, y)
    return circular_overlap_single_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r)
end
function partial(sub_ap::Subpixel{T,<:CircularAperture}, x, y) where {T}
    return circular_overlap_single_subpixel(
        x - 0.5, y - 0.5, x + 0.5, y + 0.5, sub_ap.ap.r, sub_ap.N)
end

function Base.show(io::IO, c::CircularAperture)
    print(io, "CircularAperture($(c.x), $(c.y), r=$(c.r))")
end


function bounds(c::CircularAperture)
    _xmin = c.x - c.r - 0.5
    _ymin = c.y - c.r - 0.5
    xmin = isinteger(_xmin) ? ceil(Int, _xmin) + 1 : ceil(Int, _xmin)
    ymin = isinteger(_ymin) ? ceil(Int, _ymin) + 1 : ceil(Int, _ymin)
    xmax = ceil(Int, c.x + c.r - 0.5)
    ymax = ceil(Int, c.y + c.r - 0.5)
    return (xmin, xmax, ymin, ymax)
end

#######################################################


"""
    CircularAnnulus(x, y, r_in, r_out)
    CircularAnnulus(position, r_in, r_out)

A circular annulus with inner radius `r_in` and outer radius `r_out`. 0 ≤ `r_in` ≤ `r_out`.

# Examples
```jldoctest
julia> ap = CircularAnnulus(0, 0, 5, 10)
 21×21 CircularAnnulus{Int64} with indices -10:10×-10:10:
 0.0        0.0       0.0         …  0.0         0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.00571026     0.00571026  0.0       0.0
 0.0        0.0       0.491844       0.491844    0.0       0.0
 0.0        0.170878  0.982952    …  0.982952    0.170878  0.0
 0.0        0.659735  1.0            1.0         0.659735  0.0
 0.0590655  0.975524  1.0            1.0         0.975524  0.0590655
 0.293527   1.0       1.0            1.0         1.0       0.293527
 0.445643   1.0       1.0            1.0         1.0       0.445643
 0.495832   1.0       1.0         …  1.0         1.0       0.495832
 0.445643   1.0       1.0            1.0         1.0       0.445643
 0.293527   1.0       1.0            1.0         1.0       0.293527
 0.0590655  0.975524  1.0            1.0         0.975524  0.0590655
 0.0        0.659735  1.0            1.0         0.659735  0.0
 0.0        0.170878  0.982952    …  0.982952    0.170878  0.0
 0.0        0.0       0.491844       0.491844    0.0       0.0
 0.0        0.0       0.00571026     0.00571026  0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.0            0.0         0.0       0.0
 0.0        0.0       0.0         …  0.0         0.0       0.0
```
"""
struct CircularAnnulus{T <: Number} <: AbstractAperture{T}
    x::T
    y::T
    r_in::T
    r_out::T
end


CircularAnnulus(center, r_in, r_out) = CircularAnnulus(center..., r_in, r_out)
CircularAnnulus(x, y, r_in, r_out) = CircularAnnulus(promote(x, y, r_in, r_out)...)

function Base.show(io::IO, c::CircularAnnulus)
    print(io, "CircularAnnulus($(c.x), $(c.y), r_in=$(c.r_in), r_out=$(c.r_out))")
end


function overlap(ap::CircularAnnulus, i, j)
    dist = sqrt((i - ap.x)^2 + (j - ap.y)^2)
    dr = sqrt(2) / 2 # corner-center distance of pixel
    ap.r_in - dr < dist < ap.r_out + dr || return Outside
    ap.r_in + dr < dist < ap.r_out - dr && return Inside
    return Partial
end

function partial(ap::CircularAnnulus, x, y)
    f1 = circular_overlap_single_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_out)
    f2 = circular_overlap_single_exact(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_in)
    return f1 - f2
end

function partial(sub_ap::Subpixel{T,<:CircularAnnulus}, x, y) where T
    ap = sub_ap.ap
    f1 = circular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_out, sub_ap.N)
    f2 = circular_overlap_single_subpixel(x - 0.5, y - 0.5, x + 0.5, y + 0.5, ap.r_in, sub_ap.N)
    return f1 - f2
end

function bounds(c::CircularAnnulus)
    xmin = ceil(Int, c.x - c.r_out - 0.5)
    ymin = ceil(Int, c.y - c.r_out - 0.5)
    xmax = ceil(Int, c.x + c.r_out - 0.5)
    ymax = ceil(Int, c.y + c.r_out - 0.5)
    return (xmin, xmax, ymin, ymax)
end

@testsnippet circular begin
    using Photometry.Aperture: CircularAperture, CircularAnnulus, Subpixel,
                               bounds, center, photometry
    import Random
    Random.seed!(8462852)
end

@testitem "aperture/circular: Apertures" setup=[circular] begin
    ap_circ = CircularAperture(50, 40, 10)
    @test center(ap_circ) == (50, 40)
    @test bounds(ap_circ) == (40, 60, 30, 50)
    @test size(ap_circ) == map(length, axes(ap_circ)) == (21, 21)
    @test size(ap_circ, 1) == 21
    @test all(axes(ap_circ) .== (40:60, 30:50))
    @test eachindex(ap_circ) == CartesianIndex(40, 30):CartesianIndex(60, 50)

    @test CircularAperture([50, 40], 10) == ap_circ

    ap_ann = CircularAnnulus(50, 40, 5, 10)
    @test center(ap_ann) == (50, 40)
    @test bounds(ap_ann) == (40, 60, 30, 50)
    @test size(ap_ann) == (21, 21)
    @test all(axes(ap_ann) .== (40:60, 30:50))
    @test eachindex(ap_ann) == CartesianIndex(40, 30):CartesianIndex(60, 50)

    @test CircularAnnulus([50, 40], 5, 10) == ap_ann
end

@testitem "aperture/circular: broadcasting" setup=[circular] begin
    ap = CircularAperture(3, 3, 2.5)
    data = ones(5, 7)
    weighted = ap .* data
    @test weighted == data .* ap # commutative
    # TODO: Investigate float precision issue in REPL vs. run from script
    @test sum(weighted) ≈ sum(ap) == photometry(ap, data).aperture_sum
    @test all(iszero, weighted[:, 6:7])
end

@testitem "aperture/circular: CircularAperture" setup=[circular] begin
    c0 = CircularAperture(0, 0, 0)

    @test sprint(show, c0) == "CircularAperture(0, 0, r=0)"

    c1 = CircularAperture(0, 0, 1)

    @test sprint(show, c1) == "CircularAperture(0, 0, r=1)"

    sub_c1 = Subpixel(c1, 5)
    @test sprint(show, sub_c1) == "Subpixel(CircularAperture(0, 0, r=1), 5)"
end

@testitem "aperture/circular: CircularAnnulus" setup=[circular] begin
    c0 = CircularAnnulus(0, 0, 0, 0)

    @test sprint(show, c0) == "CircularAnnulus(0, 0, r_in=0, r_out=0)"

    c1 = CircularAnnulus(0, 0, 0, 1)

    @test sprint(show, c1) == "CircularAnnulus(0, 0, r_in=0, r_out=1)"
end
