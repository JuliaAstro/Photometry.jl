#=
Part of this work is derived from astropy/photutils. The relevant derivations
are considered under a BSD 3-clause license. =#

module Aperture

using TypedTables
using Transducers
using TestItems

export photometry,
       Subpixel,
       CircularAperture,
       CircularAnnulus,
       EllipticalAperture,
       EllipticalAnnulus,
       RectangularAperture,
       RectangularAnnulus

"""
    AbstractAperture{T} <: AbstractMatrix{T}

The abstract super-type for Apertures.

Apertures can be thought of as a cutout or stamp of a geometric shape with
shading applied. For example, a circular aperture with a diameter of 3 pixels
will require a 5x5 pixel grid (when perfectly on-grid) to represent.

```jldoctest ap1
julia> ap = CircularAperture(3, 3, 2.5)
5×5 CircularAperture{Float64} with indices 1:5×1:5:
 0.136857  0.769325  0.983232  0.769325  0.136857
 0.769325  1.0       1.0       1.0       0.769325
 0.983232  1.0       1.0       1.0       0.983232
 0.769325  1.0       1.0       1.0       0.769325
 0.136857  0.769325  0.983232  0.769325  0.136857
```

This is a useful way of thinking about apertures: if we have some data,
we can weight the data with the aperture.

```jldoctest ap1
julia> data = fill(2, 5, 5);

julia> idxs = map(intersect, axes(ap), axes(data)) |> CartesianIndices;

julia> weighted_cutout = data[idxs] .* ap[idxs]
5×5 Matrix{Float64}:
 0.273713  1.53865  1.96646  1.53865  0.273713
 1.53865   2.0      2.0      2.0      1.53865
 1.96646   2.0      2.0      2.0      1.96646
 1.53865   2.0      2.0      2.0      1.53865
 0.273713  1.53865  1.96646  1.53865  0.273713
```

Performing aperture photometry is merely summing the weighted cutout shown above.

```jldoctest ap1
julia> flux = sum(weighted_cutout)
39.269908169872416

julia> flux ≈ (π * 2.5^2) * 2 # area of circle times intensity of 2
true
```

What's interesting about the implementation of apertures, though,
is they are lazy. This means there is no stored matrix of aperture values;
rather, they are calculated on the fly as needed.

```jldoctest ap1
julia> axes(ap)
(1:5, 1:5)

julia> ap[-10, -10] # out-of-bounds, but calculated on the fly
0.0

julia> ap .* ones(5, 7) # broadcasts to eachindex(data), regardless of ap bound
5×7 Matrix{Float64}:
 0.136857  0.769325  0.983232  0.769325  0.136857  0.0  0.0
 0.769325  1.0       1.0       1.0       0.769325  0.0  0.0
 0.983232  1.0       1.0       1.0       0.983232  0.0  0.0
 0.769325  1.0       1.0       1.0       0.769325  0.0  0.0
 0.136857  0.769325  0.983232  0.769325  0.136857  0.0  0.0
```

This allows extremely efficient computation of aperture photometry from small to
medium sized apertures.

```julia-repl
julia> using BenchmarkTools

julia> @btime sum(idx -> \$ap[idx] * \$data[idx], \$idxs)
  1.097 μs (0 allocations: 0 bytes)
39.26990816987243
```

This is essentially the full implementation of [`photometry`](@ref), save for
the packing of additional information into a tabular form.
"""
abstract type AbstractAperture{T}  <: AbstractMatrix{T} end

"""
    bounds(::AbstractAperture)

Return the (`xlow`, `xhigh`, `ylow`, `yhigh`) bounds for a given Aperture.
"""
bounds(::AbstractAperture)

center(ap::AbstractAperture) = ap.x, ap.y # greedy

"""
    size(::AbstractAperture)

Return (`ny`, `nx`) of the aperture.
"""
function Base.size(ap::AbstractAperture)
    xmin, xmax, ymin, ymax = bounds(ap)
    return xmax - xmin + 1, ymax - ymin + 1
end

function Base.axes(ap::AbstractAperture)
    xmin, xmax, ymin, ymax = bounds(ap)
    return xmin:xmax, ymin:ymax
end

"""
    Subpixel(ap, N=1) <: AbstractAperture

Use a subpixel quadrature approximation for pixel shading instead of exact
geometric methods.

For any pixel laying on the border of `ap`, this alters the shading algorithm by
breaking the border pixel up into `(N, N)` subpixels. The shading value is the
fraction of these subpixels within the geometric border of `ap`.

Using a subpixel shading method is sometimes faster than exact methods at the
cost of accuracy. For [`CircularAperture`](@ref) the subpixel method is only
faster than the exact method for `N` ~ 7. for [`EllipticalAperture`](@ref) the
cutoff is `N` ~ 12, and for [`RectangularAperture`](@ref) the cutoff is `N` ~ 20.


# Examples

```jldoctest
julia> ap = CircularAperture(3, 3, 2.5)
5×5 CircularAperture{Float64} with indices 1:5×1:5:
 0.136857  0.769325  0.983232  0.769325  0.136857
 0.769325  1.0       1.0       1.0       0.769325
 0.983232  1.0       1.0       1.0       0.983232
 0.769325  1.0       1.0       1.0       0.769325
 0.136857  0.769325  0.983232  0.769325  0.136857

julia> sub_ap = Subpixel(ap, 5)
5×5 Subpixel{Float64, CircularAperture{Float64}} with indices 1:5×1:5:
 0.12  0.76  1.0  0.76  0.12
 0.76  1.0   1.0  1.0   0.76
 1.0   1.0   1.0  1.0   1.0
 0.76  1.0   1.0  1.0   0.76
 0.12  0.76  1.0  0.76  0.12
```

!!! note
    `photutils` offers a `center` shading method which is equivalent to using
    the `Subpixel` method with 1 subpixel. To avoid unneccessary namespace
    cluttering, we simply instruct users to use `Subpixel(ap)` instead.
"""
struct Subpixel{T,AP<:AbstractAperture{T}} <: AbstractAperture{T}
    ap::AP
    N::Int
end

function Base.show(io::IO, c::Subpixel)
    print(io, "Subpixel(")
    print(io, c.ap)
    print(io, ", $(c.N))")
end

center(sp_ap::Subpixel) = center(sp_ap.ap)
bounds(sp_ap::Subpixel) = bounds(sp_ap.ap)
overlap(sp_ap::Subpixel, i, j) = overlap(sp_ap.ap, i, j)

Subpixel(ap::AbstractAperture) = Subpixel(ap, 1)

@enum OverlapFlag Inside Outside Partial

function Base.getindex(ap::AbstractAperture{T}, idx::Vararg{Int,2}) where T
    i, j = idx
    flag = overlap(ap, i, j)
    # TODO: revisit a better way to handle subpixel apertures
    #flag === Outside && return zero(T)
    flag === Inside && return one(T)
    cx, cy = center(ap)
    return partial(ap, i - cx, j - cy)
end

# This bypasses checking aperture axes for broadcasting
Broadcast.combine_axes(ap::AbstractAperture, arrs...) = Broadcast.combine_axes(arrs...)
Broadcast.combine_axes(arr, ap::AbstractAperture) = axes(arr)

###########

"""
    photometry(::AbstractAperture, data::AbstractMatrix, [error]; [f = sum])
    photometry(::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, [error]; [f = sum])

Perform aperture photometry on `data` given aperture(s). If `error` (the
pixel-wise standard deviation) is provided, will calculate sum error. If a list
of apertures is provided the output will be a `TypedTables.Table`, otherwise a
`NamedTuple`. An optional function `f` can be passed to return additional statistics
within each aperture. This can be useful for, e.g., computing the PSF of each source. By default, just the sum within each aperture is returned.

!!! tip
    This code is automatically multi-threaded. To take advantage of this please
    make sure `JULIA_NUM_THREADS` is set before starting your runtime.
"""
function photometry(ap::AbstractAperture, data::AbstractMatrix, error; f = sum)
    cx, cy = center(ap)
    meta = (xcenter = cx, ycenter = cy)
    idxs = map(intersect, axes(ap), axes(data), axes(error))
    if any(isempty, idxs)
        if f == sum
            return (meta..., aperture_sum = 0.0, aperture_sum_err = NaN)
        else
            return (meta..., aperture_sum = 0.0, aperture_sum_err = NaN, aperture_f = 0.0)
        end
    end
    img_ap = CartesianIndices(idxs) |> Map(idx -> ap[idx] * data[idx])
    img_ap_var = CartesianIndices(idxs) |> Map(idx -> ap[idx] * error[idx]^2)

    aperture_sum = sum(img_ap)
    aperture_sum_var = sum(img_ap_var)
    aperture_sum_err = sqrt(aperture_sum_var)

    if f == sum
        return (; meta..., aperture_sum, aperture_sum_err)
    else
        aperture_f = f(img_ap)
        # Note: not supported if f returns a Tuple
        #aperture_f_var = f(img_ap_var)
        #aperture_f_err = sqrt(aperture_f_var)
        return (; meta..., aperture_sum, aperture_sum_err, aperture_f)
    end
end

function photometry(ap::AbstractAperture, data::AbstractMatrix; f = sum)
    cx, cy = center(ap)
    meta = (xcenter = cx, ycenter = cy)
    idxs = map(intersect, axes(ap), axes(data))
    if any(isempty, idxs)
        if f == sum
            return (meta..., aperture_sum = 0.0)
        else
            return (meta..., aperture_sum = 0.0, aperture_f = 0.0)
        end
    end

    img_ap = CartesianIndices(idxs) |> Map(idx -> ap[idx] * data[idx])
    aperture_sum = sum(img_ap)

    if f == sum
        return (; meta..., aperture_sum)
    else
        aperture_f = f(img_ap)
        return (; meta..., aperture_sum, aperture_f)
    end
end

function photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix, error; f = sum)
    rows = tcollect(aps |> Map(ap -> photometry(ap, data, error; f)))
    return Table(rows)
end

function photometry(aps::AbstractVector{<:AbstractAperture}, data::AbstractMatrix; f = sum)
    rows = tcollect(aps |> Map(ap -> photometry(ap, data; f)))
    return Table(rows)
end

@testsnippet photometry begin
    using Photometry.Aperture: CircularAperture,
                               CircularAnnulus,
                               EllipticalAperture,
                               EllipticalAnnulus,
                               RectangularAperture,
                               RectangularAnnulus,
                               Subpixel,
                               photometry
    const APERTURES = [
        CircularAperture,
        CircularAperture,
        CircularAperture,
        CircularAnnulus,
        CircularAnnulus,
        CircularAnnulus,
        EllipticalAperture,
        EllipticalAperture,
        EllipticalAperture,
        EllipticalAnnulus,
        EllipticalAnnulus,
        EllipticalAnnulus,
        RectangularAperture,
        RectangularAperture,
        RectangularAperture,
        RectangularAnnulus,
        RectangularAnnulus,
        RectangularAnnulus,
    ]

    const PARAMS = [
        # CircularAperture
        (0.3),
        (1),
        (3),
        # CircularAnnulus,
        (0.3, 0.5),
        (0.5, 1.0),
        (3, 5),
        # EllipticalAperture
        (0.3, 0.3, 0),
        (1, 1, 0),
        (3, 3, 0),
        # EllipticalAnnulus
        (0.3, 0.5, 0.5, 0),
        (0.5, 1.0, 0.5, 0),
        (3, 5, 4, 0),
        # RectangularAperture
        (0.3, 0.5, 0),
        (0.5, 1.0, 0),
        (3, 5, 0),
        # RectangularAnnulus
        (0.3, 0.5, 1.0, 0),
        (0.5, 1.0, 1.0, 0),
        (3, 5, 4, 0),
    ]

    # Some helpers for testing
    area(ap::CircularAperture) = π * ap.r^2
    area(ap::CircularAnnulus) = π * (ap.r_out^2 - ap.r_in^2)
    area(ap::EllipticalAperture) = π * ap.a * ap.b
    area(ap::EllipticalAnnulus) = π * ap.a_out * ap.b_out - π * ap.a_in * ap.b_in
    area(ap::RectangularAperture) = ap.w * ap.h
    area(ap::RectangularAnnulus) = ap.w_out * ap.h_out - ap.w_in * ap.h_in
end

@testitem "aperture/Aperture: outside" setup=[photometry] begin
    @testset "outside - $AP" for (AP, params) in zip(APERTURES, PARAMS)
        data = ones(10, 10)
        aperture = AP(-60, 60, params...)
        @test photometry(aperture, data).aperture_sum ≈ 0
    end
end

@testitem "aperture/Aperture: inside zeros" setup=[photometry] begin
    @testset "inside zeros - $AP" for (AP, params) in zip(APERTURES, PARAMS)
        data = zeros(40, 40)
        aperture = AP(20.0, 20.0, params...)

        table_cent = photometry(Subpixel(aperture), data)
        table_sub = photometry(Subpixel(aperture, 10), data)
        table_ex = photometry(aperture, data)


        @test table_ex.aperture_sum ≈ 0
        @test table_sub.aperture_sum ≈ 0
        @test table_cent.aperture_sum ≈ 0
    end
end

@testitem "aperture/Aperture: inside ones" setup=[photometry] begin
    @testset "inside ones - $AP" for (AP, params) in zip(APERTURES, PARAMS)
        data = ones(40, 40)
        aperture = AP(20.0, 20.0, params...)

        table_cent = photometry(Subpixel(aperture), data)
        table_sub = photometry(Subpixel(aperture, 10), data)
        table_ex = photometry(aperture, data)

        true_flux = area(aperture)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ table_ex.aperture_sum atol = 0.1

        if any(>(1), params)
            @test table_cent.aperture_sum ≤ table_ex.aperture_sum
        end
    end
end

@testitem "aperture/Aperture: interface" setup=[photometry] begin
    data = zeros(40, 40)
    err = zeros(40, 40)
    aperture = CircularAperture(20.0, 20.0, 5.0)

    f = maximum
    t1 = photometry(aperture, data)
    t1_f = photometry(aperture, data; f)
    t2 = photometry(aperture, data, err)
    t2_f = photometry(aperture, data, err; f)

    # 1.0 compat (no hasproperty function)
    hasfunc = VERSION < v"1.1" ? haskey : hasproperty

    @test !hasfunc(t1, :aperture_sum_err)
    @test !hasfunc(t1_f, :aperture_sum_err)
    @test t2.aperture_sum_err == 0
    @test t2_f.aperture_sum_err == 0
    @test propertynames(t1) == (:xcenter, :ycenter, :aperture_sum)
    @test propertynames(t1_f) == (:xcenter, :ycenter, :aperture_sum, :aperture_f)
    @test propertynames(t2) == (:xcenter, :ycenter, :aperture_sum, :aperture_sum_err)
    @test propertynames(t2_f) == (:xcenter, :ycenter, :aperture_sum, :aperture_sum_err, :aperture_f)

    apertures = CircularAperture.(20, 20, [1, 2, 3])
    t1 = photometry(apertures, data)
    t1_f = photometry(apertures, data; f)
    t2 = photometry(apertures, data, err)
    t2_f = photometry(apertures, data, err; f)

    @test !hasfunc(t1, :aperture_sum_err)
    @test !hasfunc(t1_f, :aperture_sum_err)
    @test t2.aperture_sum_err == zeros(3)
    @test t2_f.aperture_sum_err == zeros(3)
    @test propertynames(t1) == (:xcenter, :ycenter, :aperture_sum)
    @test propertynames(t1_f) == (:xcenter, :ycenter, :aperture_sum, :aperture_f)
    @test propertynames(t2) == (:xcenter, :ycenter, :aperture_sum, :aperture_sum_err)
    @test propertynames(t2_f) == (:xcenter, :ycenter, :aperture_sum, :aperture_sum_err, :aperture_f)
end

@testitem "aperture/Aperture: type stability" setup=[photometry] begin
    @testset "type stability - $AP" for (AP, params) in zip(APERTURES, PARAMS)
        data = zeros(40, 40)
        err = zeros(40, 40)
        aperture = AP(20.0, 20.0, params...)

        @inferred photometry(Subpixel(aperture), data)
        @inferred photometry(Subpixel(aperture, 10), data)
        @inferred photometry(aperture, data)

        @inferred photometry(Subpixel(aperture), data, err)
        @inferred photometry(Subpixel(aperture, 10), data, err)
        @inferred photometry(aperture, data, err)
    end
end

@testitem "aperture/Aperture: photometry - circular" setup=[photometry] begin
    function test_aperture(data, aperture)
        error = ones(size(data))

        table_cent = photometry(Subpixel(aperture), data, error)
        table_sub = photometry(Subpixel(aperture, 10), data, error)
        table_ex = photometry(aperture, data, error)

        true_flux = area(aperture)
        true_err = sqrt(true_flux)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ table_ex.aperture_sum rtol = 1e-3
        @test table_cent.aperture_sum < table_ex.aperture_sum

        @test table_ex.aperture_sum_err ≈ true_err
        @test table_sub.aperture_sum_err ≈ table_ex.aperture_sum_err rtol = 1e-3
        @test table_cent.aperture_sum_err < table_ex.aperture_sum_err
    end

    @testset "errors - CircularAperture" begin
        data = ones(40, 40)
        aperture = CircularAperture(20, 20, 10)
        test_aperture(data, aperture)
    end

    @testset "errors - CircularAnnulus" begin
        data = ones(40, 40)
        aperture = CircularAnnulus(20, 20, 8, 10)
        test_aperture(data, aperture)
    end

    @testset "partial overlap" begin
        data = ones(20, 20)
        error = ones(size(data))
        positions = [10.5 10.5; 1 1; 1 20; 20 1; 20 20]
        apertures = [CircularAperture(positions[i, :], 5) for i in axes(positions, 1)]

        table = photometry(apertures, data, error)
        @test table.aperture_sum[1] ≈ 25π
        @test all(table.aperture_sum[2:end] .< 25π)
    end
end # photometry - circular

@testitem "aperture/Aperture: photometry - elliptical" setup=[photometry] begin
    function test_elliptical_aperture(data, aperture)
        error = ones(size(data))

        table_cent = photometry(Subpixel(aperture), data, error)
        table_sub = photometry(Subpixel(aperture, 128), data, error)
        table_ex = photometry(aperture, data, error)

        true_flux = area(aperture)
        true_err = sqrt(true_flux)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ true_flux rtol = 1e-3
        @test table_cent.aperture_sum <= table_sub.aperture_sum

        @test table_ex.aperture_sum_err ≈ true_err
        @test table_sub.aperture_sum_err ≈ true_err rtol = 1e-3
        @test table_cent.aperture_sum_err <= true_err
    end

    @testset "errors - EllipticalAperture" begin
        data = ones(40, 40)
        aperture = EllipticalAperture(20, 20, 10, 10, 0)
        test_elliptical_aperture(data, aperture)

    end

    @testset "errors - EllipticalAnnulus" begin
        data = ones(40, 40)
        aperture = EllipticalAnnulus(20, 20, 8, 10, 10, 0)
        test_elliptical_aperture(data, aperture)
    end

    @testset "partial overlap elliptical aperture" begin
        data = ones(20, 20)
        error = ones(size(data))
        positions = [10.5 10.5; 1 1; 1 20; 20 1; 20 20]
        apertures = [EllipticalAperture(positions[i, :], 5, 5) for i in axes(positions, 1)]

        table = photometry(Subpixel.(apertures, 128), data, error)
        @test table.aperture_sum[1] ≈ 25π rtol = 1e-3
        @test all(table.aperture_sum[2:end] .< 25π)
    end
end # photometry - elliptical

@testitem "aperture/Aperture: photometry - rectangular" setup=[photometry] begin
    function test_aperture(data, aperture)
        error = ones(size(data))

        table_cent = photometry(Subpixel(aperture), data, error)
        table_sub = photometry(Subpixel(aperture, 10), data, error)
        table_ex = photometry(aperture, data, error)

        true_flux = area(aperture)
        true_err = sqrt(true_flux)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ true_flux rtol = 1e-2
        @test table_cent.aperture_sum <= table_sub.aperture_sum

        @test table_ex.aperture_sum_err ≈ true_err
        @test table_sub.aperture_sum_err ≈ true_err rtol = 1e-2
        @test table_cent.aperture_sum_err <= true_err
    end

    @testset "errors - RectangularAperture" begin
        data = ones(40, 40)
        aperture = RectangularAperture(20, 20, 10, 5, 0)
        test_aperture(data, aperture)
    end

    @testset "errors - RectangularAnnulus" begin
        data = ones(40, 40)
        aperture = RectangularAnnulus(20, 20, 8, 10, 4, 0)
        test_aperture(data, aperture)
    end

    @testset "partial overlap" begin
        data = ones(20, 20)
        error = ones(size(data))
        positions = [10.5 10.5; 1 1; 1 20; 20 1; 20 20]
        apertures = [RectangularAperture(positions[i, :], 10, 10, 0) for i in axes(positions, 1)]

        table = photometry(Subpixel.(apertures, 64), data, error)
        @test table.aperture_sum[1] ≈ 100 rtol = 1e-2
        @test all(table.aperture_sum[2:end] .< 100)
    end
end # photometry - circular

include("circular.jl")
include("elliptical.jl")
include("rectangle.jl")
include("overlap.jl")
include("plotting.jl")

end
