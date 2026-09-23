using Photometry.Aperture:
    CircularAperture,
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

@testset "aperture/Aperture: outside" begin
    @testset "outside - $AP" for (AP, params) in zip(APERTURES, PARAMS)
        data = ones(10, 10)
        aperture = AP(-60, 60, params...)
        @test photometry(aperture, data).aperture_sum ≈ 0
    end
end

@testset "aperture/Aperture: inside zeros" begin
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

@testset "aperture/Aperture: inside ones" begin
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

@testset "aperture/Aperture: interface" begin
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

    # Test that the custom function gets an object with matrix properties
    matrix_f = img_ap -> (
        size = size(img_ap),
        axes = axes(img_ap),
        first = img_ap[begin, begin],
        sum = sum(img_ap),
    )
    t1_matrix_f = photometry(aperture, data; f = matrix_f)
    t2_matrix_f = photometry(aperture, data, err; f = matrix_f)

    @test t1_matrix_f.aperture_f.axes == (Base.OneTo(t1_matrix_f.aperture_f.size[1]), Base.OneTo(t1_matrix_f.aperture_f.size[2]))
    @test t1_matrix_f.aperture_f.sum == t1_matrix_f.aperture_sum
    @test t2_matrix_f.aperture_f.axes == t1_matrix_f.aperture_f.axes
    @test t2_matrix_f.aperture_f.sum == t2_matrix_f.aperture_sum
end

@testset "aperture/Aperture: type stability" begin
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

@testset "aperture/Aperture: photometry - circular" begin
    function test_aperture(data, aperture)
        error = ones(size(data))

        table_cent = photometry(Subpixel(aperture), data, error)
        table_sub = photometry(Subpixel(aperture, 10), data, error)
        table_ex = photometry(aperture, data, error)

        true_flux = area(aperture)
        true_err = sqrt(true_flux)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ table_ex.aperture_sum rtol = 1.0e-3
        @test table_cent.aperture_sum < table_ex.aperture_sum

        @test table_ex.aperture_sum_err ≈ true_err
        @test table_sub.aperture_sum_err ≈ table_ex.aperture_sum_err rtol = 1.0e-3
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

@testset "aperture/Aperture: photometry - elliptical" begin
    function test_elliptical_aperture(data, aperture)
        error = ones(size(data))

        table_cent = photometry(Subpixel(aperture), data, error)
        table_sub = photometry(Subpixel(aperture, 128), data, error)
        table_ex = photometry(aperture, data, error)

        true_flux = area(aperture)
        true_err = sqrt(true_flux)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ true_flux rtol = 1.0e-3
        @test table_cent.aperture_sum <= table_sub.aperture_sum

        @test table_ex.aperture_sum_err ≈ true_err
        @test table_sub.aperture_sum_err ≈ true_err rtol = 1.0e-3
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
        @test table.aperture_sum[1] ≈ 25π rtol = 1.0e-3
        @test all(table.aperture_sum[2:end] .< 25π)
    end
end # photometry - elliptical

@testset "aperture/Aperture: photometry - rectangular" begin
    function test_aperture(data, aperture)
        error = ones(size(data))

        table_cent = photometry(Subpixel(aperture), data, error)
        table_sub = photometry(Subpixel(aperture, 10), data, error)
        table_ex = photometry(aperture, data, error)

        true_flux = area(aperture)
        true_err = sqrt(true_flux)

        @test table_ex.aperture_sum ≈ true_flux
        @test table_sub.aperture_sum ≈ true_flux rtol = 1.0e-2
        @test table_cent.aperture_sum <= table_sub.aperture_sum

        @test table_ex.aperture_sum_err ≈ true_err
        @test table_sub.aperture_sum_err ≈ true_err rtol = 1.0e-2
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
        @test table.aperture_sum[1] ≈ 100 rtol = 1.0e-2
        @test all(table.aperture_sum[2:end] .< 100)
    end
end # photometry - rectangular

@testset "aperture/Aperture: empty aperture list" begin
    # a blank frame in a detection -> photometry pipeline yields no apertures;
    # the result is an empty table with the usual columns rather than an error
    data = ones(10, 10)
    err = ones(10, 10)
    aps = CircularAperture{Float64}[]

    for (table, names) in (
            (photometry(aps, data), (:xcenter, :ycenter, :aperture_sum)),
            (photometry(aps, data; f = maximum), (:xcenter, :ycenter, :aperture_sum, :aperture_f)),
            (photometry(aps, data, err), (:xcenter, :ycenter, :aperture_sum, :aperture_sum_err)),
            (photometry(aps, data, err; f = maximum), (:xcenter, :ycenter, :aperture_sum, :aperture_sum_err, :aperture_f)),
            (photometry(aps, data; f = collect), (:xcenter, :ycenter, :aperture_sum, :aperture_f)),
        )
        @test length(table) == 0
        @test propertynames(table) == names
    end
    @test eltype(photometry(aps, data).aperture_sum) == Float64
end

@testset "aperture/Aperture: integer apertures" begin
    # aperture weights are floating point even when the aperture parameters and
    # the data are integers, e.g. positions straight from `extract_sources`
    data = fill(7, 30, 30)
    ap = CircularAperture(10, 12, 3)
    cutout = photometry(ap, data; f = collect).aperture_f
    @test eltype(cutout) == Float64
    @test sum(cutout) ≈ photometry(ap, data).aperture_sum ≈ 7 * 9π
    @test photometry([ap], data; f = collect).aperture_f[1] == cutout
end

@testset "aperture/Aperture: broadcasting" begin
    ap = CircularAperture(3, 3, 2.5) # axes 1:5 × 1:5
    data = ones(5, 7)

    # with an array the array's axes win, whatever the aperture's bounds
    @test axes(ap .* data) == axes(data)
    @test ap .* data == [ap[i, j] for i in 1:5, j in 1:7]
    @test data .* ap == ap .* data
    @test (+).(data, ap, ap) == 1 .+ 2 .* (ap .* data)
    @test (@inferred Base.Broadcast.combine_axes(ap, data)) == axes(data)

    # alone or with scalars, the aperture's own axes apply
    @test axes(sqrt.(ap)) == axes(ap)
    @test sqrt.(ap) == sqrt.(collect(ap))
    @test axes(ap .* 2) == axes(ap)
    @test ap .* 2 == 2 .* ap == 2 .* collect(ap)
    @test (@inferred Base.Broadcast.combine_axes(ap, 2)) == axes(ap)

    # two apertures broadcast like any two arrays
    same = CircularAperture(3, 3, 2.4) # also 1:5 × 1:5
    @test axes(ap .* same) == axes(ap)
    @test ap .* same == collect(ap) .* collect(same)
    @test (+).(ap, 1, same) == collect(ap) .+ 1 .+ collect(same)
    @test_throws DimensionMismatch ap .* CircularAperture(3, 3, 1.5) # 2:4 × 2:4
    @test data .* ap .* same == (ap .* data) .* same
end
