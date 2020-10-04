APERTURES = [
    CircularAperture,
    CircularAnnulus,
    EllipticalAperture,
    EllipticalAnnulus,
    RectangularAperture,
    RectangularAnnulus
]

PARAMS = [
    (3),
    (3, 5),
    (3, 3, 0),
    (3, 5, 4, 0),
    (3, 5, 0),
    (3, 5, 4, 0)
]

###########################
# Some helpers for testing
area(ap::CircularAperture) = π * ap.r^2
area(ap::CircularAnnulus) = π * (ap.r_out^2 - ap.r_in^2)
area(ap::EllipticalAperture) = π * ap.a * ap.b
area(ap::EllipticalAnnulus) = π * ap.a_out * ap.b_out - π * ap.a_in * ap.b_in
area(ap::RectangularAperture) = ap.w * ap.h
area(ap::RectangularAnnulus) = ap.w_out * ap.h_out - ap.w_in * ap.h_in

@testset "outside - $AP" for (AP, params) in zip(APERTURES, PARAMS)
    data = ones(10, 10)
    aperture = AP(-60, 60, params...)
    @test photometry(aperture, data).aperture_sum ≈ 0
end

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

@testset "inside ones - $AP" for (AP, params) in zip(APERTURES, PARAMS)
    data = ones(40, 40)
    aperture = AP(20.0, 20.0, params...)

    table_cent = photometry(Subpixel(aperture), data)
    table_sub = photometry(Subpixel(aperture, 10), data)
    table_ex = photometry(aperture, data)

    true_flux = area(aperture)

    @test table_ex.aperture_sum ≈ true_flux
    @test table_sub.aperture_sum ≈ table_ex.aperture_sum atol = 0.1
    @test table_cent.aperture_sum ≤ table_ex.aperture_sum

end



@testset "interface" begin
    data = zeros(40, 40)
    err = zeros(40, 40)
    aperture = CircularAperture(20.0, 20.0, 5.0)

    t1 = photometry(aperture, data)
    t2 = photometry(aperture, data, err)

    # 1.0 compat (no hasproperty function)
    hasfunc = VERSION < v"1.1" ? haskey : hasproperty

    @test !hasfunc(t1, :aperture_sum_err)
    @test t2.aperture_sum_err == 0

    apertures = CircularAperture.(20, 20, [1, 2, 3])
    t1 = photometry(apertures, data)
    t2 = photometry(apertures, data, err)

    @test !hasfunc(t1, :aperture_sum_err)
    @test t2.aperture_sum_err == zeros(3)
end

# @testset "type stability - $AP" for (AP, params) in zip(APERTURES, PARAMS)
#     data = zeros(40, 40)
#     err = zeros(40, 40)
#     aperture = AP(20.0, 20.0, params...)

#     @inferred photometry(Subpixel(aperture), data)
#     @inferred photometry(Subpixel(aperture, 10), data)
#     @inferred photometry(aperture, data)

#     @inferred photometry(Subpixel(aperture), data, err)
#     @inferred photometry(Subpixel(aperture, 10), data, err)
#     @inferred photometry(aperture, data, err)
# end

@testset "photometry - circular" begin
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

@testset "photometry - elliptical" begin
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

end # photometry elliptical

@testset "photometry - rectangular" begin
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
