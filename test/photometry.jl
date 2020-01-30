APERTURES = [
    CircularAperture,
    CircularAnnulus
]

PARAMS = [
    (3),
    (3, 5)
]

###########################
# Some helpers for testing
area(c::CircularAperture) = π * c.r^2
area(c::CircularAnnulus) = π * (c.r_out^2 - c.r_in^2)


@testset "outside - $AP" for (AP, params) in zip(APERTURES, PARAMS)
    data = ones(10, 10)
    aperture = AP(-60, 60, params...)
    @test aperture_photometry(aperture, data).aperture_sum ≈ 0
end

@testset "inside zeros - $AP" for (AP, params) in zip(APERTURES, PARAMS)
    data = zeros(40, 40)
    aperture = AP(20.0, 20.0, params...)

    table_cent = aperture_photometry(aperture, data, method = :center)
    table_sub = aperture_photometry(aperture, data, method = (:subpixel, 10))
    table_ex = aperture_photometry(aperture, data, method = :exact)


    @test table_ex.aperture_sum ≈ 0
    @test table_sub.aperture_sum ≈ 0
    @test table_cent.aperture_sum ≈ 0
    
end

@testset "inside ones - $AP" for (AP, params) in zip(APERTURES, PARAMS)
    data = ones(40, 40)
    aperture = AP(20.0, 20.0, params...)

    table_cent = aperture_photometry(aperture, data, method = :center)
    table_sub = aperture_photometry(aperture, data, method = (:subpixel, 10))
    table_ex = aperture_photometry(aperture, data, method = :exact)

    true_flux = area(aperture)

    @test table_ex.aperture_sum ≈ true_flux
    @test table_sub.aperture_sum ≈ table_ex.aperture_sum atol = 0.1
    @test table_cent.aperture_sum < table_ex.aperture_sum
    
end
