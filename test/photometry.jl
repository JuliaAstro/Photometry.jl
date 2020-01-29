@testset "Trivial" for aperture in [CircularAperture(50, 50, 10), CircularAnnulus(50, 50, 5, 10)]
    data = zeros(100, 100)
    res = aperture_photometry(aperture, data)
    @test res.aperture_sum ≈ 0

    data = ones(100, 100)
    res = aperture_photometry(aperture, data)
    expected = aperture isa CircularAperture ? 100π : 75π
    @test_broken res.aperture_sum ≈ expected
end
