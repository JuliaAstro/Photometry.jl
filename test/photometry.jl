@testset "Apertures" begin  
    ap_circ = CircularAperture(50, 40, 10)
    @test Photometry.edges(ap_circ) == (-10.5, 10.5, -10.5, 10.5)
    @test Photometry.bbox(ap_circ) == (40, 60, 30, 50)
    @test size(ap_circ) == (21, 21)

    ap_ann = CircularAnnulus(50, 40, 5, 10)
    @test Photometry.edges(ap_ann) == (-10.5, 10.5, -10.5, 10.5)
    @test Photometry.bbox(ap_ann) == (40, 60, 30, 50)
    @test size(ap_ann) == (21, 21)


end

@testset "ones and zeros" begin
    ap_circ = CircularAperture(51, 51, 10)
    ap_ann = CircularAnnulus(51, 51, 5, 10)

    data = zeros(100, 100)
    res = aperture_photometry(ap_circ, data)
    @test res.aperture_sum ≈ 0
    res = aperture_photometry(ap_ann, data)
    @test res.aperture_sum ≈ 0

    data = ones(100, 100)
    res = aperture_photometry(ap_circ, data)
    @test res.aperture_sum ≈ 100π
    res = aperture_photometry(ap_ann, data)
    @test res.aperture_sum ≈ 75π
    
end
