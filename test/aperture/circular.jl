using Photometry.Aperture: CircularAperture, bounds, center

@testset "Apertures" begin
    ap_circ = CircularAperture(50, 40, 10)
    @test center(ap_circ) == (50, 40)
    @test bounds(ap_circ) == (40, 60, 30, 50)
end
