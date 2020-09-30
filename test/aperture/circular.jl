using Photometry.Aperture: bounds

@testset "Apertures" begin
    ap_circ = CircularAperture(50, 40, 10)
    @test bounds(ap_circ) == (40, 60, 30, 50)
    @test size(ap_circ) == (21, 21)

    @test CircularAperture([50, 40], 10) == ap_circ

    ap_ann = CircularAnnulus(50, 40, 5, 10)
    @test bounds(ap_ann) == (40, 60, 30, 50)
    @test size(ap_ann) == (21, 21)

    @test CircularAnnulus([50, 40], 5, 10) == ap_ann
end

@testset "Circle Aperture" begin
    c0 = CircularAperture(0, 0, 0)

    @test sprint(show, c0) == "CircularAperture(0, 0, r=0)"

    c1 = CircularAperture(0, 0, 1)

    @test sprint(show, c1) == "CircularAperture(0, 0, r=1)"
end

@testset "Circle Annulus" begin
    c0 = CircularAnnulus(0, 0, 0, 0)

    @test sprint(show, c0) == "CircularAnnulus(0, 0, r_in=0, r_out=0)"

    c1 = CircularAnnulus(0, 0, 0, 1)

    @test sprint(show, c1) == "CircularAnnulus(0, 0, r_in=0, r_out=1)"
end
