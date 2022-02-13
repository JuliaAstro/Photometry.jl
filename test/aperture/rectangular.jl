using Photometry.Aperture: bounds

@testset "Apertures" begin
    ap_rect = RectangularAperture(50, 40, 10, 10, 0)
    @test center(ap_rect) == (50, 40)
    @test bounds(ap_rect) == (45, 55, 35, 45)
    @test size(ap_rect) == (11, 11)
    @test size(ap_rect, 1) == 11
    @test RectangularAperture([50, 40], 10, 10, 0) == ap_rect

    ap_ann = RectangularAnnulus(50, 40, 5, 10, 10, 0)
    @test bounds(ap_ann) == (45, 55, 35, 45)
    @test size(ap_ann) == (11, 11)
    @test RectangularAnnulus([50, 40], 5, 10, 10, 0) == ap_ann
end

@testset "Rectangle Aperture" begin
    ap0 = RectangularAperture(0, 0, 0, 0, 0)

    @test sprint(show, ap0) == "RectangularAperture(0, 0, w=0, h=0, θ=0°)"

    ap1 = RectangularAperture(0, 0, 1, 1, 0)

    @test sprint(show, ap1) == "RectangularAperture(0, 0, w=1, h=1, θ=0°)"
end

@testset "Rectangle Annulus" begin

    ap1 = RectangularAnnulus(0, 0, 1, 1, 1, 0)
    @test center(ap1) == (0, 0)
    @test sprint(show, ap1) == "RectangularAnnulus(0.0, 0.0, w_in=1.0, w_out=1.0, h_in=1.0, h_out=1.0, θ=0.0°)"
end
