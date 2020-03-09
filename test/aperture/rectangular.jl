using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin
    ap_rect = RectangularAperture(50, 40, 10, 10, 0)
    @test edges(ap_rect) == (-5.5, 5.5, -5.5, 5.5)
    @test bbox(ap_rect) == (45, 55, 35, 45)
    @test size(ap_rect) == (11, 11)

    ap_ann = RectangularAnnulus(50, 40, 5, 10, 10, 0)
    @test edges(ap_ann) == (-5.5, 5.5, -5.5, 5.5)
    @test bbox(ap_ann) == (45, 55, 35, 45)
    @test size(ap_ann) == (11, 11)
end

@testset "Rectangle Aperture" begin
    ap0 = RectangularAperture(0, 0, 0, 0, 0)

    @test sprint(show, ap0) == "RectangularAperture(0, 0, w=0, h=0, θ=0°)"

    ap1 = RectangularAperture(0, 0, 1, 1, 0)

    @test sprint(show, ap1) == "RectangularAperture(0, 0, w=1, h=1, θ=0°)"

    @test_throws ErrorException RectangularAperture(0, 0, -12, 2, 0)
    @test_throws ErrorException RectangularAperture(0, 0, 5, -2, 0)
end

@testset "Rectangle Annulus" begin

    ap1 = RectangularAnnulus(0, 0, 1, 1, 1, 0)

    @test sprint(show, ap1) == "RectangularAnnulus(0.0, 0.0, w_in=1.0, w_out=1.0, h_in=1.0, h_out=1.0, θ=0.0°)"

    @test_throws ErrorException RectangularAnnulus(0, 0, 0, 0, 0, 0)
    @test_throws ErrorException RectangularAnnulus(0, 0, -12, 2, 3, 0)
    @test_throws ErrorException RectangularAnnulus(0, 0, 4, 2, 3, 0)
    @test_throws ErrorException RectangularAnnulus(0, 0, 5, -2, 4, 0)
    @test_throws ErrorException RectangularAnnulus(0, 0, 5, 8, -4, 0)
end
