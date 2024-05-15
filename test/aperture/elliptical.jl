using Photometry.Aperture: bounds, center, oblique_coefficients

@testset "Apertures" begin
    ap_ellipse = EllipticalAperture(0, 0, 20, 10, 0)
    @test center(ap_ellipse) == (0, 0)
    @test bounds(ap_ellipse) == (-20, 20, -10, 10)
    @test size(ap_ellipse) == (41, 21)
    @test size(ap_ellipse, 1) == 41
    @test all(axes(ap_ellipse) .== (-20:20, -10:10))
    @test eachindex(ap_ellipse) == CartesianIndex(-20, -10):CartesianIndex(20, 10)

    @test EllipticalAperture([0, 0], 20, 10, 0) == ap_ellipse

    ap_ellipse = EllipticalAperture(0, 0, 2, 1, 45)
    @test bounds(ap_ellipse) == (-1, 1, -1, 1)
end

@testset "Elliptical Aperture" begin
    e = EllipticalAperture(0, 0, 20, 10, 0)
    @test sprint(show, e) == "EllipticalAperture(0, 0, a=20, b=10, θ=0°)"
end

@testset "oblique_coefficients" begin
    @test all(oblique_coefficients(2, 2, 0) .≈ (0.25, 0.25, 0.0))
    @test all(oblique_coefficients(2, 2, 90) .≈ (0.25, 0.25, 0.0))
    @test all(oblique_coefficients(2, 1, 30) .≈ (7 / 16, 13 / 16, -6sqrt(3) / 16))
end

@testset "Elliptical Annulus" begin
    e0 = EllipticalAnnulus(0, 0, 8, 16, 4, 45)
    @test sprint(show, e0) == "EllipticalAnnulus(0.0, 0.0, a_in=8.0, a_out=16.0, b_in=2.0, b_out=4.0, θ=45.0°)"
    @test EllipticalAnnulus([0, 0], 8, 16, 4, 45) == EllipticalAnnulus(0, 0, 8, 16, 4, 45)
end

@testset "Elliptical Annulus bounding box" begin
    e = EllipticalAnnulus(0, 0, 8, 16, 4, 0)
    @test center(e) == (0, 0)
    @test bounds(e) == (-16, 16, -4, 4)
    @test size(e) == (33, 9)
    @test all(axes(e) .== (-16:16, -4:4))
    @test eachindex(e) == CartesianIndex(-16, -4):CartesianIndex(16, 4)
end

@testset "regression circular ellipse" begin
    # some weird bug where centered on-grid past a certain size (and rotated) would fail
    e = EllipticalAperture(5.5, 5.5, 4, 4, 20)
    @test sum(e) ≈ photometry(e, ones(9, 9)).aperture_sum # just a test that no errors occur
end
