using Photometry.Aperture: edges,
                           bbox,
                           oblique_coefficients

@testset "Apertures" begin

    ap_elipse = EllipticalAperture(0, 0, 20, 10, 0)
    @test bbox(ap_elipse) == (-20, 20, -10, 10)

    ap_elipse = EllipticalAperture(0, 0, 2, 1, 45)
    xmin, xmax, ymin, ymax = bbox(ap_elipse)

    @test xmin == -2
    @test xmax == 2
    @test ymax == 2
    @test ymin == -2

end

@testset "Elliptical Aperture" begin

    e = EllipticalAperture(0, 0, 20, 10, 0)
    @test sprint(show, e) == "EllipticalAperture(0, 0, a=20, b=10, θ=0°)"

    e = EllipticalAperture(0, 0, 10, 10, 0)
    @test mask(e, method = :center) == mask(e, method = (:subpixel, 1))

    # test modding of angle
    @test EllipticalAperture(0, 0, 3, 4, 380).theta == 20

    @test_throws ErrorException EllipticalAperture(0, 0, -2, 4, 0)
    @test_throws ErrorException EllipticalAperture(0, 0, 2, -3, 0)

end

@testset "oblique_coefficients" begin
    @test all(oblique_coefficients(2, 2, 0) .≈ (0.25, 0.25, 0.0))
    @test all(oblique_coefficients(2, 2, 90) .≈ (0.25, 0.25, 0.0))
    @test all(oblique_coefficients(2, 1, 30) .≈ (7 / 16, 13 / 16, -6sqrt(3) / 16))
end


@testset "Elliptical Annulus" begin

    e0 = EllipticalAnnulus(0, 0, 16, 4, 45, 2)
    @test sprint(show, e0) == "EllipticalAnnulus(0, 0, a=16, b=4, θ=45°, factor=2)"

    e = EllipticalAnnulus(0, 0, 10, 10, 45, 2)
    @test mask(e, method = :center) == mask(e, method = (:subpixel, 1))
    @test_throws ErrorException mask(e, method = :exact)

    @test EllipticalAnnulus(0, 0, 4, 3, 380, 2).theta == 20

end

@testset "Elliptical Annulus bounding box" begin

    e = EllipticalAnnulus(0, 0, 16, 4, 0, 2)
    @test all(bbox(e) == (-32, 32, -8, 8))

end
