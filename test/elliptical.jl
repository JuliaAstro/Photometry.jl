using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin

    ap_elipse = EllipticalAperture(0, 0, 20, 10, 0)
    @test bbox(ap_elipse) == (-20, 20, -10, 10)

    ap_elipse = EllipticalAperture(0,0,2,1,45)
    xmin, xmax, ymin, ymax = bbox(ap_elipse)

    @test xmin ≈ -1.5811 atol=1e-4
    @test xmax ≈ 1.5811 atol=1e-4
    @test ymax ≈ 1.5811 atol=1e-4
    @test ymin ≈ -1.5811 atol=1e-4

end

@testset "Elliptical Aperture" begin

    e0 = EllipticalAperture(0, 0, 20, 10 , 0)
    @test sprint(show, e0) == "EllipticalAperture(0, 0, a=20, b=10, theta=0°)"

end
