using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin

    ap_elipse = EllipticalAperture(0, 0, 20, 10, 0)
    @test bbox(ap_elipse) == (-20, 20, -10, 10)

    ap_elipse = EllipticalAperture(0,0,2,1,pi/4)
    xmin, xmax, ymin, ymax = bbox(ap_elipse)

    @test (-1.5811 - xmin) < 1e-4
    @test (1.5811 - xmax) < 1e-4
    @test (-1.5811 - ymin) < 1e-4
    @test (1.5811 - ymax) < 1e-4

end

@testset "Elliptical Aperture" begin

    e0 = EllipticalAperture(0, 0, 20, 10 , 0)
    @test sprint(show, e0) == "EllipticalAperture(0.0, 0.0, a=20.0, b=10.0, theta=0.0)"

end
