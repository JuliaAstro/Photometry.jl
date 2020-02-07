using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin

    ap_elipse = EllipticalAperture(0, 0, 20, 10, 0)
    @test bbox(ap_elipse) == (-20, 20, -10, 10)

    ap_elipse = EllipticalAperture(0,0,2,1,45)
    xmin, xmax, ymin, ymax = bbox(ap_elipse)

    @test xmin == -2
    @test xmax == 2
    @test ymax == 2
    @test ymin == -2

end

@testset "Elliptical Aperture" begin

    e0 = EllipticalAperture(0, 0, 20, 10 , 0)
    @test sprint(show, e0) == "EllipticalAperture(0, 0, a=20, b=10, theta=0°)"

end
