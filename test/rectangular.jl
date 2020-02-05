using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin

    ap_rectangle = RectangularAperture(0, 0, 10, 5, 50)
    @test bbox(ap_rectangle) == (-5, 5, -3, 3)

    ap_rectangle_annulus = RectangularAnnulus(0, 0, 5, 10, 6, 12, 50)
    @test bbox(ap_rectangle) == (-5, 5, -6, 6)

    
end

@testset "Rectangular Aperture" begin

    r1 = RectangularAperture(0, 0, 15, 17 , 0)
    @test sprint(show, r1) == "RectangularAperture(0.0, 0.0, w=15.0, h=17.0, theta=0.0)"

end

@testset "Rectangular Annulus" begin

    r2 = RectangularAnnulus(5, 2, 8, 9 , 0)
    @test sprint(show, r2) == "RectangularAperture(5.0, 2.0, w=8.0, h=9.0, theta=0.0)"

end