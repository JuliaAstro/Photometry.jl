using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin

    ap_rectangle = RectangularAperture(0, 0, 10, 5, 50)
    @test bbox(ap_rectangle) == (-5, 5, -3, 3)

    ap_rectangle_annulus = RectangularAnnulus(0, 0, 15, 20, 19, 21, 50)
    @test bbox(ap_rectangle_annulus) == (-10, 10, -11, 11)

    
end

@testset "Rectangular Aperture" begin

    r1 = RectangularAperture(0, 0, 15, 17 , 0)
    @test sprint(show, r1) == "RectangularAperture(0, 0, w=15, h=17, theta=0)"

end

@testset "Rectangular Annulus" begin

    r2 = RectangularAnnulus(5, 2, 8, 9 ,15, 17, 0)
    @test sprint(show, r2) == "RectangularAnnulus(5, 2, w_in=8, w_out=9, h_in=15, h_out=17, theta=0)"

end