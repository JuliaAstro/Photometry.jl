using Photometry.Aperture: edges,
                           bbox

@testset "Apertures" begin

    ap_rectangle = RectangularAperture(0, 0, 10, 5, 50)
    @test edges(ap_rectangle) == (-5, 5, -3, 3)
    @test bbox(ap_rectangle) == (-5, 5, -3, 3)
    @test size(ap_rectangle) == (7, 11)

    ap_rectangle_annulus = RectangularAnnulus(0, 0, 15, 20, 19, 21, 50)
    @test edges(ap_rectangle_annulus) == (-10, 10, -11, 11)    
    @test bbox(ap_rectangle_annulus) == (-10, 10, -11, 11)
    @test size(ap_rectangle_annulus) == (23, 21)
end

@testset "Rectangular Aperture" begin

    r0 = RectangularAperture(0, 0, 15, 17 , 0)
    @test sprint(show, r0) == "RectangularAperture(0, 0, w=15, h=17, theta=0)"

    r1 = RectangularAperture(2, 10, 15, 19 , 50)
    @test sprint(show, r1) == "RectangularAperture(2, 10, w=15, h=19, theta=50)"

    r = RectangularAperture(0, 0, 5, 10, 10)
    @test mask(r, method = :center) == mask(r, method = (:subpixel, 1))

end

@testset "Rectangular Annulus" begin

    r0 = RectangularAnnulus(5, 2, 8, 9 ,15, 17, 0)
    @test sprint(show, r0) == "RectangularAnnulus(5, 2, w_in=8, w_out=9, h_in=15, h_out=17, theta=0)"

    r1 = RectangularAnnulus(0, 0, 1, 7 ,20, 18, 55)
    @test sprint(show, r1) == "RectangularAnnulus(0, 0, w_in=1, w_out=7, h_in=20, h_out=18, theta=55)"

    r = RectangularAnnulus(6, 3, 15, 14, 20, 26, 8)
    @test mask(r, method = :center) == mask(r, method = (:subpixel, 1))

end