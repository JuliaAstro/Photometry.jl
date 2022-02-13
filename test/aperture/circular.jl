using Photometry.Aperture: bounds, center

@testset "Apertures" begin
    ap_circ = CircularAperture(50, 40, 10)
    @test center(ap_circ) == (50, 40)
    @test bounds(ap_circ) == (40, 60, 30, 50)
    @test size(ap_circ) == map(length, axes(ap_circ)) == (21, 21)
    @test size(ap_circ, 1) == 21
    @test all(axes(ap_circ) .== (40:60, 30:50))
    @test eachindex(ap_circ) == CartesianIndex(40, 30):CartesianIndex(60, 50)

    @test CircularAperture([50, 40], 10) == ap_circ

    ap_ann = CircularAnnulus(50, 40, 5, 10)
    @test center(ap_ann) == (50, 40)
    @test bounds(ap_ann) == (40, 60, 30, 50)
    @test size(ap_ann) == (21, 21)
    @test all(axes(ap_ann) .== (40:60, 30:50))
    @test eachindex(ap_ann) == CartesianIndex(40, 30):CartesianIndex(60, 50)

    @test CircularAnnulus([50, 40], 5, 10) == ap_ann
end

@testset "broadcasting" begin
    ap = CircularAperture(3, 3, 2.5)
    data = ones(5, 7)
    weighted = ap .* data
    @test weighted == data .* ap # commutative
    @test sum(weighted) == sum(ap) == photometry(ap, data).aperture_sum
    @test all(iszero, weighted[:, 6:7])

end

@testset "Circle Aperture" begin
    c0 = CircularAperture(0, 0, 0)

    @test sprint(show, c0) == "CircularAperture(0, 0, r=0)"

    c1 = CircularAperture(0, 0, 1)

    @test sprint(show, c1) == "CircularAperture(0, 0, r=1)"

    sub_c1 = Subpixel(c1, 5)
    @test sprint(show, sub_c1) == "Subpixel(CircularAperture(0, 0, r=1), 5)"
end

@testset "Circle Annulus" begin
    c0 = CircularAnnulus(0, 0, 0, 0)

    @test sprint(show, c0) == "CircularAnnulus(0, 0, r_in=0, r_out=0)"

    c1 = CircularAnnulus(0, 0, 0, 1)

    @test sprint(show, c1) == "CircularAnnulus(0, 0, r_in=0, r_out=1)"
end
