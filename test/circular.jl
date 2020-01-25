@testset "Circle Aperture" begin
    c0 = CircularAperture((0, 0), 0)

    @test area(c0) == 0

    @test sprint(show, c0) == "CircularAperture(0, 0, r=0)"

    c1 = CircularAperture((0, 0), 1)

    @test area(c1) ≈ π

    @test sprint(show, c1) == "CircularAperture(0, 0, r=1)"

    # cbig = CircularAperture{BigFloat}((0, 1), 2)

    # @test area(cbig) ≈ BigFloat(4) * π

    c = CircularAperture((0, 0), 10)
    @test mask(c, method = :center) == mask(c, method = 1)

end

@testset "Circle Annulus" begin
    c0 = CircularAnnulus((0, 0), 0, 0)

    @test area(c0) == 0

    @test sprint(show, c0) == "CircularAnnulus(0, 0, r_in=0, r_out=0)"

    c1 = CircularAnnulus((0, 0), 0, 1)

    @test area(c1) ≈ π

    @test sprint(show, c1) == "CircularAnnulus(0, 0, r_in=0, r_out=1)"

    # cbig = CircularAnnulus{BigFloat}((0, 1), 1, 2)

    # @test area(cbig) ≈ BigFloat(3) * π

    c = CircularAnnulus((0, 0), 5, 10)
    @test mask(c, method = :center) == mask(c, method = 1)

end