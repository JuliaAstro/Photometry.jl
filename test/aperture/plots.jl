using RecipesBase: apply_recipe
using Statistics

APERTURES = [
    CircularAperture(3, 3, 3),
    CircularAnnulus(3, 3, 2, 4),
    EllipticalAperture(3, 3, 4, 2, 45),
    EllipticalAnnulus(3, 3, 3, 4, 2, -26),
    RectangularAperture(3, 3, 3, 4, 15),
    RectangularAnnulus(3, 3, 3, 4, 2, -5)
]

@testset "Aperture Plots - $(typeof(ap))" for ap in APERTURES
    rec = apply_recipe(Dict{Symbol,Any}(), ap)
    for i in length(rec)
        @test getfield(rec[i], 1) == Dict{Symbol,Any}(
            :seriescolor => :red,
            :label => "",
            :seriestype => :path,
            :aspect_ratio => :equal)

        # test to make sure our position is correct (should be +0.5 the given (x,y))
        x, y = rec[i].args
        @test mean(x[1:end - 1]) ≈ 3.5
        @test mean(y[1:end - 1]) ≈ 3.5
    end
end
