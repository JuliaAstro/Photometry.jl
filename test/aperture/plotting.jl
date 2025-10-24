using Photometry.Aperture:
    CircularAperture,
    CircularAnnulus,
    EllipticalAperture,
    EllipticalAnnulus,
    RectangularAperture,
    RectangularAnnulus

using RecipesBase: apply_recipe

using Statistics

const APERTURES = [
    CircularAperture(3, 3, 3),
    CircularAnnulus(3, 3, 2, 4),
    EllipticalAperture(3, 3, 4, 2, 45),
    EllipticalAnnulus(3, 3, 3, 4, 2, -26),
    RectangularAperture(3, 3, 3, 4, 15),
    RectangularAnnulus(3, 3, 3, 4, 2, -5),
]

@testset "aperture/plotting: Aperture Plots" begin
    @testset "Aperture Plots - $(typeof(ap))" for ap in APERTURES
        recipes = apply_recipe(Dict{Symbol,Any}(), ap)
        for rec in recipes
            @test getfield(rec, 1) == Dict{Symbol,Any}(
                :seriescolor => :red,
                :fillcolor => nothing,
                :linecolor => :match,
                :label => "",
                :seriestype => :shape,
                :aspect_ratio => :equal)

            # test to make sure our position is correct (should be +0.5 the given (x,y))
            points = rec.args[1]
            x_tot = y_tot = 0
            for point in points
                x_tot += point[1]
                y_tot += point[2]
            end
            N = length(points)
            x_mean = x_tot / N
            y_mean = y_tot / N
            @test x_mean ≈ 3.5 atol = x_mean / sqrt(N)
            @test y_mean ≈ 3.5 atol = y_mean / sqrt(N)
        end
    end
end
