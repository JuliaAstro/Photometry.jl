using RecipesBase

@recipe function f(c::CircularAperture, npoints=1000)
    seriestype := :path
    aspect_ratio --> :equal
    label --> ""
    seriescolor --> :red

    t = range(0, 2pi, length=npoints)
    x = c.x .+ c.r .* sin.(t)
    y = c.y .+ c.r .* cos.(t)


    x, y
end

@recipe function f(c::CircularAnnulus, npoints=1000)
    t = range(0, 2pi, length=npoints)
    seriescolor --> :red
    seriestype := :path
    aspect_ratio --> :equal
    label --> ""
    
    # outer ring
    @series begin
        x = c.x .+ c.r_out .* sin.(t)
        y = c.y .+ c.r_out .* cos.(t)
        x, y
    end

    # inner ring
    @series begin
        x = c.x .+ c.r_in .* sin.(t)
        y = c.y .+ c.r_in .* cos.(t)
        x, y
    end


end
