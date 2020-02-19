using RecipesBase

@recipe function f(c::CircularAperture, npoints = 1000)
    seriestype := :path
    aspect_ratio --> :equal
    seriescolor --> :red
    label --> ""
    t = range(0, 2pi, length = npoints)
    
    x = @. c.x + c.r * sin(t) + 0.5
    y = @. c.y + c.r * cos(t) + 0.5


    x, y
end

@recipe function f(c::CircularAnnulus, npoints = 1000)
    t = range(0, 2pi, length = npoints)
    seriestype := :path
    seriescolor --> :red
    aspect_ratio --> :equal
    label --> ""
    
    # outer ring
    @series begin
        x = @. c.x + c.r_out * cos(t) + 0.5
        y = @. c.y + c.r_out * sin(t) + 0.5
        x, y
    end

    # inner ring
    @series begin
        x = @. c.x + c.r_in * cos(t) + 0.5
        y = @. c.y + c.r_in * sin(t) + 0.5
        x, y
    end


end


@recipe function f(e::EllipticalAperture, npoints = 1000)
    seriestype := :path
    seriescolor --> :red
    aspect_ratio --> :equal
    label --> ""

    t = range(0, 2pi, length = npoints)
    x = @. e.x + e.a * cos(t) * cosd(e.theta) - e.b * sin(t) * sind(e.theta) + 0.5
    y = @. e.y + e.a * cos(t) * sind(e.theta) + e.b * sin(t) * cosd(e.theta) + 0.5

    x, y
end
