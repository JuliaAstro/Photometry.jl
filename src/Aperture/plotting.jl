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

@recipe function f(e::EllipticalAnnulus, npoints = 1000)
    t = range(0, 2pi, length = npoints)
    seriestype := :path
    seriescolor --> :red
    aspect_ratio --> :equal
    label --> ""

    # outer ring
    @series begin
        x = @. e.x + e.a_out * cos(t) * cosd(e.theta) - e.b_out * sin(t) * sind(e.theta) + 0.5
        y = @. e.y + e.a_out * cos(t) * sind(e.theta) + e.b_out * sin(t) * cosd(e.theta) + 0.5
        x, y
    end

    # inner ring
    @series begin
        x = @. e.x + e.a_in * cos(t) * cosd(e.theta) - e.b_in * sin(t) * sind(e.theta) + 0.5
        y = @. e.y + e.a_in * cos(t) * sind(e.theta) + e.b_in * sin(t) * cosd(e.theta) + 0.5
        x, y
    end
end


@recipe function f(ap::RectangularAperture, npoints = 1000)
    seriestype := :path
    seriescolor --> :red
    aspect_ratio --> :equal
    label --> ""

    t = range(0, 2π, length = npoints)
    w2 = ap.w / 2
    h2 = ap.h / 2
    
    x = @. w2 * (abs(cos(t)) * cos(t) + abs(sin(t)) * sin(t))
    y = @. h2 * (abs(cos(t)) * cos(t) - abs(sin(t)) * sin(t))

    sinth, costh = sincos(deg2rad(ap.theta))
    u = @. ap.x + x * costh - y * sinth + 0.5
    v = @. ap.y + x * sinth + y * costh + 0.5

    u, v
end

@recipe function f(ap::RectangularAnnulus, npoints = 1000)
    seriestype := :path
    seriescolor --> :red
    aspect_ratio --> :equal
    label --> ""

    t = range(0, 2π, length = npoints)

    # outer box
    @series begin
        w2 = ap.w_out / 2
        h2 = ap.h_out / 2
        
        x = @. w2 * (abs(cos(t)) * cos(t) + abs(sin(t)) * sin(t))
        y = @. h2 * (abs(cos(t)) * cos(t) - abs(sin(t)) * sin(t))
    
        sinth, costh = sincos(deg2rad(ap.theta))
        u = @. ap.x + x * costh - y * sinth + 0.5
        v = @. ap.y + x * sinth + y * costh + 0.5

        u, v
    end

    # inner box
    @series begin
        w2 = ap.w_in / 2
        h2 = ap.h_in / 2
        
        x = @. w2 * (abs(cos(t)) * cos(t) + abs(sin(t)) * sin(t))
        y = @. h2 * (abs(cos(t)) * cos(t) - abs(sin(t)) * sin(t))
    
        sinth, costh = sincos(deg2rad(ap.theta))
        u = @. ap.x + x * costh - y * sinth + 0.5
        v = @. ap.y + x * sinth + y * costh + 0.5

        u, v
    end
end
