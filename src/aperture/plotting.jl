using RecipesBase

function circle(x, y, r, θ)
    u = x + r * cos(θ) + 0.5
    v = y + r * sin(θ) + 0.5
    return u, v
end

circle(ap::CircularAperture, θ) = circle(ap.x, ap.y, ap.r, θ)

function ellipse(x, y, a, b, θ, t=0)
    u = x + a * cos(θ) * cosd(t) - b * sin(θ) * sind(t) + 0.5
    v = y + a * cos(θ) * sind(t) + b * sin(θ) * cosd(t) + 0.5
    return u, v
end

ellipse(ap::EllipticalAperture, θ) = ellipse(ap.x, ap.y, ap.a, ap.b, θ, ap.theta)

@recipe function f(c::CircularAperture, npoints = 100)
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    t = range(0, 2pi, length = npoints)
    circle.((c,), t)
end

@recipe function f(c::CircularAnnulus, npoints = 100)
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    t = range(0, 2pi, length = npoints)
    # outer ring
    @series begin
        circle.(c.x, c.y, c.r_out, t)
    end

    # inner ring
    @series begin
        circle.(c.x, c.y, c.r_in, t)
    end
end


@recipe function f(e::EllipticalAperture, npoints = 100)
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    t = range(0, 2pi, length = npoints)
    ellipse.((e,), t)
end

@recipe function f(e::EllipticalAnnulus, npoints = 100)
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    t = range(0, 2pi, length = npoints)

    # outer ring
    @series begin
        ellipse.(e.x, e.y, e.a_out, e.b_out, t, e.theta)
    end

    # inner ring
    @series begin
        ellipse.(e.x, e.y, e.a_in, e.b_in, t, e.theta)
    end
end

function rectangle(x, y, w, h, θ, t=0)
    half_w, half_h = w / 2, h / 2
    sinth, costh = sincos(θ)
    u = half_w * (abs(costh) * costh + abs(sinth) * sinth)
    v = half_h * (abs(costh) * costh - abs(sinth) * sinth)
    sint, cost = sincos(deg2rad(t))
    j = x + u * cost - v * sint + 0.5
    k = y + u * sint + v * cost + 0.5
    return j, k
end

rectangle(ap::RectangularAperture, θ) = rectangle(ap.x, ap.y, ap.w, ap.h, θ, ap.theta)


@recipe function f(ap::RectangularAperture, npoints = 100)
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    t = range(0, 2π, length = npoints)
    return rectangle.((ap,), t)
end

@recipe function f(ap::RectangularAnnulus, npoints = 100)
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    t = range(0, 2π, length = npoints)

    # outer box
    @series begin
        rectangle.(ap.x, ap.y, ap.w_out, ap.h_out, t, ap.theta)
    end

    # inner box
    @series begin
        rectangle.(ap.x, ap.y, ap.w_in, ap.h_in, t, ap.theta)
    end
end

@recipe function f(aps::AbstractVector{<:AbstractAperture})
    seriestype := :shape
    aspect_ratio --> :equal
    seriescolor --> :red
    linecolor --> :match
    fillcolor --> nothing
    label --> ""

    for ap in aps
        @series begin ap end
    end
end
