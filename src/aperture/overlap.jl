#=
Part of this work is derived from astropy/photutils and kbarbary/sep. The relevant derivations
are considered under a BSD 3-clause license. =#
using LazySets

function circular_overlap(xmin, xmax, ymin, ymax, nx, ny, r; method = :exact)
    R = float(typeof(xmin))
    out = fill(zero(R), ny, nx)

    # width of each element
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # radius of one pixel
    pixel_radius = 0.5sqrt(dx^2 + dy^2)

    # bounding box
    bxmin = -r - 0.5dx
    bxmax =  r + 0.5dx
    bymin = -r - 0.5dy
    bymax =  r + 0.5dy

    @inbounds for i in axes(out, 2)
        # lower end of pixel
        pxmin = xmin + (i - 1) * dx
        pxcen = pxmin + 0.5dx
        # upper end of pixel
        pxmax = pxmin + dx

        if pxmax > bxmin && pxmin < bxmax
            for j in axes(out, 1)
                pymin = ymin + (j - 1) * dy
                pycen = pymin + 0.5dy
                pymax = pymin + dy
                if pymax > bymin && pymin < bymax

                    # distance from circle to pixel
                    d = sqrt(pxcen^2 + pycen^2)

                    # fully within radius
                    if d < r - pixel_radius
                        out[j, i] = 1
                    # partially within radius
                    elseif d < r + pixel_radius
                        if method === :exact
                            out[j, i] = circular_overlap_single_exact(pxmin, pymin, pxmax, pymax, r) / (dx * dy)
                        elseif method === :center
                            out[j, i] =  circular_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, r, 1)
                        elseif method[1] === :subpixel
                            out[j, i] =  circular_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, r, method[2])
                        end
                    end
                end
            end
        end

    end
    return out
end

function circular_overlap_single_subpixel(xmin, ymin, xmax, ymax, r, subpixels)
    frac = 0
    dx = (xmax - xmin) / subpixels
    dy = (ymax - ymin) / subpixels
    r2 = r^2

    x = xmin - 0.5dx
    for i in 1:subpixels
        x += dx
        y = ymin - 0.5dy
        for j in 1:subpixels
            y += dy
            if x^2 + y^2 < r2
                frac += 1
            end
        end
    end
    return frac / subpixels^2
end

"""Area of overlap between a rectangle and a circle"""
function circular_overlap_single_exact(xmin, ymin, xmax, ymax, r)
    R = float(typeof(xmin))
    r ≤ 0 && return zero(R)
    if 0 ≤ xmin
        0 ≤ ymin && return circular_overlap_core(xmin, ymin, xmax, ymax, r)
        0 ≥ ymax && return circular_overlap_core(-ymax, xmin, -ymin, xmax, r)
        return (circular_overlap_single_exact(xmin, ymin, xmax, 0, r) +
                circular_overlap_single_exact(xmin, 0, xmax, ymax, r))
    elseif 0 ≥ xmax
        0 ≤ ymin && return circular_overlap_core(-xmax, ymin, -xmin, ymax, r)
        0 ≥ ymax && return circular_overlap_core(-xmax, -ymax, -xmin, -ymin, r)
        return (circular_overlap_single_exact(xmin, ymin, xmax, 0, r) +
                circular_overlap_single_exact(xmin, 0, xmax, ymax, r))
    else
        0 ≤ ymin && return (circular_overlap_single_exact(xmin, ymin, 0, ymax, r) +
                            circular_overlap_single_exact(0, ymin, xmax, ymax, r))
        0 ≥ ymax && return (circular_overlap_single_exact(xmin, ymin, 0, ymax, r) +
                            circular_overlap_single_exact(0, ymin, xmax, ymax, r))
        return (circular_overlap_single_exact(xmin, ymin, 0, 0, r) +
                circular_overlap_single_exact(0, ymin, xmax, 0, r) +
                circular_overlap_single_exact(xmin, 0, 0, ymax, r) +
                circular_overlap_single_exact(0, 0, xmax, ymax, r))
    end
end

"""Core of circular overlap routine"""
function circular_overlap_core(xmin, ymin, xmax, ymax, r)
    R = float(typeof(xmin))
    xmin^2 + ymin^2 > r^2 && return zero(R)
    xmax^2 + ymax^2 < r^2 && return R((xmax - xmin) * (ymax - ymin))

    d1 = sqrt(xmax^2 + ymin^2)
    d2 = sqrt(xmin^2 + ymax^2)
    if d1 < r && d2 < r
        x1, y1 = sqrt(r^2 - ymax^2), ymax
        x2, y2 = xmax, sqrt(r^2 - xmax^2)
        area = ((xmax - xmin) * (ymax - ymin) -
                area_triangle(x1, y1, x2, y2, xmax, ymax) +
                area_arc(x1, y1, x2, y2, r))
    elseif d1 < r
        x1, y1 = xmin, sqrt(r^2 - xmin^2)
        x2, y2 = xmax, sqrt(r^2 - xmax^2)
        area = (area_arc(x1, y1, x2, y2, r) +
                area_triangle(x1, y1, x1, ymin, xmax, ymin) +
                area_triangle(x1, y1, x2, ymin, x2, y2))
    elseif d2 < r
        x1, y1 = sqrt(r^2 - ymin^2), ymin
        x2, y2 = sqrt(r^2 - ymax^2), ymax
        area = (area_arc(x1, y1, x2, y2, r) +
                area_triangle(x1, y1, xmin, y1, xmin, ymax) +
                area_triangle(x1, y1, xmin, y2, x2, y2))
    else
        x1, y1 = sqrt(r^2 - ymin^2), ymin
        x2, y2 = xmin, sqrt(r^2 - xmin^2)
        area = (area_arc(x1, y1, x2, y2, r) +
                area_triangle(x1, y1, x2, y2, xmin, ymin))
    end

    return R(area)
end

####################################
# Some geometric helpers

"""Area of a triangle defined by three vertices"""
area_triangle(x0, y0, x1, y1, x2, y2) = 0.5abs(x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1))


"""
Area of a circular segment above a chord between two points with circle radius `r`

[Reference](http://mathworld.wolfram.com/CircularSegment.html)
"""
function area_arc(x0, y0, x1, y1, r)
    a = sqrt((x1 - x0)^2 + (y1 - y0)^2)
    θ = 2asin(0.5a / r)
    return r^2 * (θ - sin(θ)) / 2
end

# is (x, y) contained within triangle specified by three points.
function inside_triangle(x, y, x1, y1, x2, y2, x3, y3)
    c = (y1 > y) ⊻ (y2 > y) && x < (x2 - x1) * (y - y1) / (y2 - y1) + x1
    c += (y2 > y) ⊻ (y3 > y) &&  x < (x3 - x2) * (y - y2) / (y3 - y2) + x2
    c += (y3 > y) ⊻ (y1 > y) && x < (x1 - x3) * (y - y3) / (y1 - y3) + x3

    # true if c is odd
    return c % 2 == 1
end

####################################
# Elliptical routines

"""
    inside_ellipse(x, y, h, k, cxx, cyy, cxy)

x: x coordinate of the test point
y: y coordinate of the test point
h: x coordinate of the center of ellipse
k: y coordinate of the center of ellipse
cxx, cyy, cxy: coefficients of equation of ellipse

Utility function to find whether a point is inside ellipse or not.

If point inside ellipse: Returns true else returns false

General equation of ellipse:
    cxx * (x - h)^2 + cxy * (x - h) * (y - k) + cyy * (y - k)^2 = 1
"""
inside_ellipse(x, y, h, k, cxx, cyy, cxy) = cxx * (x - h)^2 + cxy * (x - h) * (y - k) + cyy * (y - k)^2  - 1 < 0

function elliptical_overlap(xmin, xmax, ymin, ymax, nx, ny, a, b, theta; method = :exact)
    R = float(typeof(xmin))
    out = fill(zero(R), ny, nx)

    # width of each element
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # radius of one pixel
    pixel_radius = 0.5sqrt(dx^2 + dy^2)

    # bounding box
    bxmin, bxmax, bymin, bymax = bbox(EllipticalAperture(0, 0, a, b, theta))

    bxmin -= 0.5dx
    bxmax += 0.5dx
    bymin -= 0.5dy
    bymax += 0.5dy

    cxx, cyy, cxy = oblique_coefficients(a, b, theta)

    @inbounds for i in axes(out, 2)
        # lower end of pixel
        pxmin = xmin + (i - 1) * dx
        pxcen = pxmin + 0.5dx
        # upper end of pixel
        pxmax = pxmin + dx

        if pxmax > bxmin && pxmin < bxmax
            for j in axes(out, 1)
                pymin = ymin + (j - 1) * dy
                pycen = pymin + 0.5dy
                pymax = pymin + dy
                if pymax > bymin && pymin < bymax

                    # 4 flags for four different ends of the pixel
                    # each flag tells wether the point is inside the ellipse or not
                    # sample strcture of a pixel is shown below
                    #  2---4
                    #  |   |
                    #  |   |
                    #  1---3
                    flags = [
                        inside_ellipse(pxmin, pymin, 0, 0, cxx, cyy, cxy),
                        inside_ellipse(pxmin, pymax, 0, 0, cxx, cyy, cxy),
                        inside_ellipse(pxmax, pymin, 0, 0, cxx, cyy, cxy),
                        inside_ellipse(pxmax, pymax, 0, 0, cxx, cyy, cxy)
                    ]

                    # fully within radius
                    if all(flags)
                        out[j, i] = 1
                    # partially within radius
                    elseif any(flags)
                        if method === :exact
                            out[j, i] = elliptical_overlap_exact(pxmin, pymin, pxmax, pymax, a, b, theta)
                        elseif method === :center
                            out[j, i] =  elliptical_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, cxx, cyy, cxy, 1)
                        elseif method[1] === :subpixel
                            out[j, i] =  elliptical_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, cxx, cyy, cxy, method[2])
                        end
                    end
                end
            end
        end
    end
    return out
end

elliptical_overlap(xmin, xmax, ymin, ymax, nx, ny, parameters::AbstractVector; method) = elliptical_overlap(xmin, xmax, ymin, ymax, nx, ny, parameters... ; method = :center)


function elliptical_overlap_single_subpixel(xmin, ymin, xmax, ymax, cxx, cyy, cxy, subpixels)
    frac = 0
    dx = (xmax - xmin) / subpixels
    dy = (ymax - ymin) / subpixels

    x = xmin - 0.5dx
    for i in 1:subpixels
        x += dx
        y = ymin - 0.5dy
        for j in 1:subpixels
            y += dy
            if inside_ellipse(x, y, 0, 0, cxx, cyy, cxy)
                frac += 1
            end
        end
    end
    return frac / subpixels^2
end

# overlap between triangle and unit circle
function triangle_unitcircle_overlap(x1, y1, x2, y2, x3, y3)
    # distances
    d1 = x1^2 + y1^2
    d2 = x2^2 + y2^2
    d3 = x3^2 + y3^2

    # order by distances
    ds = [d1, d2, d3]
    order = sortperm(ds)
    ds = ds[order]
    x1, x2, x3 = [x1, x2, x3][order]
    y1, y2, y3 = [y1, y2, y3][order]

    # which are inside circle
    inside = ds .< 1

    # on circle
    on = isapprox.(ds, 1, atol = 1e-10)

    # triangle is completely inside circle
    if inside[3] || on[3]
        return area_triangle(x1, y1, x2, y2, x3, y3)

    # if vertex 1 or 2 are on the edge, then dot product with 3 to determine intersection
    elseif inside[2] || on[2]
        intersect13 = !on[1] || x1 * (x3 - x1) + y1 * (y3 - y1) < 0
        intersect23 = !on[2] || x2 * (x3 - x2) + y2 * (y3 - y2) < 0
        if intersect13 && intersect23
            point1 = circle_segment_single2(x1, y1, x3, y3)
            point2 = circle_segment_single2(x2, y2, x3, y3)

            return (area_triangle(x1, y1, x2, y2, point1...) +
                    area_triangle(x2, y2, point1..., point2...) +
                    area_arc(point1..., point2..., 1))
        elseif intersect13
            point1 = circle_segment_single2(x1, y1, x3, y3)

            return (area_triangle(x1, y1, x2, y2, point1...) +
                    area_arc(x2, y2, point1..., 1))
        elseif intersect23
            point2 = circle_segment_single2(x2, y2, x3, y3)

            return (area_triangle(x1, y1, x2, y2, point2...) +
                    area_arc(x2, y2, point2..., 1))
        else
            return area_arc(x1, y1, x2, y2, 1)
        end
    elseif inside[1]
        point1, point2 = circle_segment(x2, y2, x3, y3)
        point3 = circle_segment_single2(x1, y1, x2, y2)
        point4 = circle_segment_single2(x1, y1, x3, y3)

        if point1[1] > 1 # no intersection
            # check if (x1, y2) and origin are on opposite sides of segment
            if (((0 - point3[2]) * (point4[1] - point3[1]) > (point4[2] - point3[2]) * (0 - point3[1])) !=
                ((y1 - point3[2]) * (point4[1] - point3[1]) > (point4[2] - point3[2]) * (x1 - point3[1])))
                return (area_triangle(x1, y1, point3..., point4...) + π - area_arc(point3..., point4..., 1))
            else
                return area_triangle(x1, y1, point3..., point4...) + area_arc(point3..., point4..., 1)
            end
        else

              # ensure that point1 is the point closest to (x2, y2)
              if (((point2[1] - x2) * (point2[1] - x2) + (point2[2] - y2) * (point2[2] - y2)) <
                  ((point1[1] - x2) * (point1[1] - x2) + (point1[2] - y2) * (point1[2] - y2)))
                point1, point2 = point2, point1
              end

              return (area_triangle(x1, y1, point3..., point1...) +
                      area_triangle(x1, y1, point1..., point2...) +
                      area_triangle(x1, y1, point2..., point4...) +
                      area_arc(point1..., point3..., 1) +
                      area_arc(point2..., point4..., 1))

        end
    else
        point1, point2 = circle_segment(x1, y1, x2, y2)
        point3, point4 = circle_segment(x2, y2, x3, y3)
        point5, point6 = circle_segment(x3, y3, x1, y1)

        if point1[1] ≤ 1
            xp = (point1[1] + point2[1]) / 2
            yp = (point1[2] + point2[2]) / 2
            return (triangle_unitcircle_overlap(x1, y1, x3, y3, xp, yp) +
                    triangle_unitcircle_overlap(x2, y2, x3, y3, xp, yp))
        elseif point3[1] ≤ 1
            xp = (point3[1] + point4[1]) / 2
            yp = (point3[2] + point4[2]) / 2
            return (triangle_unitcircle_overlap(x3, y3, x1, y1, xp, yp) +
                    triangle_unitcircle_overlap(x2, y2, x1, y1, xp, yp))
        elseif point5[1] ≤ 1
            xp = (point5[1] + point6[1]) / 2
            yp = (point5[2] + point6[2]) / 2
            return (triangle_unitcircle_overlap(x1, y1, x3, y3, xp, yp) +
                    triangle_unitcircle_overlap(x3, y3, x2, y2, xp, yp))
        else
            return inside_triangle(0, 0, x1, y1, x2, y2, x3, y3) ? π : 0
        end
    end
end

# intersection of a segment with the unit cirlce
function circle_segment(x1, y1, x2, y2)
    point1, point2 = circle_line(x1, y1, x2, y2)

    if ((point1[1] > x1 && point1[1] > x2) || (point1[1] < x1 && point1[1] < x2) ||
        (point1[2] > y1 && point1[2] > y2) || (point1[2] < y1 && point1[2] < y2))
        point1 = (2, 2)
    end

    if ((point2[1] > x1 && point2[1] > x2) || (point2[1] < x1 && point2[1] < x2) ||
        (point2[2] > y1 && point2[2] > y2) || (point2[2] < y1 && point2[2] < y2))
        point2 = (2, 2)
    end

    return point1[1] > 1 && point2[1] < 2 ?  (point1, point2) : (point2, point1)
end

# closest intersection of a line with the unit cirlce
function circle_segment_single2(x1, y1, x2, y2)
    point1, point2 = circle_line(x1, y1, x2, y2)
    dx1 = abs(point1[1] - x2)
    dy1 = abs(point1[2] - y2)
    dx2 = abs(point2[1] - x2)
    dy2 = abs(point2[2] - y2)

    if dx1 > dy1
        return dx1 > dx2 ? point2 : point1
    else
        return dy1 > dy2 ? point2 : point1
    end
end

# intersection of a line defined by two points with a unit cirlce
function circle_line(x1, y1, x2, y2)
    dx = x2 - x1
    dy = y2 - y1

    if dx ≈ 0 && dy ≈ 0
        return (2, 2), (2, 2)
    elseif abs(dx) > abs(dy)
        # find slope and intercept
        a = dy / dx
        b = y1 - a * x1

        # find quadratic determinant
        δ2 = 1 + a^2 - b^2
        if δ2 > 0 # real solutions
            δ = sqrt(δ2)
            _x1 = (-a * b - δ) / (1 + a^2)
            _y1 = a * _x1 + b
            _x2 = (-a * b + δ) / (1 + a^2)
            _y2 = a * _x2 + b
            return (_x1, _y1), (_x2, _y2)
        else
            return (2, 2), (2, 2)
        end
    else
        # find slope and intercept
        a = dx / dy
        b = x1 - a * y1

        # find quadratic determinant
        δ2 = 1 + a^2 - b^2
        if δ2 > 0 # real solutions
            δ = sqrt(δ2)
            _y1 = (-a * b - δ) / (1 + a^2)
            _x1 = a * _y1 + b
            _y2 = (-a * b + δ) / (1 + a^2)
            _x2 = a * _y2 + b
            return (_x1, _y1), (_x2, _y2)
        else
            return (2, 2), (2, 2)
        end
    end
end

function elliptical_overlap_exact(xmin, ymin, xmax, ymax, a, b, θ)
    sint, cost = sincos(deg2rad(-θ))

    scale = a * b

    # reproject ellipse
    x1 = (xmin * cost - ymin * sint) / a
    y1 = (xmin * sint + ymin * cost) / b
    x2 = (xmax * cost - ymin * sint) / a
    y2 = (xmax * sint + ymin * cost) / b
    x3 = (xmax * cost - ymax * sint) / a
    y3 = (xmax * sint + ymax * cost) / b
    x4 = (xmin * cost - ymax * sint) / a
    y4 = (xmin * sint + ymax * cost) / b

    return scale * (triangle_unitcircle_overlap(x1, y1, x2, y2, x3, y3) +
                    triangle_unitcircle_overlap(x1, y1, x4, y4, x3, y3))
end

####################################
# Rectangular routines

function rectangular_overlap(xmin, xmax, ymin, ymax, nx, ny, w, h, θ; method = :exact)
    R = float(typeof(xmin))
    out = fill(zero(R), ny, nx)

    # width of each element
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # bounding box
    bxmin, bxmax, bymin, bymax = bbox(RectangularAperture(0, 0, w, h, θ))

    bxmin -= 0.5dx
    bxmax += 0.5dx
    bymin -= 0.5dy
    bymax += 0.5dy

    @inbounds for i in axes(out, 2)
        # lower end of pixel
        pxmin = xmin + (i - 1) * dx
        # upper end of pixel
        pxmax = pxmin + dx
        if pxmax > bxmin && pxmin < bxmax
            for j in axes(out, 1)
                pymin = ymin + (j - 1) * dy
                pymax = pymin + dy

                if bymin < pymax && pymin < bymax
                    if method === :exact
                        out[j, i] = rectangular_overlap_exact(pxmin, pymin, pxmax, pymax, w, h, θ)
                    elseif method === :center
                        out[j, i] = rectangular_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, w, h, θ, 1)
                    elseif method[1] === :subpixel
                        out[j, i] = rectangular_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, w, h, θ, method[2])
                    end
                end
            end
        end
    end

    return out
end

function rectangular_overlap_single_subpixel(x0, y0, x1, y1, w, h, θ, subpixels)
    dx = (x1 - x0) / subpixels
    dy = (y1 - y0) / subpixels

    frac = 0
    x = x0 - 0.5dx
    for i in 1:subpixels
        x += dx
        y = y0 - 0.5dy
        for j in 1:subpixels
            y += dy
            if inside_rectangle(x, y, w, h, θ)
                frac += 1
            end
        end
    end

    return frac / subpixels^2

end

# see https://math.stackexchange.com/questions/69099/equation-of-a-rectangle
"""intersection with rectangular using implicit Lamé curve"""
function inside_rectangle(x, y, w, h, θ)
    # transform into frame of rectangle
    sinth, costh = sincos(deg2rad(-θ))
    u = x * costh - y * sinth
    v = x * sinth + y * costh

    return abs(u) < w / 2 && abs(v) < h / 2
end

function rectangular_overlap_exact(xmin, ymin, xmax, ymax, w, h, θ)
    sint, cost = sincos(deg2rad(θ))
    R = [cost -sint; sint cost]
    aper = R * Hyperrectangle(zeros(2), [w / 2, h / 2])
    dy = ymax - ymin
    dx = xmax - xmin
    pix = Hyperrectangle([dx / 2 + xmin, dy / 2 + ymin], [dx / 2, dy / 2])
    return intersection_area(pix, aper)
end

# contrib: https://nbviewer.jupyter.org/github/mforets/escritoire/blob/master/2020/Week11/2D_intersection.ipynb
# the idea is to use the intersection for HPolygon, which is faster than the fallback because the
# normal vectors to the half-spaces are sorted in counter--clockwise fashion
function intersection_area(X::AbstractHyperrectangle{N},
    Y::LinearMap{N,<:AbstractHyperrectangle{N}}) where {N}
    X_clist = X |> constraints_list
    X_poly = HPolygon(X_clist, sort_constraints = true, prune = false, check_boundedness = false)

    Y_clist = linear_map(matrix(Y), set(Y)) |> constraints_list
    Y_poly = HPolygon(Y_clist, sort_constraints = true, prune = false, check_boundedness = false)

    return intersection(X_poly, Y_poly) |> area
end
