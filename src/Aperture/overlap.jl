
export Exact,
       Center,
       Subpixel

abstract type OverlapMethod end

function circular_overlap(xmin, xmax, ymin, ymax, nx, ny, r; method = :exact)
    out = fill(0.0, nx, ny)

    # width of each element
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # radius of one pixel
    pixel_radius = 0.5sqrt(dx^2 + dy^2)

    # bounding box
    bxmin = -r - 0.5dx
    bxmax = r + 0.5dx
    bymin = -r - 0.5dy
    bymax = r + 0.5dy

    for i in 1:nx
        # lower end of pixel
        pxmin = xmin + (i - 1) * dx
        pxcen = pxmin + 0.5dx
        # upper end of pixel
        pxmax = pxmin + dx

        if pxmax > bxmin && pxmin < bxmax
            for j in 1:ny
                pymin = ymin + (j - 1) * dy
                pycen = pymin + 0.5dy
                pymax = pymin + dy
                if pymax > bymin && pymin < bymax

                    # distance from circle to pixel
                    d = sqrt(pxcen^2 + pycen^2)

                    # fully within radius
                    if d < r - pixel_radius
                        @inbounds out[j, i] = 1
                    # partially within radius
                    elseif d < r + pixel_radius
                        if method === :exact
                            @inbounds out[j, i] = circular_overlap_single_exact(pxmin, pymin, pxmax, pymax, r) / (dx * dy)
                        elseif method === :center
                            @inbounds out[j, i] =  circular_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, r, 1)
                        elseif method[1] === :subpixel
                            @inbounds out[j, i] =  circular_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, r, method[2])
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
    r ≤ 0 && return 0
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
    xmin^2 + ymin^2 > r^2 && return 0
    xmax^2 + ymax^2 < r^2 && return (xmax - xmin) * (ymax - ymin)

    d1 = sqrt(xmax^2 + ymin^2)
    d2 = sqrt(xmin^2 + ymax^2)
    if d1 < r && d2 < r
        x1, y1 = sqrt(r^2 - ymax^2), ymax
        x2, y2 = xmax, sqrt(r^2 - xmax^2)
        return ((xmax - xmin) * (ymax - ymin) -
                area_triangle(x1, y1, x2, y2, xmax, ymax) +
                area_arc(x1, y1, x2, y2, r))
    elseif d1 < r
        x1, y1 = xmin, sqrt(r^2 - xmin^2)
        x2, y2 = xmax, sqrt(r^2 - xmax^2)
        return (area_arc(x1, y1, x2, y2, r) +
                area_triangle(x1, y1, x1, ymin, xmax, ymin) +
                area_triangle(x1, y1, x2, ymin, x2, y2))
    elseif d2 < r
        x1, y1 = sqrt(r^2 - ymin^2), ymin
        x2, y2 = sqrt(r^2 - ymax^2), ymax
        return (area_arc(x1, y1, x2, y2, r) +
                area_triangle(x1, y1, xmin, y1, xmin, ymax) +
                area_triangle(x1, y1, xmin, y2, x2, y2))
    else
        x1, y1 = sqrt(r^2 - ymin^2), ymin
        x2, y2 = xmin, sqrt(r^2 - xmin^2)
        return (area_arc(x1, y1, x2, y2, r) +
                area_triangle(x1, y1, x2, y2, xmin, ymin))
    end
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
    return 0.5r^2 * (θ - sin(θ))
end

"""
point_completely_inside_ellipse(x, y, h, k, cxx, cyy, cxy)

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
point_completely_inside_ellipse(x, y, h, k, cxx, cyy, cxy) = cxx * (x - h)^2 + cxy *(x - h)*(y - k) + cyy * (y - k)^2  - 1 <= 0

function elliptical_overlap(xmin, xmax, ymin, ymax, nx, ny, a, b, theta; method = :center)
    out = fill(0.0, nx, ny)

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

    for i in 1:nx
        # lower end of pixel
        pxmin = xmin + (i - 1) * dx
        pxcen = pxmin + 0.5dx
        # upper end of pixel
        pxmax = pxmin + dx

        if pxmax > bxmin && pxmin < bxmax
            for j in 1:ny
                pymin = ymin + (j - 1) * dy
                pycen = pymin + 0.5dy
                pymax = pymin + dy
                if pymax > bymin && pymin < bymax

                    #4 flags for four different ends of the pixel
                    #each flag tells wether the point is inside the ellipse or not
                    #sample strcture of a pixel is shown below
                    #  2---4
                    #  |   |
                    #  |   |
                    #  1---3

                    flag1 = point_completely_inside_ellipse(pxmin, pymin, 0, 0, cxx, cyy, cxy)
                    flag2 = point_completely_inside_ellipse(pxmin, pymax, 0, 0, cxx, cyy, cxy)
                    flag3 = point_completely_inside_ellipse(pxmax, pymin, 0, 0, cxx, cyy, cxy)
                    flag4 = point_completely_inside_ellipse(pxmax, pymax, 0, 0, cxx, cyy, cxy)

                    # fully within radius
                    if flag1 && flag2 && flag3 && flag4
                        @inbounds out[j, i] = 1
                    # partially within radius
                    elseif flag1 || flag2 || flag3 || flag4
                        if method === :exact
                            error("The exact method for elliptical overlap has not yet been implemented.")
                        elseif method === :center
                            @inbounds out[j, i] =  elliptical_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, cxx, cyy, cxy, 1)
                        elseif method[1] === :subpixel
                            @inbounds out[j, i] =  elliptical_overlap_single_subpixel(pxmin, pymin, pxmax, pymax, cxx, cyy, cxy, method[2])
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
            if point_completely_inside_ellipse(x, y, 0, 0, cxx, cyy, cxy)
                frac += 1
            end
        end
    end
    return frac / subpixels^2
end
