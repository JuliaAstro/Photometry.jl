using Photometry.Aperture: circular_overlap,
                           circular_overlap_core,
                           circular_overlap_single_exact,
                           circular_overlap_single_subpixel,
                           area_arc,
                           area_triangle,
                           inside_ellipse,
                           inside_triangle,
                           circle_line,
                           circle_segment,
                           circle_segment_single2,
                           triangle_unitcircle_overlap,
                           elliptical_overlap,
                           elliptical_overlap_exact,
                           elliptical_overlap_single_subpixel,
                           rectangular_overlap,
                           rectangular_overlap_exact,
                           rectangular_overlap_single_subpixel

@testset "overlap - circular" begin

    @testset "circular overlap" for grid_size in [50, 500, 1000], circ_size in (0.2, 0.4, 0.8), method in [:exact, :center, (:subpixel, 2), (:subpixel, 5), (:subpixel, 10)]

        g = circular_overlap(-1, -1, 1, 1, grid_size, grid_size, circ_size, method = method)
        any(g .> 0) &&  @test maximum(g) ≈ 1.0
        any(g .< 1) &&  @test minimum(g) ≈ 0.0
    end

    @testset "overlap core (r = $r)" for r in 1:10
        @test circular_overlap_core(0, 0, 0, 0, r) ≈ 0
        @test circular_overlap_single_exact(0, 0, r, r, r) ≈ r^2 * π / 4
        @test circular_overlap_core(0, 0, r, r, r) ≈ r^2 * π / 4
    end

    @testset "overlap single exact" begin
        @test circular_overlap_single_exact(0, 0, 0, 0, 1) ≈ 0
        @test circular_overlap_single_exact(0, 0, 20, 20, 15) ≈ 176.714586764429
    end

    @testset "overlap subpixel" begin
        @test circular_overlap_single_subpixel(0, 0, 0, 0, 1, 1) ≈ 1

        @test circular_overlap_single_subpixel(0, 0, 20, 20, 15, 2) ≈ 0.25
        @test circular_overlap_single_subpixel(0, 0, 20, 20, 15, 5) ≈ 0.44
        @test circular_overlap_single_subpixel(0, 0, 20, 20, 15, 10) ≈ 0.43
        @test circular_overlap_single_subpixel(0, 0, 20, 20, 15, 100) ≈ 0.4423
    end

    @testset "overlap grid" begin
        @test circular_overlap(0, 0, 0, 0, 1, 1, 1) == ones(1, 1)

        @test circular_overlap(0, 20, 0, 20, 5, 5, 15) ≈ [1 1 1 0.70506901 0;
    1 1 1 0.42240835 0;
    1 1 0.74542733 0.02213982 0;
    0.70506901 0.42240835 0.02213982 0 0;
    0 0 0 0 0]

        @test circular_overlap(-10, 10, -10, 10, 5, 5, 4) ≈ [0 0         0        0        0;
    0 0.07878669 0.45661148 0.07878669 0;
    0 0.45661148 1         0.45661148 0;
    0 0.07878669 0.45661148 0.07878669 0;
    0 0         0         0        0]


    end

end # circles

@testset "overlap - elliptical" begin

    @testset "position with respect to ellipse" begin
        @test !inside_ellipse(5, 3, 0, 0, 1 / 16, 1 / 32, 0)
        @test inside_ellipse(0, 0, 0, 0, 5, 6.2, 0)
        @test inside_ellipse(1, 2, 0, 0, 1 / 37, 1 / 36, -1 / 80)
    end

    @testset "elliptical overlap" for grid_size in [50, 500, 1000], ellipse_size in ([0.2,0.2,0], [0.4, 0.4, 0], [0.8, 0.8, 0]), method in [:center, (:subpixel, 2), (:subpixel, 5), (:subpixel, 10)]

        g = elliptical_overlap(-1, -1, 1, 1, grid_size, grid_size, ellipse_size, method = method)
        any(g .> 0) &&  @test maximum(g) ≈ 1.0
        any(g .< 1) &&  @test minimum(g) ≈ 0.0
    end

    @testset "overlap subpixel (elliptical apperture)" begin
        @test elliptical_overlap_single_subpixel(0, 0, 0, 0, 1, 1, 0, 1) ≈ 1

        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 2) ≈ 0.25
        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 5) ≈ 0.44
        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 10) ≈ 0.43
        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 100) ≈ 0.4423
    end

    @testset "circle line" for line in readlines(joinpath(@__DIR__, "..", "data", "circle_line.csv"))
        tokens = split(line, ',')
        numbers = parse.(Float64, tokens)
        x1, y1, x2, y2 = numbers[1:4]
        point1, point2 = circle_line(x1, y1, x2, y2)
        @test point1[1] ≈ numbers[5] atol = 1e-6
        @test point1[2] ≈ numbers[6] atol = 1e-6
        @test point2[1] ≈ numbers[7] atol = 1e-6
        @test point2[2] ≈ numbers[8] atol = 1e-6
    end

    @testset "circle segment" for line in readlines(joinpath(@__DIR__, "..", "data", "circle_segment.csv"))
        tokens = split(line, ',')
        numbers = parse.(Float64, tokens)
        x1, y1, x2, y2 = numbers[1:4]
        point1, point2 = circle_segment(x1, y1, x2, y2)
        @test point1[1] ≈ numbers[5] atol = 1e-6
        @test point1[2] ≈ numbers[6] atol = 1e-6
        @test point2[1] ≈ numbers[7] atol = 1e-6
        @test point2[2] ≈ numbers[8] atol = 1e-6
    end

    @testset "circle segment single" for line in readlines(joinpath(@__DIR__, "..", "data", "circle_segment_single.csv"))
        tokens = split(line, ',')
        numbers = parse.(Float64, tokens)
        x1, y1, x2, y2 = numbers[1:4]
        point1 = circle_segment_single2(x1, y1, x2, y2)
        @test point1[1] ≈ numbers[5] atol = 1e-6
        @test point1[2] ≈ numbers[6] atol = 1e-6
    end

    @testset "triangle unitcircle overlap" for line in readlines(joinpath(@__DIR__, "..", "data", "triangle_unitcircle_overlap.csv"))
        tokens = split(line, ',')
        numbers = parse.(Float64, tokens)
        x1, y1, x2, y2, x3, y3 = numbers[1:6]
        @test triangle_unitcircle_overlap(x1, y1, x2, y2, x3, y3) ≈ numbers[end] atol = 1e-6
    end

    @testset "elliptical overlap" begin
        @test elliptical_overlap_exact(0.5, 2.5, 1.5, 3.5, 3.0, 3.0, 0) ≈ 0.311725 atol = 1e-6
        @test elliptical_overlap_exact(0, 2, 1, 3, 3.0, 3.0, 0) ≈ 0.943480 atol = 1e-6
    end

end # overlap elliptical 

@testset "overlap - rectangular" begin

end

@testset "overlap - utils" begin

    @testset "inside triangle" begin
        t1 = [0, 0, 0, 1, 1, 0]
        @test inside_triangle(0, 0, t1...)
        @test inside_triangle(0.1, 0.1, t1...)
        @test inside_triangle(0.25, 0.25, t1...)
    end


    @testset "area triangle" begin
        @test area_triangle(0, 0, 0, 0, 0, 0) ≈ 0
        @test area_triangle(0, 0, 1, 0, 0, 1) ≈ 0.5
        @test area_triangle(1, 1, 1, 4, 6, 1) ≈ 7.5

    end

    @testset "area arc (r = $r)" for r in 1:10
        @test area_arc(0, 0, 0, 0, r) ≈ 0
        @test area_arc(0, r, r, 0, r) ≈ r^2 / 2 * (π / 2 - 1)
    end

end # overlap utils