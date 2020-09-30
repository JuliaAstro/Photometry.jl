using Photometry.Aperture: circular_overlap_core,
                           circular_overlap_single_exact,
                           circular_overlap_single_subpixel,
                           area_arc,
                           area_triangle,
                           inside_ellipse,
                           inside_triangle,
                           inside_rectangle,
                           circle_line,
                           circle_segment,
                           circle_segment_single2,
                           triangle_unitcircle_overlap,
                           elliptical_overlap_exact,
                           elliptical_overlap_single_subpixel,
                           rectangular_overlap_exact,
                           rectangular_overlap_single_subpixel

@testset "overlap - circular" begin

    @testset "circular overlap" for r in 0:10:100, use_subpixel in [true, false]

        ap = CircularAperture(0, 0, r)
        if use_subpixel
            ap = Subpixel(ap, 10)
        end
        @test all(p -> 0 ≤ p ≤ 1, collect(ap))
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

    @testset "type stability" begin
        r = 5
        @inferred circular_overlap_core(0, 0, r, r, r)
        @inferred circular_overlap_single_exact(0, 0, r, r, r)
        @inferred circular_overlap_single_subpixel(0, 0, 20, 20, 15, 5)
    end

end # circles

@testset "overlap - elliptical" begin

    @testset "elliptical overlap" for r in 0:10:100, use_subpixel in [true, false]
        ap = EllipticalAperture(0, 0, r, r)
        if use_subpixel
            ap = Subpixel(ap, 10)
        end
        @test all(p -> 0 ≤ p ≤ 1, collect(ap))
    end

    @testset "overlap subpixel (elliptical apperture)" begin
        @test elliptical_overlap_single_subpixel(0, 0, 0, 0, 1, 1, 0, 1) ≈ 1

        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 2) ≈ 0.25
        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 5) ≈ 0.44
        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 10) ≈ 0.43
        @test elliptical_overlap_single_subpixel(0, 0, 20, 20, 1 / 225, 1 / 225, 0, 100) ≈ 0.4423
    end

    @testset "circle line" begin
        for line in readlines(joinpath(@__DIR__, "..", "data", "circle_line.csv"))
            tokens = split(line, ',')
            numbers = parse.(Float64, tokens)
            x1, y1, x2, y2 = numbers[1:4]
            point1, point2 = circle_line(x1, y1, x2, y2)
            @test point1[1] ≈ numbers[5] atol = 1e-6
            @test point1[2] ≈ numbers[6] atol = 1e-6
            @test point2[1] ≈ numbers[7] atol = 1e-6
            @test point2[2] ≈ numbers[8] atol = 1e-6
        end
    end

    @testset "circle segment" begin
        for line in readlines(joinpath(@__DIR__, "..", "data", "circle_segment.csv"))
            tokens = split(line, ',')
            numbers = parse.(Float64, tokens)
            x1, y1, x2, y2 = numbers[1:4]
            point1, point2 = circle_segment(x1, y1, x2, y2)
            @test point1[1] ≈ numbers[5] atol = 1e-6
            @test point1[2] ≈ numbers[6] atol = 1e-6
            @test point2[1] ≈ numbers[7] atol = 1e-6
            @test point2[2] ≈ numbers[8] atol = 1e-6
        end
    end

    @testset "circle segment single" begin
        for line in readlines(joinpath(@__DIR__, "..", "data", "circle_segment_single.csv"))
            tokens = split(line, ',')
            numbers = parse.(Float64, tokens)
            x1, y1, x2, y2 = numbers[1:4]
            point1 = circle_segment_single2(x1, y1, x2, y2)
            @test point1[1] ≈ numbers[5] atol = 1e-6
            @test point1[2] ≈ numbers[6] atol = 1e-6
        end
    end

    @testset "triangle unitcircle overlap" begin
        for line in readlines(joinpath(@__DIR__, "..", "data", "triangle_unitcircle_overlap.csv"))
            tokens = split(line, ',')
            numbers = parse.(Float64, tokens)
            x1, y1, x2, y2, x3, y3 = numbers[1:6]
            @test triangle_unitcircle_overlap(x1, y1, x2, y2, x3, y3) ≈ numbers[end] atol = 1e-6
        end
    end
    @testset "elliptical overlap" begin
        @test elliptical_overlap_exact(0.5, 2.5, 1.5, 3.5, 3.0, 3.0, 0) ≈ 0.311725 atol = 1e-6
        @test elliptical_overlap_exact(0, 2, 1, 3, 3.0, 3.0, 0) ≈ 0.943480 atol = 1e-6
    end
end # overlap elliptical

@testset "overlap - rectangular" begin
    @testset "exact" begin
        @test rectangular_overlap_exact(0, 0, 1, 1, 1, 1, 0) ≈ 0.25
    end

    @testset "type stability" begin
        # @inferred rectangular_overlap_exact(0, 0, 1, 1, 2, 2, 0) # TODO area in LazySets.jl not type stable
        @inferred rectangular_overlap_single_subpixel(0, 0, 1, 1, 2, 2, 0, 5)
    end
end # overlap rectangular

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

    @testset "inside ellipse" begin
        @test inside_ellipse(0, 0, 0, 0, 1, 1, 1)
        @test !inside_ellipse(10, 10, 5, 5, 1, 1, 1)
        @test !inside_ellipse(5, 3, 0, 0, 1 / 16, 1 / 32, 0)
        @test inside_ellipse(0, 0, 0, 0, 5, 6.2, 0)
        @test inside_ellipse(1, 2, 0, 0, 1 / 37, 1 / 36, -1 / 80)
    end

    @testset "inside rectangle" begin
        @test inside_rectangle(0, 0, 3, 4, 0)
        @test !inside_rectangle(10, 10, 3, 4, 0)
    end

    @testset "type stability" begin
        @inferred area_arc(0, 0, 0, 0, 10)
        @inferred area_triangle(0, 0, 1, 0, 0, 1)
        @inferred inside_triangle(0, 0, 0, 0, 0, 1, 1, 0)
        @inferred inside_ellipse(0, 0, 5, 5, 1, 1, 1)
        @inferred inside_rectangle(0, 0, 3, 4, 0)
    end

end # overlap utils
