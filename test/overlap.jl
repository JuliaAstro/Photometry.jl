@testset "circular overlap" for grid_size in [50, 500, 1000], circ_size in (0.2, 0.4, 0.8), method in [:exact, :center, (:subpixel, 2), (:subpixel, 5), (:subpixel, 10)]

    g = Photometry.circular_overlap(-1, -1, 1, 1, grid_size, grid_size, circ_size, method = method)
    any(g .> 0) &&  @test maximum(g) ≈ 1.0
    any(g .< 1) &&  @test minimum(g) ≈ 0.0
end

@testset "area triangle" begin
    @test Photometry.area_triangle(0, 0, 0, 0, 0, 0) ≈ 0

    @test Photometry.area_triangle(0, 0, 1, 0, 0, 1) ≈ 1 / 2

    @test Photometry.area_triangle(1, 1, 1, 4, 6, 1) ≈ 7.5

end

@testset "area arc (r = $r)" for r in 1:10
    @test Photometry.area_arc(0, 0, 0, 0, r) ≈ 0
    @test Photometry.area_arc(0, r, r, 0, r) ≈ r^2 / 2 * (π / 2 - 1)
end

@testset "overlap core (r = $r)" for r in 1:10
    @test Photometry.circular_overlap_core(0, 0, 0, 0, r) ≈ 0
    @test Photometry.circular_overlap_single_exact(0, 0, r, r, r) ≈ r^2 * π / 4
    @test Photometry.circular_overlap_core(0, 0, r, r, r) ≈ r^2 * π / 4
end

@testset "overlap single exact" begin
    @test Photometry.circular_overlap_single_exact(0, 0, 0, 0, 1) ≈ 0
    @test Photometry.circular_overlap_single_exact(0, 0, 20, 20, 15) ≈ 176.714586764429
end

@testset "overlap subpixel" begin
    @test Photometry.circular_overlap_single_subpixel(0, 0, 0, 0, 1, 1) ≈ 1

    @test Photometry.circular_overlap_single_subpixel(0, 0, 20, 20, 15, 2) ≈ 0.25
    @test Photometry.circular_overlap_single_subpixel(0, 0, 20, 20, 15, 5) ≈ 0.44
    @test Photometry.circular_overlap_single_subpixel(0, 0, 20, 20, 15, 10) ≈ 0.43
    @test Photometry.circular_overlap_single_subpixel(0, 0, 20, 20, 15, 100) ≈ 0.4423
end

@testset "overlap grid" begin
    @test Photometry.circular_overlap(0, 0, 0, 0, 1, 1, 1) == ones(1, 1)

    @test Photometry.circular_overlap(0, 20, 0, 20, 5, 5, 15) ≈ [1 1 1 0.70506901 0;
    1 1 1 0.42240835 0;
    1 1 0.74542733 0.02213982 0;
    0.70506901 0.42240835 0.02213982 0 0;
    0 0 0 0 0]

    @test Photometry.circular_overlap(-10, 10, -10, 10, 5, 5, 4) ≈ [0 0         0        0        0;
    0 0.07878669 0.45661148 0.07878669 0;
    0 0.45661148 1         0.45661148 0;
    0 0.07878669 0.45661148 0.07878669 0;
    0 0         0         0        0]


end
