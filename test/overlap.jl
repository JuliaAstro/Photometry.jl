@testset "circular overlap" for grid_size in [50, 500, 1000], circ_size in (0.2, 0.4, 0.8), method in [:exact, :center, (:subpixel, 2), (:subpixel, 5), (:subpixel, 10)]

    g = Photometry.circular_overlap(-1, -1, 1, 1, grid_size, grid_size, circ_size, method = method)
    @test maximum(g) ≈ 1.0
    @test minimum(g) ≈ 0.0
end

@testset "area triangle" begin
    @test Photometry.area_triangle(0, 0, 0, 0, 0, 0) ≈ 0

    @test Photometry.area_triangle(0, 0, 1, 0, 0, 1) ≈ 1 / 2

    @test Photometry.area_triangle(1, 1, 1, 4, 6, 1) ≈ 7.5

end

@testset "area arc" begin
    @test Photometry.area_arc(0, 0, 0, 0, 1) ≈ 0

    for r in 1:10
        @test Photometry.area_arc(0, r, r, 0, r) ≈ r^2 / 2 * (π / 2 - 1)
    end
end

@testset "overlap core" begin
    @test Photometry.circular_overlap_core(0, 0, 0, 0, 1) ≈ 0
    

    for r in 1:10
        @test Photometry.circular_overlap_single_exact(0, 0, r, r, r) ≈ r^2 * π / 4
        @test Photometry.circular_overlap_core(0, 0, r, r, r) ≈ r^2 * π / 4
    end
end
