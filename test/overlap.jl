@testset "circular overlap" for grid_size in [50, 500, 1000], circ_size in (0.2, 0.4, 0.8), method in [:exact, 1, 5, 10]

    g = AperturePhotometry.circular_overlap(-1, -1, 1, 1, grid_size, grid_size, circ_size, method = method)
    @test maximum(g) ≈ 1.0
    @test minimum(g) ≈ 0.0
end
