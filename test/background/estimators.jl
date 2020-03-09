function test_ones(estimator)
    data = ones(10, 10)

    @test estimate_background(estimator, data) ≈ 1.0
    @test estimate_background(estimator, data, dims = 1) ≈ ones(1, 10)
    @test estimate_background(estimator, data, dims = 2) ≈ ones(10)
end

function test_zeros(estimator)
    data = zeros(10, 10)

    @test estimate_background(estimator, data) ≈ 0.0
    @test estimate_background(estimator, data, dims = 1) ≈ zeros(1, 10)
    @test estimate_background(estimator, data, dims = 2) ≈ zeros(10)
end

@testset "$E"  for E in [Mean, Median, Mode, MMMBackground, SourceExtractor]
    @test estimate_background(E, ones(10, 10)) == estimate_background(E(), ones(10, 10))
    test_ones(E)
    test_zeros(E)
end

@testset "mode" begin
    x = [1,2,3,4,5,6,5,4,3,4,34,3,43,43,3,3,3,3,1]
    @test estimate_background(Mode, x) == 3
end

@testset "sigma clipping" begin
    x = [1, 2, 3]
    @test sigma_clip(x, 1, 1) ≈ [1.0, 2.0, 3.0] rtol = 1e-4
    @test sigma_clip(x, 1) == sigma_clip(x, 1, 1)

    y = [1, 2, 3, 4, 5, 6]
    @test sigma_clip(y, 1, 1) ≈ [1.62917, 2.0, 3.0, 4.0, 5.0, 5.37083] rtol = 1e-4

    # using different center
    @test sigma_clip(y, 1, center = 4) ≈ [2.1291713066130296, 2.1291713066130296, 3, 4, 5, 5.87082869338697]
    # using different std
    @test sigma_clip(y, 1, std = 1) ≈ [2.5, 2.5, 3, 4, 4.5, 4.5]
end
