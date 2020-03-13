
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