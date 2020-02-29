using Photometry.Background: sigma_clip

@testset "sigma clipping" begin
    x = [1, 2, 3]
    @test sigma_clip(x, 1, 1) ≈ [1.0, 2.0, 3.0] rtol=1e-4
    @test sigma_clip(x, 1) ≈ [1.0, 2.0, 3.0] rtol=1e-4

    y = [1, 2, 3, 4, 5, 6]
    @test sigma_clip(y,1,1) ≈ [1.62917, 2.0, 3.0, 4.0, 5.0, 5.37083] rtol=1e-4
end
