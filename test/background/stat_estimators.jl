using Photometry.Background: sigma_clip

@testset "sigma clipping" begin
    x = [1.0, 2.0, 3.0]
    @test sigma_clip(x, 1, 1) ≈ [1.0, 2.0, 3.0] rtol=1e-4
    @test sigma_clip(x, 1) ≈ [1.0, 2.0, 3.0] rtol=1e-4

    y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    @test sigma_clip(y,1,1,iterations=20) ≈ [2.77002, 2.77002, 3.0, 4.0, 4.22998, 4.22998] rtol=1e-4
    @test sigma_clip(y,1,1) ≈ [2.30234, 2.30234, 3.0, 4.0, 4.69766, 4.69766] rtol=1e-4
end
