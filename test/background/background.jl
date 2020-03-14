
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

@testset "functional stuff" begin
    data = ones(100, 100)

    @test_throws ErrorException estimate_background(data, (4, 4), edge_method = :yeet)
    @test all(estimate_background(data, 20) .== estimate_background(data, (20, 20)))
end

@testset "flat background - $B, $S" for B in [Mean, Median, MMM, BiweightLocation, SourceExtractor], S in [StdRMS, MADStdRMS, BiweightScaleRMS]
    data = ones(100, 100)

    bk, rms = estimate_background(data, B(), S())
    @test bk ≈ 1
    @test rms ≈ 0

    bk, rms = estimate_background(data, (25, 25), B(), S())
    @test median(bk) ≈ 1
    @test median(rms) ≈ 0
end