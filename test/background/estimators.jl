using Photometry.Background:
    MMMBackground,
    SourceExtractorBackground,
    BiweightLocationBackground,
    StdRMS,
    MADStdRMS,
    BiweightScaleRMS

@testset "background/estimators: Trivial estimator" begin
    @testset "Trivial $E"  for E in [MMMBackground, SourceExtractorBackground, BiweightLocationBackground]
        estimator = E()

        @test estimator(ones(1)) == 1

        data = ones(10, 10)

        @test estimator(data) ≈ 1.0
        @test estimator(data, dims = 1) ≈ ones(1, 10)
        @test estimator(data, dims = 2) ≈ ones(10)

        data = zeros(10, 10)

        @test estimator(data) ≈ 0.0
        @test estimator(data, dims = 1) ≈ zeros(1, 10)
        @test estimator(data, dims = 2) ≈ zeros(10)

        data = randn(100, 100)
    end
end

@testset "background/estimators: SourceExtractorBackground" begin
    # test skewed distribution
    data = float(collect(1:100))
    data[71:end] .= 1e7

    @test SourceExtractorBackground()(data) ≈ median(data)
end

@testset "background/estimators: BiweightLocationBackground" begin
    b = BiweightLocationBackground()
    @test b([1, 3, 5, 500, 2]) ≈ 2.745 atol = 1e-3
end

###############################################################################
# RMS Estimators

@testset "background/estimators: Trivial RMS estimator" begin
    @testset "Trivial $E"  for E in [StdRMS, MADStdRMS, BiweightScaleRMS]
        estimator = E()

        @test estimator(ones(1)) == 0

        data = ones(10, 10)
        @test estimator(data) ≈ 0.0
        @test estimator(data, dims = 1) ≈ zeros(1, 10)
        @test estimator(data, dims = 2) ≈ zeros(10)


        data = randn(10000, 10000)
        @test estimator(data) ≈ 1 atol = 3e-2
    end
end

@testset "background/estimators: StdRMS" begin
    s = StdRMS()
    data = randn(100)
    @test s(data) == std(data, corrected = false)
end

@testset "background/estimators: MADStdRMS" begin
    s = MADStdRMS()
    data = randn(100)
    @test s(data) == mad(data, normalize = true)
end

@testset "background/estimators: BiweightScaleRMS" begin
    s = BiweightScaleRMS()
    @test s([1, 3, 5, 500, 2]) ≈ 1.70992562072
end
