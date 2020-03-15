using StatsBase: mad

###############################################################################
# Location Estimators

const LOC_EST = [MeanBackground, MedianBackground, ModeBackground, MMMBackground, SourceExtractorBackground, BiweightLocationBackground]

@testset "Trivial $E"  for E in LOC_EST
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

    E <: ModeBackground || @test estimator(data) ≈ 0.0 atol = 1e-2
end

@testset "MeanBackground" begin
    data = rand(100, 100)
    @test MeanBackground()(data) == mean(data)
end

@testset "MedianBackground" begin
    data = rand(100, 100)
    @test MedianBackground()(data) == median(data)
end

@testset "ModeBackground" begin
    x = [1,2,3,4,5,6,5,4,3,4,34,3,43,43,3,3,3,3,1]
    @test ModeBackground()(x) == 3
end

@testset "SourceExtractorBackground" begin
    # test skewed distribution
    data = float(collect(1:100))
    data[71:end] .= 1e7

    @test SourceExtractorBackground()(data) ≈ median(data)
end

@testset "BiweightLocationBackground" begin
    b = BiweightLocationBackground()
    @test b([1, 3, 5, 500, 2]) ≈ 2.745 atol = 1e-3
end

###############################################################################
# RMS Estimators


@testset "Trivial $E"  for E in [StdRMS, MADStdRMS, BiweightScaleRMS]
    estimator = E()

    @test estimator(ones(1)) == 0

    data = ones(10, 10)
    @test estimator(data) ≈ 0.0
    @test estimator(data, dims = 1) ≈ zeros(1, 10)
    @test estimator(data, dims = 2) ≈ zeros(10)


    data = randn(100, 100)
    @test estimator(data) ≈ 1 atol = 2e-2
end

@testset "StdRMS" begin
    s = StdRMS()
    data = randn(100)
    @test s(data) == std(data, corrected = false)
end

@testset "MADStdRMS" begin
    s = MADStdRMS()
    data = randn(100)
    @test s(data) == mad(data, normalize = true)
end

@testset "BiweightScaleRMS" begin
    s = BiweightScaleRMS()
    @test s([1, 3, 5, 500, 2]) ≈ 1.70992562072
end
