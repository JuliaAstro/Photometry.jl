using StatsBase: mad

###############################################################################
# Location Estimators

const LOC_EST = [Mean, Median, Mode, MMM, SourceExtractor, BiweightLocation]

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

    E <: Mode || @test estimator(data) ≈ 0.0 atol = 1e-2
end

@testset "Mean" begin
    data = rand(100, 100)
    @test Mean()(data) == mean(data)
end

@testset "Median" begin
    data = rand(100, 100)
    @test Median()(data) == median(data)
end

@testset "Mode" begin
    x = [1,2,3,4,5,6,5,4,3,4,34,3,43,43,3,3,3,3,1]
    @test Mode()(x) == 3
end

@testset "SourceExtractor" begin
    # test skewed distribution
    data = float(collect(1:100))
    data[71:end] .= 1e7

    @test SourceExtractor()(data) ≈ median(data)
end

@testset "BiweightLocation" begin
    b = BiweightLocation()
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

@testset "Std" begin
    s = StdRMS()
    data = randn(100)
    @test s(data) == std(data, corrected = false)
end

@testset "MADStd" begin
    s = MADStdRMS()
    data = randn(100)
    @test s(data) == mad(data, normalize = true)
end

@testset "BiweightScale" begin
    s = BiweightScaleRMS()
    @test s([1, 3, 5, 500, 2]) ≈ 1.70992562072
end
