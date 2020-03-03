using Photometry
using Test

@testset "Aperture Photometry" begin
    include("aperture/overlap.jl")
    include("aperture/circular.jl")
    include("aperture/photometry.jl")
    include("aperture/elliptical.jl")
    include("aperture/plots.jl")
end

@testset "Background Estimation" begin
    # mesh fitting not supported yet
    @test_throws ErrorException estimate_background(Mean, ones(1, 1), (1, 1), (1, 1))
    include("background/simple.jl")
end
