using Photometry
using Test
using Random

Random.seed!(8462852)

@testset "Aperture Photometry" begin
    include("aperture/overlap.jl")
    include("aperture/circular.jl")
    include("aperture/photometry.jl")
    include("aperture/elliptical.jl")
    include("aperture/rectangular.jl")
    include("aperture/plots.jl")
end

@testset "Background Estimation" begin
    include("background/estimators.jl")
    # include("background/")
end
