using Photometry
using Test

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
