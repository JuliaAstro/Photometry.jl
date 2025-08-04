using Photometry
using Random
using Statistics
using TestItemRunner

Random.seed!(8462852)

@run_package_tests

#@testset "Aperture Photometry" begin
#    include("aperture/overlap.jl")
#    include("aperture/circular.jl")
#    include("aperture/photometry.jl")
#    include("aperture/elliptical.jl")
#    include("aperture/rectangular.jl")
#    include("aperture/plots.jl")
#end

#@testset "Background Estimation" begin
#    include("background/background.jl")
#    include("background/estimators.jl")
#    include("background/interpolators.jl")
#end
#
#@testset "Source Detection" begin
#    include("detection/detection.jl")
#end
