using AperturePhotometry
using Test

@testset "AperturePhotometry.jl" begin
    include("overlap.jl")
    include("circular.jl")
end
