using Photometry
using Test

@testset "Photometry.jl" begin
    include("overlap.jl")
    include("circular.jl")
end
