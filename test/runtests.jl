using Photometry
using Random
using Statistics
using ParallelTestRunner

const init_code = quote
    import StatsBase: median, mean, std, mad

    const DATA_DIR = joinpath(@__DIR__, "data")
end


runtests(Photometry, ["--verbose"]; init_code)
