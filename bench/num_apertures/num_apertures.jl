using Photometry
using BenchmarkTools
using CSV
using DataFrames
using ProgressLogging
using Random

rng = Random.seed!(11256)

data = randn(rng, 512, 512) .+ 10

rows = []
using Statistics
@progress for N in [1, 10, 50, 100, 200, 400, 500, 1000, 2000]
    xs, ys = size(data) .* (rand(rng, N), rand(rng, N))
    apers = CircularAperture.(xs, ys, 10)
    time = @belapsed photometry($apers, $data)
    push!(rows, (N=N, time=time))
end

results = DataFrame(rows)
outfile = "julia_num_apertures_nt=$(Threads.nthreads()).csv"
@info "saving to $outfile"
CSV.write(joinpath(@__DIR__, outfile), results)

nothing
