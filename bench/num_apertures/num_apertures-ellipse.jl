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
    apers = EllipticalAperture.(xs, ys, 10, 10, 20)
    time = @belapsed photometry($apers, $data)
    push!(rows, (nt=Threads.nthreads(), N=N, time=time))
end

path = joinpath(@__DIR__, "julia_num_apertures-ellipse.csv")
results = CSV.File(path)

rows_to_update = @. results[:nt] == Threads.nthreads()
results[rows_to_update] = DataFrame(rows)

CSV.write(path, results)

nothing
