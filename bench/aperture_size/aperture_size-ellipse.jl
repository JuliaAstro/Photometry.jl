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
@progress for r in 1:5:200
    ap = EllipticalAperture(256.5, 256.5, r, r, 20)
    time = @belapsed photometry($ap, $data)
    push!(rows, (nt=Threads.nthreads(), r=r, time=time))
end

path = joinpath(@__DIR__, "julia_aperture_size-ellipse.csv")
results = CSV.File(path)

rows_to_update = @. results[:nt] == Threads.nthreads()
results[rows_to_update] = DataFrame(rows)

CSV.write(path, results)
