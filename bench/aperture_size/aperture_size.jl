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
    ap = CircularAperture(256.5, 256.5, r)
    time = @belapsed photometry($ap, $data)
    push!(rows, (r=r, time=time))
end

results = DataFrame(rows)
outfile = "julia_aperture_size_nt=$(Threads.nthreads()).csv"
@info "saving to $outfile"
CSV.write(joinpath(@__DIR__, outfile), results)

nothing
