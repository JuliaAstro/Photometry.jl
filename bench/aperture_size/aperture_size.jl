using Photometry
using Chairmarks
using CSV
using DataFramesMeta
using ProgressMeter
using Random

rng = Random.seed!(11256)

data = randn(rng, 512, 512) .+ 10

rows = @showprogress map(1:5:200) do r
    ap = CircularAperture(256.5, 256.5, r)
    bench = @b photometry($ap, $data)
    (nt=Threads.nthreads(), r=r, time=bench.time)
end

path = joinpath(@__DIR__, "julia_aperture_size.csv")
results = CSV.read(path, DataFrame)

rows_to_update = @. results.nt == Threads.nthreads()
results[rows_to_update, :] .= DataFrame(rows)

CSV.write(path, results)

@info :Updated path

rows = @showprogress map(1:5:200) do r
    ap = EllipticalAperture(256.5, 256.5, r, r, 20)
    bench = @b photometry($ap, $data)
    (nt=Threads.nthreads(), r=r, time=bench.time)
end

path = joinpath(@__DIR__, "julia_aperture_size-ellipse.csv")
results = CSV.read(path, DataFrame)

rows_to_update = @. results.nt == Threads.nthreads()
results[rows_to_update, :] .= DataFrame(rows)

CSV.write(path, results)

@info :Updated path
