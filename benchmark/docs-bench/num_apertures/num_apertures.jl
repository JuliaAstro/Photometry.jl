using Photometry
using Chairmarks
using CSV
using DataFramesMeta
using ProgressMeter
using Random
using PythonCall

rng = Random.seed!(11256)

data = randn(rng, 512, 512) .+ 10

rows = @showprogress map((1, 10, 50, 100, 200, 400, 500, 1000, 2000)) do N
    xs, ys = size(data) .* (rand(rng, N), rand(rng, N))
    apers = CircularAperture.(xs, ys, 10)
    bench = @b photometry($apers, $data)
    (nt=Threads.nthreads(), N=N, time=bench.time)
end

path = joinpath(@__DIR__, "julia_num_apertures.csv")
results = CSV.read(path, DataFrame)

rows_to_update = @. results.nt == Threads.nthreads()
results[rows_to_update, :] .= DataFrame(rows)

CSV.write(path, results)

@info :Updated path

rows = @showprogress map((1, 10, 50, 100, 200, 400, 500, 1000, 2000)) do N
    xs, ys = size(data) .* (rand(rng, N), rand(rng, N))
    apers = EllipticalAperture.(xs, ys, 10, 10, 20)
    bench = @b photometry($apers, $data)
    (nt=Threads.nthreads(), N=N, time=bench.time)
end

path = joinpath(@__DIR__, "julia_num_apertures-ellipse.csv")
results = CSV.read(path, DataFrame)

rows_to_update = @. results.nt == Threads.nthreads()
results[rows_to_update, :] .= DataFrame(rows)

CSV.write(path, results)

@info :Updated path
