using Photometry
using BenchmarkTools
using CSV
using DataFrames
using ProgressLogging

data = randn(512, 512) .+ 10

rows = []

@progress for N in [1, 10, 50, 100, 500, 1000]
    aper = CircularAperture(256, 256, 3)
    apers = repeat([aper], N)
    time = @belapsed aperture_photometry($apers, $data)
    push!(rows, (N=N, time=time))
end

results = DataFrame(rows)
CSV.write(joinpath(@__DIR__, "julia_circle_apertures.csv"), results)

nothing
