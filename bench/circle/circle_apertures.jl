using Photometry
using BenchmarkTools
using CSV
using DataFrames
using ProgressLogging

data = randn(512, 512) .+ 10

rows = []

@progress for N in [1, 10, 50, 100, 200, 400, 500, 1000, 2000]
    apers = CircularAperture.(256, 256, 10 .* rand(N))
    time = @belapsed aperture_photometry($apers, $data)
    push!(rows, (N=N, time=time))
end

results = DataFrame(rows)
outfile = "julia_circle_apertures_nt=$(Threads.nthreads()).csv"
@info "saving to $outfile"
CSV.write(joinpath(@__DIR__, outfile), results)

nothing
