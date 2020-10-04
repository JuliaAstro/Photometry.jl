using StatsPlots
using CSV

jt_s = CSV.read(joinpath(@__DIR__, "julia_aperture_size_nt=1.csv"))
jt_4 = CSV.read(joinpath(@__DIR__, "julia_aperture_size_nt=4.csv"))
jt = CSV.read(joinpath(@__DIR__, "julia_aperture_size_nt=8.csv"))
pt = CSV.read(joinpath(@__DIR__, "python_aperture_size.csv"))

plot(markerstrokealpha=0, legend=:bottomright)

@df pt plot!(:r, :time, c=2, label="photutils")
@df jt_s plot!(:r, :time, markerstrokealpha=0, c=1, label="Photometry.jl - 1 thread", yscale=:log10)
@df jt_4 plot!(:r, :time, c=3,  label="Photometry.jl - 4 threads")
@df jt plot!(:r, :time, c=3, ls=:dash, label="Photometry.jl - 8 threads")
title!("size(512, 512) - CircularAperture - exact")
ylabel!("time [s]")
xlabel!("aperture radius [px]")

savefig(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "aperture_size_benchmark.png"))

nothing
