using StatsPlots
using CSV

jt_s = CSV.read(joinpath(@__DIR__, "julia_circle_apertures_nt=1.csv"))
jt_4 = CSV.read(joinpath(@__DIR__, "julia_circle_apertures_nt=4.csv"))
jt = CSV.read(joinpath(@__DIR__, "julia_circle_apertures_nt=8.csv"))
pt = CSV.read(joinpath(@__DIR__, "python_circle_apertures.csv"))

plot(markerstrokealpha=0, xscale=:log10, legend=:topleft)

@df pt plot!(:N, :time, c=2, label="photutils")
@df jt_s plot!(:N, :time, markerstrokealpha=0, c=1, label="Photometry.jl - 1 thread", yscale=:log10)
@df jt_4 plot!(:N, :time, c=3,  label="Photometry.jl - 4 threads")
@df jt plot!(:N, :time, c=3, ls=:dash, label="Photometry.jl - 8 threads")
title!("size(512, 512) - CircularAperture - exact")
ylabel!("time [s]")
xlabel!("number of apertures")

savefig(joinpath(@__DIR__, "circle_apertures_benchmark.png"))
