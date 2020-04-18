using StatsPlots
using CSV

jt_s = CSV.read(joinpath(@__DIR__, "julia_circle_apertures_single.csv"))
jt = CSV.read(joinpath(@__DIR__, "julia_circle_apertures.csv"))
pt = CSV.read(joinpath(@__DIR__, "python_circle_apertures.csv"))

@df jt_s plot(:N, :time, markerstrokealpha=0, c=1, label="Photometry.jl - 1 thread", legend=:topleft)
@df jt plot!(:N, :time, c=1, ls=:dash, label="Photometry.jl - 8 threads")
@df pt plot!(:N, :time, c=2, label="photutils")
title!("(512,512) - CircularAperture - exact")
ylabel!("time [s]")
xlabel!("number of apertures")

savefig(joinpath(@__DIR__, "apertures.png"))
