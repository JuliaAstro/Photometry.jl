using StatsPlots
using CSV
using DataFrames

jt = CSV.read(joinpath(@__DIR__, "julia_num_apertures.csv")) |> DataFrame
pt = CSV.read(joinpath(@__DIR__, "python_num_apertures.csv")) |> DataFrame

jt_ell = CSV.read(joinpath(@__DIR__, "julia_num_apertures-ellipse.csv")) |> DataFrame
pt_ell = CSV.read(joinpath(@__DIR__, "python_num_apertures-ellipse.csv")) |> DataFrame

plot(markerstrokealpha=0, xscale=:log10, yscale=:log10, link=:y, legend=:topleft, layout=2, size=(800, 300), dpi=100, bottom_margin=5StatsPlots.mm)

@df pt plot!(:N, :time, c=2, label="photutils", sp=1)
jt_group = groupby(jt, :nt)
@df jt_group[(nt=1,)] plot!(:N, :time, markerstrokealpha=0, c=1, label="Photometry.jl - 1 thread", sp=1)
@df jt_group[(nt=4,)] plot!(:N, :time, c=3,  label="Photometry.jl - 4 threads", sp=1)
@df jt_group[(nt=8,)] plot!(:N, :time, c=3, ls=:dash, label="Photometry.jl - 8 threads", sp=1)
title!("CircularAperture(r=10)", sp=1)
ylabel!("time [s]", sp=1)
xlabel!("number of apertures", sp=1)

@df pt_ell plot!(:N, :time, c=2, label="", sp=2)
jt_ell_group = groupby(jt, :nt)
@df jt_ell_group[(nt=1,)] plot!(:N, :time, markerstrokealpha=0, c=1, label="", sp=2)
@df jt_ell_group[(nt=4,)] plot!(:N, :time, c=3,  label="", sp=2)
@df jt_ell_group[(nt=8,)] plot!(:N, :time, c=3, ls=:dash, label="", sp=2)
title!("EllipticalAperture(a=10, b=10)", sp=2)
xlabel!("number of apertures", sp=2)

savefig(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "num_apertures_benchmark.png"))

nothing
