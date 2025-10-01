using StatsPlots
using CSV
using DataFramesMeta

jt = CSV.read(joinpath(@__DIR__, "julia_aperture_size.csv"), DataFrame)
pt = CSV.read(joinpath(@__DIR__, "python_aperture_size.csv"), DataFrame)

jt_ell = CSV.read(joinpath(@__DIR__, "julia_aperture_size-ellipse.csv"), DataFrame)
pt_ell = CSV.read(joinpath(@__DIR__, "python_aperture_size-ellipse.csv"), DataFrame)

plot(xticks=0:50:200, yticks=logrange(1e-6, 1e-2, 5), markerstrokealpha=0, yscale=:log10, link=:y, legend=:bottomright, layout=2, size=(800, 300), dpi=100, bottom_margin=5StatsPlots.mm)

@df pt plot!(:r, :time, c=2, label="", sp=1)
jt_groups = groupby(jt, :nt)
@df jt_groups[(nt=1,)] plot!(:r, :time, markerstrokealpha=0, c=1, label="", sp=1)
@df jt_groups[(nt=4,)] plot!(:r, :time, c=3,  label="", sp=1)
@df jt_groups[(nt=8,)] plot!(:r, :time, c=3, ls=:dash, label="", sp=1)
title!("CircularAperture", sp=1)
ylabel!("time [s]", sp=1)
xlabel!("aperture radius [px]", sp=1)


@df pt_ell plot!(:r, :time, c=2, label="photutils", sp=2)
jt_ell_groups = groupby(jt_ell, :nt)
@df jt_ell_groups[(nt=1,)] plot!(:r, :time, markerstrokealpha=0, c=1, label="Photometry.jl - 1 thread", sp=2)
@df jt_ell_groups[(nt=4,)]  plot!(:r, :time, c=3,  label="Photometry.jl - 4 threads", sp=2)
@df jt_ell_groups[(nt=8,)]  plot!(:r, :time, c=3, ls=:dash, label="Photometry.jl - 8 threads", sp=2)
title!("EllipticalAperture", sp=2)
xlabel!("aperture radius [px]", sp=2)

savefig(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "aperture_size_benchmark.png"))
