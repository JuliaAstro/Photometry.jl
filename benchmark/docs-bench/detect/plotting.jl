using StatsPlots
using CSV
using DataFramesMeta

jt = CSV.read(joinpath(@__DIR__, "julia_extract_sources.csv"), DataFrame)
pt = CSV.read(joinpath(@__DIR__, "python_extract_sources.csv"), DataFrame)

plot(
    xticks = 3:8:51,
    yticks = logrange(1e-6, 1e-1, 6),
    ylims  = (exp10(-6.5), exp10(-0.5)),
    markerstrokealpha = 0,
    yscale  = :log10,
    legend  = :bottomright,
    size    = (500, 350),
    dpi     = 100,
    bottom_margin = 5StatsPlots.mm,
)

@df pt plot!(:box_size, :median_s, c=2, label="photutils (find_peaks)")

jt_groups = groupby(jt, :nt)
@df jt_groups[(nt=1,)] plot!(:box_size, :median_s, c=1, label="Photometry.jl - 1 thread")
@df jt_groups[(nt=4,)] plot!(:box_size, :median_s, c=3, label="Photometry.jl - 4 threads")
@df jt_groups[(nt=8,)] plot!(:box_size, :median_s, c=3, ls=:dash, label="Photometry.jl - 8 threads")

title!("PeakMesh source detection")
ylabel!("time [s]")
xlabel!("box size [px]")

savefig(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "extract_sources_benchmark.png"))
