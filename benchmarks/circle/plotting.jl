using StatsPlots
using CSV

jt_s = CSV.read(joinpath(@__DIR__, "julia_circle_apertures_nt=1.csv"))
jt = CSV.read(joinpath(@__DIR__, "julia_circle_apertures_nt=8.csv"))
pt = CSV.read(joinpath(@__DIR__, "python_circle_apertures.csv"))

plot(markerstrokealpha=0, xscale=:log10, legend=:topleft)#, layout=(1, 2), size=(1000, 400))

@df jt_s plot!(:N, :time, markerstrokealpha=0, c=1, label="Photometry.jl - 1 thread", yscale=:log10, subplot=1)
@df jt plot!(:N, :time, c=1, ls=:dash, label="Photometry.jl - 8 threads", subplot=1)
@df pt plot!(:N, :time, c=2, label="photutils", subplot=1)
title!("(512,512) - CircularAperture - exact", subplot=1)
ylabel!("time [s]", subplot=1)
xlabel!("number of apertures", subplot=1)


# plot!(jt_s.N, jt_s.time ./ pt.time, c=1, label="Photometry.jl - 1 thread", yscale=:log10, subplot=2)
# plot!(jt.N, jt.time ./ pt.time, c=1, ls=:dash, label="Photometry.jl - 8 threads", subplot=2)
# hline!([1], c=:black, label="photutils", subplot=2)

# ylabel!("relative time (Photometry.jl / photutils)", subplot=2)
# xlabel!("number of apertures", subplot=2)

savefig(joinpath(@__DIR__, "circle_apertures.png"))
