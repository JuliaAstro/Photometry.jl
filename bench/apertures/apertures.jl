using Photometry
using BenchmarkTools
using CSV
using DataFrames
using ProgressLogging

AP = [
    CircularAperture,
    CircularAnnulus,
    EllipticalAperture,
    EllipticalAnnulus,
    RectangularAperture,
    RectangularAnnulus,
]
PARAMS = [
    (3),
    (3, 5),
    (3, 3, 0),
    (3, 5, 4, 0),
    (3, 5, 0),
    (3, 5, 4, 0)
]

data = randn(512, 512) .+ 10

rows = []

@progress for N in [1, 10, 50, 100], method in [:exact, :center, (:subpixel, 5)], (A, par) in zip(AP, PARAMS)
    aper = A(256, 256, par...)
    apers = repeat([aper], N)
    func(a, d) = aperture_photometry(a, d, method=method)
    time = @belapsed $func($apers, $data)
    m = method isa Tuple ? "$(String(method[1]))-$(method[2])" : String(method)
    push!(rows, (N=N, aperture=A, method=m, time=time))
end

results = DataFrame(rows)
CSV.write(joinpath(@__DIR__, "julia_apertures.csv"), results)
