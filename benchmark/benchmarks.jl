using BenchmarkTools
using StableRNGs
using Photometry

const SUITE = BenchmarkGroup()
const IMG_SIZE = 512
const RNG = StableRNG(1234)
const DATA = randn(RNG, IMG_SIZE, IMG_SIZE) .+ 10
const ERR = fill(1.0, size(DATA))

SUITE["circular_aperture"] = BenchmarkGroup()

cx = cy = (IMG_SIZE + 1) / 2
for r in range(1, 100; step = 11)
    circ_ap = CircularAperture(cx, cy, r)
    SUITE["circular_aperture"][lpad(string(r), 3)] = @benchmarkable photometry($circ_ap, $DATA)
    SUITE["circular_aperture"][lpad(string(r, " + error"), 11)] = @benchmarkable photometry($circ_ap, $DATA, $ERR)
end

# TODO: Consider merging this with `bench/`
## Aperuter sizes
#let
#    SUITE["aperture_size"] = BenchmarkGroup()
#    SUITE["aperture_size"]["CircularAperture"] = BenchmarkGroup()
#    SUITE["aperture_size"]["EllipticalAperture"] = BenchmarkGroup()
#
#    cx = cy = (IMG_SIZE + 1) / 2
#    for r in 1:5:200
#        circ_ap = CircularAperture(cx, cy, r)
#        ell_ap = EllipticalAperture(cx, cy, r, r, 20)
#        SUITE["aperture_size"]["CircularAperture"][r] = @benchmarkable photometry($circ_ap, $DATA)
#        SUITE["aperture_size"]["EllipticalAperture"][r] = @benchmarkable photometry($ell_ap, $DATA)
#    end
#end
#
## Number of apertures
#let
#    SUITE["num_apertures"] = BenchmarkGroup()
#    SUITE["num_apertures"]["CircularAperture"] = BenchmarkGroup()
#    SUITE["num_apertures"]["EllipticalAperture"] = BenchmarkGroup()
#
#    size_x, size_y = size(DATA)
#    for N in (1, 10, 50, 100, 200, 400, 500, 1000, 2000)
#        xs = size_x * rand(RNG, N)
#        ys = size_y * rand(RNG, N)
#        circ_apers = CircularAperture.(xs, ys, 10)
#        ell_apers = EllipticalAperture.(xs, ys, 10, 10, 20)
#        SUITE["num_apertures"]["CircularAperture"][N] = @benchmarkable photometry($circ_apers, $DATA)
#        SUITE["num_apertures"]["EllipticalAperture"][N] = @benchmarkable photometry($ell_apers, $DATA)
#    end
#end
#
## Source detection
#let
#    SUITE["detect"] = BenchmarkGroup()
#    SUITE["detect"]["PeakMesh"] = BenchmarkGroup()
#
#    nsigma = 3.0
#    for box_size in 3:2:51
#        alg = PeakMesh(; box_size, nsigma)
#        SUITE["detect"]["PeakMesh"][box_size] = @benchmarkable extract_sources($alg, $DATA)
#    end
#end
