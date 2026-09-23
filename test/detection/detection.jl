using Photometry.Detection:
    PeakMesh,
    extract_sources
using Photometry.Aperture:
    CircularAperture,
    photometry
import Photometry

import Random
Random.seed!(8462852)

# arrays with non-1-based axes, without adding a test dependency
const OffsetArrays = Photometry.Detection.ImageFiltering.OffsetArrays

@testset "detection/Detection: peak finding" begin
    @testset "peak finding - $P" for P in [PeakMesh()]
        data = randn(100, 100)
        idxs = [1, 700, 1524]
        fake_peaks = randn(length(idxs)) .+ 10
        data[idxs] .= fake_peaks

        table = extract_sources(P, data)
        @test table.value[1:3] == sort(fake_peaks, rev = true)
    end
end

@testset "detection/Detection: Peak Mesh" begin
    @test PeakMesh(box_size = 3) == PeakMesh(box_size = (3, 3))
    # box sizes are converted to Int, so any integer type can be passed on
    @test PeakMesh(box_size = Int32(3)) == PeakMesh(box_size = (3, 3))
    @test PeakMesh(box_size = (Int32(3), 5)).box_size === (3, 5)
    @test_throws ArgumentError PeakMesh(box_size = 0)
    @test_throws ArgumentError PeakMesh(box_size = (3, -1))
end

@testset "detection/Detection: position convention" begin
    # x indexes the first array axis and y the second, matching the coordinate
    # convention of the aperture types
    data = zeros(20, 30)
    data[10, 3] = 3.0
    data[4, 27] = 5.0

    table = extract_sources(PeakMesh(), data)
    @test table.x == [4, 10]
    @test table.y == [27, 3]
    @test table.value == [5.0, 3.0]

    # detected positions feed directly into apertures: the flux is exactly
    # where detection reported it
    fluxes = photometry(CircularAperture.(table.x, table.y, 2.0), data).aperture_sum
    @test fluxes ≈ table.value

    # sort is a keyword; with sort = false rows come in array (column-major) order
    unsorted = extract_sources(PeakMesh(), data; sort = false)
    @test unsorted.value == [3.0, 5.0]

    # the error map is looked up with the same x/y convention as the returned
    # positions. A non-uniform map on a square image catches a transposed but
    # in-bounds lookup, which a uniform map or a non-square image cannot.
    square = zeros(30, 30)
    square[10, 3] = 3.0
    square[4, 27] = 5.0
    errs = zeros(size(square))
    errs[4, 27] = 2.0 # threshold 8 at the bright peak, 0 at its transpose
    filtered = extract_sources(PeakMesh(nsigma = 4.0), square, errs)
    @test filtered.x == [10]
    @test filtered.y == [3]
    @test filtered.value == [3.0]

    # an error map with other axes is rejected rather than read at the wrong pixel
    @test_throws DimensionMismatch extract_sources(PeakMesh(), data, ones(30, 20))
    @test_throws DimensionMismatch extract_sources(PeakMesh(), data, 1.0)
end

@testset "detection/Detection: offset axes" begin
    # the default error map must follow the axes of `data`, e.g. the output of
    # `imfilter(img, kernel, Inner())`
    data = zeros(20, 30)
    data[10, 3] = 3.0
    data[4, 27] = 5.0
    odata = OffsetArrays.OffsetArray(data, 100:119, -5:24)

    table = extract_sources(PeakMesh(), odata)
    @test table.x == [103, 109]
    @test table.y == [21, -3]
    @test table.value == [5.0, 3.0]
    @test extract_sources(PeakMesh(), odata, nothing) == table
    @test extract_sources(PeakMesh(), odata, zero(odata)) == table
    @test_throws DimensionMismatch extract_sources(PeakMesh(), odata, zeros(20, 30))
end

@testset "detection/Detection: empty result" begin
    for T in (Float64, Int, Real)
        table = extract_sources(PeakMesh(), Matrix{T}(zeros(5, 5)))
        @test length(table) == 0
        @test propertynames(table) == (:x, :y, :value)
        @test eltype(table.value) == T
    end
end

@testset "detection/Detection: interface" begin
    @testset "interface - $P" for P in [PeakMesh()]
        data = randn(100, 100)
        idxs = [1, 700, 1524]
        fake_peaks = randn(length(idxs)) .+ 10
        data[idxs] .= fake_peaks

        table = extract_sources(P, data)
        table2 = extract_sources(P, data, zeros(100, 100))
        @test table == table2
    end
end
