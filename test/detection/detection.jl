using Photometry.Detection:
    PeakMesh,
    extract_sources
using Photometry.Aperture:
    CircularAperture,
    photometry

import Random
Random.seed!(8462852)

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
