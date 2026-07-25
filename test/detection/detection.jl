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
    data[10, 3] = 5.0
    data[4, 27] = 3.0

    table = extract_sources(PeakMesh(), data)
    @test table.x == [10, 4]
    @test table.y == [3, 27]
    @test table.value == [5.0, 3.0]

    # detected positions feed directly into apertures: the flux is exactly
    # where detection reported it
    fluxes = photometry(CircularAperture.(table.x, table.y, 2.0), data).aperture_sum
    @test fluxes ≈ table.value

    # threshold filtering on a non-square image (regression: the filter used
    # to index the error array transposed, which errored or silently compared
    # against the wrong pixel)
    errs = ones(size(data))
    filtered = extract_sources(PeakMesh(nsigma = 4.0), data, errs)
    @test filtered.x == [10]
    @test filtered.y == [3]
    @test filtered.value == [5.0]
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
