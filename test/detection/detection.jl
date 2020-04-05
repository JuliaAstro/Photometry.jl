@testset "peak finding - $P" for P in [PeakMesh()]
    data = randn(100, 100)
    idxs = [1, 700, 1524]
    fake_peaks = randn(length(idxs)) .+ 10
    data[idxs] .= fake_peaks

    table = extract_sources(P, data)
    @test table.value[1:3] == sort(fake_peaks, rev = true)
end

@testset "Peak Mesh" begin
    @test PeakMesh(box_size = 3) == PeakMesh(box_size = (3, 3))
end

@testset "interface - $P" for P in [PeakMesh()]
    data = randn(100, 100)
    idxs = [1, 700, 1524]
    fake_peaks = randn(length(idxs)) .+ 10
    data[idxs] .= fake_peaks

    table = extract_sources(P, data)
    table2 = extract_sources(P, data, zeros(100, 100))
    @test table == table2
end
