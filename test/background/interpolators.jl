@testset "IDWInterpolator" begin
    coordinates = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0;
                    1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]

    weights = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0]
    value = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]

    alg = IDWInterpolator(coordinates,value, weights, 8)
    positions = [9.9 1.0; 10.0 1.0]

    @test alg(positions) == IDWInterpolator(coordinates, value)(positions)
    @test_throws ArgumentError alg(positions, n_neighbors = 11)
    @test_throws ErrorException IDWInterpolator(coordinates, value, nothing, -5)
end
