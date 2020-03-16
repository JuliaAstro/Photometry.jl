using Photometry.Background: ShepherdInterpolator

@testset "ShepherdInterpolator" begin
    coordinates = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0;
                        1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]
    value = [1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]
    positions = [9.9 1.0; 10.0 1.0]

    output = ShepherdInterpolator(coordinates, value, positions)
    @test output[1] ≈ 9.57581 rtol = 1e-5
    @test output[2] ≈ 1.0 rtol = 1e-5
end
