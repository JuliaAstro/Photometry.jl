@testset "IDWInterpolator" begin
    @test IDWInterpolator((2,3), n_neighbors = 4)(ones(3, 2)) == ones(6,6)
    @test IDWInterpolator(2)(ones(3,4)) == IDWInterpolator((2,2))(ones(3,4))
    @test_throws ErrorException IDWInterpolator(2)(ones(3,2))
end
