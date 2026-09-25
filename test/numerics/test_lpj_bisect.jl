using NeuralCrop
using Test

@testset "LPJ-compatible bisection" begin
    root, iterations = lpj_bisect(
        x -> x * x - 2.0f0,
        0.0f0,
        2.0f0;
        y_accuracy = 1.0f-5,
        max_iterations = 100,
    )

    @test root ≈ sqrt(2.0f0) atol = 1.0f-5
    @test iterations < 100
end
