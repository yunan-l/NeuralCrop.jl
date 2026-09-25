using NeuralCrop
using Test

@testset "LPJmL equilibrium evaporation" begin
    for T in (Float32, Float64)
        # LPJmL's petpar only clips negative values. This deliberately
        # high-energy forcing guards against reintroducing a 15 mm d⁻¹ cap.
        evaporation = NeuralCrop.compute_equilibrium_evaporation(
            T(40), T(1_000), T(500), T(0.2), T(14), T(86_400),
        )
        @test evaporation > T(15)
    end
end
