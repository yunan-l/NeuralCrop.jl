using NeuralCrop
using Enzyme
using Test

function _mutating_scalar_objective(theta, state, scale)
    state[1] = theta[1] * scale
    state[2] = theta[2] * theta[2]
    return state[1] + state[2]
end

@testset "Enzyme adapter directional and coordinate gradients" begin
    theta = Float32[2, 3]
    state_factory() = (zeros(Float32, 2), zeros(Float32, 2))

    directional = enzyme_forward_directional(
        _mutating_scalar_objective,
        theta,
        Float32[1, -1],
        state_factory,
        2.0f0,
    )
    @test directional ≈ -4.0f0

    gradient = enzyme_forward_gradient(
        _mutating_scalar_objective,
        theta,
        state_factory,
        2.0f0,
    )
    @test gradient ≈ Float32[2, 6]

    tangent = enzyme_zero_tangent((values = Float32[1, 2],))
    @test tangent.values == Float32[0, 0]
end
