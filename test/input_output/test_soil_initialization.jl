using NeuralCrop
using Test

@testset "Generic soil initial state" begin
    selection = (cell_ids = Int32[101],)
    targets = (
        selection,
        soil_organic_carbon = reshape(Float32[100, 200, 300, 400, 500], 5, 1),
        total_nitrogen = reshape(Float32[10, 20, 30, 40, 50], 5, 1),
    )
    soil = (
        selection,
        saturation = fill(0.45f0, 5, 1),
        sand = Float32[0.4],
        clay = Float32[0.2],
        layer_depth = Float32[200, 300, 500, 1000, 1000],
    )

    state = soil_initial_state(targets, soil)
    @test state.fastc + state.slowc == targets.soil_organic_carbon
    @test state.fastc == 0.4f0 .* targets.soil_organic_carbon
    @test state.litc == zeros(Float32, 3, 1)
    @test state.fastn + state.slown + 0.02f0 .* state.slown ≈ targets.total_nitrogen
    @test all(state.swc .> 0)
    @test all(state.swc .< soil.saturation .* reshape(soil.layer_depth, :, 1))
end
