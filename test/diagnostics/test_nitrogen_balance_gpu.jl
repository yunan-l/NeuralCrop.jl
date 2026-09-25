using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA daily nitrogen-balance diagnostics" begin
    cell_size = 2
    crop = init_crop(cell_size, CuArray)
    soil = init_soil(cell_size, soilparams.soildepth, CuArray)
    state = test_model_state(crop, soil)
    balance = init_nitrogen_balance(1, cell_size, CuArray)

    crop.state.nitrogen.total .= 1.0f0
    soil.nitrogen.nitrate .= 0.2f0
    soil.nitrogen.ammonium .= 0.1f0
    soil.nitrogen.litter .= 0.3f0
    soil.nitrogen.fast .= 0.4f0
    soil.nitrogen.slow .= 0.5f0

    NeuralCrop.record_nitrogen_balance_start!(balance, 1, state, state)

    # Root uptake only transfers N within the tracked plant-soil system.
    @views soil.nitrogen.nitrate[1, :] .-= 0.25f0
    crop.state.nitrogen.total .+= 0.25f0
    crop.fluxes.nitrogen.uptake .= 0.25f0

    NeuralCrop.record_nitrogen_balance_end!(balance, 1, state, state)

    @test all(Array(balance.root_uptake) .== 0.25f0)
    @test maximum(abs, Array(balance.residual)) <= 2.0f-6
    @test maximum(abs, Array(balance.relative_residual)) <= 1.0f-6
end
