using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL-compatible potential and actual LAI" begin
    cells = 4096
    crop = init_crop(Float32, cells, CuArray)
    pet = init_pet(Float32, cells, CuArray)
    state = test_model_state(crop; pet)

    crop.state.phenology.is_growing .= Int32(1)
    crop.state.phenology.growing_days .= Int32(20)
    crop.state.phenology.senescence .= true
    crop.state.canopy.lai .= 0.1f0
    crop.state.canopy.lai_previous_potential .= 0.1f0
    crop.state.canopy.lai_npp_deficit .= 0.3f0
    crop.state.carbon.biomass .= 10.0f0
    crop.state.carbon.leaf .= 1.0f0
    crop.state.carbon.root .= 2.0f0
    crop.state.carbon.pool .= 7.0f0
    crop.state.nitrogen.sufficiency .= 1.0f0
    crop.state.water.sufficiency .= 1.0f0

    carbon_allocation!(cft1, state)
    synchronize()
    @test all(==(0.1f0), Array(crop.state.canopy.lai))
    @test all(==(0.3f0), Array(crop.state.canopy.lai_npp_deficit))
    @test all(iszero, Array(crop.auxiliary.canopy.actual_lai))

    pet.par .= 20.0f0
    apar_crop!(cft1, state, pet)
    pet.eeq .= 2.0f0
    rain = CUDA.fill(5.0f0, cells)
    interception!(state, cft1, pet.eeq, rain)
    synchronize()
    @test all(iszero, Array(crop.auxiliary.canopy.fpar))
    @test all(iszero, Array(crop.auxiliary.canopy.apar))
    @test all(iszero, Array(crop.auxiliary.canopy.canopy_wet))
    @test all(iszero, Array(crop.fluxes.water.interception))

    # Exercise the non-senescence branch that uses LPJmL's `lai000` state.
    crop.state.phenology.senescence .= false
    crop.state.canopy.lai .= 1.0f0
    crop.state.canopy.lai_previous_potential .= 1.0f0
    crop.auxiliary.canopy.flaimax .= 0.8f0
    crop.state.water.sufficiency .= 0.6f0
    crop.state.nitrogen.sufficiency .= 0.9f0
    lai_crop!(state, cft1)
    synchronize()
    expected = 1.0f0 + (0.8f0 * cft1.laimax - 1.0f0) * 0.6f0
    @test all(≈(expected), Array(crop.state.canopy.lai))
end
