using NeuralCrop
using Test

@testset "LPJ CO2 units and minimum conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= Float32[100, 150, 250, 400, 400]
    pedotransfer!(state)

    crop.state.phenology.is_growing .= 1
    crop.state.carbon.root .= 1000.0f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1] = 1.0f0
    crop.auxiliary.canopy.fpar .= 0.8f0
    crop.auxiliary.canopy.canopy_wet .= 0.0f0
    pet.eeq .= 1.0f0
    pet.daylength .= 12.0f0

    adtmm = Float32[5.0]
    co2_pa = Float32[40.0]
    expected_gp = 1.6f0 * adtmm[1] /
                  (co2_pa[1] * 1.0f-5 * (1.0f0 - 0.8f0) * 12.0f0 * 3600.0f0) +
                  cft1.gmin * crop.auxiliary.canopy.fpar[1]

    transpiration!(adtmm, cft1, state, pet, state, co2_pa)

    @test crop.auxiliary.canopy.canopy_conductance[1] ≈ expected_gp rtol = 1.0f-6
    @test crop.auxiliary.canopy.canopy_conductance[1] > cft1.gmin * crop.auxiliary.canopy.fpar[1]
end

@testset "Inactive crop has zero conductance" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    state = test_model_state(crop, soil; pet)

    crop.state.phenology.is_growing .= 0
    crop.auxiliary.canopy.fpar .= 0.8f0
    pet.daylength .= 12.0f0
    pet.eeq .= 1.0f0

    transpiration!(Float32[5.0], cft1, state, pet, state, Float32[40.0])

    @test crop.auxiliary.canopy.canopy_conductance[1] == 0.0f0
end
