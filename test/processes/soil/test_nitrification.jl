using NeuralCrop
using Test

@testset "LPJmL-style nitrification" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.properties.ph .= 6.5f0
    soil.nitrogen.ammonium[1, 1] = 0.8f0
    soil.water.relative_content[1, 1] = 0.6f0
    soil.water.holding_capacity_storage[1, 1] = 100.0f0
    soil.water.saturation_storage[1, 1] = 100.0f0
    soil.thermal.temperature[1, 1] = 20.0f0
    ammonium_before = soil.nitrogen.ammonium[1, 1]
    nitrate_before = soil.nitrogen.nitrate[1, 1]

    NeuralCrop.launch_1D!(
        NeuralCrop.nitrify_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        soil.nitrogen.nitrate,
        soil.water.relative_content,
        soil.water.holding_capacity_storage,
        soil.water.wilting_storage,
        soil.water.wilting_ice_fraction,
        soil.water.free_water,
        soil.water.saturation_storage,
        soil.thermal.temperature,
        soil.nitrogen.nitrification,
        soil.nitrogen.n2o_nitrification,
        (; lpjmlparams, soil_layers = 5),
    )

    gross = soil.nitrogen.nitrification[1, 1]
    n2o = soil.nitrogen.n2o_nitrification[1, 1]
    @test 0.0f0 < gross <= ammonium_before
    @test n2o ≈ lpjmlparams.k_2 * gross atol = 1.0f-7
    @test ammonium_before - soil.nitrogen.ammonium[1, 1] ≈ gross atol = 1.0f-7
    @test soil.nitrogen.nitrate[1, 1] - nitrate_before ≈ gross - n2o atol = 1.0f-7
end
