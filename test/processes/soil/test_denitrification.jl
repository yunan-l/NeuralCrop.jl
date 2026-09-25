using NeuralCrop
using Test

@testset "LPJmL-style denitrification" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.nitrogen.nitrate[1, 1] = 0.8f0
    soil.carbon.fast[1, 1] = 10.0f0
    soil.water.relative_content[1, 1] = 0.95f0
    soil.water.holding_capacity_storage[1, 1] = 100.0f0
    soil.water.saturation_storage[1, 1] = 100.0f0
    soil.thermal.temperature[1, 1] = 20.0f0
    nitrate_before = soil.nitrogen.nitrate[1, 1]

    NeuralCrop.launch_1D!(
        NeuralCrop.denitrify_kernel!,
        soil.properties.ph,
        soil.carbon.fast,
        soil.carbon.slow,
        soil.water.relative_content,
        soil.water.holding_capacity_storage,
        soil.water.wilting_storage,
        soil.water.wilting_ice_fraction,
        soil.water.free_water,
        soil.water.saturation_storage,
        soil.thermal.temperature,
        soil.nitrogen.nitrate,
        soil.nitrogen.denitrification,
        soil.nitrogen.n2o_denitrification,
        soil.nitrogen.n2_denitrification,
        (; lpjmlparams, soil_layers = 5),
    )

    gross = soil.nitrogen.denitrification[1, 1]
    n2o = soil.nitrogen.n2o_denitrification[1, 1]
    n2 = soil.nitrogen.n2_denitrification[1, 1]
    @test 0.0f0 < gross <= nitrate_before
    @test nitrate_before - soil.nitrogen.nitrate[1, 1] ≈ gross atol = 1.0f-7
    @test n2o ≈ lpjmlparams.n2o_denit_frac * gross atol = 1.0f-7
    @test n2 + n2o ≈ gross atol = 1.0f-7
end
