using NeuralCrop
using Test

@testset "LPJmL-style mineralization and immobilization" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.nitrogen.decomposed_litter .= 0.0f0
    soil.carbon.decomposed_litter .= 0.0f0
    soil.nitrogen.decomposed_litter[1, 1] = 0.1f0
    soil.carbon.decomposed_litter[1, 1] = 3.0f0
    soil.decomposition.shift_fast[1, 1] = 1.0f0
    soil.decomposition.shift_slow[1, 1] = 1.0f0
    soil.nitrogen.ammonium .= 0.2f0
    soil.nitrogen.nitrate .= 0.1f0

    mineral_before = sum(soil.nitrogen.ammonium .+ soil.nitrogen.nitrate)
    organic_before = sum(soil.nitrogen.fast .+ soil.nitrogen.slow)
    NeuralCrop.launch_custom!(
        NeuralCrop.mineralize_immobilize_kernel!,
        soil.carbon.decomposed_litter,
        1,
        soil.nitrogen.decomposed_litter,
        soil.nitrogen.decomposed_fast,
        soil.nitrogen.decomposed_slow,
        soil.decomposition.shift_fast,
        soil.decomposition.shift_slow,
        soil.nitrogen.ammonium,
        soil.nitrogen.nitrate,
        soil.nitrogen.fast,
        soil.nitrogen.slow,
        soil.properties.layer_depth,
        soil.nitrogen.mineralization,
        soil.nitrogen.immobilization,
        (; lpjmlparams, soil_layers = 5),
    )

    gross_mineralization = sum(soil.nitrogen.mineralization)
    immobilization = sum(soil.nitrogen.immobilization)
    mineral_after = sum(soil.nitrogen.ammonium .+ soil.nitrogen.nitrate)
    organic_after = sum(soil.nitrogen.fast .+ soil.nitrogen.slow)

    @test gross_mineralization ≈ 0.05f0 atol = 1.0f-6
    @test immobilization > 0.0f0
    @test mineral_after - mineral_before ≈
          gross_mineralization - immobilization atol = 1.0f-6
    @test organic_after - organic_before ≈ immobilization atol = 1.0f-6
    @test minimum(soil.nitrogen.ammonium) >= 0.0f0
    @test minimum(soil.nitrogen.nitrate) >= 0.0f0
end
