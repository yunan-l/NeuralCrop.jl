using NeuralCrop
using Test

@testset "LPJmL tillage–hydraulic coupling" begin
    cells = 2
    settled = init_soil(cells, soilparams.soildepth, identity)
    crop = init_crop(cells, identity)
    settled_state = test_model_state(crop, settled)
    settled.properties.sand_fraction .= 0.4f0
    settled.properties.clay_fraction .= 0.2f0
    settled.carbon.fast .= 60.0f0
    settled.carbon.slow .= 300.0f0
    settled.water.saturation_fraction .= 0.45f0
    settled.water.storage .= reshape(Float32[40, 60, 100, 200, 200], 5, 1)
    tilled = deepcopy(settled)
    tilled_state = test_model_state(crop, tilled)

    crop.events.sowing .= Int32[1, 0]
    tillage_hydraulics!(tilled_state)
    expected_density_factor = 1.0f0 - (1.0f0 - 0.667f0) * 0.9f0
    @test tilled.management.tillage_density_factor ≈
        reshape(Float32[expected_density_factor, 1.0f0], 1, :)

    water_before = vec(sum(
        tilled.water.storage .+ tilled.water.ice_storage; dims = 1,
    ))
    pedotransfer!(settled_state)
    pedotransfer!(tilled_state)

    @test tilled.water.saturation_fraction[1, 1] >
        settled.water.saturation_fraction[1, 1]
    @test tilled.water.field_capacity[1, 1] > settled.water.field_capacity[1, 1]
    @test tilled.water.saturated_conductivity[1, 1] >
        settled.water.saturated_conductivity[1, 1]
    @test tilled.water.saturation_fraction[:, 2] ≈
        settled.water.saturation_fraction[:, 2]
    @test tilled.water.saturation_fraction[2:end, 1] ≈
        settled.water.saturation_fraction[2:end, 1]
    @test vec(sum(tilled.water.storage .+ tilled.water.ice_storage; dims = 1)) ≈
        water_before

    density_before_rain = tilled.management.tillage_density_factor[1, 1]
    soil_infiltration!(tilled_state, tilled_state, Float32[20, 0])
    top_infiltration = tilled.water.influx[1, 1]
    settling_index = 0.2f0 * top_infiltration * (
        1.0f0 + 2.0f0 * 40.0f0 / (40.0f0 + exp(8.597f0 - 0.075f0 * 40.0f0))
    ) / 0.2f0^0.6f0
    settling_fraction = settling_index /
        (settling_index + exp(3.92f0 - 0.0226f0 * settling_index))
    expected_settled_factor = density_before_rain +
        settling_fraction * (1.0f0 - density_before_rain)
    @test tilled.management.tillage_density_factor[1, 1] ≈ expected_settled_factor
    @test density_before_rain < tilled.management.tillage_density_factor[1, 1] <= 1.0f0
    @test tilled.management.tillage_density_factor[1, 2] == 1.0f0
end
