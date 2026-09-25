using NeuralCrop
using Test

@testset "LPJmL-style ammonia volatilization" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.properties.ph .= 7.0f0
    soil.nitrogen.ammonium[1, 1] = 1.0f0
    air_temperature = Float32[20.0]
    # Use low, non-saturating wind speeds so the monotonic wind response is
    # observable instead of both cases being capped by the available NH4.
    wind_speed = Float32[0.01]
    ammonium_before = soil.nitrogen.ammonium[1, 1]

    NeuralCrop.launch_1D!(
        NeuralCrop.volatilization_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        air_temperature,
        wind_speed,
        soil.properties.layer_depth,
        soil.nitrogen.volatilization,
        lpjmlparams,
    )

    flux = soil.nitrogen.volatilization[1]
    @test 0.0f0 < flux <= ammonium_before
    @test ammonium_before - soil.nitrogen.ammonium[1, 1] ≈ flux atol = 1.0f-7
    @test all(soil.nitrogen.ammonium .>= 0.0f0)

    soil.nitrogen.ammonium[1, 1] = ammonium_before
    high_wind = Float32[0.02]
    NeuralCrop.launch_1D!(
        NeuralCrop.volatilization_kernel!,
        soil.properties.ph,
        soil.nitrogen.ammonium,
        air_temperature,
        high_wind,
        soil.properties.layer_depth,
        soil.nitrogen.volatilization,
        lpjmlparams,
    )
    @test soil.nitrogen.volatilization[1] > flux
end
