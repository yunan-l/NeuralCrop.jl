using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

function run_soil_cn_gpu_fixture(device)
    cells = 2
    soil = init_soil(Float32, cells, Float32.(soilparams.soildepth), device)
    state = test_model_state(soil)
    soil.carbon.litter .= device(Float32[30 18; 20 12; 10 6])
    soil.nitrogen.litter .= soil.carbon.litter ./ 12.0f0
    soil.carbon.fast .= 12.0f0
    soil.carbon.slow .= 24.0f0
    soil.nitrogen.fast .= 1.2f0
    soil.nitrogen.slow .= 1.2f0
    soil.nitrogen.ammonium .= 0.4f0
    soil.nitrogen.nitrate .= 0.6f0
    soil.carbon.litter_response .= 0.08f0
    soil.nitrogen.litter_response .= 0.08f0
    soil.decomposition.shift_fast .= 0.0f0
    soil.decomposition.shift_slow .= 0.0f0
    for shift in (
        soil.decomposition.shift_fast, soil.decomposition.shift_slow,
    )
        @views shift[1, :] .= 1.0f0
    end
    soil.thermal.temperature .= 10.0f0
    soil.surface_litter.temperature .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.water.holding_capacity_storage .= 60.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.relative_content .= 0.5f0
    soil.properties.ph .= 6.5f0

    soil_cn_decomposition!(state)
    post_crop_nitrogen_losses!(
        state;
        air_temperature = device(fill(20.0f0, cells)),
        wind_speed = device(fill(2.0f0, cells)),
    )
    return soil
end

@testset "CUDA coupled soil C-N decomposition" begin
    cpu = run_soil_cn_gpu_fixture(identity)
    gpu = run_soil_cn_gpu_fixture(CuArray)

    fields = (
        (:carbon_litter, cpu.carbon.litter, gpu.carbon.litter),
        (:carbon_fast, cpu.carbon.fast, gpu.carbon.fast),
        (:carbon_slow, cpu.carbon.slow, gpu.carbon.slow),
        (:respiration, cpu.carbon.heterotrophic_respiration, gpu.carbon.heterotrophic_respiration),
        (:nitrogen_litter, cpu.nitrogen.litter, gpu.nitrogen.litter),
        (:nitrogen_fast, cpu.nitrogen.fast, gpu.nitrogen.fast),
        (:nitrogen_slow, cpu.nitrogen.slow, gpu.nitrogen.slow),
        (:ammonium, cpu.nitrogen.ammonium, gpu.nitrogen.ammonium),
        (:nitrate, cpu.nitrogen.nitrate, gpu.nitrogen.nitrate),
        (:mineralization, cpu.nitrogen.mineralization, gpu.nitrogen.mineralization),
        (:immobilization, cpu.nitrogen.immobilization, gpu.nitrogen.immobilization),
        (:nitrification, cpu.nitrogen.nitrification, gpu.nitrogen.nitrification),
        (:denitrification, cpu.nitrogen.denitrification, gpu.nitrogen.denitrification),
        (:volatilization, cpu.nitrogen.volatilization, gpu.nitrogen.volatilization),
    )
    for (name, expected, actual) in fields
        @test Array(actual) ≈ expected rtol = 2.0f-5 atol = 5.0f-6
    end
end
