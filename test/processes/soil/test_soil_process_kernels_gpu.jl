using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA fused soil process kernels" begin
    cells = 1024
    depths = Float32[200, 300, 500, 700, 1000]
    reference = init_soil(cells, depths, identity)
    gpu = init_soil(cells, depths, CuArray)
    crop_reference = init_crop(cells, identity)
    crop_gpu = init_crop(cells, CuArray)
    state_reference = test_model_state(crop_reference, reference)
    state_gpu = test_model_state(crop_gpu, gpu)

    sand = reshape(Float32.(range(0.15, 0.75; length = cells)), 1, :)
    clay = reshape(Float32.(range(0.45, 0.10; length = cells)), 1, :)
    fast = reshape(Float32.(range(40, 220; length = 5 * cells)), 5, cells)
    slow = reshape(Float32.(range(200, 900; length = 5 * cells)), 5, cells)
    saturation = fill(0.45f0, 5, cells)
    storage = fill(80.0f0, 5, cells)
    density_factor = reshape(Float32.(range(0.7, 1.0; length = cells)), 1, :)
    for soil in (reference,)
        soil.properties.sand_fraction .= sand
        soil.properties.clay_fraction .= clay
        soil.carbon.fast .= fast
        soil.carbon.slow .= slow
        soil.water.saturation_fraction .= saturation
        soil.water.storage .= storage
        soil.management.tillage_density_factor .= density_factor
    end
    gpu.properties.sand_fraction .= CuArray(sand)
    gpu.properties.clay_fraction .= CuArray(clay)
    gpu.carbon.fast .= CuArray(fast)
    gpu.carbon.slow .= CuArray(slow)
    gpu.water.saturation_fraction .= CuArray(saturation)
    gpu.water.storage .= CuArray(storage)
    gpu.management.tillage_density_factor .= CuArray(density_factor)

    crop_reference.events.sowing .= Int32.(isodd.(1:cells))
    crop_gpu.events.sowing .= CuArray(Int32.(isodd.(1:cells)))
    NeuralCrop.tillage_hydraulics!(state_reference)
    tillage_hydraulics!(state_gpu)
    synchronize()
    @test Array(gpu.management.tillage_density_factor) ≈
        reference.management.tillage_density_factor

    NeuralCrop.pedotransfer!(state_reference)
    pedotransfer!(state_gpu)
    synchronize()
    for field in (:field_capacity, :saturation_fraction, :beta, :saturated_conductivity)
        @test Array(getproperty(gpu.water, field)) ≈
            getproperty(reference.water, field) rtol = 8.0f-6 atol = 8.0f-6
    end
    for soil in (reference,)
        soil.water.saturation_storage .= 180.0f0
        soil.water.holding_capacity_storage .= 100.0f0
        soil.water.wilting_storage .= 30.0f0
        soil.water.relative_content .= 0.55f0
        soil.water.free_water .= 0.0f0
        soil.water.ice_storage .= 0.0f0
        soil.water.wilting_ice_fraction .= 0.0f0
        soil.water.available_ice_storage .= 0.0f0
        soil.water.free_ice_storage .= 0.0f0
        soil.thermal.temperature .= 15.0f0
        soil.surface_litter.temperature .= 12.0f0
        soil.surface_litter.water_capacity .= 4.0f0
        soil.surface_litter.water_storage .= 2.0f0
        soil.carbon.litter .= 12.0f0
        soil.nitrogen.litter .= 0.5f0
        soil.carbon.fast .= 60.0f0
        soil.carbon.slow .= 300.0f0
        soil.nitrogen.fast .= 5.0f0
        soil.nitrogen.slow .= 25.0f0
        soil.carbon.litter_response .= Float32[0.0012, 0.0008, 0.0005]
        soil.nitrogen.litter_response .= Float32[0.0012, 0.0008, 0.0005]
        soil.decomposition.shift_fast .= 0.2f0
        soil.decomposition.shift_slow .= 0.2f0
        soil.nitrogen.ammonium .= 0.4f0
        soil.nitrogen.nitrate .= 0.8f0
        soil.properties.ph .= 6.5f0
    end
    for (field, value) in (
        (:saturation_storage, 180.0f0), (:holding_capacity_storage, 100.0f0),
        (:wilting_storage, 30.0f0), (:relative_content, 0.55f0),
        (:free_water, 0.0f0), (:ice_storage, 0.0f0),
        (:wilting_ice_fraction, 0.0f0), (:available_ice_storage, 0.0f0),
        (:free_ice_storage, 0.0f0),
    )
        getproperty(gpu.water, field) .= value
    end
    gpu.thermal.temperature .= 15.0f0
    gpu.surface_litter.temperature .= 12.0f0
    gpu.surface_litter.water_capacity .= 4.0f0
    gpu.surface_litter.water_storage .= 2.0f0
    gpu.carbon.litter .= 12.0f0
    gpu.nitrogen.litter .= 0.5f0
    gpu.carbon.fast .= 60.0f0
    gpu.carbon.slow .= 300.0f0
    gpu.nitrogen.fast .= 5.0f0
    gpu.nitrogen.slow .= 25.0f0
    gpu.carbon.litter_response .= CuArray(Float32[0.0012, 0.0008, 0.0005])
    gpu.nitrogen.litter_response .= CuArray(Float32[0.0012, 0.0008, 0.0005])
    gpu.decomposition.shift_fast .= 0.2f0
    gpu.decomposition.shift_slow .= 0.2f0
    gpu.nitrogen.ammonium .= 0.4f0
    gpu.nitrogen.nitrate .= 0.8f0
    gpu.properties.ph .= 6.5f0

    NeuralCrop.soil_carbon!(state_reference, state_reference)
    soil_carbon!(state_gpu, state_gpu)
    air = fill(15.0f0, cells)
    wind = fill(2.0f0, cells)
    NeuralCrop.soil_nitrogen!(
        state_reference, state_reference; air_temperature = air, wind_speed = wind,
    )
    air_gpu = CuArray(air)
    wind_gpu = CuArray(wind)
    soil_nitrogen!(
        state_gpu, state_gpu; air_temperature = air_gpu, wind_speed = wind_gpu,
    )
    synchronize()
    @test Array(gpu.decomposition.response) ≈
        reference.decomposition.response rtol = 1.0f-5 atol = 5.0f-6
    @test Array(gpu.decomposition.litter_response) ≈
        reference.decomposition.litter_response rtol = 1.0f-5 atol = 5.0f-6
    @test Array(gpu.carbon.decomposed_litter) ≈
        reference.carbon.decomposed_litter rtol = 1.0f-5 atol = 5.0f-6
    @test Array(gpu.nitrogen.decomposed_litter) ≈
        reference.nitrogen.decomposed_litter rtol = 1.0f-5 atol = 5.0f-6
    @test Array(gpu.carbon.heterotrophic_respiration) ≈
        reference.carbon.heterotrophic_respiration rtol = 1.0f-5 atol = 5.0f-6
    @test Array(gpu.nitrogen.nitrate) ≈ reference.nitrogen.nitrate rtol = 1.0f-5 atol = 5.0f-6
    @test Array(gpu.nitrogen.ammonium) ≈ reference.nitrogen.ammonium rtol = 1.0f-5 atol = 5.0f-6

    soil_carbon!(state_gpu, state_gpu)
    soil_nitrogen!(state_gpu, state_gpu; air_temperature = air_gpu, wind_speed = wind_gpu)
    synchronize()
    carbon_bytes = CUDA.@allocated begin
        soil_carbon!(state_gpu, state_gpu)
        synchronize()
    end
    nitrogen_bytes = CUDA.@allocated begin
        soil_nitrogen!(state_gpu, state_gpu; air_temperature = air_gpu, wind_speed = wind_gpu)
        synchronize()
    end

    gpu.management.tillage_fraction .= CuArray(Float32[
        0.05 0 0
        0.95 1 0
        0 0 1
    ])
    crop_gpu.events.sowing .= Int32(1)
    litter_tillage!(state_gpu, state_gpu)
    litter_bioturbation!(state_gpu)
    synchronize()
    litter_bytes = CUDA.@allocated begin
        litter_tillage!(state_gpu, state_gpu)
        litter_bioturbation!(state_gpu)
        synchronize()
    end

    crop_gpu.fluxes.water.transpiration_layer .= 0.01f0
    gpu.water.evaporation .= 0.01f0
    soil_evapotranspiration!(state_gpu, state_gpu)
    synchronize()
    water_update_bytes = CUDA.@allocated begin
        soil_evapotranspiration!(state_gpu, state_gpu)
        synchronize()
    end
    # Measure this mutating kernel last: a second pedotransfer call uses the
    # previously updated saturation fraction and therefore changes the model
    # state. Running it before the equivalence checks would no longer compare
    # the same one-step CPU and GPU trajectories.
    pedotransfer_bytes = CUDA.@allocated begin
        pedotransfer!(state_gpu)
        synchronize()
    end
    @test pedotransfer_bytes == 0
    @test carbon_bytes == 0
    @test nitrogen_bytes == 0
    @test litter_bytes == 0
    @test water_update_bytes == 0
end
