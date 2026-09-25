using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA allocation-free crop process kernels" begin
    cells = 4096
    latitude_cpu = Float32.(range(-70, 70; length = cells))
    temperature_cpu = Float32.(range(-10, 40; length = cells))
    longwave_cpu = Float32.(range(-100, 20; length = cells))
    shortwave_cpu = Float32.(range(0, 350; length = cells))
    apar_cpu = Float32.(range(0, 30; length = cells))
    daylength_cpu = Float32.(range(6, 18; length = cells))
    stress_cpu = Float32.(range(0, 1; length = cells))

    pet_reference = init_pet(cells, identity)
    pet_gpu = init_pet(cells, CuArray)
    albedo_cpu = fill(0.2f0, cells)
    pet_reference.albedo .= albedo_cpu
    pet_gpu.albedo .= CuArray(albedo_cpu)
    NeuralCrop.petpar!(
        pet_reference, 172, latitude_cpu, temperature_cpu, longwave_cpu, shortwave_cpu,
    )
    petpar!(
        pet_gpu, 172, CuArray(latitude_cpu), CuArray(temperature_cpu),
        CuArray(longwave_cpu), CuArray(shortwave_cpu),
    )
    synchronize()
    @test Array(pet_gpu.daylength) ≈ pet_reference.daylength rtol = 4.0f-6
    @test Array(pet_gpu.par) ≈ pet_reference.par rtol = 4.0f-6
    @test Array(pet_gpu.eeq) ≈ pet_reference.eeq rtol = 5.0f-6 atol = 3.0f-7

    crop_reference = init_crop(cells, identity)
    crop_gpu = init_crop(cells, CuArray)
    state_reference = test_model_state(crop_reference)
    state_gpu = test_model_state(crop_gpu)
    crop_reference.auxiliary.photosynthesis.temperature_stress .= stress_cpu
    crop_gpu.auxiliary.photosynthesis.temperature_stress .= CuArray(stress_cpu)
    c3_gross_destination = crop_gpu.fluxes.carbon.gross_assimilation
    c3_vcmax_destination = crop_gpu.auxiliary.photosynthesis.vcmax
    NeuralCrop.photosynthesis_C3!(
        cft1, state_reference, apar_cpu, daylength_cpu,
        temperature_cpu, Float32[40]; comp_vcmax = true,
    )
    apar_gpu = CuArray(apar_cpu)
    daylength_gpu = CuArray(daylength_cpu)
    temperature_gpu = CuArray(temperature_cpu)
    co2_gpu = CuArray(Float32[40])
    photosynthesis_C3!(
        cft1, state_gpu, apar_gpu, daylength_gpu,
        temperature_gpu, co2_gpu; comp_vcmax = true,
    )
    synchronize()
    @test crop_gpu.fluxes.carbon.gross_assimilation === c3_gross_destination
    @test crop_gpu.auxiliary.photosynthesis.vcmax === c3_vcmax_destination
    @test Array(crop_gpu.auxiliary.photosynthesis.vcmax) ≈
        crop_reference.auxiliary.photosynthesis.vcmax rtol = 5.0f-6 atol = 3.0f-7
    @test Array(crop_gpu.fluxes.carbon.gross_assimilation) ≈
        crop_reference.fluxes.carbon.gross_assimilation rtol = 5.0f-6 atol = 3.0f-7

    crop_c4_reference = init_crop(cells, identity)
    crop_c4_gpu = init_crop(cells, CuArray)
    state_c4_reference = test_model_state(crop_c4_reference)
    state_c4_gpu = test_model_state(crop_c4_gpu)
    crop_c4_reference.auxiliary.photosynthesis.temperature_stress .= stress_cpu
    crop_c4_gpu.auxiliary.photosynthesis.temperature_stress .= CuArray(stress_cpu)
    NeuralCrop.photosynthesis_C4!(
        cft3, state_c4_reference, apar_cpu, daylength_cpu,
        temperature_cpu; comp_vcmax = true,
    )
    photosynthesis_C4!(
        cft3, state_c4_gpu, apar_gpu, daylength_gpu,
        temperature_gpu; comp_vcmax = true,
    )
    synchronize()
    @test Array(crop_c4_gpu.auxiliary.photosynthesis.vcmax) ≈
        crop_c4_reference.auxiliary.photosynthesis.vcmax rtol = 5.0f-6 atol = 3.0f-7
    @test Array(crop_c4_gpu.fluxes.carbon.gross_assimilation) ≈
        crop_c4_reference.fluxes.carbon.gross_assimilation rtol = 5.0f-6 atol = 3.0f-7

    root = Float32.(range(1, 30; length = cells))
    storage = Float32.(range(0, 20; length = cells))
    pool = Float32.(range(2, 12; length = cells))
    growing = fill(Int32(1), cells)
    gross = Float32.(range(0, 18; length = cells))
    leaf = Float32.(range(0, 2; length = cells))
    soil_temperature_cpu = reshape(Float32.(range(-10, 25; length = cells)), 1, :)
    for (field, values) in ((:root, root), (:storage, storage), (:pool, pool))
        getproperty(crop_reference.state.carbon, field) .= values
        getproperty(crop_gpu.state.carbon, field) .= CuArray(values)
    end
    crop_reference.state.phenology.is_growing .= growing
    crop_gpu.state.phenology.is_growing .= CuArray(growing)
    NeuralCrop.respiration!(
        state_reference, cft1, temperature_cpu, soil_temperature_cpu, gross, leaf,
    )
    gross_gpu = CuArray(gross)
    leaf_gpu = CuArray(leaf)
    soil_temperature_gpu = CuArray(soil_temperature_cpu)
    respiration_destination = crop_gpu.fluxes.carbon.respiration
    respiration!(
        state_gpu, cft1, temperature_gpu, soil_temperature_gpu, gross_gpu, leaf_gpu,
    )
    synchronize()
    @test crop_gpu.fluxes.carbon.respiration === respiration_destination
    @test Array(crop_gpu.fluxes.carbon.respiration) ≈
        crop_reference.fluxes.carbon.respiration rtol = 4.0f-6 atol = 3.0f-7

    # Warmed kernels must not request new device memory; all results are
    # written into arrays allocated by init_crop/init_pet.
    pet_device_bytes = CUDA.@allocated begin
        petpar!(
            pet_gpu, 172, CuArray(latitude_cpu), CuArray(temperature_cpu),
            CuArray(longwave_cpu), CuArray(shortwave_cpu),
        )
        synchronize()
    end
    # Input construction above is intentionally counted, so separately check
    # the steady-state call with already resident inputs.
    latitude_gpu = CuArray(latitude_cpu)
    longwave_gpu = CuArray(longwave_cpu)
    shortwave_gpu = CuArray(shortwave_cpu)
    petpar!(pet_gpu, 172, latitude_gpu, temperature_gpu, longwave_gpu, shortwave_gpu)
    synchronize()
    steady_pet_device_bytes = CUDA.@allocated begin
        petpar!(pet_gpu, 172, latitude_gpu, temperature_gpu, longwave_gpu, shortwave_gpu)
        synchronize()
    end
    steady_photo_device_bytes = CUDA.@allocated begin
        photosynthesis_C3!(
            cft1, state_gpu, apar_gpu, daylength_gpu,
            temperature_gpu, co2_gpu; comp_vcmax = true,
        )
        synchronize()
    end
    steady_c4_device_bytes = CUDA.@allocated begin
        photosynthesis_C4!(
            cft3, state_c4_gpu, apar_gpu, daylength_gpu,
            temperature_gpu; comp_vcmax = true,
        )
        synchronize()
    end
    steady_respiration_device_bytes = CUDA.@allocated begin
        respiration!(
            state_gpu, cft1, temperature_gpu, soil_temperature_gpu, gross_gpu, leaf_gpu,
        )
        synchronize()
    end

    @test pet_device_bytes > 0
    @test steady_pet_device_bytes == 0
    @test steady_photo_device_bytes == 0
    @test steady_c4_device_bytes == 0
    @test steady_respiration_device_bytes == 0

    pet_seconds = @elapsed begin
        for _ in 1:50
            petpar!(pet_gpu, 172, latitude_gpu, temperature_gpu, longwave_gpu, shortwave_gpu)
        end
        synchronize()
    end
    photo_seconds = @elapsed begin
        for _ in 1:50
            photosynthesis_C3!(
                cft1, state_gpu, apar_gpu, daylength_gpu,
                temperature_gpu, co2_gpu; comp_vcmax = true,
            )
        end
        synchronize()
    end
    c4_seconds = @elapsed begin
        for _ in 1:50
            photosynthesis_C4!(
                cft3, state_c4_gpu, apar_gpu, daylength_gpu,
                temperature_gpu; comp_vcmax = true,
            )
        end
        synchronize()
    end
    respiration_seconds = @elapsed begin
        for _ in 1:50
            respiration!(
                state_gpu, cft1, temperature_gpu, soil_temperature_gpu,
                gross_gpu, leaf_gpu,
            )
        end
        synchronize()
    end
    benchmark = (
        cells = cells,
        kernel_device_bytes = (
            pet = steady_pet_device_bytes,
            c3 = steady_photo_device_bytes,
            c4 = steady_c4_device_bytes,
            respiration = steady_respiration_device_bytes,
        ),
        kernel_seconds_50_calls = (
            pet = pet_seconds,
            c3 = photo_seconds,
            c4 = c4_seconds,
            respiration = respiration_seconds,
        ),
    )
    @info "CUDA crop-kernel benchmark" benchmark
end
