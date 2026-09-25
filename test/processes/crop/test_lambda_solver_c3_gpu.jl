using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA C3 lambda solver" begin
    crop = init_crop(2, CuArray)
    photos = crop.auxiliary.photosynthesis
    pet = init_pet(2, CuArray)
    state = test_model_state(crop; pet)

    target_lambda = 0.5f0
    vcmax = 2.0f0
    tstress = 1.0f0
    # Match the full simulation: one annual CO2 value shared by all cells.
    co2_cpu = Float32[40.0]
    temp_cpu = Float32[20.0, 20.0]
    apar = 1.0f6
    daylength = 12.0f0
    fpar = 0.8f0

    target_adtmm = NeuralCrop.c3_adtmm_scalar(
        target_lambda, vcmax, tstress, cft1.b, co2_cpu[1], temp_cpu[1], apar, daylength,
    )
    fac = target_adtmm / (1.0f0 - target_lambda)
    gpd = fac * 1.6f0 / (co2_cpu[1] * 1.0f-5)
    target_conductance = gpd / (daylength * 3600.0f0) + cft1.gmin * fpar

    photos.vcmax .= vcmax
    photos.temperature_stress .= tstress
    crop.auxiliary.canopy.apar .= apar
    crop.auxiliary.canopy.fpar .= fpar
    crop.auxiliary.canopy.canopy_conductance .= target_conductance
    pet.daylength .= daylength

    solve_lambda_c3!(
        cft1, state, pet, CuArray(temp_cpu), CuArray(co2_cpu),
    )

    lambda_cpu = Array(photos.lambda)
    @test lambda_cpu[1] ≈ target_lambda atol = 2.0f-3
    @test lambda_cpu[2] ≈ target_lambda atol = 2.0f-3

    crop.auxiliary.canopy.canopy_conductance .= CuArray(Float32[target_conductance, 0.0])
    solve_lambda_c3!(
        cft1, state, pet, CuArray(temp_cpu), CuArray(co2_cpu),
    )
    @test Array(photos.lambda)[2] == 0.0f0
end
