using NeuralCrop
using Test

@testset "Backend-compatible C3 lambda kernel" begin
    crop = init_crop(2, identity)
    photos = crop.auxiliary.photosynthesis
    pet = init_pet(2, identity)
    state = test_model_state(crop; pet)

    target_lambda = 0.5f0
    vcmax = 2.0f0
    tstress = 1.0f0
    # Annual CO2 is a shared scalar buffer in the full simulation.
    co2 = Float32[40.0]
    temp = Float32[20.0, 20.0]
    apar = 1.0f6
    daylength = 12.0f0
    fpar = 0.8f0

    target_adtmm = NeuralCrop.c3_adtmm_scalar(
        target_lambda, vcmax, tstress, cft1.b, co2[1], temp[1], apar, daylength,
    )
    fac = target_adtmm / (1.0f0 - target_lambda)
    gpd = fac * 1.6f0 / (co2[1] * 1.0f-5)
    target_conductance = gpd / (daylength * 3600.0f0) + cft1.gmin * fpar

    photos.vcmax .= vcmax
    photos.temperature_stress .= tstress
    crop.auxiliary.canopy.apar .= apar
    crop.auxiliary.canopy.fpar .= fpar
    crop.auxiliary.canopy.canopy_conductance .= target_conductance
    pet.daylength .= daylength

    solve_lambda_c3!(cft1, state, pet, temp, co2)

    @test photos.lambda[1] ≈ target_lambda atol = 2.0f-3
    @test photos.lambda[2] ≈ target_lambda atol = 2.0f-3

    crop.auxiliary.canopy.canopy_conductance[2] = 0.0f0
    solve_lambda_c3!(cft1, state, pet, temp, co2)
    @test photos.lambda[2] == 0.0f0
end

@testset "LPJ-compatible C3 lambda solver" begin
    vcmax = 2.0f0
    tstress = 1.0f0
    b = cft1.b
    co2 = 40.0f0
    temp = 20.0f0
    apar = 1.0f6
    daylength = 12.0f0
    target_lambda = 0.5f0

    target_adtmm = NeuralCrop.c3_adtmm_scalar(
        target_lambda, vcmax, tstress, b, co2, temp, apar, daylength,
    )
    @test target_adtmm > 0.0f0

    crop = init_crop(1, identity)
    photos = crop.auxiliary.photosynthesis
    state = test_model_state(crop)
    photos.temperature_stress .= tstress
    photos.lambda .= target_lambda
    photos.vcmax .= vcmax
    photosynthesis_C3!(
        cft1,
        state,
        Float32[apar],
        Float32[daylength],
        Float32[temp],
        Float32[co2];
        comp_vcmax = false,
    )
    @test target_adtmm ≈ crop.fluxes.carbon.water_limited_assimilation[1] atol = 1.0f-6

    fac = target_adtmm / (1.0f0 - target_lambda)

    lambda, iterations, residual = solve_lambda_c3_lpj(
        fac, vcmax, tstress, b, co2, temp, apar, daylength,
    )

    @test 0.02f0 <= lambda <= 0.85f0
    @test lambda ≈ target_lambda atol = 2.0f-3
    @test abs(residual) < 1.0f-3
    @test iterations <= 30
end
