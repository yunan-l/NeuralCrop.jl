using NeuralCrop
using Test

@testset "Photosynthesis scalar formulas" begin
    T = Float32
    gross = NeuralCrop.compute_co_limited_assimilation(T(3), T(5), T(0.9), T(12))
    @test gross > zero(T)
    @test isfinite(gross)

    net, daily_net = NeuralCrop.compute_net_assimilation(T(5), T(1), T(12))
    @test net == daily_net == T(4.5)
    @test NeuralCrop.compute_water_limited_assimilation(daily_net, T(12), T(20), T(101325)) > zero(T)

    net, daily_net = NeuralCrop.compute_net_assimilation(T(1), T(1), T(24))
    @test net == zero(T)
    @test NeuralCrop.compute_water_limited_assimilation(daily_net, T(12), T(20), T(101325)) == zero(T)
end

@testset "C3 photosynthesis CPU smoke test" begin
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    photos = crop.auxiliary.photosynthesis
    photos.temperature_stress .= 1.0f0

    photosynthesis_C3!(
        cft1,
        state,
        Float32[10.0],
        Float32[12.0],
        Float32[20.0],
        Float32[40.0];
        comp_vcmax = true,
    )

    @test photos.lambda == Float32[0.8]
    @test all(isfinite, photos.vcmax)
    @test all(isfinite, crop.fluxes.carbon.gross_assimilation)
    @test all(isfinite, crop.fluxes.carbon.leaf_respiration)
    @test all(isfinite, crop.fluxes.carbon.net_assimilation)
    @test all(isfinite, crop.fluxes.carbon.water_limited_assimilation)
    @test all(photos.vcmax .>= 0.0f0)
    @test all(crop.fluxes.carbon.gross_assimilation .>= 0.0f0)
    @test all(crop.fluxes.carbon.net_assimilation .>= 0.0f0)

    photos.temperature_stress .= 0.0f0
    photosynthesis_C3!(
        cft1,
        state,
        Float32[10.0],
        Float32[12.0],
        Float32[20.0],
        Float32[40.0];
        comp_vcmax = true,
    )
    @test all(iszero, photos.vcmax)
    @test all(iszero, crop.fluxes.carbon.gross_assimilation)
    @test all(iszero, crop.fluxes.carbon.leaf_respiration)
    @test all(iszero, crop.fluxes.carbon.net_assimilation)
    @test all(iszero, crop.fluxes.carbon.water_limited_assimilation)

    photos.temperature_stress .= 1.0f0
    photos.lambda .= 0.8f0
    photos.vcmax .= 1.0f0
    photosynthesis_C3!(
        cft1,
        state,
        Float32[0.0],
        Float32[12.0],
        Float32[20.0],
        Float32[40.0];
        comp_vcmax = false,
    )
    @test all(iszero, crop.fluxes.carbon.net_assimilation)
    @test all(iszero, crop.fluxes.carbon.water_limited_assimilation)
end

@testset "C4 photosynthesis CPU smoke test" begin
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    photos = crop.auxiliary.photosynthesis
    photos.temperature_stress .= 1.0f0

    photosynthesis_C4!(
        cft3,
        state,
        Float32[10.0],
        Float32[12.0],
        Float32[25.0];
        comp_vcmax = true,
    )

    @test photos.lambda == Float32[0.8]
    @test all(isfinite, photos.vcmax)
    @test all(isfinite, crop.fluxes.carbon.gross_assimilation)
    @test all(isfinite, crop.fluxes.carbon.leaf_respiration)
    @test all(isfinite, crop.fluxes.carbon.net_assimilation)
    @test all(isfinite, crop.fluxes.carbon.water_limited_assimilation)
    @test all(photos.vcmax .>= 0.0f0)
    @test all(crop.fluxes.carbon.gross_assimilation .>= 0.0f0)
    @test all(crop.fluxes.carbon.net_assimilation .>= 0.0f0)

    photos.temperature_stress .= 0.0f0
    photosynthesis_C4!(
        cft3,
        state,
        Float32[10.0],
        Float32[12.0],
        Float32[25.0];
        comp_vcmax = true,
    )
    @test all(iszero, photos.vcmax)
    @test all(iszero, crop.fluxes.carbon.gross_assimilation)
    @test all(iszero, crop.fluxes.carbon.leaf_respiration)
    @test all(iszero, crop.fluxes.carbon.net_assimilation)
    @test all(iszero, crop.fluxes.carbon.water_limited_assimilation)
end
