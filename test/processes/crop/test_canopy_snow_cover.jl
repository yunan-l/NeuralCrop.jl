using NeuralCrop
using Test

@testset "LPJmL crop snow cover suppresses absorbed PAR" begin
    for T in (Float32, Float64), cft in (cft1, cft3)
        crop = init_crop(T, 2, identity)
        pet = init_pet(T, 2, identity)
        state = test_model_state(crop; pet)
        crop.state.canopy.lai .= T(2)
        pet.par .= T(20)
        snow_height = T[0, 0.1]

        if cft === cft3
            apar_crop_maize!(cft, state, pet, snow_height)
        else
            apar_crop!(cft, state, pet, snow_height)
        end

        @test crop.auxiliary.canopy.fpar[1] > zero(T)
        @test crop.auxiliary.canopy.apar[1] > zero(T)
        @test crop.auxiliary.canopy.fpar[2] == zero(T)
        @test crop.auxiliary.canopy.apar[2] == zero(T)
    end
end

@testset "LPJmL surface albedo mixes crop, litter, soil, and snow" begin
    cells = 5
    crop = init_crop(cells, identity)
    soil = init_soil(cells, soilparams.soildepth, identity)
    pet = init_pet(cells, identity)
    state = test_model_state(crop, soil; pet)
    crop.state.phenology.is_growing .= Int32[0, 1, 1, 0, 1]
    crop.state.canopy.lai .= Float32[0, 0, 0, 0, 2]

    litter_cover = Float32[0, 0, 0.5, 0, 0.25]
    soil.carbon.litter[1, :] .=
        -log.(1.0f0 .- litter_cover) ./ 6.0f-3 .* 0.42f0
    soil.snow.height .= Float32[0, 0, 0, 0.1, 0.1]
    soil.snow.fraction .= Float32[0, 0, 0, 0.5, 0.6]

    albedo!(cft1, state, state, pet)

    @test pet.albedo[1] ≈ 0.3f0
    @test pet.albedo[2] ≈ 0.3f0
    @test pet.albedo[3] ≈ 0.5f0 * cft1.albedo_litter + 0.5f0 * 0.3f0
    @test pet.albedo[4] ≈ 0.5f0 * 0.65f0 + 0.5f0 * 0.3f0

    green_fraction = 1.0f0 - exp(-cft1.lightextcoeff * 2.0f0)
    expected_snow_albedo = green_fraction * 0.65f0 +
        0.25f0 * (1.0f0 - green_fraction) * 0.65f0 +
        0.75f0 * (1.0f0 - green_fraction) * 0.6f0 * 0.65f0
    @test pet.albedo[5] ≈ expected_snow_albedo
    @test crop.auxiliary.canopy.albedo[1] == 0.0f0
    @test crop.auxiliary.canopy.albedo[5] ≈ expected_snow_albedo
end
