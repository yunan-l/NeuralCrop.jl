using NeuralCrop
using Test

@testset "LPJmL-compatible potential and actual LAI" begin
    for T in (Float32, Float64)
        crop = init_crop(T, 1, identity)
        pet = init_pet(T, 1, identity)
        state = test_model_state(crop; pet)

        crop.state.phenology.is_growing .= Int32(1)
        crop.state.phenology.growing_days .= Int32(20)
        crop.state.phenology.senescence .= true
        crop.state.canopy.lai .= T(0.1)
        crop.state.canopy.lai_npp_deficit .= T(0.3)
        crop.state.carbon.biomass .= T(10)
        crop.state.carbon.leaf .= T(1)
        crop.state.carbon.root .= T(2)
        crop.state.carbon.pool .= T(7)
        crop.state.nitrogen.sufficiency .= one(T)
        crop.state.water.sufficiency .= one(T)

        carbon_allocation!(cft1, state)

        # LPJmL retains potential phenological LAI and clips only actual LAI.
        @test only(crop.state.canopy.lai) == T(0.1)
        @test only(crop.state.canopy.lai_npp_deficit) == T(0.3)
        @test only(crop.auxiliary.canopy.actual_lai) == zero(T)

        pet.par .= T(20)
        apar_crop!(cft1, state, pet)
        @test only(crop.auxiliary.canopy.fpar) == zero(T)
        @test only(crop.auxiliary.canopy.apar) == zero(T)

        pet.eeq .= T(2)
        rain = fill(T(5), 1)
        interception!(state, cft1, pet.eeq, rain)
        @test only(crop.auxiliary.canopy.canopy_wet) == zero(T)
        @test only(crop.fluxes.water.interception) == zero(T)
    end
end

@testset "LPJmL canopy-growth water multiplier" begin
    for T in (Float32, Float64)
        crop = init_crop(T, 1, identity)
        state = test_model_state(crop)
        crop.state.phenology.is_growing .= Int32(1)
        crop.state.phenology.senescence .= false
        crop.state.canopy.lai .= T(1)
        # `lai_previous_potential` is LPJmL's `lai000`; establish a
        # self-consistent state for this isolated growth-step test.
        crop.state.canopy.lai_previous_potential .= T(1)
        crop.auxiliary.canopy.flaimax .= T(0.8)
        crop.state.water.sufficiency .= T(0.6)
        crop.state.nitrogen.sufficiency .= T(0.9)

        lai_crop!(state, cft1)

        potential_lai = T(0.8) * T(cft1.laimax)
        expected = T(1) + (potential_lai - T(1)) * T(0.6)
        @test only(crop.state.canopy.lai) ≈ expected
    end
end
