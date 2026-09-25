using NeuralCrop
using Test

@testset "Fertilizer process has no implicit unlimited-N source" begin
    crop = init_crop(1, identity)
    managed_land = init_managed_land(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil; managed_land)
    no3_before = copy(soil.nitrogen.nitrate)
    nh4_before = copy(soil.nitrogen.ammonium)

    fertilizer!(state, managed_land, state, 1)

    @test soil.nitrogen.nitrate == no3_before
    @test soil.nitrogen.ammonium == nh4_before
end

@testset "Prescribed fertilizer is split and conserved" begin
    crop = init_crop(1, identity)
    managed_land = init_managed_land(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil; managed_land)
    crop.auxiliary.calendar.sowing_date .= 1
    managed_land.fertilizer .= 10.0f0
    no3_before = sum(soil.nitrogen.nitrate)
    nh4_before = sum(soil.nitrogen.ammonium)

    fertilizer!(state, managed_land, state, 1)
    @test crop.state.nitrogen.pending_fertilizer[1] ≈ 8.0f0
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] ≈ 2.0f0

    crop.auxiliary.phenology.phu .= 1.0f0
    crop.state.phenology.husum .= 0.3f0
    fertilizer!(state, managed_land, state, 2)
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] ≈ 8.0f0

    no3_input = sum(soil.nitrogen.nitrate) - no3_before
    nh4_input = sum(soil.nitrogen.ammonium) - nh4_before
    @test no3_input ≈ 5.0f0 atol = 1.0f-6
    @test nh4_input ≈ 5.0f0 atol = 1.0f-6
    @test no3_input + nh4_input ≈ 10.0f0 atol = 1.0f-6
    @test crop.state.nitrogen.pending_fertilizer[1] == 0.0f0
end

@testset "Prescribed manure is split and conserved" begin
    crop = init_crop(1, identity)
    managed_land = init_managed_land(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil; managed_land)
    crop.auxiliary.calendar.sowing_date .= 1
    managed_land.manure .= 10.0f0
    mineral_before = sum(soil.nitrogen.ammonium)
    organic_before = sum(soil.nitrogen.litter)

    fertilizer!(state, managed_land, state, 1; manure = true)
    @test crop.state.nitrogen.pending_manure[1] ≈ 8.0f0
    @test crop.fluxes.nitrogen.prescribed_manure_input[1] ≈ 2.0f0

    crop.auxiliary.phenology.phu .= 1.0f0
    crop.state.phenology.husum .= 0.3f0
    fertilizer!(state, managed_land, state, 2; manure = true)
    @test crop.fluxes.nitrogen.prescribed_manure_input[1] ≈ 8.0f0

    mineral_input = sum(soil.nitrogen.ammonium) - mineral_before
    organic_input = sum(soil.nitrogen.litter) - organic_before
    @test mineral_input + organic_input ≈ 10.0f0 atol = 2.0f-6
    @test crop.state.nitrogen.pending_manure[1] == 0.0f0
end

@testset "Prescribed mineral fertilizer can be disabled independently" begin
    crop = init_crop(1, identity)
    managed_land = init_managed_land(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil; managed_land)
    crop.auxiliary.calendar.sowing_date .= 1
    managed_land.fertilizer .= 10.0f0
    no3_before = copy(soil.nitrogen.nitrate)
    nh4_before = copy(soil.nitrogen.ammonium)

    fertilizer!(state, managed_land, state, 1; fertilizer = false)

    @test soil.nitrogen.nitrate == no3_before
    @test soil.nitrogen.ammonium == nh4_before
    @test crop.state.nitrogen.pending_fertilizer[1] == 0.0f0
end

@testset "Manure remains independent of mineral fertilizer mode" begin
    crop = init_crop(1, identity)
    managed_land = init_managed_land(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil; managed_land)
    crop.auxiliary.calendar.sowing_date .= 1
    managed_land.fertilizer .= 10.0f0
    managed_land.manure .= 6.0f0
    no3_before = sum(soil.nitrogen.nitrate)
    litter_before = sum(soil.nitrogen.litter)

    fertilizer!(state, managed_land, state, 1; fertilizer = false, manure = true)

    @test sum(soil.nitrogen.nitrate) == no3_before
    @test sum(soil.nitrogen.litter) > litter_before
    @test crop.fluxes.nitrogen.prescribed_fertilizer_input[1] == 0.0f0
    @test crop.fluxes.nitrogen.prescribed_manure_input[1] > 0.0f0
end

@testset "LPJmL-style fertilizer modes" begin
    @test fertilizer_mode(:no) === :no
    @test fertilizer_mode("yes") === :yes
    @test fertilizer_mode(:auto) === :auto
    @test_throws ArgumentError fertilizer_mode(:invalid)
end
