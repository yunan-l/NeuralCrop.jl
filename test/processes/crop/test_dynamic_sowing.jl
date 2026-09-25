using NeuralCrop
using Test

function dynamic_sowing_fixture(cft; T = Float32, device = identity)
    crop = init_crop(T, 1, device)
    climbuf = init_climbuf(T, 1, device)
    crop.auxiliary.calendar.prescribed_sowing_date .= Int32(100)
    crop.auxiliary.phenology.winter_type .= false
    climbuf.mprec20 .= T(20)
    climbuf.mpet20 .= T(20)
    climbuf.mtemp20 .= T(15)
    return crop, climbuf
end

@testset "LPJmL climate calendars cover four target crops" begin
    # Wheat resolves a cooling-season climate candidate for the winter variety.
    wheat, wheat_climate = dynamic_sowing_fixture(cft1)
    wheat.auxiliary.phenology.winter_type .= true
    wheat.auxiliary.calendar.prescribed_sowing_date .= Int32(244)
    wheat_climate.mtemp20 .= reshape(Float32[20, 20, 18, 16, 15, 14, 13, 12, 10, 7, 5, 4], 12, 1)
    update_dynamic_sowing_calendar!(wheat_climate, cft1, wheat.auxiliary.phenology.winter_type)
    @test only(wheat_climate.climate_sowing_day) > 0
    @test only(wheat_climate.reference_sowing_day) == only(wheat_climate.climate_sowing_day)

    # Maize chooses a temperature candidate in a cold seasonal climate.
    maize, maize_climate = dynamic_sowing_fixture(cft3)
    maize_climate.mtemp20 .= reshape(Float32[2, 4, 15, 18, 20, 22, 23, 22, 18, 12, 7, 3], 12, 1)
    update_dynamic_sowing_calendar!(maize_climate, cft3, maize.auxiliary.phenology.winter_type)
    @test only(maize_climate.climate_sowing_day) > 0

    # Rice and soybean use the four-month precipitation window, which starts
    # in May for this climate.
    rainy_months = Float32[10, 10, 10, 10, 100, 100, 100, 100, 10, 10, 10, 10]
    for cft in (cft2, cft9)
        crop, climate = dynamic_sowing_fixture(cft)
        climate.mprec20 .= reshape(rainy_months, 12, 1)
        update_dynamic_sowing_calendar!(climate, cft, crop.auxiliary.phenology.winter_type)
        @test only(climate.climate_sowing_day) == 121
    end
end

@testset "Dynamic sowing applies bounded climate anomalies to prescribed dates" begin
    crop, climate = dynamic_sowing_fixture(cft1)
    crop.auxiliary.calendar.prescribed_sowing_date .= Int32(305)
    climate.reference_sowing_day .= Int32(280)
    climate.climate_sowing_day .= Int32(270)
    dynamic_sowing_date!(crop, climate, cft1, 295)
    @test only(crop.auxiliary.calendar.sowing_date) == 295

    # A large circular climate-calendar jump is capped at 30 days rather than
    # moving an autumn prescribed date to January.
    climate.climate_sowing_day .= Int32(100)
    dynamic_sowing_date!(crop, climate, cft1, 275)
    @test only(crop.auxiliary.calendar.sowing_date) == 275

    # The calendar bootstrap has no reference candidate and therefore uses the
    # prescribed spatial date exactly.
    climate.climate_sowing_day .= Int32(0)
    climate.reference_sowing_day .= Int32(0)
    dynamic_sowing_date!(crop, climate, cft1, 305)
    @test only(crop.auxiliary.calendar.sowing_date) == 305
end

@testset "Dynamic precipitation uses the four-month rainfall window" begin
    crop, climate = dynamic_sowing_fixture(cft3)
    climate.mprec20 .= reshape(
        Float32[50, 50, 50, 50, 60, 60, 60, 60, 10, 10, 10, 10], 12, 1,
    )
    climate.mpet20 .= reshape(
        Float32[20, 20, 20, 20, 10, 10, 10, 10, 20, 20, 20, 20], 12, 1,
    )
    # LPJmL selects the first four-month window with the largest sum of
    # monthly P/PET, not the largest raw precipitation total.
    @test NeuralCrop._dynamic_sowing_month_from_precipitation(
        climate.mprec20, climate.mpet20, 1, 1,
    ) == 5
end

@testset "LPJmL dynamic calendar uses persistent P/PET sowing months" begin
    crop, climate = dynamic_sowing_fixture(cft2)
    # Raw precipitation alone would select January, but LPJmL selects the
    # largest four-month sum of P/PET, which begins in May here.
    climate.mprec20 .= reshape(
        Float32[200, 200, 200, 200, 80, 80, 80, 80, 10, 10, 10, 10], 12, 1,
    )
    climate.mpet20 .= reshape(
        Float32[100, 100, 100, 100, 10, 10, 10, 10, 100, 100, 100, 100], 12, 1,
    )
    climate.mtemp20 .= 20f0
    update_dynamic_sowing_calendar!(climate, cft2, crop.auxiliary.phenology.winter_type)
    @test only(climate.seasonality_type) == NeuralCrop._DYNAMIC_PRECIPITATION
    @test only(climate.sowing_month) == 5

    # The first valid calendar becomes the reference and the current climate
    # candidate therefore maps back to the prescribed spatial date.
    climate.mprec20 .= 1f0
    climate.mpet20 .= 1f0
    dynamic_sowing_date!(crop, climate, cft2, 100)
    @test only(crop.auxiliary.calendar.sowing_date) == 100
end

@testset "Dynamic sowing retains annual precipitation climatology" begin
    climate = init_climbuf(Float32, 1, identity)
    annual_temperature = fill(15f0, 365, 1)
    annual_precipitation = fill(1f0, 365, 1)
    annual_pet = fill(2f0, 365, 1)
    annual_climbuf!(
        annual_temperature, climate, cft1;
        daily_prec = annual_precipitation, daily_pet = annual_pet,
    )
    @test vec(climate.mprec[:, 1]) == Float32[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    @test climate.mprec20 == climate.mprec
    @test vec(climate.mpet[:, 1]) == Float32[62, 56, 62, 60, 62, 60, 62, 62, 60, 62, 60, 62]
    @test climate.mpet20 == climate.mpet
end

@testset "Sowing-mode configuration remains explicit" begin
    initial = lifecycle_initial_data(Float32)
    prescribed = initialize_simulation(cft1, initial; days = 1, diagnostics = false)
    dynamic = initialize_simulation(
        cft1, initial; days = 1, diagnostics = false, sowing_mode = :dynamic_sdate,
    )
    @test prescribed.config.sowing_mode === :prescribed_sdate
    @test dynamic.config.sowing_mode === :dynamic_sdate
    @test_throws ArgumentError initialize_simulation(
        cft1, initial; days = 1, diagnostics = false, sowing_mode = :unknown,
    )
    @test_throws ArgumentError initialize_simulation(
        cft4, initial; days = 1, diagnostics = false, sowing_mode = :dynamic_sdate,
    )
end

@testset "Prescribed sowing remains independent of dynamic sowing" begin
    simulation = initialize_simulation(cft1, lifecycle_initial_data(Float32); days = 2, diagnostics = false)
    climate = lifecycle_climate(Float32, 2)
    management = (
        sdate = Int32[1], phu = Float32[-700], manure = Float32[0],
        fertilizer = Float32[0], residuefrac = Float32[0.67],
    )
    run_simulation!(simulation, climate; spinup = false, management)
    @test simulation.config.sowing_mode === :prescribed_sdate
    @test event_days(simulation.output.calendar.sowing_event) == [1]
end

@testset "Daily driver consumes a dynamic sowing trigger" begin
    days = 244
    simulation = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days, diagnostics = false, sowing_mode = :dynamic_sdate,
    )
    climate = lifecycle_climate(Float32, days)
    temperature = copy(climate.temp)
    temperature[end - 1, 1] = 13f0
    temperature[end, 1] = 11f0
    climate = merge(climate, (temp = temperature,))
    management = (
        sdate = Int32[244], phu = Float32[-700], manure = Float32[0],
        fertilizer = Float32[0], residuefrac = Float32[0.67],
    )
    run_simulation!(simulation, climate; spinup = false, management)
    @test event_days(simulation.output.calendar.sowing_event) == [244]
end

@testset "Dynamic calendar closes only at a completed climate-year boundary" begin
    days = 366
    simulation = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days, diagnostics = false, sowing_mode = :dynamic_sdate,
    )
    management = (
        sdate = Int32[244], phu = Float32[-700], manure = Float32[0],
        fertilizer = Float32[0], residuefrac = Float32[0.67],
    )
    run_simulation!(simulation, lifecycle_climate(Float32, days); spinup = false, management)
    @test all(Array(simulation.climbuf.mpet20) .> -9998f0)
    @test only(Array(simulation.climbuf.seasonality_type)) !=
          NeuralCrop._DYNAMIC_CLIMATOLOGY_UNAVAILABLE
    @test 1 <= only(Array(simulation.climbuf.sowing_month)) <= 12
    @test only(Array(simulation.climbuf.climate_sowing_day)) > 0
    @test only(Array(simulation.climbuf.reference_sowing_day)) > 0
end

@testset "Dynamic calendar uses initialized winter type without annual management" begin
    simulation = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days = 366, diagnostics = false, sowing_mode = :dynamic_sdate,
    )
    run_simulation!(simulation, lifecycle_climate(Float32, 366); spinup = false)
    @test only(Array(simulation.climbuf.seasonality_type)) !=
          NeuralCrop._DYNAMIC_CLIMATOLOGY_UNAVAILABLE
end
