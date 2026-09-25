using NeuralCrop
using Test

function lifecycle_simulation(days)
    return initialize_simulation(
        cft1,
        lifecycle_initial_data(Float32);
        T = Float32,
        days = days,
        diagnostics = false,
        fertilizer = :yes,
    )
end

@testset "Annual climate blocks preserve crop events and soil memory" begin
    days = 730
    continuous = lifecycle_simulation(days)
    chunked = lifecycle_simulation(days)
    full_climate = lifecycle_climate(Float32, days)
    first_year = NamedTuple(
        field => getproperty(full_climate, field) isa AbstractMatrix ?
            getproperty(full_climate, field)[1:365, :] : getproperty(full_climate, field)
        for field in propertynames(full_climate)
    )
    second_year = NamedTuple(
        field => getproperty(full_climate, field) isa AbstractMatrix ?
            getproperty(full_climate, field)[366:730, :] : getproperty(full_climate, field)
        for field in propertynames(full_climate)
    )

    run_simulation!(continuous, full_climate; spinup = false)
    run_simulation!(chunked, [first_year, second_year]; spinup = false)

    @test chunked.simulated_days == days
    for field in (:sowing_event, :harvest_event, :harvesting_mask)
        @test getproperty(chunked.output.calendar, field) ==
            getproperty(continuous.output.calendar, field)
    end
    for field in (:fphu, :biomass, :lai, :npp)
        @test getproperty(chunked.output.crop, field) ≈
            getproperty(continuous.output.crop, field)
    end
    for (chunked_state, continuous_state) in (
        (chunked.state.prognostic.soil.water.storage, continuous.state.prognostic.soil.water.storage),
        (chunked.state.prognostic.soil.water.ice_storage, continuous.state.prognostic.soil.water.ice_storage),
        (chunked.state.prognostic.soil.thermal.temperature, continuous.state.prognostic.soil.thermal.temperature),
        (chunked.state.prognostic.soil.carbon.litter, continuous.state.prognostic.soil.carbon.litter),
        (chunked.state.prognostic.soil.carbon.fast, continuous.state.prognostic.soil.carbon.fast),
        (chunked.state.prognostic.soil.carbon.slow, continuous.state.prognostic.soil.carbon.slow),
        (chunked.state.prognostic.soil.nitrogen.litter, continuous.state.prognostic.soil.nitrogen.litter),
        (chunked.state.prognostic.soil.nitrogen.fast, continuous.state.prognostic.soil.nitrogen.fast),
        (chunked.state.prognostic.soil.nitrogen.slow, continuous.state.prognostic.soil.nitrogen.slow),
        (chunked.state.prognostic.soil.nitrogen.nitrate, continuous.state.prognostic.soil.nitrogen.nitrate),
        (chunked.state.prognostic.soil.nitrogen.ammonium, continuous.state.prognostic.soil.nitrogen.ammonium),
        (chunked.state.prognostic.soil.management.tillage_density_factor,
         continuous.state.prognostic.soil.management.tillage_density_factor),
    )
        @test chunked_state ≈ continuous_state
    end
    @test event_days(chunked.output.calendar.sowing_event) == [100, 465]
    @test event_days(chunked.output.calendar.harvest_event) == [137, 502]
end
