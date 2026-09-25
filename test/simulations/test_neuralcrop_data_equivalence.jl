using NeuralCrop
using JLD2
using Test

if !isdefined(@__MODULE__, :NeuralCropDataFixture)
    @eval module NeuralCropDataFixture
        include(joinpath(@__DIR__, "..", "..", "lib", "NeuralCropData", "src", "NeuralCropData.jl"))
    end
end
const FixtureNeuralCropData = NeuralCropDataFixture.NeuralCropData

function fixture_array_snapshot(value, path = "")
    arrays = Dict{String, Any}()
    if value isa AbstractArray
        arrays[path] = Array(value)
        return arrays
    end
    for name in fieldnames(typeof(value))
        child = getfield(value, name)
        child_path = isempty(path) ? string(name) : string(path, ".", name)
        merge!(arrays, fixture_array_snapshot(child, child_path))
    end
    return arrays
end

function test_fixture_arrays_equal(new_value, old_value)
    new_arrays = fixture_array_snapshot(new_value)
    old_arrays = fixture_array_snapshot(old_value)
    @test keys(new_arrays) == keys(old_arrays)
    for path in keys(old_arrays)
        @test new_arrays[path] == old_arrays[path]
    end
end

@testset "NeuralCropData ten-cell fixture matches the direct input contract" begin
    initial = load(joinpath(@__DIR__, "..", "..", "examples", "initial_wheat.jld2"), "initial_data")
    climate = load(joinpath(@__DIR__, "..", "..", "examples", "climate_2000_2009.jld2"), "climate")
    cells = length(initial.latitude)
    days = 365
    indices = collect(1:cells)

    fixture_simulation = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days, fertilizer = :yes,
    )
    run_simulation!(fixture_simulation, climate; end_day = days, spinup = false)

    selection = FixtureNeuralCropData.CellSelection(indices, Int32.(0:(cells - 1)))
    grid = FixtureNeuralCropData.GridIndex(
        Float32[0],
        copy(initial.latitude),
        reshape(copy(selection.cell_ids), 1, :),
        copy(selection.cell_ids),
        ones(Int32, cells),
        Int32.(indices),
    )
    soil = FixtureNeuralCropData.SoilData(
        selection,
        Int32.(initial.soilparam.soilcode),
        copy(initial.soilparam.soilph),
        copy(initial.soilparam.w_sat),
        copy(initial.soilparam.sand),
        copy(initial.soilparam.silt),
        copy(initial.soilparam.clay),
        copy(initial.soilparam.tdiff_0),
        copy(initial.soilparam.tdiff_15),
        copy(initial.soilparam.soildepth),
        (fixture = "examples/initial_wheat.jld2",),
    )
    state_fields = (:swc, :litc, :fastc, :slowc, :litn, :fastn, :slown)
    initial_state = NamedTuple{state_fields}(map(
        name -> copy(getproperty(initial.initial_state, name)), state_fields,
    ))
    initial_state = merge(initial_state, (
        c_shift_fast = copy(initial.initial_state.c_shift_fast),
        c_shift_slow = copy(initial.initial_state.c_shift_slow),
    ))
    prepared = model_initial_data(grid, soil, initial.crop, initial_state)
    adapted = initialize_simulation(
        cft1, prepared;
        T = Float32, days, fertilizer = :yes,
    )

    block_size = 73
    blocks = FixtureNeuralCropData.ClimateBlock[]
    for first_day in 1:block_size:days
        last_day = min(days, first_day + block_size - 1)
        rows = first_day:last_day
        push!(blocks, FixtureNeuralCropData.ClimateBlock(
            collect(rows),
            copy(climate.temp[rows, :]),
            copy(climate.prec[rows, :]),
            copy(climate.swdown[rows, :]),
            copy(climate.lwnet[rows, :]),
            nothing,
            nothing,
            fill(Float32(climate.co2[1]), length(rows)),
            selection,
            (fixture = "examples/climate_2000_2009.jld2",),
        ))
    end
    run_simulation!(adapted, FixtureNeuralCropData.climate_forcings(blocks); spinup = false)

    @test adapted.simulated_days == fixture_simulation.simulated_days == days
    test_fixture_arrays_equal(adapted.state.prognostic, fixture_simulation.state.prognostic)
    test_fixture_arrays_equal(adapted.state.fluxes, fixture_simulation.state.fluxes)
    test_fixture_arrays_equal(adapted.state.auxiliary, fixture_simulation.state.auxiliary)
    test_fixture_arrays_equal(adapted.state.events, fixture_simulation.state.events)
    test_fixture_arrays_equal(adapted.state.inputs.crop, fixture_simulation.state.inputs.crop)
    test_fixture_arrays_equal(adapted.state.inputs.soil, fixture_simulation.state.inputs.soil)
    test_fixture_arrays_equal(adapted.state.inputs.management, fixture_simulation.state.inputs.management)
    test_fixture_arrays_equal(adapted.state.output, fixture_simulation.state.output)
    test_fixture_arrays_equal(adapted.diagnostics, fixture_simulation.diagnostics)
    @test adapted.daily_weather.temp == fixture_simulation.daily_weather.temp
    @test adapted.daily_weather.prec == fixture_simulation.daily_weather.prec
    @test adapted.daily_weather.swr == fixture_simulation.daily_weather.swr
    @test adapted.daily_weather.lwr == fixture_simulation.daily_weather.lwr
    @test adapted.daily_weather.annual_co2 == fixture_simulation.daily_weather.annual_co2
end
