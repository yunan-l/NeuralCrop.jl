using JLD2
using NCDatasets

function stream_fixture_climate(climate, indices)
    return ClimateDataLoader(climate, indices, identity; T = Float32)
end

function split_stream_climate(climate, block_days, days = size(climate.temp, 1))
    return [(
        temp = copy(climate.temp[rows, :]),
        prec = copy(climate.prec[rows, :]),
        sw = copy(climate.sw[rows, :]),
        lw = copy(climate.lw[rows, :]),
        co2 = ndims(climate.co2) == 1 ? fill(climate.co2[1], length(rows)) :
            copy(climate.co2[rows, :]),
        co2_daily = true,
    ) for first_day in 1:block_days:days
      for rows in (first_day:min(days, first_day + block_days - 1),)]
end

@testset "streamed ten-cell daily output equals in-memory output" begin
    example_dir = joinpath(@__DIR__, "..", "..", "examples")
    initial = JLD2.load(joinpath(example_dir, "initial_wheat.jld2"), "initial_data")
    raw_climate = JLD2.load(joinpath(example_dir, "climate_2000_2009.jld2"), "climate")
    cells = length(initial.latitude)
    indices = collect(1:cells)
    days = 365
    climate = stream_fixture_climate(raw_climate, indices)

    baseline = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days, diagnostics = false, fertilizer = :yes,
    )
    run_simulation!(baseline, climate; end_day = days, spinup = false)

    chunks = OutputChunk[]
    stream = OutputStream(
        [
            OutputVariable(:crop, :gpp),
            OutputVariable(:crop, :npp),
            OutputVariable(:crop, :lai),
            OutputVariable(:soil, :ecosystem_respiration),
            OutputVariable(:soil, :heterotrophic_respiration),
            OutputVariable(:soil, :evapotranspiration),
            OutputVariable(:crop, :yield),
        ];
        frequency = :daily,
        writer = chunk -> push!(chunks, chunk),
        cell_ids = indices,
    )
    streamed = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days, diagnostics = false, fertilizer = :yes,
    )
    blocks = split_stream_climate(climate, 73, days)
    run_simulation!(streamed, blocks; spinup = false, output_stream = stream)

    daily_chunks = filter(chunk -> chunk.frequency === :daily, chunks)
    annual_chunks = filter(chunk -> chunk.frequency === :annual, chunks)
    @test reduce(vcat, [chunk.values[:crop_gpp] for chunk in daily_chunks]) ==
        Array(baseline.output.crop.gpp)
    @test reduce(vcat, [chunk.values[:crop_npp] for chunk in daily_chunks]) ==
        Array(baseline.output.crop.npp)
    @test reduce(vcat, [chunk.values[:crop_lai] for chunk in daily_chunks]) ==
        Array(baseline.output.crop.lai)
    @test reduce(vcat, [chunk.values[:soil_ecosystem_respiration] for chunk in daily_chunks]) ==
        Array(baseline.output.soil.ecosystem_respiration)
    @test reduce(vcat, [chunk.values[:soil_heterotrophic_respiration] for chunk in daily_chunks]) ==
        Array(baseline.output.soil.heterotrophic_respiration)
    @test reduce(vcat, [chunk.values[:soil_evapotranspiration] for chunk in daily_chunks]) ==
        Array(baseline.output.soil.evapotranspiration)
    @test reduce(vcat, [chunk.values[:crop_yield] for chunk in annual_chunks]) ==
        Array(baseline.output.crop.yield)
    @test all(length(chunk.time) <= 73 for chunk in daily_chunks)
    @test all(Set(keys(chunk.values)) == Set((
        :crop_gpp, :crop_npp, :crop_lai,
        :soil_ecosystem_respiration, :soil_heterotrophic_respiration,
        :soil_evapotranspiration,
    ))
              for chunk in daily_chunks)
    @test NeuralCrop._output_timeseries_empty(streamed.output)

    diagnostic_simulation = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days = 1, diagnostics = true, fertilizer = :no,
    )
    diagnostic_stream = OutputStream(
        [OutputVariable(:crop, :npp)]; cell_ids = indices,
    )
    @test_throws ArgumentError run_simulation!(
        diagnostic_simulation, [blocks[1]]; spinup = false,
        output_stream = diagnostic_stream,
    )
end

@testset "monthly and annual stream aggregation" begin
    monthly_chunks = OutputChunk[]
    monthly = OutputStream(
        [
            OutputVariable(:crop, :npp; reduction = :sum),
            OutputVariable(:crop, :lai; reduction = :mean),
        ];
        frequency = :monthly,
        writer = chunk -> push!(monthly_chunks, chunk),
        cell_ids = 1:2,
    )
    output = init_output(Float32, 2, identity)
    NeuralCrop.prepare_output_block!(output, 30, 0)
    output.crop.npp .= 2
    output.crop.lai .= 4
    consume_output!(monthly, output, 1)
    clear_output_timeseries!(output)
    NeuralCrop.prepare_output_block!(output, 30, 0)
    output.crop.npp .= 2
    output.crop.lai .= 4
    consume_output!(monthly, output, 31)
    finish_output_stream!(monthly, 60)
    @test getfield.(monthly_chunks, :time) == [[31], [59], [60]]
    @test monthly_chunks[1].values[:crop_npp] == fill(62.0f0, 1, 2)
    @test monthly_chunks[1].values[:crop_lai] == fill(4.0f0, 1, 2)
    @test monthly_chunks[2].values[:crop_npp] == fill(56.0f0, 1, 2)
    @test monthly_chunks[3].values[:crop_npp] == fill(2.0f0, 1, 2)

    annual_chunks = OutputChunk[]
    annual = OutputStream(
        [
            OutputVariable(:crop, :npp; reduction = :sum),
            OutputVariable(:crop, :yield),
            OutputVariable(:crop, :season_gpp),
            OutputVariable(:crop, :season_lai_days),
            OutputVariable(:crop, :season_length),
            OutputVariable(:crop, :season_water_deficit),
            OutputVariable(:crop, :season_evapotranspiration),
            OutputVariable(:crop, :harvest_aboveground_carbon),
        ];
        frequency = :annual,
        writer = chunk -> push!(annual_chunks, chunk),
        cell_ids = 1:2,
    )
    clear_output_timeseries!(output)
    NeuralCrop.prepare_output_block!(output, 365, 1)
    output.crop.npp .= 1
    output.crop.yield .= reshape(Float32[10, 20], 1, :)
    output.crop.season_gpp .= reshape(Float32[11, 21], 1, :)
    output.crop.season_lai_days .= reshape(Float32[12, 22], 1, :)
    output.crop.season_length .= reshape(Float32[13, 23], 1, :)
    output.crop.season_water_deficit .= reshape(Float32[14, 24], 1, :)
    output.crop.season_evapotranspiration .= reshape(Float32[15, 25], 1, :)
    output.crop.harvest_aboveground_carbon .= reshape(Float32[16, 26], 1, :)
    consume_output!(annual, output, 1)
    finish_output_stream!(annual, 365)
    @test length(annual_chunks) == 1
    @test annual_chunks[1].time == [365]
    @test annual_chunks[1].values[:crop_npp] == fill(365.0f0, 1, 2)
    @test annual_chunks[1].values[:crop_yield] == reshape(Float32[10, 20], 1, :)
    @test annual_chunks[1].values[:crop_season_gpp] == reshape(Float32[11, 21], 1, :)
    @test annual_chunks[1].values[:crop_season_lai_days] == reshape(Float32[12, 22], 1, :)
    @test annual_chunks[1].values[:crop_season_length] == reshape(Float32[13, 23], 1, :)
    @test annual_chunks[1].values[:crop_season_water_deficit] == reshape(Float32[14, 24], 1, :)
    @test annual_chunks[1].values[:crop_season_evapotranspiration] == reshape(Float32[15, 25], 1, :)
    @test annual_chunks[1].values[:crop_harvest_aboveground_carbon] ==
        reshape(Float32[16, 26], 1, :)
end

@testset "stream block writers" begin
    chunk = OutputChunk(
        1, :daily, [1, 2], Int32[7, 9],
        Dict{Symbol, Any}(:crop_npp => Float32[1 2; 3 4]),
    )
    mktempdir() do directory
        jld_path = JLD2BlockWriter(directory; prefix = "jld")(chunk)
        @test isfile(jld_path)
        @test JLD2.load(jld_path, "values")[:crop_npp] == chunk.values[:crop_npp]

        nc_path = NetCDFBlockWriter(directory; prefix = "nc")(chunk)
        @test isfile(nc_path)
        NCDataset(nc_path, "r") do dataset
            @test dataset["time"][:] == Int32[1, 2]
            @test dataset["cell_id"][:] == Int32[7, 9]
            @test dataset["crop_npp"][:, :] == chunk.values[:crop_npp]
            @test dataset["crop_npp"].attrib["units"] == "gC m-2 day-1"
            @test dataset["crop_npp"].attrib["long_name"] == "Net primary production"
        end
    end
end


@testset "stream output reuses equal-sized block buffers" begin
    output = init_output(Float32, 2, identity)
    NeuralCrop.prepare_output_block!(output, 31, 0; reuse = true)
    gpp = output.crop.gpp
    output.crop.gpp .= 3
    rows = NeuralCrop.prepare_output_block!(output, 31, 0; reuse = true)
    @test rows.first_daily_row == 1
    @test output.crop.gpp === gpp
    @test all(iszero, output.crop.gpp)
    NeuralCrop.prepare_output_block!(output, 30, 0; reuse = true)
    @test size(output.crop.gpp) == (30, 2)
end

@testset "stream output can begin on a later simulation year" begin
    chunks = OutputChunk[]
    stream = OutputStream(
        [OutputVariable(:crop, :gpp; reduction = :sum), OutputVariable(:crop, :yield)];
        frequency = :annual,
        writer = chunk -> push!(chunks, chunk),
        cell_ids = 1:2,
        first_output_day = 366,
    )
    output = init_output(Float32, 2, identity)
    NeuralCrop.prepare_output_block!(output, 730, 2)
    output.crop.gpp .= 1
    output.crop.yield .= Float32[10 20; 30 40]

    consume_output!(stream, output, 1; rows = 730)
    finish_output_stream!(stream, 730)

    @test length(chunks) == 1
    @test chunks[1].time == [730]
    @test chunks[1].values[:crop_gpp] == fill(365.0f0, 1, 2)
    @test chunks[1].values[:crop_yield] == reshape(Float32[30, 40], 1, :)
    @test_throws ArgumentError OutputStream(
        [OutputVariable(:crop, :yield)]; cell_ids = 1:2, first_output_day = 365,
    )
end
