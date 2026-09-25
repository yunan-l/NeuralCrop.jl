using NeuralCropData
using Dates
using Test

include("fixtures/fixture_data.jl")
using .FixtureData

include("test_hwsd.jl")
include("test_prepare_global_cft_subset.jl")
include("test_climate_deposition.jl")

@testset "NeuralCropData" begin
    mktempdir() do directory
        paths = FixtureData.write_fixture(directory)
        catalog = load_catalog(paths.catalog_path)

        # Keep legacy LPJmL source dimensions at the file-I/O boundary only.
        @test NeuralCropData._canonical_dimension(:pft) === :cft

        @test dataset(catalog, :grid).path == paths.grid_path
        @test dataset(catalog, :landuse).management_bands.rainfed == Int32[1, 2]
        @test dataset(catalog, :landuse).management_bands.irrigated == Int32[3, 4]
        @test cft_index(catalog.cfts, 20) == 2
        @test cft_index(catalog.cfts, "crop_b") == 2
        @test cft_name(catalog.cfts, 10) == "crop_a"

        grid = read_grid(dataset(catalog, :grid))
        @test grid.cell_ids == Int32[0, 1, 2, 3]
        @test grid.longitude_indices == Int32[1, 3, 1, 2]
        @test grid.latitude_indices == Int32[2, 1, 1, 2]

        spatial = Float32[12 10; -1 13; 11 -1]
        compact = compact_spatial(spatial, grid, 1, 2)
        @test compact == Float32[10, 11, 12, 13]
        @test expand_to_grid(compact, grid; fill_value = -1.0f0) == spatial

        subset = select_cells(grid, [1, 3])
        @test subset.cell_ids == Int32[0, 2]
        @test compact_spatial(spatial, grid, 1, 2; selection = subset) == Float32[10, 12]

        climate_reader = climate_blocks(
            catalog, grid; co2_path = paths.co2_path, block_days = 2,
        )
        @test length(climate_reader) == 3
        @test climate_days(climate_reader) == 6
        climate = [block for block in climate_reader]
        @test size.(getfield.(climate, :temperature)) == fill((2, 4), 3)
        @test climate[1].co2 == Float32[369.5, 369.5]
        @test climate[2].co2 == Float32[369.5, 371.0]
        @test climate[1].provenance.calendar == "unspecified"
        @test climate[1].provenance.model_units == (
            temp = "degC", prec = "mm/day", lwnet = "W/m2",
            swdown = "W/m2", co2 = "ppm",
            no3_deposition = nothing, nh4_deposition = nothing,
        )
        @test isnothing(climate[1].no3_deposition)
        @test isnothing(climate[1].nh4_deposition)
        @test climate_forcing(climate[1]).co2_daily
        @test climate_forcing(climate[1]).backend_neutral
        @test !hasproperty(climate_forcing(climate[1]), :provenance)
        forcings = climate_forcings(climate_reader)
        @test forcings isa AbstractVector{NamedTuple}
        @test length(forcings) == length(climate_reader)
        @test forcings[2].temp == climate[2].temperature
        @test forcings[2].co2 == climate[2].co2
        prefetched = prefetch_climate_forcings(climate_reader)
        @test prefetched[1] == forcings[1]
        @test prefetched[2] == forcings[2]
        close(prefetched)
        @test_throws ArgumentError prefetched[3]

        eager_temp = read_compact_variable(
            dataset(catalog, :temp), grid; order = (:time, :cell), T = Float32,
        )
        @test reduce(vcat, getfield.(climate, :temperature)) == eager_temp.values
        @test reduce(vcat, getfield.(climate, :precipitation)) ==
            read_compact_variable(
                dataset(catalog, :prec), grid; order = (:time, :cell), T = Float32,
            ).values
        year_2001 = climate_blocks(
            catalog, grid; co2_path = paths.co2_path, start_year = 2001,
            end_year = 2001, block_days = 31,
        )
        @test length(year_2001) == 1
        @test read_climate_block(year_2001, 1).co2 == fill(371.0f0, 3)

        landuse = read_management(catalog, :landuse, grid, 20; years = 2000:2001)
        @test size(landuse.values) == (2, 4)
        @test landuse.values == Float32[0 0.2 0 0.3; 0.4 0 0 0.1]
        irrigated_landuse = read_management(catalog, :landuse, grid, 20; years = 2000:2001, irrigated = true)
        @test irrigated_landuse.values == landuse.values .+ 0.1f0
        extended_landuse = read_management(
            catalog, :landuse, grid, 20; simulation_years = 1998:2003,
        )
        @test extended_landuse.time == collect(1998:2003)
        @test extended_landuse.values == landuse.values[[1, 1, 1, 2, 2, 2], :]
        crop_mask = build_crop_mask(grid, landuse.values)
        @test crop_mask.selection.cell_ids == Int32[0, 1, 3]
        @test crop_mask.fraction == Float32[0 0.2 0.3; 0.4 0 0.1]
        @test crop_mask.active == Bool[0 1 1; 1 0 1]

        rainfed_patches = build_patch_domain(grid, landuse.values, 20)
        irrigated_patches = build_patch_domain(
            grid, irrigated_landuse.values, 20; irrigated = true,
        )
        patches = combine_patch_domains([rainfed_patches, irrigated_patches])
        @test rainfed_patches.cell_ids == Int32[0, 1, 3]
        @test rainfed_patches.landfrac == Float32[0.4, 0.2, 0.3]
        @test patches.patch_ids == Int32.(1:7)
        @test patches.cell_ids == Int32[0, 1, 3, 0, 1, 2, 3]
        @test patches.cft_ids == fill(Int32(20), 7)
        @test patches.irrigated == BitVector([false, false, false, true, true, true, true])
        @test patches.landfrac == Float32[0.4, 0.2, 0.3, 0.5, 0.3, 0.1, 0.4]

        residue = read_management(
            catalog,
            :residue_fraction,
            grid,
            20;
            years = 2000:2001,
            selection = crop_mask.selection,
        )
        @test size(residue.values) == (2, 3)
        @test all(0 .<= residue.values .<= 1)

        single_year = (; (
            name => read_management(
                catalog,
                name,
                grid,
                20;
                years = [2000],
                selection = crop_mask.selection,
                active = crop_mask.active[1:1, :],
            ) for name in (:sowing_date, :phu, :manure, :fertilizer, :residue_fraction)
        )...)
        crop = crop_inputs(; single_year...)
        @test crop.sdate == Int32[0, 100, 100]
        @test crop.phu == Float32[0, 1200, 1200]
        @test crop.fertilizer == Float32[0, 4, 6]
        automatic_crop = crop_inputs(
            sowing_date = single_year.sowing_date,
            phu = single_year.phu,
            residue_fraction = single_year.residue_fraction,
            fertilizer_mode = :auto,
            manure_enabled = false,
        )
        @test automatic_crop.fertilizer == zeros(Float32, 3)
        @test automatic_crop.manure == zeros(Float32, 3)

        transient_active = read_management(
            catalog, :landuse, grid, 20;
            simulation_years = 1999:2002, selection = crop_mask.selection,
        ).values .> 0
        transient = (; (
            name => read_management(
                catalog,
                name,
                grid,
                20;
                simulation_years = 1999:2002,
                selection = crop_mask.selection,
                active = transient_active,
            ) for name in (:sowing_date, :phu, :manure, :fertilizer, :residue_fraction)
        )...)
        schedule = management_schedule(; transient...)
        @test schedule.years == Int32.(1999:2002)
        @test schedule.phu == transient.phu.values
        @test schedule.sdate == Int32.(transient.sowing_date.values)
        @test schedule.fertilizer == transient.fertilizer.values

        soil = read_soil_data(catalog, grid)
        @test soil.soilcode == Int32[1, 6, 9, 14]
        @test soil.ph == Float32[6, 7, 8, 9]
        @test soil.sand == Float32[0.22, 0.58, 0.58, 0.99]
        @test size(soil.saturation) == (5, 4)
        @test soil_properties(soil).soilph === soil.ph

        selected_soil = read_soil_data(
            catalog, grid; selection = crop_mask.selection,
        )
        baseline_selection = CellSelection(1:10, 0:9)
        baseline_codes = Int32[6, 7, 9, 9, 9, 9, 9, 9, 9, 9]
        baseline_ph = Float32[6.5, 7, 7, 7, 5.5, 5.5, 5.5, 5.5, 7, 5.5]
        baseline = soil_data_from_values(baseline_codes, baseline_ph, baseline_selection)
        @test baseline.sand == Float32[0.58, 0.43, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58]
        @test baseline.saturation[:, 1] == fill(0.404f0, 5)

        @test_throws ArgumentError validate_management(:landuse, Float32[1.1 0.0])
        @test_throws ArgumentError validate_management(:sowing_date, Float32[0 100]; active = Bool[1 1])
        @test validate_management(:phu, Float32[-1200 1200]) == Float32[-1200 1200]
        @test_throws ArgumentError validate_management(:phu, Float32[0 1200]; active = Bool[1 1])
        @test_throws DimensionMismatch read_management(
            DatasetSpec(paths.management_path, "landfrac"; units = "%"),
            :landuse,
            grid,
            catalog.cfts,
            20,
        )
    end

    example = load_catalog(joinpath(@__DIR__, "..", "config", "catalog.example.toml"))
    @test example.cfts.ids == Int32.(1:12)
    @test cft_name(example.cfts, 4) == "tropical cereals"
    @test cft_name(example.cfts, 9) == "oil crops soybean"
    @test dataset(example, :landuse).management_bands.irrigated == Int32.(17:28)
    @test dataset(example, :sowing_date).management_bands.irrigated == Int32.(13:24)
    @test dataset(example, :residue_fraction).management_bands.irrigated == Int32.(1:12)

    leap_dates = collect(Date(2000, 1, 1):Day(1):Date(2001, 12, 31))
    leap_indices, calendar = NeuralCropData._normalize_calendar_indices(
        leap_dates, collect(eachindex(leap_dates)), "standard",
    )
    @test length(leap_indices) == 730
    @test calendar == "noleap"
    @test all(leap_dates[index] != Date(2000, 2, 29) for index in leap_indices)
    @test_throws ArgumentError NeuralCropData._normalize_calendar_indices(
        leap_dates, collect(1:10), "standard",
    )
    @test_throws ArgumentError NeuralCropData._normalize_calendar_indices(
        leap_dates, collect(eachindex(leap_dates)), "360_day",
    )

    kelvin, kelvin_units = NeuralCropData._normalize_climate_units(
        :temp, Float32[273.15;;], "K",
    )
    celsius, celsius_units = NeuralCropData._normalize_climate_units(
        :temp, Float32[15;;], "celsius",
    )
    precipitation, precipitation_units = NeuralCropData._normalize_climate_units(
        :prec, Float32[1.0f0 / 86400;;], "kg m-2 s-1",
    )
    daily_precipitation, daily_precipitation_units =
        NeuralCropData._normalize_climate_units(
            :prec, Float32[2;;], "kg/m2/day",
        )
    @test kelvin ≈ Float32[0;;] atol = 1.0f-5
    @test kelvin_units == "degC"
    @test celsius == Float32[15;;]
    @test celsius_units == "degC"
    @test precipitation ≈ Float32[1;;] atol = 1.0f-6
    @test precipitation_units == "mm/day"
    @test daily_precipitation == Float32[2;;]
    @test daily_precipitation_units == "mm/day"
    @test_throws ArgumentError NeuralCropData._normalize_climate_units(
        :temp, Float32[273.15;;], "fahrenheit",
    )
end
