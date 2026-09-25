using NeuralCrop
using Test

@testset "CPU state initialization" begin
    cell_size = 2
    crop = @inferred init_crop(cell_size, identity)
    managed_land = @inferred init_managed_land(cell_size, identity)
    soil = @inferred init_soil(cell_size, soilparams.soildepth, identity)
    weather = @inferred init_weather(cell_size, identity)
    pet = @inferred init_pet(cell_size, identity)
    climbuf = @inferred init_climbuf(cell_size, identity)
    output = @inferred init_output(cell_size, identity)

    @test propertynames(crop) == (:state, :fluxes, :auxiliary, :events, :workspace)
    @test propertynames(crop.state) == (:phenology, :canopy, :carbon, :nitrogen, :water)
    @test propertynames(crop.fluxes) == (:carbon, :nitrogen, :water)
    @test propertynames(crop.auxiliary) ==
          (:phenology, :calendar, :root, :canopy, :photosynthesis, :stress)
    @test propertynames(crop.events) == (:sowing, :harvest)
    @test size(crop.state.carbon.leaf) == (cell_size,)
    @test isempty(fieldnames(typeof(crop.workspace)))
    @test size(crop.fluxes.water.transpiration_layer) == (5, cell_size)
    @test length(crop.auxiliary.calendar.sowing_date) == cell_size
    @test length(crop.auxiliary.calendar.prescribed_sowing_date) == cell_size
    @test length(managed_land.latitude) == cell_size
    @test length(crop.fluxes.carbon.gross_assimilation) == cell_size
    @test length(pet.daylength) == cell_size
    @test size(climbuf.temp) == (31, cell_size)
    @test eltype(crop.state.canopy.lai) == Float32
    @test eltype(crop.state.phenology.is_growing) == Int32
    @test all(iszero, crop.state.phenology.is_growing)
    @test all(iszero, crop.state.canopy.lai)
    @test all(iszero, crop.fluxes.carbon.gross_assimilation)
    @test crop.state.phenology isa CropPhenology
    @test crop.state.canopy isa CropCanopyState
    @test crop.auxiliary.root isa CropRootAuxiliary
    @test crop.state.carbon isa CropCarbonState
    @test crop.state.nitrogen isa CropNitrogenState
    @test crop.state.water isa CropWaterState
    @test crop.auxiliary.calendar isa CropCalendarAuxiliary
    @test crop.auxiliary.photosynthesis isa CropPhotosynthesisAuxiliary
    @test soil.properties isa SoilProperties
    @test soil.water isa SoilWater
    @test soil.thermal isa SoilThermal
    @test length(soil.thermal.initialized) == cell_size
    @test all(.!soil.thermal.initialized)
    @test size(soil.thermal.enthalpy) == (5, cell_size)
    @test size(soil.thermal.frozen_fraction) == (5, cell_size)
    @test soil.carbon isa SoilCarbon
    @test soil.nitrogen isa SoilNitrogen
    @test soil.decomposition isa SoilDecomposition
    @test size(soil.decomposition.litter_response) == (3, cell_size)
    @test size(soil.decomposition.shift_fast) == (5, cell_size)
    @test size(soil.decomposition.shift_slow) == (5, cell_size)
    @test !hasproperty(soil.carbon, :shift_fast)
    @test !hasproperty(soil.nitrogen, :shift_fast)
    @test size(soil.decomposition.layer_scratch_1) == (5, cell_size)
    @test length(soil.decomposition.surface_scratch_1) == cell_size
    @test soil.management isa SoilManagement
    @test size(soil.management.tillage_density_factor) == (1, cell_size)
    @test all(isone, soil.management.tillage_density_factor)
    @test soil.surface_litter isa SoilSurfaceLitter
    @test soil.snow isa SoilSnow
    @test size(soil.water.storage) == (5, cell_size)
    @test size(soil.water.ice_storage) == (5, cell_size)
    @test size(soil.water.wilting_ice_fraction) == (5, cell_size)
    @test size(soil.water.available_ice_storage) == (5, cell_size)
    @test size(soil.water.free_ice_storage) == (5, cell_size)
    @test size(soil.carbon.litter) == (3, cell_size)
    @test size(soil.carbon.litter_to_fast) == (5, cell_size)
    @test size(soil.carbon.litter_to_slow) == (5, cell_size)
    @test size(soil.nitrogen.litter_to_fast) == (5, cell_size)
    @test size(soil.nitrogen.litter_to_slow) == (5, cell_size)
    @test length(soil.surface_litter.water_storage) == cell_size
    @test length(weather.temp) == cell_size
    @test output.crop isa CropOutput
    @test output.soil isa SoilOutput
    @test output.climate isa ClimateOutput
    @test output.calendar isa CalendarOutput
    @test output.annual isa AnnualOutputAccumulator

    rows = NeuralCrop.prepare_output_block!(output, 3, 1)
    @test rows == (first_daily_row = 1, first_annual_row = 1)
    @test size(output.crop.npp) == (3, cell_size)
    @test size(output.crop.yield) == (1, cell_size)
    @test size(output.calendar.harvest_event) == (3, cell_size)
    @test size(output.calendar.harvest_date) == (1, cell_size)
end

@testset "Selectable floating-point state precision" begin
    cells = 2
    crop = init_crop(Float64, cells, identity)
    managed_land = init_managed_land(Float64, cells, identity)
    soil = init_soil(Float64, cells, soilparams.soildepth, identity)
    weather = init_weather(Float64, cells, identity)
    pet = init_pet(Float64, cells, identity)
    climbuf = init_climbuf(Float64, cells, identity)
    output = init_output(Float64, cells, identity)
    water = init_water_balance(2, cells, identity; T = Float64)
    nitrogen = init_nitrogen_balance(2, cells, identity; T = Float64)
    carbon = init_carbon_balance(2, cells, identity; T = Float64)
    thermal = init_thermal_balance(2, cells, identity; T = Float64)

    @test eltype(crop.state.canopy.lai) == Float64
    @test eltype(crop.state.carbon.leaf) == Float64
    @test eltype(crop.auxiliary.photosynthesis.lambda) == Float64
    @test eltype(managed_land.latitude) == Float64
    @test eltype(soil.properties.layer_depth) == Float64
    @test eltype(soil.water.storage) == Float64
    @test eltype(soil.thermal.enthalpy) == Float64
    @test eltype(soil.carbon.fast) == Float64
    @test eltype(soil.nitrogen.nitrate) == Float64
    @test eltype(soil.management.tillage_density_factor) == Float64
    @test eltype(weather.temp) == Float64
    @test eltype(weather.annual_co2) == Float64
    @test eltype(pet.eeq) == Float64
    @test eltype(climbuf.atemp) == Float64
    @test eltype(output.crop.npp) == Float64
    @test eltype(water.residual) == Float64
    @test eltype(nitrogen.residual) == Float64
    @test eltype(carbon.residual) == Float64
    @test eltype(thermal.energy_residual) == Float64
    @test eltype(crop.state.phenology.is_growing) == Int32
    @test eltype(crop.auxiliary.calendar.sowing_date) == Int32
    @test eltype(crop.auxiliary.calendar.prescribed_sowing_date) == Int32

    cft64 = convert_precision(Float64, cft1)
    parameters64 = ModelParameters(Float64)
    soilparams64 = convert_precision(Float64, soilparams)
    @test cft64 isa CFTParameters{Float64, Int32}
    @test parameters64 isa ModelParameters{Float64}
    @test parameters64.lpjml isa LPJmLParams{Float64}
    @test parameters64.photosynthesis isa PhotoParams{Float64}
    @test parameters64.snow isa SnowParams{Float64}
    @test parameters64.soil_thermal isa SoilThermalParams{Float64}
    @test parameters64.soil_decomposition isa SoilDecompParams{Float64}
    @test soilparams64 isa SoilParams{Float64}

    layers = 5
    initial_data = (
        latitude = Float32[45, 46],
        soilparams = (
            ph = fill(6.5f0, cells),
            w_sat = fill(0.45f0, layers, cells),
            sand = fill(0.4f0, 1, cells),
            clay = fill(0.2f0, 1, cells),
            tdiff_0 = fill(0.2f0, cells),
            tdiff_15 = fill(0.5f0, cells),
            soildepth = Float32[200, 300, 500, 1000, 1000],
        ),
        ModelState = (
            crop = (
                sdate = Int32[100, 101],
                phu = Float32[-543, 620],
                manure = zeros(Float32, cells),
                fertilizer = fill(20.0f0, cells),
                residuefrac = fill(0.5f0, cells),
            ),
            u0 = (
                swc = fill(100.0f0, layers, cells),
                litc = fill(1.0f0, 3, cells),
                fastc = fill(10.0f0, layers, cells),
                slowc = fill(100.0f0, layers, cells),
                litn = fill(0.1f0, 3, cells),
                fastn = fill(1.0f0, layers, cells),
                slown = fill(10.0f0, layers, cells),
            ),
        ),
    )
    _, crop64, pet64, soil64, managed64, weather64, output64 = init_states!(
        cft1, initial_data, cells, identity; T = Float64,
    )
    @test eltype(crop64.state.canopy.lai) == Float64
    @test crop64.auxiliary.phenology.phu == Float64[543, 620]
    @test crop64.auxiliary.phenology.winter_type == Bool[true, false]
    @test eltype(pet64.eeq) == Float64
    @test eltype(soil64.water.storage) == Float64
    @test eltype(managed64.latitude) == Float64
    @test eltype(weather64.temp) == Float64
    @test eltype(output64.crop.npp) == Float64
end

@testset "Soil c_shift initialization strategies" begin
    cells = 2
    soil = init_soil(cells, soilparams.soildepth, identity)
    model_state = (u0 = nothing,)

    NeuralCrop.initialize_soil_c_shift!(soil, model_state, :default)
    expected = Float32[0.55, 0.1125, 0.1125, 0.1125, 0.1125]
    @test soil.decomposition.shift_fast == repeat(expected, 1, cells)
    @test soil.decomposition.shift_slow == repeat(expected, 1, cells)
    @test vec(sum(soil.decomposition.shift_fast; dims = 1)) ≈ ones(Float32, cells)

    restart_fast = repeat(Float32[0.4, 0.25, 0.15, 0.1, 0.1], 1, cells)
    restart_slow = repeat(Float32[0.5, 0.2, 0.15, 0.1, 0.05], 1, cells)
    restart_state = (c_shift_fast = restart_fast, c_shift_slow = restart_slow)
    NeuralCrop.initialize_soil_c_shift!(soil, restart_state, :restart)
    @test soil.decomposition.shift_fast == restart_fast
    @test soil.decomposition.shift_slow == restart_slow

    @test_throws ArgumentError NeuralCrop.initialize_soil_c_shift!(
        soil, model_state, :restart,
    )
    @test_throws ArgumentError NeuralCrop.initialize_soil_c_shift!(
        soil, model_state, :unknown,
    )
end

@testset "Initial data loader makes c_shift optional" begin
    cells = 2
    layers = 5
    u0 = (
        swc = fill(100.0f0, layers, cells),
        litc = fill(1.0f0, 3, cells),
        fastc = fill(10.0f0, layers, cells),
        slowc = fill(100.0f0, layers, cells),
        litn = fill(0.1f0, 3, cells),
        fastn = fill(1.0f0, layers, cells),
        slown = fill(10.0f0, layers, cells),
    )
    initial_without_shift = u0
    data_without_shift = (
        latitude = Float32[45, 46],
        crop = (
            sdate = Int32[100, 101],
            phu = Float32[1200, 1250],
            manure = zeros(Float32, cells),
            fertilizer = fill(50.0f0, cells),
            residuefrac = fill(0.3f0, cells),
        ),
        soilparam = (
            soilph = fill(6.5f0, cells),
            w_sat = fill(0.45f0, layers, cells),
            sand = fill(0.4f0, cells),
            clay = fill(0.2f0, cells),
            tdiff_0 = fill(0.2f0, cells),
            tdiff_15 = fill(0.5f0, cells),
            soildepth = Float32[200, 300, 500, 700, 1300],
        ),
        initial_state = initial_without_shift,
    )

    loaded = InitialDataLoader(data_without_shift, [1, 2], identity)
    @test !hasproperty(loaded.ModelState, :c_shift_fast)
    @test !hasproperty(loaded.ModelState, :c_shift_slow)

    loaded64 = InitialDataLoader(
        data_without_shift, [1, 2], identity; T = Float64,
    )
    @test eltype(loaded64.latitude) == Float64
    @test eltype(loaded64.ModelState.crop.phu) == Float64
    @test eltype(loaded64.ModelState.crop.sdate) == Int32
    @test eltype(loaded64.ModelState.u0.swc) == Float64

    fast_shift = repeat(Float32[0.4, 0.25, 0.15, 0.1, 0.1], 1, cells)
    slow_shift = repeat(Float32[0.5, 0.2, 0.15, 0.1, 0.05], 1, cells)
    data_with_shift = merge(data_without_shift, (
        initial_state = merge(initial_without_shift, (
            c_shift_fast = fast_shift,
            c_shift_slow = slow_shift,
        )),
    ))
    restart = InitialDataLoader(
        data_with_shift, [1, 2], identity; load_c_shift_restart = true,
    )
    @test restart.ModelState.c_shift_fast == fast_shift
    @test restart.ModelState.c_shift_slow == slow_shift

    simulation = initialize_simulation(
        cft1, data_with_shift;
        indices = [1, 2], T = Float32, days = 1, diagnostics = false,
    )
    @test simulation.state.inputs.soil.decomposition.shift_fast == fast_shift
    @test simulation.state.inputs.soil.decomposition.shift_slow == slow_shift
end

@testset "Mineral-N initialization strategies" begin
    slow_n = reshape(Float32[100, 200, 300, 400, 500], 5, 1)
    u0 = (
        soil_NO3 = fill(9000.0f0, 5, 1),
        soil_NH4 = fill(8000.0f0, 5, 1),
    )
    soil = init_soil(1, soilparams.soildepth, identity)
    soil.nitrogen.slow .= slow_n

    NeuralCrop.initialize_soil_mineral_nitrogen!(soil, u0, :restart)
    @test soil.nitrogen.nitrate == u0.soil_NO3
    @test soil.nitrogen.ammonium == u0.soil_NH4

    NeuralCrop.initialize_soil_mineral_nitrogen!(soil, u0, :from_slow_organic_nitrogen)
    @test soil.nitrogen.nitrate ≈ slow_n ./ 100
    @test soil.nitrogen.ammonium ≈ slow_n ./ 100
    @test_throws ArgumentError NeuralCrop.initialize_soil_mineral_nitrogen!(
        soil, u0, :unknown,
    )
    @test_throws ArgumentError NeuralCrop.initialize_soil_mineral_nitrogen!(
        soil, (slown = slow_n,), :restart,
    )
end
