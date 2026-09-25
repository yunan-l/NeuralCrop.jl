using NeuralCrop

if !isdefined(@__MODULE__, :NeuralCropData)
    include(joinpath(@__DIR__, "..", "lib", "NeuralCropData", "src", "NeuralCropData.jl"))
end
import .NeuralCropData: CFTRegistry, DATA_SCHEMA_VERSION, DatasetCatalog, DatasetSpec,
    HWSD_CN_PREPROCESSING_VERSION, SoilCNTargets, build_crop_mask, climate_forcings,
    crop_inputs, dataset, read_grid, read_management, read_soil_cn_targets,
    read_soil_data, select_cells

length(ARGS) == 2 || error(
    "usage: julia --project=. scripts/global_wheat_smoke.jl " *
    "INPUT_DATA_DIR HWSD_PROFILE_DIR",
)

input_directory, hwsd_directory = abspath.(ARGS)
management_directory = joinpath(input_directory, "management")
climate_directory = joinpath(input_directory, "climate")
soil_directory = joinpath(input_directory, "soil")
grid_path = joinpath(soil_directory, "grid.nc")
soilcode_path = joinpath(soil_directory, "soil_30arcmin_13_types.nc")
soilph_path = joinpath(soil_directory, "soil_pH30arcmin.nc")

registry = CFTRegistry([1], ["temperate cereals"])
single_cft(path, variable; units = "") = DatasetSpec(
    joinpath(management_directory, path), variable; units, cft_ids = [1],
)
catalog = DatasetCatalog(
    Dict{Symbol, DatasetSpec}(
        :grid => DatasetSpec(grid_path, "cellid"),
        :soilcode => DatasetSpec(soilcode_path, "soilcode"),
        :soilph => DatasetSpec(soilph_path, "soilph"),
        :landuse => single_cft("landuse_wheat_rainfed.nc", "landfrac"),
        :sowing_date => single_cft("sdate_wheat_rainfed.nc", "sdate"),
        :phu => single_cft("phu_wheat_rainfed.nc", "phusum"),
        :fertilizer => single_cft("fertilizer_wheat_rainfed.nc", "fertilizer"),
        :manure => single_cft("manure_wheat_rainfed.nc", "manure"),
        :residue_fraction => single_cft(
            "residue_wheat_rainfed.nc", "residuefrac",
        ),
        :temp => DatasetSpec(joinpath(climate_directory, "temp_2015_2016.nc"), "temp"),
        :prec => DatasetSpec(joinpath(climate_directory, "prec_2015_2016.nc"), "prec"),
        :lwnet => DatasetSpec(joinpath(climate_directory, "lwnet_2015_2016.nc"), "lwnet"),
        :swdown => DatasetSpec(joinpath(climate_directory, "swdown_2015_2016.nc"), "swdown"),
        :co2 => DatasetSpec(joinpath(climate_directory, "co2_2015_2016.txt"), "co2"),
    ),
    registry,
)

grid = read_grid(dataset(catalog, :grid); T = Float32)
landuse = read_management(catalog, :landuse, grid, 1; T = Float32)
crop_mask = build_crop_mask(grid, landuse.values)
length(crop_mask.selection.cell_ids) >= 10 || error(
    "rainfed wheat land use contains fewer than ten active cells",
)
selection = select_cells(grid, crop_mask.selection.compact_indices[1:10])
active = trues(1, length(selection.cell_ids))

sowing_date = read_management(
    catalog, :sowing_date, grid, 1; selection, active, T = Float32,
)
phu = read_management(catalog, :phu, grid, 1; selection, active, T = Float32)
fertilizer = read_management(
    catalog, :fertilizer, grid, 1; selection, T = Float32,
)
manure = read_management(catalog, :manure, grid, 1; selection, T = Float32)
residue_fraction = read_management(
    catalog, :residue_fraction, grid, 1; selection, T = Float32,
)
crop = crop_inputs(
    ; sowing_date, phu, fertilizer, manure, residue_fraction,
    fertilizer_mode = :yes, manure_enabled = true,
)

soil = read_soil_data(catalog, grid; selection)
profiles = [
    read_soil_cn_targets(joinpath(hwsd_directory, "cell_$(cell_id).nc"))
    for cell_id in selection.cell_ids
]
for (profile, cell_id) in zip(profiles, selection.cell_ids)
    profile.selection.cell_ids == [cell_id] || error(
        "HWSD profile for cell $cell_id contains $(profile.selection.cell_ids)",
    )
    profile.layer_bounds == profiles[1].layer_bounds || error(
        "HWSD layer bounds differ for cell $cell_id",
    )
end
targets = SoilCNTargets(
    selection,
    profiles[1].layer_bounds,
    reduce(hcat, (profile.soil_organic_carbon for profile in profiles)),
    reduce(hcat, (profile.total_nitrogen for profile in profiles)),
    reduce(hcat, (profile.coverage for profile in profiles)),
    BitMatrix(reduce(hcat, (profile.uncertain for profile in profiles))),
    (
        schema_version = DATA_SCHEMA_VERSION,
        preprocessing_version = HWSD_CN_PREPROCESSING_VERSION,
        source_version = "HWSD v2.01",
        deep_rule = :extend_deepest_density,
        source_files = ["cell_$(cell_id).nc" for cell_id in selection.cell_ids],
    ),
)
initial_state = soil_initial_state(targets, soil)
initial_data = model_initial_data(grid, soil, crop, initial_state)

reader = climate_blocks(
    catalog, grid; selection, block_days = 365, T = Float32,
)
simulation = initialize_simulation(
    cft1,
    initial_data;
    days = climate_days(reader),
    cell_ids = selection.cell_ids,
    T = Float32,
    diagnostics = true,
    irrigation = false,
    manure = true,
    fertilizer = :yes,
    with_tillage = true,
)
warmup_report = agricultural_warmup!(simulation, climate_forcings(reader))
warmup_drift = agricultural_warmup_drift(warmup_report)
run_simulation!(simulation, climate_forcings(reader); spinup = false)

summary = simulation_summary(simulation)
all(isfinite, simulation.output.crop.gpp) || error("GPP contains non-finite values")
all(isfinite, simulation.output.crop.npp) || error("NPP contains non-finite values")
all(isfinite, simulation.output.crop.biomass) || error("biomass contains non-finite values")
all(>=(0), simulation.state.prognostic.soil.water.storage) ||
    error("soil water contains negative values")
for pool in values(simulation.state.prognostic.soil.carbon)
    all(>=(0), pool) || error("soil carbon contains negative values")
end
for pool in values(simulation.state.prognostic.soil.nitrogen)
    all(>=(0), pool) || error("soil nitrogen contains negative values")
end

println("cell_ids = ", selection.cell_ids)
println("landfrac = ", crop_mask.fraction[1, 1:10])
println("climate_days = ", climate_days(reader))
println("warmup_years = ", warmup_report.years)
println("warmup_drift = ", warmup_drift)
println("summary = ", summary)
