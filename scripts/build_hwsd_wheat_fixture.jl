using NeuralCrop
using JLD2

if !isdefined(@__MODULE__, :NeuralCropData)
    include(joinpath(@__DIR__, "..", "lib", "NeuralCropData", "src", "NeuralCropData.jl"))
end
import .NeuralCropData: CFTRegistry, DATA_SCHEMA_VERSION, DatasetCatalog, DatasetSpec,
    HWSD_CN_PREPROCESSING_VERSION, SoilCNTargets, SoilPoolAllocation, crop_inputs,
    dataset, read_grid, read_management, read_soil_cn_targets, read_soil_data,
    read_soil_pool_allocation, select_cells

# Keep this fixture aligned with the ten daily forcing columns in
# `notebooks/climate_2000_2009.jld2`.  These are `(longitude, latitude)`
# indices in the canonical grid used by the original notebook fixture.
const WHEAT_FIXTURE_GRID_INDICES = (
    (222, 3), (223, 3), (215, 4), (216, 4), (217, 4),
    (219, 4), (222, 4), (223, 4), (214, 5), (215, 5),
)

length(ARGS) in (4, 5) || error(
    "usage: julia --project=. scripts/build_hwsd_wheat_fixture.jl " *
    "INPUT_DATA_DIR HWSD_PROFILE_DIR ALLOCATION_PATH EXAMPLES_OUTPUT [NOTEBOOK_OUTPUT]",
)

input_directory, hwsd_directory, allocation_path, examples_output = abspath.(ARGS[1:4])
notebook_output = length(ARGS) == 5 ? abspath(ARGS[5]) : nothing
management_directory = joinpath(input_directory, "management")
soil_directory = joinpath(input_directory, "soil")

single_cft(path, variable) = DatasetSpec(
    joinpath(management_directory, path), variable; cft_ids = [1],
)
catalog = DatasetCatalog(
    Dict{Symbol, DatasetSpec}(
        :grid => DatasetSpec(joinpath(soil_directory, "grid.nc"), "cellid"),
        :soilcode => DatasetSpec(joinpath(soil_directory, "soil_30arcmin_13_types.nc"), "soilcode"),
        :soilph => DatasetSpec(joinpath(soil_directory, "soil_pH30arcmin.nc"), "soilph"),
        :sowing_date => single_cft("sdate_wheat_rainfed.nc", "sdate"),
        :phu => single_cft("phu_wheat_rainfed.nc", "phusum"),
        :fertilizer => single_cft("fertilizer_wheat_rainfed.nc", "fertilizer"),
        :manure => single_cft("manure_wheat_rainfed.nc", "manure"),
        :residue_fraction => single_cft("residue_wheat_rainfed.nc", "residuefrac"),
    ),
    CFTRegistry([1], ["temperate cereals"]),
)

grid = read_grid(dataset(catalog, :grid); T = Float32)
compact_index = Dict(cell_id => index for (index, cell_id) in pairs(grid.cell_ids))
fixture_cell_ids = Int32[
    grid.cellid[longitude_index, latitude_index]
    for (longitude_index, latitude_index) in WHEAT_FIXTURE_GRID_INDICES
]
all(haskey(compact_index, cell_id) for cell_id in fixture_cell_ids) ||
    error("fixture cell IDs are not present in the canonical grid")
# NeuralCropData reads compact selections in canonical cell-ID order.  The
# notebook climate columns retain their historical longitude/latitude order,
# so restore that order only at the final fixture boundary.
selection = select_cells(grid, sort([compact_index[cell_id] for cell_id in fixture_cell_ids]))
fixture_order = [findfirst(==(cell_id), selection.cell_ids) for cell_id in fixture_cell_ids]
active = trues(1, length(selection.cell_ids))

sowing_date = read_management(catalog, :sowing_date, grid, 1; selection, active, T = Float32)
phu = read_management(catalog, :phu, grid, 1; selection, active, T = Float32)
fertilizer = read_management(catalog, :fertilizer, grid, 1; selection, T = Float32)
manure = read_management(catalog, :manure, grid, 1; selection, T = Float32)
residue_fraction = read_management(catalog, :residue_fraction, grid, 1; selection, T = Float32)
crop = crop_inputs(
    ; sowing_date, phu, fertilizer, manure, residue_fraction,
    fertilizer_mode = :yes, manure_enabled = true,
)
soil = read_soil_data(catalog, grid; selection)

profiles = [
    read_soil_cn_targets(joinpath(hwsd_directory, "cell_$(cell_id).nc"))
    for cell_id in selection.cell_ids
]
all(profile -> profile.layer_bounds == profiles[1].layer_bounds, profiles) ||
    error("HWSD fixture profiles use inconsistent layer bounds")
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
        source_files = ["cell_$(cell_id).nc" for cell_id in selection.cell_ids],
    ),
)

full_allocation = read_soil_pool_allocation(allocation_path; T = Float32)
allocation_index = Dict(cell_id => index for (index, cell_id) in pairs(full_allocation.selection.cell_ids))
all(haskey(allocation_index, cell_id) for cell_id in selection.cell_ids) ||
    error("calibrated allocation does not cover every fixture cell")
columns = [allocation_index[cell_id] for cell_id in selection.cell_ids]
allocation = SoilPoolAllocation(
    selection,
    full_allocation.fast_carbon_fraction[:, columns],
    full_allocation.fast_nitrogen_fraction[:, columns],
    full_allocation.c_shift_fast[:, columns],
    full_allocation.c_shift_slow[:, columns];
    cft_id = full_allocation.cft_id,
    irrigated = full_allocation.irrigated,
    provenance = full_allocation.provenance,
)
allocation.cft_id == 1 || error("fixture requires the CFT 1 allocation")
allocation.irrigated && error("fixture requires the rainfed allocation")

initial_state = soil_initial_state(targets, soil; allocation)
canonical_data = model_initial_data(grid, soil, crop, initial_state)

reorder_vector(values) = values[fixture_order]
reorder_layer_cell(values) = values[:, fixture_order]
initial_data = merge(canonical_data, (
    coords = copy(fixture_cell_ids),
    latitude = reorder_vector(canonical_data.latitude),
    crop = (
        sdate = reorder_vector(canonical_data.crop.sdate),
        phu = reorder_vector(canonical_data.crop.phu),
        manure = reorder_vector(canonical_data.crop.manure),
        fertilizer = reorder_vector(canonical_data.crop.fertilizer),
        residuefrac = reorder_vector(canonical_data.crop.residuefrac),
    ),
    soilparam = (
        soilcode = reorder_vector(canonical_data.soilparam.soilcode),
        soilph = reorder_vector(canonical_data.soilparam.soilph),
        w_sat = reorder_layer_cell(canonical_data.soilparam.w_sat),
        sand = reorder_vector(canonical_data.soilparam.sand),
        silt = reorder_vector(canonical_data.soilparam.silt),
        clay = reorder_vector(canonical_data.soilparam.clay),
        tdiff_0 = reorder_vector(canonical_data.soilparam.tdiff_0),
        tdiff_15 = reorder_vector(canonical_data.soilparam.tdiff_15),
        soildepth = canonical_data.soilparam.soildepth,
    ),
    initial_state = (
        swc = reorder_layer_cell(canonical_data.initial_state.swc),
        litc = reorder_layer_cell(canonical_data.initial_state.litc),
        fastc = reorder_layer_cell(canonical_data.initial_state.fastc),
        slowc = reorder_layer_cell(canonical_data.initial_state.slowc),
        litn = reorder_layer_cell(canonical_data.initial_state.litn),
        fastn = reorder_layer_cell(canonical_data.initial_state.fastn),
        slown = reorder_layer_cell(canonical_data.initial_state.slown),
        c_shift_fast = reorder_layer_cell(canonical_data.initial_state.c_shift_fast),
        c_shift_slow = reorder_layer_cell(canonical_data.initial_state.c_shift_slow),
    ),
    c_shift_fast = reorder_layer_cell(canonical_data.c_shift_fast),
    c_shift_slow = reorder_layer_cell(canonical_data.c_shift_slow),
))

for output in filter(!isnothing, (examples_output, notebook_output))
    mkpath(dirname(output))
    jldsave(output; initial_data)
    println("wrote ", output)
end
