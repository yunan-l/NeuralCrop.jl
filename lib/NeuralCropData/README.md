# NeuralCropData.jl

For bounded server-side extraction of all 12 crops in rainfed and irrigated
form, see `scripts/prepare_global_cft_subset.jl` and
`config/global_cft_subset.example.toml`. The complete workflow is documented
in the NeuralCrop manual under **Global CFT input data**.

`NeuralCropData.jl` is the backend-neutral input-data layer for NeuralCrop. It
aligns external NetCDF data to the canonical NeuralCrop `cellid` grid, constructs
land-use/CFT selections, and loads compact soil and management arrays.

The package intentionally does not depend on CUDA or run scientific processes.
NeuralCrop.jl remains responsible for device transfer, state initialization, and
simulation.

The implementation roadmap is maintained in
[`docs/neuralcrop_data_roadmap.md`](../../docs/neuralcrop_data_roadmap.md).

```julia
using NeuralCrop
using NeuralCropData

catalog = load_catalog("catalog.toml")
grid = read_grid(dataset(catalog, :grid))

landuse = read_management(catalog, :landuse, grid, 1; years = 2000:2019)
crop_mask = build_crop_mask(grid, landuse.values)
soil = read_soil_data(catalog, grid; selection = crop_mask.selection)

# `crop` comes from crop_inputs(...) using single-time management readers.
hwsd_targets = read_soil_cn_targets("hwsd_cn_targets.nc")
# NeuralCrop, not NeuralCropData, applies its initialization policy:
initial_state = soil_initial_state(hwsd_targets, soil)
initial = model_initial_data(grid, soil, crop, initial_state)

reader = climate_blocks(
    catalog,
    grid;
    selection = crop_mask.selection,
    start_year = 2000,
    end_year = 2019,
    block_days = 31,
)

simulation = initialize_simulation(
    cft1,
    initial;
    days = climate_days(reader),
    device = identity,
)
run_simulation!(simulation, climate_forcings(reader); spinup = false)
```

Before allocating a large run, estimate the final model payload and
forcing-buffer peak from the selected-cell count:

```julia
memory = estimate_memory(
    length(reader.selection.compact_indices), climate_days(reader);
    T = Float32,
    diagnostics = false,
    block_days = reader.block_days,
    backend = :accelerator,
    prefetch = false,
)
memory.recommended_host_peak_gib
memory.recommended_device_peak_gib
```

After initialization, `estimate_memory(simulation, reader)` replaces the
persistent-state formula with the bytes of the actual allocated arrays.

Setting `prefetch=true` only reserves one additional host forcing block in the
estimate. It does not enable asynchronous reading.

`climate_forcings(reader)` remains lazy: each block is read only when
`run_simulation!` requests it, then NeuralCrop converts it to the simulation
precision and transfers it to the configured backend. Crop, soil, climate
buffers, diagnostics, and checkpoint state remain continuous between blocks.

`crop_mask.selection` has a stable, `cellid`-sorted order for the full run.
`crop_mask.active` retains the year-specific land-use state for switching crop
processes on and off without reallocating soil state.

The catalog uses the canonical 12-crop registry and declares separate band
positions for each file. Pass `irrigated=true` to select its irrigated band.
Mineral fertilizer follows the LPJmL-style `:no`, `:yes`, and `:auto` modes;
only `:yes` requires fertilizer input. Manure remains an independent switch and
requires manure input when enabled. Tillage is a model configuration switch,
not a data input.

Climate blocks read only the requested daily rows from `temp`, `prec`, `lwnet`,
and `swdown`. Annual global CO₂ from the two-column text file is matched by
calendar year and emitted as a small daily vector, so arbitrary block boundaries
remain correct.

The model calendar is explicitly 365-day. Standard/Gregorian inputs have
February 29 removed, while `noleap` and `365_day` inputs are retained. A
`360_day` calendar or an incomplete year is rejected. Climate values are
normalized before handoff to °C, mm/day, and W/m²; provenance retains both the
source metadata and these canonical model units.

## HWSD carbon and nitrogen targets

The HWSD preprocessing API converts the official seven HWSD v2.x layers to
NeuralCrop's five soil layers without assigning fast, slow, or litter pools.
Processing is tile-based: `hwsd_tile_mapping` returns compact target indices
and spherical pixel areas, and `accumulate_soil_cn!` adds one source tile to a
bounded-memory accumulator. `finish_soil_cn` produces SOC and total-N targets,
coverage and uncertainty maps, and reproducible provenance.

HWSD ends at 2 m. By default, `remap_hwsd_layers` extends the stock density of
the 1.5–2 m layer through NeuralCrop's 2–3 m layer and marks every extrapolated
value uncertain. Pass `deep_rule=:missing` when extrapolation is unsuitable.
The generated targets remain source totals rather than runtime pools.

For a real-data check after downloading the official `HWSD2.bil` and
`HWSD2.mdb`, install `mdbtools` and run:

```bash
julia --project=. lib/NeuralCropData/scripts/extract_hwsd_cell.jl \
  /path/to/HWSD2 /path/to/grid.nc 172 197 /path/to/hwsd_cn_test.nc 0
```

The longitude and latitude arguments are one-based indices into the canonical
`grid.nc`; the script reads their coordinates instead of reconstructing the
grid from an assumed origin. It reads all 3600 source pixels within that 0.5°
cell, mixes soil components by HWSD `SHARE`, aggregates stocks by spherical
pixel area, and writes the five-layer target file with coverage and uncertainty
flags. Ocean, water, and glacier pixels are excluded from the soil-area mean;
shallow components contribute zero below their HWSD `ROOT_DEPTH` class.
If the target cell lacks a complete profile, the extractor uses the nearest
complete 0.5° HWSD cell within the configured search radius. Filled layers are
marked uncertain, and donor coordinates, distance, and the original minimum
layer coverage are stored in the NetCDF provenance attributes. The command
fails if no donor is found; it never silently accepts an incomplete profile.
