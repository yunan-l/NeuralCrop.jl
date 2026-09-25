module NeuralCropData

using Dates
using NCDatasets
using TOML

include("contracts.jl")
include("catalog.jl")
include("grid.jl")
include("netcdf.jl")
include("masks.jl")
include("soil.jl")
include("hwsd.jl")
include("management.jl")
include("climate.jl")

export DATA_SCHEMA_VERSION
export DataProvenance, DatasetSpec, DatasetCatalog, CFTRegistry, ManagementBands
export GridIndex, CellSelection, PatchDomain, CompactVariable, TimeCellData, CropMask
export CO2Series, ClimateBlock, ClimateBlockReader, ClimateForcingReader
export PrefetchedClimateForcingReader
export SoilLookup, SoilData, DEFAULT_SOIL_LOOKUP_VERSION
export SoilCNTargets, SoilPoolAllocation, SoilCNAggregator
export HWSD_LAYER_BOUNDS, NEURALCROP_SOIL_LAYER_BOUNDS, HWSD_CN_PREPROCESSING_VERSION
export load_catalog, dataset, cft_index, cft_name
export read_grid, all_cells, select_cells, compact_spatial, expand_to_grid
export read_compact_variable, read_static_cell
export build_crop_mask, build_patch_domain, combine_patch_domains
export default_soil_lookup, soil_data_from_values, read_soil_data, soil_properties
export hwsd_layer_stocks, mix_hwsd_components, remap_hwsd_layers
export hwsd_tile_mapping
export init_soil_cn_aggregator, accumulate_soil_cn!, finish_soil_cn
export soil_cn_conservation
export preprocess_hwsd_cn, write_soil_cn_targets, read_soil_cn_targets
export write_soil_pool_allocation, read_soil_pool_allocation
export read_management, validate_management, crop_inputs, management_schedule
export read_co2_series, climate_blocks, climate_days, read_climate_block
export climate_forcing, climate_forcings
export prefetch_climate_forcings

end
