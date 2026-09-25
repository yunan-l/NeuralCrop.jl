module NeuralCrop

# Write your package code here.
# NUMERICS
using Statistics
using Random

# GPU PARALLEL
# import KernelAbstractions: @kernel, @index, @inbounds # get_backend, synchronize
using KernelAbstractions # GPU/CPU parallelization
using Adapt

# INPUT OUTPUT
using NCDatasets, Dates
using SHA
import JLD2: @load, @save, jldsave, load

# PARAMETER HANDLING
import Parameters: @with_kw, @unpack
import MuladdMacro: @muladd

# STRUCTURES
export LPJmLParams, CFTParameters, PhotoParams, ModelParameters, PetPar, Output
export DailyWeather, ClimBuf, CO2
export Crop, CropState, CropFluxes, CropAuxiliary, CropWorkspace
export crop_restart_payload
export CropPhenology, CropCarbonState, CropNitrogenState, CropWaterState
export CropCalendarAuxiliary, CropCanopyState, CropCanopyAuxiliary, CropPhotosynthesisAuxiliary
export CropCarbonFluxes, CropNitrogenFluxes, CropWaterFluxes, CropEvents
export CropStressAuxiliary, CropRootAuxiliary, ManagedLand
export ProcessModules, ModelState, model_state
export AbstractPhotosynthesisPathway, C3Pathway, C4Pathway
export ActiveLandDomain, ExecutionContext, SimulationConfiguration
export HostArchitecture, AcceleratorArchitecture
export VariableSpec, state_schema, validate_state_schema, output_variable_spec
export float_type, array_device, architecture_name
export SoilParams, SoilDecompParams, SoilThermalParams, SnowParams, Soil
export SoilProperties, SoilWater, SoilThermal, SoilCarbon, SoilNitrogen
export SoilDecomposition, SoilManagement, SoilSurfaceLitter, SoilSnow
export CropOutput, SoilOutput, ClimateOutput, CalendarOutput, AnnualOutputAccumulator

# PARAMETERS (CFTs)
export lpjmlparams, photoparams, soilparams, soil_decomp_params, soil_thermal_params, snowparams
export cft1, cft2, cft3, cft4, cft5, cft6, cft7, cft8, cft9, cft10, cft11, cft12
export CFT_NAMES, CFTS, crop_cft
export convert_precision
export FERTILIZER_MODES, fertilizer_mode

# INITIALIZATION
export init_states!, init_climbuf, init_crop, init_pet, init_soil, init_output
export init_weather, init_managed_land
export WaterBalance, init_water_balance
export NitrogenBalance, init_nitrogen_balance
export CarbonBalance, init_carbon_balance
export ThermalBalance, init_thermal_balance

# CLIMATE
export annual_climbuf!, daily_climbuf!, monthlyprec!, infil_perc!, spin_up_climbuf!, update_climbuf!, record_potential_evaporation!, readclimate!, snow!

# PHYSICS FUNCTIONS
# RADIATION
export albedo!, petpar!, apar_crop!, apar_crop_maize!

# CROP
export photosynthesis_C3!, photosynthesis_C4!, carbon_allocation!, respiration!
export photosynthesis!, solve_lambda!
export phenology_crop!, lai_crop!, cultivate!, dynamic_sowing_date!, update_dynamic_sowing_calendar!, harvest_crop!, fertilizer!
export transpiration!, interception!
export crop_carbon!, crop_nitrogen!, ndemand_crop!, nuptake_crop!
export limit_vcmax_by_nitrogen!
export root_distribution, temp_stress
export lpj_bisect, solve_lambda_c3_lpj, solve_lambda_c4_lpj
export solve_lambda_c3!, solve_lambda_c4!

# SOIL
export apply_percolation_enthalpy!, soil_temperature!
export pedotransfer!, soil_carbon!
export evaporation!, soil_infiltration!, soil_evapotranspiration!
export soil_nitrogen!, nitrogen_transform!, soil_cn_decomposition!, nitrogen_deposition!, post_crop_nitrogen_losses!
export soil_decomp_response!
export accumulate_c_shift_response!, equilibrated_c_shift!
export update_surface_litter_properties!, surface_litter_interception!
export litter_tillage!, tillage_hydraulics!, litter_bioturbation!

# UNITS
export deg2rad, ppm2Pa, ppm2bar, hour2day, hour2sec, degCtoK

# DATA
export InitialDataLoader, ClimateDataLoader, DataLoader, DataLoader_winter_wheat
export field_capacity_water, soil_initial_state, model_initial_data
export write_output_nc
export OutputVariable, OutputChunk, OutputStream
export JLD2BlockWriter, NetCDFBlockWriter
export consume_output!, finish_output_stream!, clear_output_timeseries!

# DAILY CROP SIMULATIONS
export daily_crop_C3!, daily_crop_C4!
export CropSimulation, initialize_simulation, agricultural_warmup!, agricultural_warmup_drift
export transition_day!, run_simulation!, simulation_summary
export save_checkpoint, restore_checkpoint!
export estimate_memory

# OPTIONAL AUTOMATIC DIFFERENTIATION
export ADSeasonContext, NeuralSeasonContext, ManagementAdaptationContext
export enzyme_forward_directional, enzyme_forward_gradient, enzyme_zero_tangent
export enzyme_prepare_daily_state!, enzyme_daily_transition_objective
export enzyme_seasonal_loss
export enzyme_seasonal_gradient_blockwise
export enzyme_seasonal_soil_loss
export enzyme_seasonal_soil_gradient_blockwise
export enzyme_seasonal_joint_loss
export enzyme_seasonal_joint_gradient_blockwise
export enzyme_neural_seasonal_loss, enzyme_neural_seasonal_gradient_blockwise
export enzyme_neural_seasonal_loss_buffer!, enzyme_neural_predictions
export enzyme_management_yield_loss, enzyme_management_yield_split_loss
export enzyme_joint_adaptation_yield_loss

# NEURALCROP HYBRID COMPONENTS
export MLPLayout, NeuralCropLayout, NeuralComponents
export neural_parameter_count, neural_trainable_parameter_count
export neural_parameter_ranges, initialize_neural_parameters
export neural_gpp, neural_gpp_multiplier, neural_lambda, neural_vcmax
export neural_allocation_fractions
export neural_decomposition_response, neural_crop_respiration_multiplier
export neural_snow_melt, neural_transpiration_multiplier
export neural_soil_evaporation
export neural_gpp!, neural_crop_carbon!, neural_transpiration!


# process-based crop model
# Parameters
include("parameters/default_params.jl")
include("parameters/cft.jl")
include("utils/kernel_constant.jl")

# Flat-parameter MLPs used by the opt-in NeuralCrop hybrid. Keeping this file
# independent of a neural-network framework makes the daily transition directly
# differentiable by Enzyme and keeps global forward runs lightweight.
include("neural_networks/components.jl")

# Numerics
include("numerics/lpj_bisect.jl")

# Backend adaptation
include("utils/device.jl")

# Initialization
include("processes/initialization/climate/climate.jl")
include("processes/initialization/management/managed_land.jl")
include("processes/initialization/output/output.jl")
include("processes/initialization/crop/phenology.jl")
include("processes/initialization/crop/canopy.jl")
include("processes/initialization/crop/carbon.jl")
include("processes/initialization/crop/nitrogen.jl")
include("processes/initialization/crop/water.jl")
include("processes/initialization/crop/calendar.jl")
include("processes/initialization/crop/photosynthesis.jl")
include("processes/initialization/crop/crop.jl")
include("processes/initialization/soil/properties.jl")
include("processes/initialization/soil/water.jl")
include("processes/initialization/soil/thermal.jl")
include("processes/initialization/soil/carbon.jl")
include("processes/initialization/soil/nitrogen.jl")
include("processes/initialization/soil/decomposition.jl")
include("processes/initialization/soil/management.jl")
include("processes/initialization/soil/surface_litter.jl")
include("processes/initialization/soil/snow.jl")
include("processes/initialization/soil/soil.jl")
include("processes/initialization/init_states.jl")
include("simulations/model_runtime.jl")
include("simulations/runtime_contracts.jl")

# Diagnostics
include("diagnostics/water_balance.jl")
include("diagnostics/nitrogen_balance.jl")
include("diagnostics/carbon_balance.jl")
include("diagnostics/thermal_balance.jl")

# Climate
include("processes/climate/climbuf.jl")
include("processes/climate/temp_stress.jl")
include("processes/climate/spinup_climbuf.jl")
include("processes/climate/readclimate.jl")
include("processes/climate/snow.jl")

# Crop
include("processes/crop/cultivate.jl")
include("processes/crop/dynamic_sowing.jl")
include("processes/crop/phenology.jl")
include("processes/crop/photosynthesis.jl")
include("processes/crop/lambda_solver.jl")
include("processes/crop/carbon_allocation.jl")
include("processes/crop/crop_carbon.jl")
include("processes/crop/lai_crop.jl")
include("processes/crop/radiation.jl")
include("processes/crop/albedo.jl")
include("processes/crop/respiration.jl")
include("processes/crop/interception.jl")
include("processes/crop/transpiration.jl")
include("processes/crop/nitrogen_allocation.jl")
include("processes/crop/nitrogen_demand.jl")
include("processes/crop/nitrogen_uptake.jl")
include("processes/crop/nitrogen_vcmax_limit.jl")
include("processes/crop/fertilizer.jl")
include("processes/crop/harvesting.jl")

# Soil
include("processes/soil/water_ice_pools.jl")
include("processes/soil/tillage.jl")
include("processes/soil/pedotransfer.jl")
include("processes/soil/evaporation.jl")
include("processes/soil/soil_temp.jl")
include("processes/soil/surface_litter.jl")
include("processes/soil/litter_routing.jl")
include("processes/soil/nitrogen_transform.jl")
include("processes/soil/infil_perc.jl")
include("processes/soil/soil_water.jl")
include("processes/soil/soil_carbon.jl")
include("processes/soil/soil_nitrogen.jl")
include("processes/soil/nitrogen_deposition.jl")
include("processes/soil/soil_response.jl")

# Input and output
include("input_output/climate_data_loader.jl")
include("input_output/initial_data_loader.jl")
include("input_output/soil_initialization.jl")
include("input_output/write_output_nc.jl")
include("input_output/stream_output.jl")

# Utilities
include("utils/kernel_launch.jl")
include("utils/conversions.jl")
include("utils/tools.jl")

# Inference-only hybrid process replacements used by ordinary CPU/GPU forward
# simulations. Enzyme-specific differentiation remains in the optional
# extension.
include("neural_networks/runtime.jl")

# Daily crop simulations
include("simulations/daily_crop.jl")
include("simulations/simulation_api.jl")
include("simulations/agricultural_warmup.jl")
include("simulations/memory_estimate.jl")

# The Enzyme implementation is an optional package extension. This file only
# defines the dependency-free season contract and fallback API.
include("ad/adapter.jl")

end
