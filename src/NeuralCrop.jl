module NeuralCrop

# Write your package code here.

# NUMERICS
using DiffEqFlux, OrdinaryDiffEq, SciMLSensitivity, Statistics, LinearAlgebra, StatsBase
# import SciMLBase: solve, AbstractDEProblem
using SciMLBase

# GPU PARALLEL
# import KernelAbstractions: @kernel, @index, @inbounds # get_backend, synchronize
using KernelAbstractions ##  GPU/CPU parallelization
using Lux, CUDA, LuxCUDA, Adapt

# TRAINING
using ComponentArrays, Zygote, Optimisers, Optimization

# INPUT OUTPUT
using DataFrames, NCDatasets, Random, Plots, XLSX, Dates, Printf,  ParameterSchedulers, ProgressMeter
import JLD2: @load, @save

# PARAMETER HANDLING
import Parameters: @with_kw, @unpack
import MuladdMacro: @muladd

# STRUCTURES
export LPJmLParam, PftParameters, PhotoPar, SoilPar, Photos, PetPar, ClimBuf, Crop, Calendar, Managed_land, Soil, Output

# PARAMETERS (PFTs)
export lpjmlparam, photopar, soilpar, cft1, k1, cft2, k2, cft3, k3, cft4, k4

# INITIALIZATION
export init_structs!, init_climbuf, init_crop, init_pet, init_soil, init_data_norm, init_output

# CLIMBUF
export annual_climbuf!, daily_climbuf!, infil_perc!, spin_up_climbuf!, update_climbuf!, readclimate!         

# PHYSICS FUNCTIONS
### RADIATION
export albedo!, petpar!, apar_crop!, apar_crop_maize!

### CROP
export photosynthesis_C3!, photosynthesis_C4!, carbon_allocation!, respiration!
export phenology_crop!, lai_crop!, lai_deficit!, cultivate!, harvest_crop!, fertilizer!
export transpiration!, interception!
export crop_nitrogen!, crop_nitrogen_old!, ndemand_crop!, nuptake_crop!
export root_distribution, temp_stress
export crop_carbon!, crop_carbon_node!, crop_carbon_hybrid!, hybrid_photos_C3!, hybrid_photos_C4!

### SOIL
export soiltemp_lag!
export pedotransfer!, update_lit_tillage!, update_lit_winter_wheat!,  soil_carbon!
export evaporation!, soil_water!
export nitrogen_transform!, soil_nitrogen!, update_litc_tillage!, update_litn_tillage!

# OUTPUT
export output_training!, output_predict!, output_finetune!, output_yield!

# UNITS
export deg2rad, ppm2Pa, ppm2bar, hour2day, hour2sec, degCtoK
export min_max_norm, z_score_norm # OnePoint_one_dimension, OnePoint_dimensions, z_score_one_dimension, min_max_one_dimension, divide_data_one_dimension, divide_data_dimensions

# DATA
export InitilDataLoader, ClimateDataLoader, DataLoader, DataLoader_winter_wheat

# NEURAL NETWORK
export NODE, MLP, solve, SciMLEuler, SciMLEuler_litc, SciMLEuler_soilc, neural_gpp, neural_lambda, neural_vmax, neural_stoc, neural_allocation, hybrid_litc, hybrid_soilc, hybrid_litn, hybrid_soiln,
       neural_moisture, get_mlp, get_node
export train_loop_winter_wheat_rollout!, train_loop_rollout!, loss_crop_rollout!
export daily_crop_C3_training!

# DAILY CROP SIMULATIONS
export daily_crop_C3!, daily_crop_C4!

# PLOT
 export load_nc_file_one_dimension, load_nc_file_dimensions, plot_loss_curve


# process-based crop model
### Variables
include("physics/variables/define_structs.jl")
include("physics/variables/default_param.jl")
include("physics/variables/units.jl")
include("physics/variables/init_var.jl")
include("physics/variables/init_struct.jl")
include("physics/variables/callback.jl")
include("physics/variables/output.jl")
include("physics/variables/DataLoader.jl")

### Climate
include("physics/climate/climbuf.jl")
include("physics/climate/temp_stress.jl")
include("physics/climate/spinup_climbuf.jl")
include("physics/climate/readclimate.jl")

### Crop
include("physics/crop/cultivate.jl")
include("physics/crop/phenology.jl")
include("physics/crop/photosynthesis.jl")
include("physics/crop/carbon_allocation.jl")
include("physics/crop/crop_carbon.jl")
include("physics/crop/lai_crop.jl")
include("physics/crop/radiation.jl")
include("physics/crop/albedo.jl")
include("physics/crop/respiration.jl")
include("physics/crop/interception.jl")
include("physics/crop/transpiration.jl")
include("physics/crop/nitrogen_allocation.jl")
include("physics/crop/nitrogen_demand.jl")
include("physics/crop/nitrogen_uptake.jl")
include("physics/crop/fertilizer.jl")
include("physics/crop/harvesting.jl")

### Soil
include("physics/soil/pedotransfer.jl")
include("physics/soil/evaporation.jl")
include("physics/soil/soil_temp.jl")
include("physics/soil/nitrogen_transform.jl")
include("physics/soil/infil_perc.jl")
include("physics/soil/soil_water.jl")
include("physics/soil/soil_carbon.jl")
include("physics/soil/soil_nitrogen.jl")

### Hybrid
include("hybrid/crop_carbon.jl")
include("hybrid/photosynthesis.jl")
include("hybrid/soil_carbon.jl")
include("hybrid/soil_nitrogen.jl")
include("hybrid/soil_water.jl")

### Neural network
include("neural_network/define_net_struct.jl")
include("neural_network/init_net.jl")
include("neural_network/neural_emulator.jl")
include("neural_network/solver.jl")
include("neural_network/loss.jl")
include("neural_network/training_loop.jl")

### Training
include("training/daily_crop_C3_training.jl")

### Utilities
include("utilities/data_loader.jl")
include("utilities/data_norm.jl")
include("utilities/utils.jl")
include("utilities/visualization.jl")

### Simulations
include("simulations/daily_crop_C3.jl")
include("simulations/daily_crop_C4.jl")

end