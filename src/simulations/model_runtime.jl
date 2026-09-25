abstract type AbstractPhotosynthesisPathway end
struct C3Pathway <: AbstractPhotosynthesisPathway end
struct C4Pathway <: AbstractPhotosynthesisPathway end

"""Immutable process configuration. Numerical arrays live in `ModelState`."""
struct ProcessModules{C, G, P <: AbstractPhotosynthesisPathway}
    crop::C
    global_parameters::G
    pathway::P
end

function ProcessModules(crop, global_parameters)
    pathway = crop.path == 1 ? C3Pathway() : crop.path == 2 ? C4Pathway() : throw(
        ArgumentError("unsupported photosynthetic pathway $(crop.path)"),
    )
    return ProcessModules(crop, global_parameters, pathway)
end

pathway_value(::C3Pathway) = Val(:C3)
pathway_value(::C4Pathway) = Val(:C4)

"""
Numerical variables grouped by lifecycle, independently of the crop and soil
process hierarchy. The leaves are the original backend arrays; constructing
this container does not copy model data.
"""
struct ModelState{P, F, A, I, E, W, O}
    prognostic::P
    fluxes::F
    auxiliary::A
    inputs::I
    events::E
    workspace::W
    output::O
end

# Process kernels consume only the canonical lifecycle state. Legacy direct
# `Crop`/`Soil` access is confined to test fixtures while those tests migrate.
crop_prognostic(x::ModelState) = x.prognostic.crop
crop_fluxes(x::ModelState) = x.fluxes.crop
crop_events(x::ModelState) = x.events.crop
crop_canopy_auxiliary(x::ModelState) = x.auxiliary.crop.canopy
crop_photosynthesis_auxiliary(x::ModelState) = x.auxiliary.crop.photosynthesis
crop_stress_auxiliary(x::ModelState) = x.auxiliary.crop.stress
crop_phenology_auxiliary(x::ModelState) = x.auxiliary.crop.phenology
crop_phenology_input(x::ModelState) = x.inputs.crop.phenology
crop_calendar_input(x::ModelState) = x.inputs.crop.calendar
crop_root_auxiliary(x::ModelState) = x.auxiliary.crop.root
crop_root_input(x::ModelState) = x.inputs.crop.root

soil_properties(x::ModelState) = x.inputs.soil.properties
for (selector, lifecycle_group) in (
    (:soil_water_prognostic, :(x.prognostic.soil.water)),
    (:soil_water_fluxes, :(x.fluxes.soil.water)),
    (:soil_water_auxiliary, :(x.auxiliary.soil.water)),
    (:soil_thermal_prognostic, :(x.prognostic.soil.thermal)),
    (:soil_thermal_fluxes, :(x.fluxes.soil.thermal)),
    (:soil_thermal_input, :(x.inputs.soil.thermal)),
    (:soil_carbon_prognostic, :(x.prognostic.soil.carbon)),
    (:soil_carbon_fluxes, :(x.fluxes.soil.carbon)),
    (:soil_carbon_auxiliary, :(x.auxiliary.soil.carbon)),
    (:soil_nitrogen_prognostic, :(x.prognostic.soil.nitrogen)),
    (:soil_nitrogen_fluxes, :(x.fluxes.soil.nitrogen)),
    (:soil_nitrogen_auxiliary, :(x.auxiliary.soil.nitrogen)),
    (:soil_decomposition_auxiliary, :(x.auxiliary.soil.decomposition)),
    (:soil_decomposition_input, :(x.inputs.soil.decomposition)),
    (:soil_decomposition_workspace, :(x.workspace.soil.decomposition)),
    (:soil_management_prognostic, :(x.prognostic.soil.management)),
    (:soil_management_fluxes, :(x.fluxes.soil.management)),
    (:soil_management_input, :(x.inputs.soil.management)),
    (:soil_surface_litter_prognostic, :(x.prognostic.soil.surface_litter)),
    (:soil_surface_litter_fluxes, :(x.fluxes.soil.surface_litter)),
    (:soil_surface_litter_auxiliary, :(x.auxiliary.soil.surface_litter)),
    (:soil_snow_prognostic, :(x.prognostic.soil.snow)),
    (:soil_snow_fluxes, :(x.fluxes.soil.snow)),
)
    @eval begin
        $selector(x::ModelState) = $lifecycle_group
    end
end

_fields(value, names::Tuple) = NamedTuple{names}(map(name -> getproperty(value, name), names))

function _soil_lifecycle_views(soil::Soil)
    prognostic = (
        water = _fields(soil.water, (
            :storage, :ice_storage, :wilting_ice_fraction,
            :available_ice_storage, :free_ice_storage, :saturation_fraction,
        )),
        thermal = _fields(soil.thermal, (
            :temperature, :enthalpy, :frozen_fraction, :freeze_depth,
            :heat_capacity_frozen, :heat_capacity_unfrozen, :latent_heat,
            :conductivity_frozen, :conductivity_unfrozen, :water_reference,
            :initialized,
        )),
        carbon = _fields(soil.carbon, (:litter, :fast, :slow)),
        nitrogen = _fields(soil.nitrogen, (:nitrate, :ammonium, :litter, :fast, :slow)),
        management = _fields(soil.management, (:tillage_density_factor,)),
        surface_litter = _fields(soil.surface_litter, (
            :dry_matter, :depth, :cover, :water_storage, :temperature, :conductivity,
        )),
        # Height and fraction currently feed the following day's pre-snow
        # albedo calculation, so they are transition state, not diagnostics.
        snow = _fields(soil.snow, (:pack, :height, :fraction)),
    )
    fluxes = (
        water = _fields(soil.water, (
            :evaporation, :influx, :outflux, :surface_runoff, :lateral_runoff,
            :bottom_drainage, :infiltration, :percolation,
        )),
        thermal = _fields(soil.thermal, (
            :percolation_energy, :surface_energy_flux, :energy_residual,
            :untracked_water_energy_flux, :rain_energy_input, :snowmelt_energy_input,
            :lateral_runoff_energy_output, :bottom_drainage_energy_output,
            :percolation_energy_residual,
        )),
        carbon = _fields(soil.carbon, (
            :input, :decomposed_litter, :decomposed_fast, :decomposed_slow,
            :litter_to_fast, :litter_to_slow, :heterotrophic_respiration,
        )),
        nitrogen = _fields(soil.nitrogen, (
            :input, :decomposed_litter, :decomposed_fast, :decomposed_slow,
            :litter_to_fast, :litter_to_slow, :mineralization, :immobilization,
            :nitrification, :n2o_nitrification, :denitrification,
            :n2o_denitrification, :n2_denitrification, :volatilization, :leaching,
        )),
        management = _fields(soil.management, (
            :tillage_carbon, :tillage_nitrogen,
            :bioturbation_carbon, :bioturbation_nitrogen,
        )),
        surface_litter = _fields(soil.surface_litter, (:interception, :evaporation)),
        snow = _fields(soil.snow, (:melt, :sublimation, :runoff)),
    )
    auxiliary = (
        water = _fields(soil.water, (
            :relative_content, :free_water, :wilting_fraction, :wilting_storage,
            :field_capacity, :saturation_storage, :beta,
            :holding_capacity_fraction, :holding_capacity_storage,
            :saturated_conductivity,
        )),
        carbon = _fields(soil.carbon, (:litter_response,)),
        nitrogen = _fields(soil.nitrogen, (:litter_response,)),
        decomposition = _fields(soil.decomposition, (:response, :litter_response)),
        surface_litter = _fields(soil.surface_litter, (:water_capacity,)),
    )
    inputs = (
        properties = soil.properties,
        thermal = _fields(soil.thermal, (:diffusivity_0, :diffusivity_15)),
        decomposition = _fields(soil.decomposition, (:shift_fast, :shift_slow)),
        management = _fields(soil.management, (:tillage_fraction,)),
    )
    workspace = (
        decomposition = _fields(soil.decomposition, (
            :layer_scratch_1, :layer_scratch_2,
            :surface_scratch_1, :surface_scratch_2,
        )),
    )
    return (; prognostic, fluxes, auxiliary, inputs, workspace)
end


function model_state(climbuf, crop, pet, soil, managed_land, weather, output)
    soil_views = _soil_lifecycle_views(soil)
    crop_auxiliary = (
        phenology = _fields(crop.auxiliary.phenology, (:fphu,)),
        root = _fields(crop.auxiliary.root, (:zone_available_water,)),
        canopy = crop.auxiliary.canopy,
        photosynthesis = crop.auxiliary.photosynthesis,
        stress = crop.auxiliary.stress,
    )
    crop_inputs = (
        phenology = _fields(crop.auxiliary.phenology, (:phu, :winter_type)),
        calendar = _fields(crop.auxiliary.calendar, (:sowing_date, :prescribed_sowing_date)),
        root = _fields(crop.auxiliary.root, (:distribution,)),
    )
    return ModelState(
        (crop = crop.state, soil = soil_views.prognostic, climate = climbuf),
        (crop = crop.fluxes, soil = soil_views.fluxes),
        (crop = crop_auxiliary, soil = soil_views.auxiliary, pet = pet),
        (weather = weather, management = managed_land,
         crop = crop_inputs, soil = soil_views.inputs),
        (crop = crop.events,),
        (crop = crop.workspace, soil = soil_views.workspace),
        output,
    )
end
