"""Persistent crop state required to advance the next daily transition."""
mutable struct CropState{P, C, B, N, W}
    phenology::P # Integrated phenology and discrete crop modes.
    canopy::C    # Checkpointed canopy variables required by the next daily transition.
    carbon::B    # Conserved living-plant carbon pools.
    nitrogen::N  # Conserved plant-N pools and seasonal management memory.
    water::W     # Integrated seasonal crop-water demand and supply.
end

"""Current-day carbon, nitrogen, and water transfers."""
mutable struct CropFluxes{C, N, W}
    carbon::C   # Current-day crop carbon transfers.
    nitrogen::N # Current-day crop and management nitrogen transfers.
    water::W    # Current-day canopy and root-water transfers.
end

"""Current-day nitrogen and water stress diagnostics."""
mutable struct CropStressAuxiliary{A}
    nitrogen_demand_total::A # Potential whole-plant nitrogen demand (gN m⁻² day⁻¹).
    nitrogen_demand_leaf::A  # Potential leaf nitrogen demand (gN m⁻² day⁻¹).
    nitrogen_deficit::A      # Unmet crop nitrogen demand (gN m⁻² day⁻¹).
    water_deficit::A         # Water-deficit factor used by allocation (percent, 0–100).
end

"""Static root geometry plus the current-day root-zone water diagnostic."""
mutable struct CropRootAuxiliary{A, V}
    zone_available_water::A # Top-three-layer root-weighted plant-available water (mm).
    distribution::V         # Static fraction of crop roots in each soil layer (sums to 1).
end

mutable struct CropAuxiliary{H, K, R, C, P, S}
    phenology::H     # Static/derived phenological variables and process memory.
    calendar::K      # Calendar configuration.
    root::R          # Static root distribution and daily root-zone water diagnostic.
    canopy::C        # Canopy geometry, radiation, and conductance variables.
    photosynthesis::P # Photosynthetic capacities and limiting factors.
    stress::S        # Nitrogen/water demand, limitation, and root-zone diagnostics.
end

"""Implementation-only buffers, currently empty and excluded from restart/output."""
struct CropWorkspace end

"""Lifecycle-grouped crop storage; all leaves remain backend arrays (CPU/GPU SoA)."""
mutable struct Crop{S, F, A, E, W}
    state::S     # Prognostic stocks, integrals, and discrete modes retained across days.
    fluxes::F    # Transfers produced and overwritten during the current day.
    auxiliary::A # Derived variables and explicitly retained process memory.
    events::E    # One-day discrete sowing and harvest events.
    workspace::W # Implementation-only preallocated scratch storage.
end

"""
    crop_restart_payload(crop)

Return the scientific crop checkpoint boundary. Daily fluxes, events,
recomputable auxiliary fields, and implementation workspace are intentionally
excluded. Static configuration and derived geometry are reconstructed from
parameters and input data.
"""
crop_restart_payload(crop::Crop) = (
    state = crop.state,
    process_memory = (
        phenology = (
            phu = crop.auxiliary.phenology.phu,
            winter_type = crop.auxiliary.phenology.winter_type,
        ),
        calendar = (
            sowing_date = crop.auxiliary.calendar.sowing_date,
            prescribed_sowing_date = crop.auxiliary.calendar.prescribed_sowing_date,
        ),
    ),
)

"""
init_crop(cell_size, device; soil_layers=5)

Allocate the complete lifecycle-grouped crop storage on `device`. Calendar
state and photosynthesis auxiliary variables are owned by `Crop`.
"""
init_crop(cell_size::Int, device; kwargs...) =
    init_crop(Float32, cell_size, device; kwargs...)
function init_crop(::Type{T},
                   cell_size::Int,
                   device;
                   soil_layers::Int = 5) where {T <: AbstractFloat}

    float_auxiliary() = device(zeros(T, cell_size))
    crop = Crop(
        CropState(
            init_crop_phenology(T, cell_size, device),
            init_crop_canopy_state(T, cell_size, device),
            init_crop_carbon_state(T, cell_size, device),
            init_crop_nitrogen_state(T, cell_size, device),
            init_crop_water_state(T, cell_size, device),
        ),
        CropFluxes(
            init_crop_carbon_fluxes(T, cell_size, device),
            init_crop_nitrogen_fluxes(T, cell_size, device),
            init_crop_water_fluxes(T, cell_size, device; soil_layers = soil_layers),
        ),
        CropAuxiliary(
            init_crop_phenology_auxiliary(T, cell_size, device),
            init_crop_calendar_auxiliary(cell_size, device),
            CropRootAuxiliary(
                float_auxiliary(),
                device(zeros(T, soil_layers)),
            ),
            init_crop_canopy_auxiliary(T, cell_size, device),
            init_crop_photosynthesis_auxiliary(T, cell_size, device),
            CropStressAuxiliary(
                float_auxiliary(), float_auxiliary(), float_auxiliary(),
                float_auxiliary(),
            ),
        ),
        init_crop_events(cell_size, device),
        CropWorkspace(),
    )

    return crop
end
