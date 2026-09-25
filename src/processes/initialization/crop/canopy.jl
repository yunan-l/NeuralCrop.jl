"""Checkpointed canopy variables needed by the next daily transition."""
mutable struct CropCanopyState{A}
    lai::A                    # Phenological leaf-area index (m² leaf m⁻² ground).
    lai_previous_potential::A # Previous potential LAI (`lai000` in LPJmL).
    laimax_adjusted::A        # Maximum LAI retained at onset of senescence (m² m⁻²).
    lai_npp_deficit::A        # LAI not supported by available plant carbon (m² m⁻²).
end

"""Current-day canopy radiation, geometry, and conductance diagnostics."""
mutable struct CropCanopyAuxiliary{A}
    actual_lai::A           # Carbon-supported LAI after NPP-deficit correction (m² m⁻²).
    flaimax::A              # Fraction of maximum LAI reached by current-day phenology (0–1).
    albedo::A               # Effective crop-covered surface albedo (fraction, 0–1).
    fpar::A                 # Fraction of incident PAR absorbed by the crop canopy (0–1).
    apar::A                 # Absorbed photosynthetically active radiation (J m⁻² day⁻¹).
    canopy_conductance::A   # Bulk canopy conductance used by transpiration (mm s⁻¹).
    canopy_wet::A           # Fraction of daily evaporative demand used by wet-canopy evaporation (0–1).
end

function init_crop_canopy_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropCanopyState(ntuple(_ -> float_state(), 4)...)
end

function init_crop_canopy_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_auxiliary() = device(zeros(T, cell_size))
    return CropCanopyAuxiliary(ntuple(_ -> float_auxiliary(), 7)...)
end
