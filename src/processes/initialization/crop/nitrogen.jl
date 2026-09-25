"""Persistent plant nitrogen stocks and seasonal memory."""
mutable struct CropNitrogenState{A}
    total::A              # Total living crop nitrogen stock (gN m⁻²).
    leaf::A               # Leaf nitrogen stock (gN m⁻²).
    root::A               # Root nitrogen stock (gN m⁻²).
    pool::A               # Mobile/intermediate plant nitrogen pool (gN m⁻²).
    storage::A            # Harvestable storage-organ nitrogen stock (gN m⁻²).
    pending_manure::A     # Prescribed manure N awaiting a scheduled split application (gN m⁻²).
    pending_fertilizer::A # Mineral fertilizer N awaiting a scheduled split application (gN m⁻²).
    stress_sum::A         # Seasonal sum of the daily nitrogen-sufficiency factor (dimensionless day).
    sufficiency::A        # Prior-day N-sufficiency multiplier needed by next-day canopy growth (0–1).
end

"""Current-day plant and management nitrogen fluxes."""
mutable struct CropNitrogenFluxes{A}
    uptake::A                      # Total mineral-N transfer from soil to crop (gN m⁻² day⁻¹).
    auto_fertilizer::A             # Portion of uptake supplied by automatic fertilizer (gN m⁻² day⁻¹).
    seed_input::A                  # Nitrogen introduced with seed at cultivation (gN m⁻² day⁻¹).
    prescribed_manure_input::A     # Prescribed manure N applied today (gN m⁻² day⁻¹).
    prescribed_fertilizer_input::A # Prescribed mineral fertilizer N applied today (gN m⁻² day⁻¹).
    harvest_export::A              # Plant nitrogen removed at harvest (gN m⁻² day⁻¹).
end

function init_crop_nitrogen_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropNitrogenState(
        ntuple(_ -> float_state(), 8)...,
        device(ones(T, cell_size)),
    )
end

function init_crop_nitrogen_fluxes(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropNitrogenFluxes(ntuple(_ -> float_flux(), 6)...)
end
