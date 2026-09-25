"""Current-day crop photosynthetic capacities and limiting factors."""
mutable struct CropPhotosynthesisAuxiliary{A}
    potential_vcmax::A    # Potential Rubisco carboxylation capacity before N limitation (gC m⁻² day⁻¹).
    vcmax::A              # Active Rubisco carboxylation capacity (gC m⁻² day⁻¹).
    nitrogen_limitation::A # Retained fraction of potential Vcmax after N limitation (0–1).
    lambda::A             # Ratio of intercellular to ambient CO₂ partial pressure (0–1).
    temperature_stress::A # Photosynthetic temperature-response multiplier (0–1).
end

function init_crop_photosynthesis_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_auxiliary() = device(zeros(T, cell_size))
    return CropPhotosynthesisAuxiliary(ntuple(_ -> float_auxiliary(), 5)...)
end
