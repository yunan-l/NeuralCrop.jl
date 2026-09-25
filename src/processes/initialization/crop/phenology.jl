"""Integrated phenology coordinates and discrete prognostic modes."""
mutable struct CropPhenology{A, B, I}
    vdsum::A               # Accumulated effective vernalization days (day equivalent).
    husum::A               # Accumulated heat units since cultivation (°C day).
    senescence::B          # Current phenological senescence-mode flag.
    senescence_previous::B # Previous-day senescence-mode flag used at transitions.
    harvesting::B          # Current phenological harvest-readiness flag.
    harvesting_previous::B # Previous-day harvest-readiness flag used to detect harvest.
    growing_days::I        # Number of simulated days since cultivation (day).
    is_growing::I          # Active crop-presence/growth mode (0/1).
end

"""Static and current-day algebraically derived phenology variables."""
mutable struct CropPhenologyAuxiliary{A, B}
    phu::A         # Potential heat units required for maturity (°C day).
    winter_type::B # Winter-crop/vernalization requirement flag.
    fphu::A        # Current heat-unit fraction derived from `husum / phu` (0–1).
end

init_crop_phenology(cell_size::Int, device) = init_crop_phenology(Float32, cell_size, device)
function init_crop_phenology(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    bool_state(value = false) = device(fill(value, cell_size))

    return CropPhenology(
        float_state(),
        float_state(),
        bool_state(),
        bool_state(),
        bool_state(true),
        bool_state(true),
        device(zeros(Int32, cell_size)),
        device(zeros(Int32, cell_size)),
    )
end

function init_crop_phenology_auxiliary(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    return CropPhenologyAuxiliary(
        device(zeros(T, cell_size)),
        device(fill(false, cell_size)),
        device(zeros(T, cell_size)),
    )
end
