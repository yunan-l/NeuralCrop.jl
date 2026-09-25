"""Shared soil C/N decomposition responses, fixed routing configuration, and workspace."""
mutable struct SoilDecomposition{A, M}
    response::M          # Combined temperature–moisture response by soil layer (0–1+).
    litter_response::M   # Combined response by litter class and grid cell (0–1+).
    shift_fast::M        # Fixed post-spin-up vertical routing fractions for fast C/N inputs (sum to 1 by cell).
    shift_slow::M        # Fixed post-spin-up vertical routing fractions for slow C/N inputs (sum to 1 by cell).
    layer_scratch_1::M   # Preallocated layer-wise temporary buffer; no scientific state.
    layer_scratch_2::M   # Second preallocated layer-wise temporary buffer; no scientific state.
    surface_scratch_1::A # Preallocated grid-cell temporary buffer; no scientific state.
    surface_scratch_2::A # Second preallocated grid-cell temporary buffer; no scientific state.
end

init_soil_decomposition(cell_size::Int, device; kwargs...) =
    init_soil_decomposition(Float32, cell_size, device; kwargs...)
function init_soil_decomposition(::Type{T}, cell_size::Int, device;
                                 soil_layers::Int = 5,
                                 litter_layers::Int = 3) where {T <: AbstractFloat}
    layer_state() = device(zeros(T, soil_layers, cell_size))
    return SoilDecomposition(
        layer_state(),
        device(zeros(T, litter_layers, cell_size)),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        device(zeros(T, cell_size)),
        device(zeros(T, cell_size)),
    )
end
