"""Apply the sowing-day reduction in topsoil bulk density caused by tillage."""

function tillage_hydraulics!(soil, crop;
                             lpjmlparams::LPJmLParams = lpjmlparams)
    launch_1D!(
        tillage_hydraulics_kernel!,
        soil_management_prognostic(soil).tillage_density_factor,
        crop_events(crop).sowing,
        kernel_constant(
            soil_management_prognostic(soil).tillage_density_factor, lpjmlparams,
        ),
    )
    return nothing
end

function tillage_hydraulics!(state::ModelState;
                             lpjmlparams::LPJmLParams = lpjmlparams)
    density_factor = state.prognostic.soil.management.tillage_density_factor
    launch_1D!(
        tillage_hydraulics_kernel!, density_factor, state.events.crop.sowing,
        kernel_constant(density_factor, lpjmlparams),
    )
    return nothing
end

@kernel inbounds = true function tillage_hydraulics_kernel!(
    density_factor::AbstractMatrix{T},
    sowing_event::AbstractVector{S},
    fixed_parameters,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    lpjmlparams = kernel_value(fixed_parameters)
    if sowing_event[cell] != 0
        density_factor[1, cell] -=
            (density_factor[1, cell] - T(0.667)) * T(lpjmlparams.mixing_efficiency)
    end
end
