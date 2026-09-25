"""
crop_carbon!(photos, CFT, crop, pet, soil, temp, co2)

Run the daily crop carbon process chain: respiration, allocation, and phenology coupling.
"""
function crop_carbon!(crop,
                      output::Output,
                      CFT::CFTParameters,
                      air_temperature::AbstractVector{T},
                      soil_temperature::AbstractMatrix{T};
                      output_row::Union{Nothing, Integer} = nothing,
                      crop_resp_fix::Bool = false,
                      lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(
        crop, CFT, air_temperature, soil_temperature,
        crop_fluxes(crop).carbon.gross_assimilation,
        crop_fluxes(crop).carbon.leaf_respiration;
        crop_resp_fix,
        lpjmlparams = lpjmlparams,
    )

    # compute crop carbon allocation
    carbon_allocation!(CFT, crop)

    sources = (
        gpp = crop_fluxes(crop).carbon.gross_assimilation,
        npp = crop_fluxes(crop).carbon.npp,
        lambda = crop_photosynthesis_auxiliary(crop).lambda,
        potential_vcmax = crop_photosynthesis_auxiliary(crop).potential_vcmax,
        vcmax = crop_photosynthesis_auxiliary(crop).vcmax,
        nitrogen_limitation = crop_photosynthesis_auxiliary(crop).nitrogen_limitation,
        respiration = crop_fluxes(crop).carbon.respiration,
        lai = crop_canopy_auxiliary(crop).actual_lai,
        fphu = crop_phenology_auxiliary(crop).fphu,
        water_deficit = crop_stress_auxiliary(crop).water_deficit,
        biomass = crop_prognostic(crop).carbon.biomass,
    )
    for (field, source) in pairs(sources)
        if output_row === nothing
            setproperty!(
                output.crop,
                field,
                _append_output_row(getproperty(output.crop, field), source),
            )
        else
            _write_output_row!(getproperty(output.crop, field), output_row, source)
        end
    end

end
