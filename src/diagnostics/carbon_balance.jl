"""
    CarbonBalance

Optional daily carbon-budget diagnostics. All fields have shape
`(number_of_days, number_of_cells)` and use the model carbon unit (g C m⁻²).

The tracked stock is crop organ carbon plus litter, fast-soil, and slow-soil
carbon. NPP, seed carbon, and manure carbon are boundary inputs; harvest
export and heterotrophic respiration are boundary outputs. Residue transfer
is recorded for inspection but is internal to the tracked system.
"""
mutable struct CarbonBalance{M <: AbstractArray{<:AbstractFloat}}
    plant_before::M              # Plant C stock at start of day (gC m⁻²).
    plant_after::M               # Plant C stock at end of day (gC m⁻²).
    soil_before::M               # Litter plus soil C stock at start of day (gC m⁻²).
    soil_after::M                # Litter plus soil C stock at end of day (gC m⁻²).
    total_before::M              # Tracked ecosystem C stock at start of day (gC m⁻²).
    total_after::M               # Tracked ecosystem C stock at end of day (gC m⁻²).
    net_primary_production::M    # Plant NPP boundary input (gC m⁻² day⁻¹).
    seed_input::M                # Seed-carbon boundary input (gC m⁻² day⁻¹).
    manure_input::M              # Manure-carbon boundary input (gC m⁻² day⁻¹).
    residue_transfer::M          # Internal crop-to-litter transfer (gC m⁻² day⁻¹).
    harvest_export::M            # Harvested carbon leaving the system (gC m⁻² day⁻¹).
    heterotrophic_respiration::M # Total microbial respiration loss (gC m⁻² day⁻¹).
    litter_respiration::M        # Respiration attributed to litter decomposition (gC m⁻² day⁻¹).
    fast_pool_respiration::M     # Respiration attributed to fast SOC (gC m⁻² day⁻¹).
    slow_pool_respiration::M     # Respiration attributed to slow SOC (gC m⁻² day⁻¹).
    residual::M                  # Absolute daily carbon-budget closure error (gC m⁻²).
    relative_residual::M         # Carbon residual normalized by daily budget magnitude.
end

function init_carbon_balance(number_of_days::Integer,
                             number_of_cells::Integer,
                             device = identity;
                             T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))
    return CarbonBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(),
    )
end

function crop_carbon_stock(crop)
    return crop_prognostic(crop).carbon.leaf .+ crop_prognostic(crop).carbon.root .+
           crop_prognostic(crop).carbon.pool .+ crop_prognostic(crop).carbon.storage
end

function soil_carbon_stock(soil)
    return vec(sum(soil_carbon_prognostic(soil).litter; dims = 1)) .+
           vec(sum(soil_carbon_prognostic(soil).fast .+ soil_carbon_prognostic(soil).slow; dims = 1))
end

function record_carbon_balance_start!(balance::CarbonBalance,
                                      day_index::Integer,
                                      crop,
                                      soil)
    @views begin
        balance.plant_before[day_index, :] .= crop_carbon_stock(crop)
        balance.soil_before[day_index, :] .= soil_carbon_stock(soil)
        balance.total_before[day_index, :] .=
            balance.plant_before[day_index, :] .+
            balance.soil_before[day_index, :]
    end
    return nothing
end

function record_carbon_balance_after_cultivate!(balance::CarbonBalance,
                                                day_index::Integer,
                                                crop;
                                                lpjmlparams::LPJmLParams = lpjmlparams)
    @views begin
        balance.seed_input[day_index, :] .=
            max.(crop_carbon_stock(crop) .-
                 balance.plant_before[day_index, :], zero(eltype(balance.seed_input))) .*
            crop_events(crop).sowing
        balance.manure_input[day_index, :] .=
            crop_fluxes(crop).nitrogen.prescribed_manure_input .* lpjmlparams.manure_cn
    end
    return nothing
end

function record_carbon_balance_after_harvest!(balance::CarbonBalance,
                                              day_index::Integer,
                                              crop,
                                              soil,
                                              residue_fraction)
    event = crop_events(crop).harvest
    @views begin
        balance.residue_transfer[day_index, :] .=
            vec(sum(soil_carbon_fluxes(soil).input; dims = 1))
        balance.harvest_export[day_index, :] .=
            crop_fluxes(crop).carbon.harvest_export
    end
    return nothing
end

function record_carbon_balance_end!(balance::CarbonBalance,
                                    day_index::Integer,
                                    crop,
                                    soil;
                                    lpjmlparams::LPJmLParams = lpjmlparams)
    @views begin
        balance.plant_after[day_index, :] .= crop_carbon_stock(crop)
        balance.soil_after[day_index, :] .= soil_carbon_stock(soil)
        balance.total_after[day_index, :] .=
            balance.plant_after[day_index, :] .+
            balance.soil_after[day_index, :]
        balance.net_primary_production[day_index, :] .= crop_fluxes(crop).carbon.npp
        balance.heterotrophic_respiration[day_index, :] .=
            soil_carbon_fluxes(soil).heterotrophic_respiration
        balance.litter_respiration[day_index, :] .=
            vec(sum(soil_carbon_fluxes(soil).decomposed_litter; dims = 1)) .* lpjmlparams.atmfrac
        balance.fast_pool_respiration[day_index, :] .=
            vec(sum(soil_carbon_fluxes(soil).decomposed_fast; dims = 1))
        balance.slow_pool_respiration[day_index, :] .=
            vec(sum(soil_carbon_fluxes(soil).decomposed_slow; dims = 1))

        balance.residual[day_index, :] .=
            balance.total_before[day_index, :] .+
            balance.seed_input[day_index, :] .+
            balance.manure_input[day_index, :] .+
            balance.net_primary_production[day_index, :] .-
            balance.harvest_export[day_index, :] .-
            balance.heterotrophic_respiration[day_index, :] .-
            balance.total_after[day_index, :]
        balance.relative_residual[day_index, :] .=
            balance.residual[day_index, :] ./ max.(
                abs.(balance.total_before[day_index, :]) .+
                balance.seed_input[day_index, :] .+
                balance.manure_input[day_index, :] .+
                abs.(balance.net_primary_production[day_index, :]),
                eps(eltype(balance.residual)),
            )
    end
    return nothing
end
