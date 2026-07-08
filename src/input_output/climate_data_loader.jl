"""
ClimateDataLoader(climate, data_index, device)

Extract climate forcing slices for selected grid points and years.
"""
function ClimateDataLoader(climate::NamedTuple, 
                           data_index::Vector{Int},
                           device
)

    climate = (
        temp_spinup = climate.temp_spinup[:, data_index],
        temp = climate.temp[:, data_index],
        prec = climate.prec[:, data_index],
        sw = climate.swdown[:, data_index],
        lw = climate.lwnet[:, data_index],
        co2 = climate.co2,
        temp_n = climate.temp_n[:, data_index],
        tmax_n = climate.tmax_n[:, data_index],
        tmin_n = climate.tmin_n[:, data_index],
        prec_n = climate.prec_n[:, data_index],
        sw_n = climate.sw_n[:, data_index],
        lw_n = climate.lw_n[:, data_index],
        vpd_n = climate.vpd_n[:, data_index],
    ) |> device

    return climate
end
