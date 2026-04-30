function InitilDataLoader(data::NamedTuple, 
                          data_index::Vector{Int},
                          device;
                          training = false,
                          training_by_yield = false
)


    @unpack latitude, crop, soilparam, lpjml = data
    
    latitude_set = latitude[data_index] |> device
    
    crop = (
        sdate = Int32.(crop.sdate[data_index]),
        phu = crop.phu[data_index],
        manure = crop.manure[data_index],
        fertilizer = crop.fertilizer[data_index],
        residuefrac = crop.residuefrac[data_index]
    ) |> device

    soilparams = (
        ph = soilparam.soilph[data_index],
        w_sat = soilparam.w_sat[:, data_index],
        sand = reshape(soilparam.sand[data_index], (1, :)),
        clay = reshape(soilparam.clay[data_index], (1, :)),
        # silt = soilparam.silt[data_index],
        tdiff_0 = soilparam.tdiff_0[data_index],
        tdiff_15 = soilparam.tdiff_15[data_index],
        soildepth = soilparam.soildepth,
    ) |> device
      
    u0_set = (
        swc = lpjml.u0.swc[:, data_index],
        litc = lpjml.u0.litc[:, data_index],
        fastc = lpjml.u0.fastc[:, data_index],
        slowc = lpjml.u0.slowc[:, data_index],
        litn = lpjml.u0.litn[:, data_index],
        fastn = lpjml.u0.fastn[:, data_index],
        slown = lpjml.u0.slown[:, data_index],
        soil_NH4 = lpjml.u0.soil_NH4[:, data_index],
        soil_NO3 = lpjml.u0.soil_NO3[:, data_index],
    ) |> device
    
    if training && training_by_yield
        lpjml = (
            crop = crop,
            c_shift_fast = lpjml.c_shift_fast[:, data_index],
            c_shift_slow = lpjml.c_shift_slow[:, data_index],
            u0 = u0_set,
            output = lpjml.output[:, data_index],   
            output_n = lpjml.output_n[:, data_index],   
            yield = lpjml.yield[:, data_index],   
            μ = lpjml.μ[data_index],
            σ = lpjml.σ[data_index],
            gdhy_yield_n = lpjml.gdhy_yield_n[:, data_index]
        ) |> device
    elseif training
        lpjml = (
            crop = crop,
            c_shift_fast = lpjml.c_shift_fast[:, data_index],
            c_shift_slow = lpjml.c_shift_slow[:, data_index],
            u0 = u0_set,
            output = lpjml.output[:, data_index],   
            output_n = lpjml.output_n[:, data_index],   
            yield = lpjml.yield[:, data_index],   
            μ = lpjml.μ[data_index],
            σ = lpjml.σ[data_index]
        ) |> device
    else
        lpjml = (
            crop = crop,
            c_shift_fast = lpjml.c_shift_fast[:, data_index],
            c_shift_slow = lpjml.c_shift_slow[:, data_index],
            u0 = u0_set,
        ) |> device
    end

    data_set = (
        latitude = latitude_set,
        lpjml = lpjml,
        soilparams = soilparams
    )

    return data_set
end


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
        prec_n = climate.prec_n[:, data_index],
        sw_n = climate.sw_n[:, data_index],
        lw_n = climate.lw_n[:, data_index],
    ) |> device


    return climate
end



function DataLoader(data::NamedTuple, 
                    data_index::Vector{Int},
                    device
)


    @unpack latitude, crop, climate, soilparams, lpjml = data
    
    latitude_set = latitude[data_index] |> device
    
    crop = (
        sdate = Int32.(crop.sdate[:, data_index]),
        phu = crop.phu[data_index],
        manure = crop.manure[data_index],
        fertilizer = crop.fertilizer[data_index],
        residuefrac = crop.residuefrac[data_index]
    ) |> device
    
    # μ_set = (prec = climate.μ.prec[data_index], temp = climate.μ.temp[data_index], swdown = climate.μ.swdown[data_index], lwnet = climate.μ.swdown[data_index])
    # σ_set = (prec = climate.σ.prec[data_index], temp = climate.σ.temp[data_index], swdown = climate.σ.swdown[data_index], lwnet = climate.σ.lwnet[data_index])
    
    climate = (
        temp_spinup = climate.temp_spinup[:, data_index],
        temp = climate.temp[:, data_index],
        prec = climate.prec[:, data_index],
        sw = climate.swdown[:, data_index],
        lw = climate.lwnet[:, data_index],
        co2 = climate.co2,
        temp_n = climate.temp_n[:, data_index],
        prec_n = climate.prec_n[:, data_index],
        sw_n = climate.sw_n[:, data_index],
        lw_n = climate.lw_n[:, data_index],
        # μ = μ_set,
        # σ = σ_set,
    ) |> device

    soilparams = (
        ph = soilparams.soilph[data_index],
        w_sat = soilparams.w_sat[:, data_index],
        sand = reshape(soilparams.sand[data_index], (1, :)),
        clay = reshape(soilparams.clay[data_index], (1, :)),
        # silt = soilparams.silt[data_index],
        tdiff_0 = soilparams.tdiff_0[data_index],
        tdiff_15 = soilparams.tdiff_15[data_index],
        soildepth = soilparams.soildepth,
    ) |> device
    
    vegc_n = (hcat(lpjml.vegc[data_index]...) .- reshape(lpjml.μ.vegc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.vegc[:, data_index], (1, :)) - reshape(lpjml.μ.vegc[:, data_index], (1, :)))
    litc_n = (hcat(lpjml.litc[data_index]...) .- reshape(lpjml.μ.litc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.litc[:, data_index], (1, :)) - reshape(lpjml.μ.litc[:, data_index], (1, :)))
    fastc_n = (hcat(lpjml.fastc'[data_index]...) .- reshape(lpjml.μ.fastc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.fastc[:, data_index], (1, :)) - reshape(lpjml.μ.fastc[:, data_index], (1, :)))
    slowc_n = (hcat(lpjml.slowc'[data_index]...) .- reshape(lpjml.μ.slowc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.slowc[:, data_index], (1, :)) - reshape(lpjml.μ.slowc[:, data_index], (1, :)))
    swc_n = (hcat(lpjml.swc'[data_index]...).- reshape(lpjml.μ.swc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.swc[:, data_index], (1, :)) - reshape(lpjml.μ.swc[:, data_index], (1, :)))
    
    vegc_μ = lpjml.μ.vegc[:, data_index]
    vegc_σ = lpjml.σ.vegc[:, data_index]
    litc_μ = lpjml.μ.litc[:, data_index]
    litc_σ = lpjml.σ.litc[:, data_index]
    fastc_μ = lpjml.μ.fastc[:, data_index]
    fastc_σ = lpjml.σ.fastc[:, data_index]
    slowc_μ = lpjml.μ.slowc[:, data_index]
    slowc_σ = lpjml.σ.slowc[:, data_index]
    swc_μ = lpjml.μ.swc[:, data_index]
    swc_σ = lpjml.σ.swc[:, data_index]
    
    μ_set = (gpp = lpjml.μ.gpp[data_index], lambda = lpjml.μ.lambda[data_index], vmax = lpjml.μ.vmax[data_index], vegc = vegc_μ, resp = lpjml.μ.resp[data_index], swc = swc_μ, litc = litc_μ, fastc = fastc_μ, slowc = slowc_μ) |> device
    σ_set = (gpp = lpjml.σ.gpp[data_index], lambda = lpjml.σ.lambda[data_index], vmax = lpjml.σ.vmax[data_index], vegc = vegc_σ, resp = lpjml.σ.resp[data_index], swc = swc_σ, litc = litc_σ, fastc = fastc_σ, slowc = slowc_σ) |> device
    
    u0_set = (
        swc = lpjml.u0.swc[:, data_index],
        litc = lpjml.u0.litc[:, data_index],
        fastc = lpjml.u0.fastc[:, data_index],
        slowc = lpjml.u0.slowc[:, data_index],
        litn = lpjml.u0.litn[:, data_index],
        fastn = lpjml.u0.fastn[:, data_index],
        slown = lpjml.u0.slown[:, data_index],
        soil_NH4 = lpjml.u0.soil_NH4[:, data_index],
        soil_NO3 = lpjml.u0.soil_NO3[:, data_index],
    ) |> device
    
    lpjml = (
        crop = crop,
        c_shift_fast = lpjml.c_shift_fast[:, data_index],
        c_shift_slow = lpjml.c_shift_slow[:, data_index],
        lambda_n = lpjml.lambda_n[:, data_index],
        vmax_n = lpjml.vmax_n[:, data_index],
        gpp_n = lpjml.gpp_n[:, data_index],
        resp_n = lpjml.resp_n[:, data_index],
        vegc_n = vegc_n,
        litc_n = litc_n[2:end, :],
        fastc_n = fastc_n[2:end, :],
        slowc_n = slowc_n[2:end, :],
        swc_n = swc_n[2:end, :],
        μ = μ_set,
        σ = σ_set,
        u0 = u0_set,
    ) |> device

    data_set = (
        latitude = latitude_set,
        climate = climate,
        lpjml = lpjml,
        soilparams = soilparams
    )

    return data_set
end


function DataLoader_winter_wheat(data::NamedTuple, 
                                 data_index::Vector{Int},
                                 device
)

    @unpack latitude, crop, climate, soilparams, lpjml = data
    
    latitude_set = latitude[data_index] |> device
    
    crop = (
        sdate = Int32.(crop.sdate[:, data_index]),
        phu = crop.phu[data_index],
        manure = crop.manure[data_index],
        fertilizer = crop.fertilizer[data_index],
        residuefrac = crop.residuefrac[data_index],
        hdate = crop.hdate[data_index],
    ) |> device
    
    # μ_set = (prec = climate.μ.prec[data_index], temp = climate.μ.temp[data_index], swdown = climate.μ.swdown[data_index], lwnet = climate.μ.swdown[data_index])
    # σ_set = (prec = climate.σ.prec[data_index], temp = climate.σ.temp[data_index], swdown = climate.σ.swdown[data_index], lwnet = climate.σ.lwnet[data_index])
    
    climate = (
        temp_spinup = climate.temp_spinup[:, data_index],
        temp = climate.temp[:, data_index],
        prec = climate.prec[:, data_index],
        sw = climate.swdown[:, data_index],
        lw = climate.lwnet[:, data_index],
        co2 = climate.co2,
        temp_n = climate.temp_n[:, data_index],
        prec_n = climate.prec_n[:, data_index],
        sw_n = climate.sw_n[:, data_index],
        lw_n = climate.lw_n[:, data_index],
        # μ = μ_set,
        # σ = σ_set,
    ) |> device

    soilparams = (
        ph = soilparams.soilph[data_index],
        w_sat = soilparams.w_sat[:, data_index],
        sand = reshape(soilparams.sand[data_index], (1, :)),
        clay = reshape(soilparams.clay[data_index], (1, :)),
        # silt = soilparams.silt[data_index],
        tdiff_0 = soilparams.tdiff_0[data_index],
        tdiff_15 = soilparams.tdiff_15[data_index],
        soildepth = soilparams.soildepth,
    ) |> device
    
    vegc_n = (hcat(lpjml.vegc[data_index]...) .- reshape(lpjml.μ.vegc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.vegc[:, data_index], (1, :)) - reshape(lpjml.μ.vegc[:, data_index], (1, :)))
    litc_n = (hcat(lpjml.litc[data_index]...) .- reshape(lpjml.μ.litc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.litc[:, data_index], (1, :)) - reshape(lpjml.μ.litc[:, data_index], (1, :)))
    fastc_n = (hcat(lpjml.fastc'[data_index]...) .- reshape(lpjml.μ.fastc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.fastc[:, data_index], (1, :)) - reshape(lpjml.μ.fastc[:, data_index], (1, :)))
    slowc_n = (hcat(lpjml.slowc'[data_index]...) .- reshape(lpjml.μ.slowc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.slowc[:, data_index], (1, :)) - reshape(lpjml.μ.slowc[:, data_index], (1, :)))
    swc_n = (hcat(lpjml.swc'[data_index]...).- reshape(lpjml.μ.swc[:, data_index], (1, :))) ./ (reshape(lpjml.σ.swc[:, data_index], (1, :)) - reshape(lpjml.μ.swc[:, data_index], (1, :)))
    
    litch = reshape(hcat(lpjml.litch[data_index]...), (3, :))
    litnh = reshape(hcat(lpjml.litnh[data_index]...), (3, :))
    
    vegc_μ = lpjml.μ.vegc[:, data_index]
    vegc_σ = lpjml.σ.vegc[:, data_index]
    litc_μ = lpjml.μ.litc[:, data_index]
    litc_σ = lpjml.σ.litc[:, data_index]
    fastc_μ = lpjml.μ.fastc[:, data_index]
    fastc_σ = lpjml.σ.fastc[:, data_index]
    slowc_μ = lpjml.μ.slowc[:, data_index]
    slowc_σ = lpjml.σ.slowc[:, data_index]
    swc_μ = lpjml.μ.swc[:, data_index]
    swc_σ = lpjml.σ.swc[:, data_index]
    
    μ_set = (gpp = lpjml.μ.gpp[data_index], lambda = lpjml.μ.lambda[data_index], vmax = lpjml.μ.vmax[data_index], vegc = vegc_μ, resp = lpjml.μ.resp[data_index], swc = swc_μ, litc = litc_μ, fastc = fastc_μ, slowc = slowc_μ) |> device
    σ_set = (gpp = lpjml.σ.gpp[data_index], lambda = lpjml.σ.lambda[data_index], vmax = lpjml.σ.vmax[data_index], vegc = vegc_σ, resp = lpjml.σ.resp[data_index], swc = swc_σ, litc = litc_σ, fastc = fastc_σ, slowc = slowc_σ) |> device
    
    u0_set = (
        swc = lpjml.u0.swc[:, data_index],
        litc = lpjml.u0.litc[:, data_index],
        fastc = lpjml.u0.fastc[:, data_index],
        slowc = lpjml.u0.slowc[:, data_index],
        litn = lpjml.u0.litn[:, data_index],
        fastn = lpjml.u0.fastn[:, data_index],
        slown = lpjml.u0.slown[:, data_index],
        soil_NH4 = lpjml.u0.soil_NH4[:, data_index],
        soil_NO3 = lpjml.u0.soil_NO3[:, data_index],
    ) |> device
    
    lpjml = (
        crop = crop,
        c_shift_fast = lpjml.c_shift_fast[:, data_index],
        c_shift_slow = lpjml.c_shift_slow[:, data_index],
        lambda_n = lpjml.lambda_n[:, data_index],
        vmax_n = lpjml.vmax_n[:, data_index],
        gpp_n = lpjml.gpp_n[:, data_index],
        resp_n = lpjml.resp_n[:, data_index],
        vegc_n = vegc_n,
        litc_n = litc_n[2:end, :],
        litch = litch,
        litnh = litnh,
        fastc_n = fastc_n[2:end, :],
        slowc_n = slowc_n[2:end, :],
        swc_n = swc_n[2:end, :],
        μ = μ_set,
        σ = σ_set,
        u0 = u0_set,
    ) |> device

    data_set = (
        latitude = latitude_set,
        climate = climate,
        lpjml = lpjml,
        soilparams = soilparams
    )

    return data_set
end