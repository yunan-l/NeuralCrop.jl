### Hybrid modelling
function daily_crop_C4!(day_start,
                        day_end,
                        model,
                        ps, 
                        st,
                        parameters,
                        data_set,
                        cell_size,
                        climbuf, 
                        crop, 
                        crop_cal, 
                        photos, 
                        pet, 
                        soil, 
                        managed_land, 
                        output,
                        device;
                        maize = false
)

    @unpack cft, lpjmlparam, photopar, k = parameters
    @unpack latitude, climate, lpjml = data_set

    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        temp, prec, swr, lwr, temp_n, prec_n, swr_n, lwr_n, co2 = readclimate!(climate, day)

        # initial crop variables in sowing day and fertilizer
        cultivate!(crop, crop_cal, lpjml.crop.sdate, lpjmlparam, managed_land, soil, day_of_year, device)

        update_climbuf!(cft, temp, climate.temp, climbuf, day, device) # update climate buffer
        albedo!(cft, crop, pet.albedo)  # compute albedo
        petpar!(pet, day_of_year, latitude, temp, lwr, swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf, device)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, cft, temp, pet.daylength)
        
        harvest_crop!(crop_cal, crop, soil, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        if maize
            apar_crop_maize!(cft, crop, pet) # crop absorbed photosynthetic radiation
        else
            apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
        end

        temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function

        # C4 photosynthesis
        hybrid_photos_C4!(model, ps, st, lpjmlparam, cft, photopar, photos, crop, pet.daylength, soil.swc./soil.layer_depth, temp_n, temp)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, cft, lpjmlparam, temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle         
        
        # evapotranspiration
        interception!(crop, lpjmlparam, cft, pet.eeq, prec)
        pedotransfer!(soil, lpjmlparam)
        transpiration!(photos.adtmm, lpjmlparam, cft, crop, pet, soil, co2)
        evaporation!(pet.eeq, lpjmlparam, crop, soil)

        # soil carbon cycle
        soil_carbon!(model, ps, st, lpjmlparam, temp_n, swr_n, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(model, ps, st, temp_n, swr_n, crop_cal, soil)

        # soil water cycle
        soil_water!(model.swc, ps.ps_swc, st.st_swc, lpjmlparam, soil, crop, prec, swr_n, lwr_n)

        # output
        output_yield!(output, crop, day_of_year)

    end

end


### Purely process-based modelling
function daily_crop_C4!(day_start,
                        day_end,
                        parameters,
                        data_set,
                        cell_size,
                        climbuf, 
                        crop, 
                        crop_cal, 
                        photos, 
                        pet, 
                        soil, 
                        managed_land, 
                        output,
                        device;
                        maize = false
)

    @unpack cft, lpjmlparam, photopar, k = parameters
    @unpack latitude, climate, lpjml = data_set

    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        temp, prec, swr, lwr, temp_n, prec_n, swr_n, lwr_n, co2 = readclimate!(climate, day)

        # initial crop variables in sowing day and fertilizer
        cultivate!(crop, crop_cal, lpjml.crop.sdate, lpjmlparam, managed_land, soil, day_of_year, device)

        update_climbuf!(cft, temp, climate.temp, climbuf, day, device) # update climate buffer
        albedo!(cft, crop, pet.albedo)  # compute albedo
        petpar!(pet, day_of_year, latitude, temp, lwr, swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf, device)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, cft, temp, pet.daylength)
        
        harvest_crop!(crop_cal, crop, soil, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        if maize
            apar_crop_maize!(cft, crop, pet) # crop absorbed photosynthetic radiation
        else
            apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
        end

        temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function

        # C4 photosynthesis
        photosynthesis_C4!(cft, lpjmlparam, photopar, photos, crop.apar, pet.daylength, temp; comp_vmax = true)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, cft, lpjmlparam, temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle         
        
        # evapotranspiration
        interception!(crop, lpjmlparam, cft, pet.eeq, prec)
        pedotransfer!(soil, lpjmlparam)
        transpiration!(photos.adtmm, lpjmlparam, cft, crop, pet, soil, co2)
        evaporation!(pet.eeq, lpjmlparam, crop, soil)

        # soil carbon cycle
        soil_carbon!(lpjmlparam, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(crop_cal, soil)

        # soil water cycle
        soil_water!(lpjmlparam, soil, crop, prec)

        # output
        output_yield!(output, crop, day_of_year)

    end

end