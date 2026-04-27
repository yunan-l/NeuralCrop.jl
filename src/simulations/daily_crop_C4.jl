### Hybrid modelling
function daily_crop_C4!(day_start,
                        day_end,
                        model,
                        ps, 
                        st,
                        pftparameters,
                        data_set,
                        cell_size,
                        climbuf, 
                        crop, 
                        crop_cal, 
                        photos, 
                        pet, 
                        soil, 
                        managed_land,
                        dailyWeather, 
                        output,
                        device;
                        maize = false
)

    @unpack latitude, climate, lpjml = data_set

    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        readclimate!(climate, dailyWeather, day)

        # sonw
        snow!(soil, dailyWeather)

        # initial crop variables in sowing day and fertilizer
        cultivate!(crop, crop_cal, managed_land, soil, day_of_year, device)

        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day, device) # update climate buffer
        albedo!(pftparameters, crop, pet.albedo)  # compute albedo
        petpar!(pet, day_of_year, latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf, device)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)
        
        harvest_crop!(crop_cal, crop, soil, output, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        if maize
            apar_crop_maize!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        else
            apar_crop!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        end

        temp_stress(pftparameters, pet, photos, dailyWeather.temp) # temperature stress function

        # C4 photosynthesis
        hybrid_photos_C4!(model, ps, st, pftparameters, photos, crop, pet.daylength, soil.swc./soil.layer_depth, dailyWeather.temp_n, temp)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, pftparameters, dailyWeather.temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, pftparameters, soil, photos.vmax, pet.daylength, dailyWeather.temp) # nitrogen cycle         
        
        # evapotranspiration
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec)
        pedotransfer!(soil)
        transpiration!(photos.adtmm, pftparameters, crop, pet, soil, dailyWeather.annual_co2)
        evaporation!(pet.eeq, crop, soil)

        # soil carbon cycle
        soil_carbon!(model, ps, st, dailyWeather.temp_n, dailyWeather.swr_n, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(model, ps, st, dailyWeather.temp_n, dailyWeather.swr_n, crop_cal, soil)

        # soil water cycle
        soil_water!(model.swc, ps.swc, st.swc, soil, crop, dailyWeather.prec, dailyWeather.swr_n, dailyWeather.lwr_n)

        # # output
        # output_yield!(output, crop, day_of_year)

    end
end



### Purely process-based modelling
function daily_crop_C4!(day_start,
                        day_end,
                        pftparameters,
                        data_set,
                        cell_size,
                        climbuf, 
                        crop, 
                        crop_cal, 
                        photos, 
                        pet, 
                        soil, 
                        managed_land, 
                        dailyWeather,
                        output,
                        device;
                        maize = false
)

    @unpack latitude, climate, lpjml = data_set

    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        readclimate!(climate, dailyWeather, day)

        # sonw
        snow!(soil, dailyWeather)

        # initial crop variables in sowing day and fertilizer
        cultivate!(crop, crop_cal, managed_land, soil, day_of_year, device)

        update_climbuf!(pftparameters, dailyWeather.temp, climbuf, day, device) # update climate buffer
        albedo!(pftparameters, crop, pet.albedo)  # compute albedo
        petpar!(pet, day_of_year, latitude, dailyWeather.temp, dailyWeather.lwr, dailyWeather.swr) # compute crop potential evapotraspiration variables
        soiltemp_lag!(soil, climbuf, device)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same

        # compute phenology variables
        phenology_crop!(crop, climbuf.V_req, pftparameters, dailyWeather.temp, pet.daylength)
        
        harvest_crop!(crop_cal, crop, soil, output, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        if maize
            apar_crop_maize!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        else
            apar_crop!(pftparameters, crop, pet) # crop absorbed photosynthetic radiation
        end

        temp_stress(pftparameters, pet, photos, dailyWeather.temp) # temperature stress function

        # C4 photosynthesis
        photosynthesis_C4!(pftparameters, photos, crop.apar, pet.daylength, dailyWeather.temp; comp_vmax = true)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, pftparameters, dailyWeather.temp)

        # crop nitrogen allocation
        crop_nitrogen!(crop, pftparameters, soil, photos.vmax, pet.daylength, dailyWeather.temp) # nitrogen cycle         
        
        # evapotranspiration
        interception!(crop, pftparameters, pet.eeq, dailyWeather.prec)
        pedotransfer!(soil)
        transpiration!(photos.adtmm, pftparameters, crop, pet, soil, dailyWeather.annual_co2)
        evaporation!(pet.eeq, crop, soil)

        # soil carbon cycle
        soil_carbon!(crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(crop_cal, soil)

        # soil water cycle
        soil_water!(soil, crop, dailyWeather.prec)

        # # output
        # output_yield!(output, crop, day_of_year)

    end
end