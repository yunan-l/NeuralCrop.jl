function daily_crop_C3_training!(day_start,
                                 day_end,
                                 model,
                                 ps,
                                 ps_frozen,
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
                                 device
)

    @unpack cft, lpjmlparam, photopar, k = parameters
    @unpack latitude, climate, lpjml = data_set

    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        temp, prec, swr, lwr, temp_n, swr_n, lwr_n, co2 = readclimate!(climate, day)

        # initial crop variables in sowing day and fertilizer
        cultivate!(crop, crop_cal, lpjml.crop.sdate, lpjmlparam, managed_land, soil, day_of_year, device)

        Zygote.ignore() do
            update_climbuf!(cft, temp, climate.temp, climbuf, day, device) # update climate buffer
            albedo!(cft, crop, pet.albedo)  # compute albedo
            petpar!(pet, day_of_year, latitude, temp, lwr, swr) # compute crop potential evapotraspiration variables
            soiltemp_lag!(soil, climbuf, device)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same
        end

        # compute phenology variables
        Zygote.ignore() do
            phenology_crop!(crop, climbuf.V_req, cft, temp, pet.daylength)
        end
        
        harvest_crop!(crop_cal, crop, soil, output, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        Zygote.ignore() do
            apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
            temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function
        end

        # C3 photosynthesis
        photosynthesis_C3!(cft, lpjmlparam, photopar, photos, crop.apar, pet.daylength, temp, co2; comp_vmax = true)

        # crop respiration and carbon allocation
        crop_carbon_hybrid!(model.stoc, ps.ps_stoc, st.st_stoc, photos, crop, cft, lpjmlparam, temp, temp_n, soil.swc./soil.layer_depth)

        Zygote.ignore() do
            crop_nitrogen!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle
            # evapotranspiration
            interception!(crop, lpjmlparam, cft, pet.eeq, prec)
            pedotransfer!(soil, lpjmlparam)
            transpiration!(photos.adtmm, lpjmlparam, cft, crop, pet, soil, co2)
            evaporation!(pet.eeq, lpjmlparam, crop, soil)
        end

        # soil carbon cycle
        soil_carbon!(model, ps_frozen, st, lpjmlparam, temp_n, swr_n, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(model, ps_frozen, st, temp_n, swr_n, crop_cal, soil)

        # soil water cycle
        soil_water!(model.swc, ps_frozen.ps_swc, st.st_swc, lpjmlparam, soil, crop, prec, swr_n, lwr_n)

    end
        
end