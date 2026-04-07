function daily_crop_C3!(start_day,
                        end_day,
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
                        device
)

    @unpack cft, lpjmlparam, photopar, k = parameters
    @unpack latitude, climate, lpjml = data_set

    for day = start_day : end_day

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
        
        apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
        temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function

        # C3 photosynthesis
        hybrid_photos_C3!(model, ps, st, lpjmlparam, cft, photopar, photos, crop, pet.daylength, soil.swc./soil.layer_depth, temp_n, temp, co2; comp_vmax = true)

        # crop respiration and carbon allocation
        crop_carbon!(photos, crop, cft, lpjmlparam, temp)

        # crop nitrogen allocation
        # crop_nitrogen_old!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle
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
        # output_predict!(output, photos, crop, soil, day_)
        output_yield!(output, crop, day_of_year)   

    end

end


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
        hybrid_photos_C4!(model, ps, st, lpjmlparam, cft, photopar, photos, crop, pet.daylength, soil.swc./soil.layer_depth, temp_n, temp; comp_vmax = true)

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


function daily_crop_C3_training_rollout!(day_start,
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
                                         device
)

    @unpack cft, lpjmlparam, photopar, k = parameters
    @unpack latitude, climate, lpjml = data_set


    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        temp, prec, swr, lwr, temp_n, prec_n, swr_n, lwr_n, co2 = readclimate!(climate, day)

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
        
        harvest_crop!(crop_cal, crop, soil, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        Zygote.ignore() do
            apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
            temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function
        end

        # C3 photosynthesis
        hybrid_photos_C3!(model, ps, st, lpjmlparam, cft, photopar, photos, crop, pet.daylength, soil.swc./soil.layer_depth, temp_n, temp, co2)

        # crop respiration and carbon allocation
        crop_carbon!(model.vegc, ps.ps_vegc, st.st_vegc, photos, crop, cft, lpjmlparam, temp, temp_n, soil.swc./soil.layer_depth)

        Zygote.ignore() do
            crop_nitrogen!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle
            # evapotranspiration
            interception!(crop, lpjmlparam, cft, pet.eeq, prec)
            pedotransfer!(soil, lpjmlparam)
            transpiration!(photos.adtmm, lpjmlparam, cft, crop, pet, soil, co2)
            evaporation!(pet.eeq, lpjmlparam, crop, soil)
        end

        # soil carbon cycle
        soil_carbon!(model, ps, st, lpjmlparam, temp_n, swr_n, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(model, ps, st, temp_n, swr_n, crop_cal, soil)

        # soil water cycle
        soil_water!(model.swc, ps.ps_swc, st.st_swc, lpjmlparam, soil, crop, prec, swr_n, lwr_n)

        # output
        output_training!(output, photos, crop, soil, lpjml.μ, lpjml.σ)

    end

    return output
        
end

function daily_crop_winter_wheat_training_rollout!(day_start,
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
                                                   device
)

    @unpack cft, lpjmlparam, photopar, k = parameters
    @unpack latitude, climate, lpjml = data_set

    for day = day_start : day_end

        day_of_year = day % 365 != 0 ? day % 365 : 365

        temp, prec, swr, lwr, temp_n, prec_n, swr_n, lwr_n, co2 = readclimate!(climate, day)

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
        
        harvest_crop!(crop_cal, crop, soil, lpjml.crop.residuefrac, device, cell_size, day_of_year) # crop harvesting
        
        Zygote.ignore() do
            apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
            temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function
        end

        # C3 photosynthesis
        hybrid_photos_C3!(model, ps, st, lpjmlparam, cft, photopar, photos, crop, pet.daylength, soil.swc./soil.layer_depth, temp_n, temp, co2)

        # crop respiration and carbon allocation
        crop_carbon!(model.vegc, ps.ps_vegc, st.st_vegc, photos, crop, cft, lpjmlparam, temp, temp_n, soil.swc./soil.layer_depth)

        Zygote.ignore() do
            crop_nitrogen!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle
            # evapotranspiration
            interception!(crop, lpjmlparam, cft, pet.eeq, prec)
            pedotransfer!(soil, lpjmlparam)
            transpiration!(photos.adtmm, lpjmlparam, cft, crop, pet, soil, co2)
            evaporation!(pet.eeq, lpjmlparam, crop, soil)
        end

        # soil carbon cycle
        soil_carbon!(model, ps, st, lpjmlparam, temp_n, swr_n, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(model, ps, st, temp_n, swr_n, crop_cal, soil)

        # soil water cycle
        soil_water!(model.swc, ps.ps_swc, st.st_swc, lpjmlparam, soil, crop, prec, swr_n, lwr_n)

        update_lit_winter_wheat!(soil, lpjml.litch, lpjml.litnh, crop.wtype, lpjml.crop.hdate, crop_cal.hcallback, day)

        # output
        output_training!(output, photos, crop, soil, lpjml.μ, lpjml.σ)

    end

    return output
        
end

function daily_crop_C4_training_rollout!(day_start,
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

        day_ = day % 365 != 0 ? day % 365 : 365

        temp, prec, swr, lwr, temp_n, prec_n, swr_n, lwr_n, co2 = readclimate!(climate, day)

        # initial crop variables in sowing day and fertilizer
        cultivate!(crop, crop_cal, lpjml.crop.sdate, lpjmlparam, managed_land, soil, day, device)

        Zygote.ignore() do
            update_climbuf!(cft, temp, climate.temp, climbuf, day, device) # update climate buffer
            # daily_climbuf!(temp, climbuf.temp) # daily climate buf update
            albedo!(cft, crop, pet.albedo)  # compute albedo
            petpar!(pet, day_, latitude, temp, lwr, swr) # compute crop potential evapotraspiration variables
            soiltemp_lag!(soil, climbuf, device)  # compute soil temperature, using very siample linear method, now the five soil-layer temperature is same
        end

        # compute phenology variables
        Zygote.ignore() do
            phenology_crop!(crop, climbuf.V_req, cft, temp, pet.daylength)
        end
        
        harvest_crop!(crop_cal, crop, soil, lpjml.crop.residuefrac, device, cell_size, day_) # crop harvesting
        
        Zygote.ignore() do
            if maize
                apar_crop_maize!(cft, crop, pet) # crop absorbed photosynthetic radiation
            else
                apar_crop!(cft, crop, pet) # crop absorbed photosynthetic radiation
            end
            temp_stress(cft, photopar, pet, photos, k, temp) # temperature stress function
        end

        # C4 photosynthesis
        hybrid_photos_C4!(model, ps, st, lpjmlparam, cft, photopar, photos, crop, pet.daylength, soil.swc./soil.layer_depth, temp_n, temp)

        # crop respiration and carbon allocation
        crop_carbon!(model.vegc, ps.ps_vegc, st.st_vegc, photos, crop, cft, lpjmlparam, temp, temp_n, soil.swc./soil.layer_depth)
        # Zygote.ignore() do
        #     crop_carbon_old!(photos, crop, cft, lpjmlparam, temp)
        # end

        Zygote.ignore() do
            crop_nitrogen!(crop, cft, soil, lpjmlparam, photos.vmax, pet.daylength, temp) # nitrogen cycle 
            # evapotranspiration
            interception!(crop, lpjmlparam, cft, pet.eeq, prec)
            pedotransfer!(soil, lpjmlparam)
            transpiration!(photos.adtmm, lpjmlparam, cft, crop, pet, soil, co2)
            evaporation!(pet.eeq, lpjmlparam, crop, soil)
        end

        # soil carbon cycle
        soil_carbon!(model, ps, st, lpjmlparam, temp_n, swr_n, crop_cal, soil)

        # soil nitrogen cycle
        soil_nitrogen!(model, ps, st, temp_n, swr_n, crop_cal, soil)

        # soil water cycle
        soil_water!(model.swc, ps.ps_swc, st.st_swc, lpjmlparam, soil, crop, prec, swr_n, lwr_n)

        # output
        output_training!(output, photos, crop, soil, lpjml.μ, lpjml.σ)

    end

    return output
        
end