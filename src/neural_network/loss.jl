function loss_crop_rollout!(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_i, batch_size, climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
    
    @unpack lpjml = data_i
    
    output = daily_crop(day_start, day_end, nn_model, ps, st, parameters, data_i, batch_size, climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
    
    # generate crop growing mask
    growing_mask_gpp = output.growing_mask[2:end, :][day_start:day_end, :]
    growing_mask_vegc = repeat(output.growing_mask[2:end, :][day_start:day_end, :], 1, 4)[:, sort(repeat(1:size(output.growing_mask[2:end, :][day_start:day_end, :], 2), 4))] 

    if growing_mask_gpp[growing_mask_gpp .!= 0] == []
        lambda_loss = mean(abs2, ((output.lambda[2:end, :][day_start:day_end, :] .- lpjml.lambda_n[day_start:day_end, :]) .* growing_mask_gpp))
        vmax_loss = mean(abs2, ((output.vmax[2:end, :][day_start:day_end, :] .- lpjml.vmax_n[day_start:day_end, :]) .* growing_mask_gpp))
        gpp_loss = mean(abs2, ((output.gpp[2:end, :][day_start:day_end, :] .- lpjml.gpp_n[day_start:day_end, :]) .* growing_mask_gpp))
        vegc_loss = mean(abs2, ((output.vegc[2:end, :][day_start:day_end, :] .- lpjml.vegc_n[day_start:day_end, :]) .* growing_mask_vegc))
        # mass_balance = mean(abs2, ((output.carbon_sum[2:end, :][day_start:day_end, :] .- (output.biomass[2:end, :][day_start:day_end, :]))/1000.0f0 .* growing_mask_gpp))
    else
        lambda_loss = mean(abs2, ((output.lambda[2:end, :][day_start:day_end, :] .- lpjml.lambda_n[day_start:day_end, :]) .* growing_mask_gpp)[growing_mask_gpp .!= 0])
        vmax_loss = mean(abs2, ((output.vmax[2:end, :][day_start:day_end, :] .- lpjml.vmax_n[day_start:day_end, :]) .* growing_mask_gpp)[growing_mask_gpp .!= 0])
        gpp_loss = mean(abs2, ((output.gpp[2:end, :][day_start:day_end, :] .- lpjml.gpp_n[day_start:day_end, :]) .* growing_mask_gpp)[growing_mask_gpp .!= 0])
        vegc_loss = mean(abs2, ((output.vegc[2:end, :][day_start:day_end, :] .- lpjml.vegc_n[day_start:day_end, :]) .* growing_mask_vegc)[growing_mask_vegc .!= 0])
        # mass_balance = mean(abs2, ((output.carbon_sum[2:end, :][day_start:day_end, :] .- (output.biomass[2:end, :][day_start:day_end, :]))/1000.0f0 .* growing_mask_gpp)[growing_mask_gpp .!= 0])
    end

    litc_loss = mean(abs2, output.litc[2:end, :][day_start:day_end, :]  .- lpjml.litc_n[day_start:day_end, :]) 
    fastc_loss = mean(abs2, output.fastc[2:end, :][day_start:day_end, :] .- lpjml.fastc_n[day_start:day_end, :]) 
    # slowc_loss = mean(abs2, output.slowc[2:end, :][day_start:day_end, :]  .- lpjml.slowc_n[day_start:day_end, :])

    swc_loss = mean(abs2, output.swc[2:end, :][day_start:day_end, :] .- lpjml.swc_n[day_start:day_end, :])

    data_loss = lambda_loss + vmax_loss + gpp_loss + vegc_loss + litc_loss/5 + fastc_loss/5 + swc_loss/5
     
    return data_loss
end