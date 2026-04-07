function soil_carbon!(param::LPJmLParam,
                      crop_cal::Calendar,
                      soil::Soil
)

    @unpack atmfrac = param

    # compute soil carbon: litter carbon and soil carbon
    soil.decom_litc = (1.0f0 .- exp.(-soil.respose_litc / 100)) .* soil.litc
    soil.litc = soil.litc  - soil.decom_litc

    # using 'callback' to adjust litter carbon due to tillage, 'scallback' means the tillage of sowing day and 'hcallback' means the tillage of harvesting day
    soil.litc = soil.litc .* (1 .- reshape(crop_cal.scallback, (1, :))) .* (1 .- reshape(crop_cal.hcallback, (1, :))) + 
                soil.tillage_frac * soil.litc .* reshape(crop_cal.scallback, (1, :)) +
                (soil.tillage_frac * (soil.litc .+ soil.c_input)) .* reshape(crop_cal.hcallback, (1, :)) 
 
    temp_response = exp.(308.56f0 * (1.0f0 / (56.02f0 + 10) .- 1.0f0 ./ (soil.temp / 10 .+ 56.02f0)))
    moist = (soil.w .* soil.whcs + soil.wpwps + soil.w_fw) ./ (soil.wsats - soil.wpwps)
    response = temp_response .* (0.04021601f0 .- 5.00505434f0 * moist.^3 + 4.26937932f0 * moist.^2 + 0.71890122f0 * moist)
    
    soil.decom_fastc = (1.0f0 .- exp.(-soil.respose_fastc .* response / 50)) .* soil.fastc
    soil.fastc = soil.fastc + soil.c_shift_fast .* sum(soil.decom_litc, dims = 1) - soil.decom_fastc
    
    soil.decom_slowc = (1.0f0 .- exp.(-soil.respose_slowc .* response / 10)) .* soil.slowc
    soil.slowc = soil.slowc + soil.c_shift_slow .* sum(soil.decom_litc, dims = 1) - soil.decom_slowc

    soil.rh = vec(sum(soil.decom_litc, dims = 1) * atmfrac .+ sum(soil.decom_fastc, dims = 1) .+ sum(soil.decom_slowc, dims = 1))
    
end