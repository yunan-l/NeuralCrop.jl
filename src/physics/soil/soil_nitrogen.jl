function soil_nitrogen!(crop_cal::Calendar,
                        soil::Soil
)

    # compute soil carbon: litter carbon and soil carbon
    soil.decom_litn = (1.0f0 .- exp.(-soil.respose_litn / 100)) .* soil.litn
    soil.litn = soil.litn  - soil.decom_litn

    # using 'callback' to adjust litter carbon due to tillage, 'scallback' means the tillage of sowing day and 'hcallback' means the tillage of harvesting day
    soil.litn = soil.litn .* (1 .- reshape(crop_cal.scallback, (1, :))) .* (1 .- reshape(crop_cal.hcallback, (1, :))) + 
                soil.tillage_frac * soil.litn .* reshape(crop_cal.scallback, (1, :)) +
                (soil.tillage_frac * (soil.litn .+ soil.c_input)) .* reshape(crop_cal.hcallback, (1, :)) 

    temp_response = exp.(308.56f0 * (1.0f0 / (56.02f0 + 10) .- 1.0f0 ./ (soil.temp / 10 .+ 56.02f0)))
    moist = (soil.w .* soil.whcs + soil.wpwps + soil.w_fw) ./ (soil.wsats - soil.wpwps)
    response = temp_response .* (0.04021601f0 .- 5.00505434f0 * moist.^3 + 4.26937932f0 * moist.^2 + 0.71890122f0 * moist)

    soil.decom_fastn = (1.0f0 .- exp.(-soil.respose_fastn .* response / 50)) .* soil.fastn
    soil.fastn = soil.fastn + soil.n_shift_fast .* sum(soil.decom_litn, dims = 1) - soil.decom_fastn
 
    soil.decom_slown = (1.0f0 .- exp.(-soil.respose_slown .* response / 10)) .* soil.slown
    soil.slown = soil.slown + soil.n_shift_slow .* sum(soil.decom_litn, dims = 1) - soil.decom_slown
end