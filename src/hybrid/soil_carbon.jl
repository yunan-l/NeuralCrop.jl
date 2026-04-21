function soil_carbon!(nn_model,
                      ps,
                      st,
                      temp_n::AbstractArray{T},
                      sw_n::AbstractArray{T},
                      crop_cal::Calendar,
                      soil::Soil;
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    @unpack atmfrac = lpjmlparams

    # compute soil carbon: litter carbon and soil carbon
    input = vcat(reshape((soil.swc./soil.layer_depth)[1, :], (1, :)), reshape(temp_n, (1, :)))
    soil.litc, soil.decom_litc = hybrid_litc(nn_model.lit, soil.litc, ps.ps_litd, st.st_litd, input, soil.respose_litc)
    # using 'callback' to adjust litter carbon due to tillage, 'scallback' means the tillage of sowing day and 'hcallback' means the tillage of harvesting day
    update_litc_tillage!(soil, crop_cal)

    input = vcat(mean(soil.swc./soil.layer_depth, dims = 1), reshape(temp_n, (1, :)), reshape(sw_n, (1, :)))
    soil.fastc, soil.decom_fastc = hybrid_soilc(nn_model.soil, soil.fastc, ps.ps_soild, st.st_soild, input, soil.respose_fastc, soil.c_shift_fast, sum(soil.decom_litc, dims = 1))
    soil.slowc, soil.decom_slowc = hybrid_soilc(nn_model.soil, soil.slowc, ps.ps_soild, st.st_soild, input, soil.respose_slowc, soil.c_shift_slow, sum(soil.decom_litc, dims = 1))
    soil.rh = vec(sum(soil.decom_litc, dims = 1) * atmfrac .+ sum(soil.decom_fastc, dims = 1) .+ sum(soil.decom_slowc, dims = 1))
end