function soil_nitrogen!(nn_model,
                        ps,
                        st,
                        temp_n::AbstractArray{T},
                        sw_n::AbstractArray{T},
                        crop_cal::Calendar,
                        soil::Soil
)  where {T <: AbstractFloat}


    # compute soil nitrogen: litter nitrogen and soil nitrogen
    input = vcat(reshape((soil.swc./soil.layer_depth)[1, :], (1, :)), reshape(temp_n, (1, :)))
    soil.litn, soil.decom_litn = hybrid_litn(nn_model.lit, soil.litn, ps.litd, st.litd, input, soil.respose_litn)
    # using 'callback' to adjust litter nitrogen due to tillage, 'scallback' means the tillage of sowing day and 'hcallback' means the tillage of harvesting day
    update_litn_tillage!(soil, crop_cal)

    input = vcat(mean(soil.swc./soil.layer_depth, dims = 1), reshape(temp_n, (1, :)), reshape(sw_n, (1, :)))
    soil.fastn, soil.decom_fastn = hybrid_soiln(nn_model.soil, soil.fastn, ps.soild, st.soild, input, soil.respose_fastn, soil.n_shift_fast, sum(soil.decom_litn, dims = 1))
    soil.slown, soil.decom_slown = hybrid_soiln(nn_model.soil, soil.slown, ps.soild, st.soild, input, soil.respose_slown, soil.n_shift_slow, sum(soil.decom_litn, dims = 1))
end