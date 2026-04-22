function pedotransfer!(soil::Soil;
                       lpjmlparams::LPJmLParams = lpjmlparams
)

    @unpack MINERALDENS = lpjmlparams # mineral density in kg/m3

    om_layer = 2 * ((soil.fastc + soil.slowc) ./ ((1 .- soil.wsat) * MINERALDENS .* soil.layer_depth)) * 100 #calculation of soil organic matter in % 
    
    # idx = om_layer .> 8
    # om_layer[idx] .= T(8.0)
    # om_layer .= ifelse.(om_layer .> 8, 8.0, om_layer)
    om_layer .= min.(om_layer, 8.0f0)

    wpwpt = -0.024f0 * soil.sand + 0.487f0 * soil.clay .+ 0.006f0 * om_layer + 0.005f0 * (soil.sand .* om_layer) - 0.013f0 * (soil.clay .* om_layer) .+ 0.068f0 * (soil.sand .* soil.clay) .+ 0.031f0
    soil.wpwp .= wpwpt + (0.14 * wpwpt .- 0.02)
    soil.wpwps .= soil.wpwp .* soil.layer_depth
    ws33t = 0.278f0 * soil.sand + 0.034f0 * soil.clay .+ 0.022f0 * om_layer - 0.018f0 * (soil.sand .* om_layer) - 0.027f0 * (soil.clay .* om_layer) .- 0.584f0 * (soil.sand .* soil.clay) .+ 0.078f0
    ws33 = ws33t + (0.636f0 * ws33t .- 0.107f0)

    wfct = -0.251f0 * soil.sand + 0.195f0 * soil.clay .+ 0.011f0 * om_layer + 0.006f0 * (soil.sand .* om_layer) - 0.027f0 * (soil.clay .* om_layer) .+ 0.452f0 * (soil.sand .* soil.clay) .+ 0.299f0
    soil.wfc .= (wfct + (((1.283f0 * wfct) .^ 2) - 0.374f0 * wfct .- 0.015f0))

    soil.wsat .= soil.wfc + ws33 .- 0.097f0 * soil.sand .+ 0.043f0
    soil.wsats .= soil.wsat .* soil.layer_depth

    # here, we ignore the effects of tillage to soil water content at saturation.
    # if(l < NTILLLAYER)
    # {
    #     soil->wsat[l] = 1 - (1-w_sat)*soil->df_tillage[l];
    #     soil->wfc[l] = w_fc - 0.2 * (w_sat - soil->wsat[l]);
    # }
    # else
    # {
    #     soil->wsat[l] = w_sat;
    #     soil->wfc[l] = w_fc;
    # }  

    # idx = (soil.wsat - soil.wfc) .< 0.05
    # soil.wfc[idx] .= soil.wsat[idx] .- 0.05
    soil.wfc .= ifelse.((soil.wsat - soil.wfc) .< 0.05, soil.wsat .- 0.05f0, soil.wfc)

    soil.beta_soil .= -2.655f0 ./ log10.(soil.wfc ./ soil.wsat)
    soil.whc .= soil.wfc - soil.wpwp
    soil.whcs .= soil.whc .* soil.layer_depth

    # idx = (soil.swc - soil.wpwps) .> 1.0e-10
    # soil.w[idx] .= min.((soil.swc[idx] - soil.wpwps[idx]) ./ soil.whcs[idx], one(T))
    # soil.w_fw[idx] .= min.(soil.swc[idx] - soil.wpwps[idx] - soil.w[idx] .* soil.whcs[idx], soil.wsats[idx] - soil.wfc[idx] .* soil.layer_depth)
    # idx = (soil.swc - soil.wpwps) <= 1.0e-10
    # soil.w[idx] .= zero(T)
    # soil.w_fw[idx] .= zero(T)
    soil.w .= ifelse.((soil.swc - soil.wpwps) .> 1.0e-10, min.((soil.swc - soil.wpwps) ./ soil.whcs, 1.0f0), 0.0f0)
    soil.w_fw .= ifelse.((soil.swc - soil.wpwps) .> 1.0e-10, min.(soil.swc - soil.wpwps - soil.w .* soil.whcs, soil.wsats - soil.wfc .* soil.layer_depth), 0.0f0)

    # Calculation of Ks
    lambda = (log.(soil.wfc) - log.(soil.wpwp)) / (log(1500) - log(33))
    soil.Ks .= 1930 * (soil.wsat - soil.wfc) .^ (3 .- lambda)

    # update agtop_cover
    dm_sum = soil.litc[1, :]/0.42 # dm_sum+=stand->soil.litter.item[l].agtop.leaf.carbon/0.42; Accounting that C content in plant dry matter is 42%
    # idx = dm_sum .< 0
    # dm_sum[idx] .= zero(T)
    # dm_sum .= ifelse.(dm_sum .< 0, 0.0, dm_sum)
    dm_sum .= max.(dm_sum, 0.0f0)
    soil.agtop_cover .= 1 .- exp.(-6e-3 * dm_sum)
end