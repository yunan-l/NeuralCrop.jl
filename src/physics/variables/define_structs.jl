mutable struct Photos{T <: AbstractArray{ <: AbstractFloat}} 
    agd::T
    adt::T
    adtmm::T
    rd::T
    vmax::T
    lambda::T
    tstress::T
end

mutable struct DailyWeather{T <: AbstractArray{ <: AbstractFloat}, F <: AbstractFloat} 
    temp::T
    prec::T
    swr::T
    lwr::T
    temp_n::T
    swr_n::T
    lwr_n::T
    daily_co2::T
    annual_co2::F 
end
    
mutable struct PetPar{T <: AbstractArray{ <: AbstractFloat}} 
    daylength::T
    par::T
    eeq::T
    albedo::T
end

mutable struct ClimBuf{T <: AbstractArray{ <: AbstractFloat}, M <: AbstractArray{<: AbstractFloat}} 
    temp::M
    mtemp::M
    mtemp20::M
    min_temp::M
    atemp::M
    atemp_mean::T
    V_req_a::T    
    V_req::T         
end

mutable struct Crop{T <: AbstractArray{ <: AbstractFloat}, B <: AbstractArray{ <: Bool}, S <: AbstractArray{ <: Integer}, M <: AbstractArray{ <: AbstractFloat}} 
    phu::T               
    vdsum::T             
    husum::T             
    fphu::T              
    senescence::B
    senescence0::B             
    harvesting::B
    harvesting0::B          
    growingdays::S
    wtype::B
    lai::T
    flaimax::T
    laimax_adjusted::T
    lai_nppdeficit::T
    biomass::T 
    leafc::T 
    rootc::T
    poolc::T
    stoc::T
    vegc::M
    carbon_sum::T
    nitrogen::T
    leafn::T 
    rootn::T
    pooln::T
    ston::T
    ndemand_tot::T
    yield::T
    npp::T
    resp::T
    phen::T
    albedo::T
    fpar::T
    apar::T
    gp::T
    transp::T
    canopy_wet::T
    intercep::T
    trans_layer::M
    rootdist::T
    isgrowing::S
    nmanure::T
    nfertilizer::T
    vscal_sum::T
    vscal::T
    ndf::T
    wdf::T
    ndemand_leaf::T
    w_demandsum::T
    w_supplysum::T
    wscal::T
end

mutable struct Calendar{S <: AbstractArray{ <: Integer}} 
    sdate::S 
    hdate::S
    scallback::S
    hcallback::S
end

mutable struct Managed_land{T <: AbstractArray{ <: AbstractFloat}} 
    manure::T # reactive nitrogen fertilizer
    fertilizer::T # manure nitrogen fertilizer
end


mutable struct Soil{T <: AbstractArray{ <: AbstractFloat}, M <: AbstractArray{ <: AbstractFloat}} 
    c_input::M
    litc::M
    decom_litc::M
    fastc::M
    slowc::M
    decom_fastc::M
    decom_slowc::M
    swc::M           # soil moisture
    sand::M 
    clay::M
    evap::M          # soil evaporation
    w::M             # oil water as fraction of whc (fractional water holding capacity)
    w_fw::M          # free water or gravitational water (mm), absolute water content between field capacity and saturation
    wpwp::M          # relative water content at wilting point
    wpwps::M         # water at permanent wilting point in mm, depends on soildepth
    wfc::M           # relative water content at field capacity
    wsat::M          # relative water content at saturation
    wsats::M         # absolute water content at saturation (mm), wsats = wsat * soildepth
    beta_soil::M
    whc::M           # water holding capacity (fraction), whc = wfc - wpwp
    whcs::M          # absolute water holding capacity (mm), whcs = whc * soildepth
    Ks::M            # saturated hydraulic conductivity (mm/h) per layer
    agtop_cover::T   # fraction of soil coverd by ag litter
    layer_depth::T
    temp::M          # soil temperature
    tdiff_0::T
    tdiff_15::T
    NO3::M
    NH4::M
    litn::M
    decom_litn::M
    fastn::M
    slown::M
    decom_fastn::M
    decom_slown::M
    n_input::M
    ph::T
    w_influx::M
    w_outflux::M
    srunoff::T
    lrunoff::M
    outflux_f::T
    infil::T
    perc::M # percolation
    tillage_frac::M
    c_shift_fast::M
    c_shift_slow::M
    respose_litc::T
    respose_fastc::T
    respose_slowc::T
    n_shift_fast::M
    n_shift_slow::M
    respose_litn::T
    respose_fastn::T
    respose_slown::T
    rh::T
    snowpack::T
    snowheight::T
    snowfraction::T
end


mutable struct DataNorm{T <: AbstractArray{ <: AbstractFloat}, M <: AbstractArray{ <: AbstractFloat}} 
    npp::T
    leafc::T 
    rootc::T
    poolc::T
    soc::T
    vegc::M
    litc::M    
    fastc::M    
    slowc::M
    swc::M           
end


mutable struct Output{T <: AbstractArray{ <: AbstractFloat}, M <: AbstractArray{<: AbstractFloat}} 
    growing_mask::T
    gpp::T
    lambda::T
    vmax::T
    resp::T
    biomass::T
    lai::T
    stoc::T
    yield::T
    reco::T
    vegc::M 
    litc::M    
    fastc::M    
    slowc::M
    swc::M
    vegn::M
    litn::M    
    fastn::M    
    slown::M
    rh::T
    eeq::T
    et::T
    prec::T
    temp::T
    harvesting_mask::T
    fphu::T
end