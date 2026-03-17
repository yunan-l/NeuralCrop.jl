""" LPJmL Parameters"""
struct K_Soil10{T} # lower and upper coldest monthly mean temperature(deg C)
    fast::T
    slow::T
end

@with_kw struct LPJmLParam{T}
    ko25::T = 3.0e4
    kc25::T = 30.0
    theta::T = 0.9
    alphac3::T = 0.08
    alphac4::T = 0.053
    k::T = 0.0548
    r_growth::T = 0.25
    e0::T = 308.56
    temp_response::T = 56.02
    residue_frac::T = 0.95 # fraction of residues to be submerged by tillage
    k_soil10::K_Soil10{T} = K_Soil10{T}(0.03, 0.001)
    fastfrac::T = 0.98
    atmfrac::T = 0.6
    ALPHAM::T = 1.485
    GM::T = 2.41
    LAMBDA_OPT::T = 0.8
    PRIESTLEY_TAYLOR::T = 1.32 # Priestley-Taylor coefficient
    MINERALDENS::T = 2700 # mineral density in kg/m3
    soildepth_evap::T = 300.0
    p::T = 25
    k_temp::T = 0.02 # factor of temperature dependence of nitrogen demand for Rubisco activity
    T_0::T = -25.0 # parameter in N uptake temperature function
    T_m::T = 15.0 # parameter in N uptake temperature function
    T_r::T = 15.0 # parameter in N uptake temperature function
    k_max::T = 0.10 # maximum fraction of soil->NH4 assumed to be nitrified
    k_2::T = 0.01 # fraction of nitrified N lost as N20 flux
    soil_infil::T = 2.0 # default soil infiltration
    soil_infil_litter::T = 2.0 # soil infiltration intensification by litter cover
    percthres::T = 1.0
    manure_cn::T = 14.5 # CN ration of manure gC/gN
    nfert_split_frac::T = 0.2 # fraction of fertilizer input at sowing
    nmanure_nh4_frac::T = 0.666667 # fraction of NH4 in manure input
    nfert_no3_frac::T = 0.5 # fraction of NO3 in fertilizer input
end


""" PFT Parameters"""
struct Temp{T} # lower and upper coldest monthly mean temperature(deg C)
    low::T
    high::T
end

struct TempCO2{T}  # lower and upper temperature limit for co2 (deg C)
    low::T
    high::T
end

struct TempPhotos{T} # lower and upper limit of temperature optimum for photosynthesis(deg C)
    low::T
    high::T
end

struct TvEff{T} # min & max tv: lower and upper temperature threshold under which vernalization is possible (deg C)
    low::T
    high::T
end

struct TvOpt{T}  # min & max tv: lower and upper temperature threshold under which vernalization is optimal (deg C)
    low::T
    high::T
end

struct BaseTemp{T} # min & max basetemp: base temperature
    low::T
    high::T
end

struct nc_ratio{T} # N:C mass ratio
    root::T
    sto::T
    pool::T
end

struct ratio{T} # relative C:N ratios
    root::T
    sto::T
    pool::T
end

struct ncleaf{T} # relative C:N ratios
    low::T
    median::T
    high::T
end

struct K_Litter10{T} # lower and upper coldest monthly mean temperature(deg C)
    leaf::T
    root::T
end

struct PftParameters{T <: AbstractFloat, S <: Integer}
    name::S
    plant_type::S
    path::S
    temp::Temp{T}
    temp_co2::TempCO2{T}
    temp_photos::TempPhotos{T}
    tv_eff::TvEff{T}
    tv_opt::TvOpt{T}
    psens::T
    pb::T
    ps::T
    basetemp::BaseTemp{T}
    fphuc::T
    flaimaxc::T
    fphuk::T
    flaimaxk::T
    fphusen::T
    flaimaxharvest::T
    laimax::T
    laimin::T
    hlimit::S
    pvd_max::S
    b::T   #leaf respiration as fraction of vmax
    albedo_leaf::T
    albedo_litter::T
    alphaa::T
    lightextcoeff::T #light extinction coeffcient in Lambert-Beer equation
    longevity::T
    sla::T #fscanpft_crop.c, 2e-4 * 10^(2.25 - 0.4 * log10(longevity * 12)) / 0.4763
    respcoeff::T
    shapesenescencenorm::T
    fpc::T
    nc_ratio::nc_ratio{T}
    ratio::ratio{T}
    ncleaf::ncleaf{T}
    k_litter10::K_Litter10{T}
    beta_root::T
    intc::T # Interception storage parameter
    emax::T
    knstore::T
    vmax_up::T
    kNmin::T
    KNmin::T
    hiopt::T
    himin::T
end

""" Photosynthesis Parameters"""
@with_kw struct PhotoPar{T}
    po2::T = 20.9e3
    p::T = 1.0e5
    q10ko::T = 1.2
    q10kc::T = 2.1
    q10tau::T = 0.57
    tau25::T = 2600.0
    cmass::T = 12.0
    cq::T = 4.6e-6
    lambdamc4::T = 0.4
    lambdamc3::T = 0.8
    tmc3::T = 45.0
    tmc4::T = 55.0
end


""" Soil Parameters"""
sand = Float32.([0.22, 0.06, 0.52, 0.32, 0.10, 0.58, 0.43, 0.17, 0.58, 0.10, 0.82, 0.92, 0.24, 0.99])
silt = Float32.([0.20, 0.47, 0.06, 0.34, 0.56, 0.15, 0.39, 0.70, 0.32, 0.60, 0.12, 0.05, 0.28, 0.00])
clay = Float32.([0.58, 0.47, 0.42, 0.34, 0.34, 0.27, 0.18, 0.13, 0.10, 0.30, 0.06, 0.03, 0.48, 0.01])
w_sat = Float32.([0.468, 0.468, 0.406, 0.465, 0.464, 0.404, 0.439, 0.476, 0.434, 0.476, 0.421, 0.339, 0.468, 0.006])
tdiff_0 = Float32.([0.572, 0.502, 0.785, 0.650, 0.556, 0.780, 0.701, 0.637, 0.640, 0.637, 0.403, 0.201, 0.572, 4.137])
tdiff_15 = Float32.([0.571, 0.503, 0.791, 0.656, 0.557, 0.808, 0.740, 0.657, 0.713, 0.657, 0.529, 0.196, 0.571, 4.127])

# soil_data = hcat(sand, silt, clay)

soildepth = Float32.([200.0, 300.0, 500.0, 1000.0, 1000.0]) # five-layer soil depth (mm)
layerbound = Float32.([200.0, 500.0, 1000.0, 2000.0, 3000.0])
# beta_root = Float32.([0.969, 0.969, 0.969, 0.969]) # for crop, wheat, rice, mazie, soybean

@with_kw struct SoilPar{T}
    sand::Vector{T} = sand
    silt::Vector{T} = silt
    clay::Vector{T} = clay
    w_sat::Vector{T} = w_sat
    tdiff_0::Vector{T} = tdiff_0
    tdiff_15::Vector{T} = tdiff_15
    soildepth::Vector{T} = soildepth
end


struct K{T}
    k1::T
    k2::T
    k3::T
end

mutable struct Photos{T <: AbstractArray{ <: AbstractFloat}} 
    agd::T
    adt::T
    adtmm::T
    rd::T
    vmax::T
    lambda::T
    tstress::T
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
    carbon_sum::T
    biomass::T
    lai::T
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
end