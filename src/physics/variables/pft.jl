""" PFT structure"""
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

@kwdef struct PftParameters{T <: AbstractFloat, S <: Integer}
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

""" temperate cereals (wheat)"""
cft1 = PftParameters{Float32, Int32}(
    name = 1,
    plant_type = 1,
    path = 1,
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(0.0, 40.0),
    temp_photos = TempPhotos{Float32}(12.0, 17.0),
    tv_eff = TvEff{Float32}(-4.0, 17.0),
    tv_opt = TvOpt{Float32}(3.0, 10.0),
    psens = 1.0,
    pb = 8.0,
    ps = 20,
    basetemp = BaseTemp{Float32}(0.0, 0.0),
    fphuc = 0.05,
    flaimaxc = 0.05,
    fphuk = 0.45,
    flaimaxk = 0.95,
    fphusen = 0.7,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 2.0,
    hlimit = 330,
    pvd_max = 70,
    b = 0.015,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.5,
    sla = 0.0364661,
    respcoeff = 1.0,
    shapesenescencenorm = 2.0,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 0.99, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    knstore = 0.1,
    vmax_up = 5.51,
    kNmin = 0.05,
    KNmin = 1.48,
    hiopt = 0.60,
    himin = 0.20
)

""" rice """
cft2 = PftParameters{Float32, Int32}(
    name = 2,
    plant_type = 1,
    path = 1,
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(6.0, 55.0),
    temp_photos = TempPhotos{Float32}(20.0, 45.0),
    tv_eff = TvEff{Float32}(1000.0, 1000.0),
    tv_opt = TvOpt{Float32}(1000.0, 1000.0),
    psens = 1.0,
    pb = 24,
    ps = 0,
    basetemp = BaseTemp{Float32}(10.0, 10.0),
    fphuc = 0.1,
    flaimaxc = 0.05,
    fphuk = 0.5,
    flaimaxk = 0.95,
    fphusen = 0.8,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 5.0,
    hlimit = 180,
    pvd_max = 0,
    b = 0.015,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.33,
    sla = 0.0430598,
    respcoeff = 1.0,
    shapesenescencenorm = 2.0,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 1.3, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    knstore = 0.1,
    vmax_up = 5.51,
    kNmin = 0.05,
    KNmin = 1.48,
    hiopt = 0.60,
    himin = 0.25
)

""" mazie """
cft3 = PftParameters{Float32, Int32}(
    name = 3,
    plant_type = 1,
    path = 2, # C4
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(8.0, 42.0),
    temp_photos = TempPhotos{Float32}(21.0, 26.0),
    tv_eff = TvEff{Float32}(1000.0, 1000.0),
    tv_opt = TvOpt{Float32}(1000.0, 1000.0),
    psens = 1.0,
    pb = 0,
    ps = 24,
    basetemp = BaseTemp{Float32}(5.0, 15.0),
    fphuc = 0.1,
    flaimaxc = 0.05,
    fphuk = 0.5,
    flaimaxk = 0.95,
    fphusen = 0.75,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 4.0,
    hlimit = 240,
    pvd_max = 0,
    b = 0.035,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.33,
    sla = 0.0430598,
    respcoeff = 1.0,
    shapesenescencenorm = 2.0,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 0.83, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    knstore = 0.1,
    vmax_up = 5.51,
    kNmin = 0.05,
    KNmin = 1.48,
    hiopt = 0.60,
    himin = 0.30
)

""" soybean """
cft4 = PftParameters{Float32, Int32}(
    name = 4,
    plant_type = 1,
    path = 1, # C3
    temp = Temp{Float32}(-1000.0, 1000.0),
    temp_co2 = TempCO2{Float32}(5.0, 45.0),
    temp_photos = TempPhotos{Float32}(28.0, 32.0),
    tv_eff = TvEff{Float32}(1000.0, 1000.0),
    tv_opt = TvOpt{Float32}(1000.0, 1000.0),
    psens = 1.0,
    pb = 0,
    ps = 24,
    basetemp = BaseTemp{Float32}(10.0, 10.0),
    fphuc = 0.15,
    flaimaxc = 0.05,
    fphuk = 0.5,
    flaimaxk = 0.95,
    fphusen = 0.7,
    flaimaxharvest = 0.0,
    laimax = 5.0,
    laimin = 5.0,
    hlimit = 240,
    pvd_max = 70,
    b = 0.015,
    albedo_leaf = 0.18,
    albedo_litter = 0.06,
    alphaa = 1.0,
    lightextcoeff = 0.5,
    longevity = 0.66,
    sla = 0.0326332,
    respcoeff = 1.0,
    shapesenescencenorm = 0.5,
    fpc = 1.0,
    nc_ratio = nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio = ratio{Float32}(1.16, 0.42, 3),
    ncleaf = ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    k_litter10 = K_Litter10{Float32}(0.97, 0.3),
    beta_root = 0.969,
    intc = 0.01,
    emax = 8.0,
    knstore = 0.1,
    vmax_up = 5.51,
    kNmin = 0.05,
    KNmin = 1.48,
    hiopt = 0.40,
    himin = 0.10
)