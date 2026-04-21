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

""" temperate cereals (wheat)"""
cft1 = PftParameters{Float32, Int32}(
    1,
    1,
    1,
    Temp{Float32}(-1000.0, 1000.0),
    TempCO2{Float32}(0.0, 40.0),
    TempPhotos{Float32}(12.0, 17.0),
    TvEff{Float32}(-4.0, 17.0),
    TvOpt{Float32}(3.0, 10.0),
    1.0,
    8.0,
    20,
    BaseTemp{Float32}(0.0, 0.0),
    0.05,
    0.05,
    0.45,
    0.95,
    0.7,
    0.0,
    5.0,
    2.0,
    330,
    70,
    0.015,
    0.18,
    0.06,
    1.0,
    0.5,
    0.5,
    0.0364661,
    1.0,
    2.0,
    1.0,
    nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio{Float32}(1.16, 0.99, 3),
    ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    K_Litter10{Float32}(0.97, 0.3),
    0.969,
    0.01,
    8.0,
    0.1,
    5.51,
    0.05,
    1.48,
    0.60,
    0.20
)

""" rice """
cft2 = PftParameters{Float32, Int32}(
    2,
    1,
    1,
    Temp{Float32}(-1000.0, 1000.0),
    TempCO2{Float32}(6.0, 55.0),
    TempPhotos{Float32}(20.0, 45.0),
    TvEff{Float32}(1000.0, 1000.0),
    TvOpt{Float32}(1000.0, 1000.0),
    1.0,
    24,
    0,
    BaseTemp{Float32}(10.0, 10.0),
    0.1,
    0.05,
    0.5,
    0.95,
    0.8,
    0.0,
    5.0,
    5.0,
    180,
    0.0,
    0.015,
    0.18,
    0.06,
    1.0,
    0.5,
    0.33,
    0.0430598,
    1.0,
    2.0,
    1.0,
    nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio{Float32}(1.16, 1.3, 3),
    ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    K_Litter10{Float32}(0.97, 0.3),
    0.969,
    0.01,
    8.0,
    0.1,
    5.51,
    0.05,
    1.48,
    0.60,
    0.25
)

""" mazie """
cft3 = PftParameters{Float32, Int32}(
    3,
    1,
    2, #C4
    Temp{Float32}(-1000.0, 1000.0),
    TempCO2{Float32}(8.0, 42.0),
    TempPhotos{Float32}(21.0, 26.0),
    TvEff{Float32}(1000.0, 1000.0),
    TvOpt{Float32}(1000.0, 1000.0),
    1.0,
    0,
    24,
    BaseTemp{Float32}(5.0, 15.0),
    0.1,
    0.05,
    0.5,
    0.95,
    0.75,
    0.0,
    5.0,
    4.0,
    240,
    0.0,
    0.035,
    0.18,
    0.06,
    1.0,
    0.5,
    0.33,
    0.0430598,
    1.0,
    2.0,
    1.0,
    nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio{Float32}(1.16, 0.83, 3),
    ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    K_Litter10{Float32}(0.97, 0.3),
    0.969,
    0.01,
    8.0,
    0.1,
    5.51,
    0.05,
    1.48,
    0.60,
    0.30
)

""" soybean """
cft4 = PftParameters{Float32, Int32}(
    4,
    1,
    1, #C3
    Temp{Float32}(-1000.0, 1000.0),
    TempCO2{Float32}(5.0, 45.0),
    TempPhotos{Float32}(28.0, 32.0),
    TvEff{Float32}(1000.0, 1000.0),
    TvOpt{Float32}(1000.0, 1000.0),
    1.0,
    0,
    24,
    BaseTemp{Float32}(10.0, 10.0),
    0.15,
    0.05,
    0.5,
    0.95,
    0.7,
    0.0,
    5.0,
    5.0,
    240,
    70,
    0.015,
    0.18,
    0.06,
    1.0,
    0.5,
    0.66,
    0.0326332,
    1.0,
    0.5,
    1.0,
    nc_ratio{Float32}(1/25.0, 1/25.0, 1/25.0),
    ratio{Float32}(1.16, 0.42, 3),
    ncleaf{Float32}(1/58.8, 1/25.0, 1/14.3),
    K_Litter10{Float32}(0.97, 0.3),
    0.969,
    0.01,
    8.0,
    0.1,
    5.51,
    0.05,
    1.48,
    0.40,
    0.10
)