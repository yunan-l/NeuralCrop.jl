"""
CFTParameters{T,S}

Crop functional type parameter bundle for one crop type.
Contains phenology, photosynthesis, allocation, and nutrient traits.
"""
struct Temp{T} # lower and upper coldest monthly mean temperature(deg C)
    low::T  # Lower coldest-month temperature limit (°C).
    high::T # Upper coldest-month temperature limit (°C).
end

struct TempCO2{T}  # lower and upper temperature limit for co2 (deg C)
    low::T  # Lower temperature limit of CO₂ response (°C).
    high::T # Upper temperature limit of CO₂ response (°C).
end

struct TempPhotos{T} # lower and upper limit of temperature optimum for photosynthesis(deg C)
    low::T  # Lower optimum temperature for photosynthesis (°C).
    high::T # Upper optimum temperature for photosynthesis (°C).
end

struct TvEff{T} # min & max tv: lower and upper temperature threshold under which vernalization is possible (deg C)
    low::T  # Lower temperature allowing effective vernalization (°C).
    high::T # Upper temperature allowing effective vernalization (°C).
end

struct TvOpt{T}  # min & max tv: lower and upper temperature threshold under which vernalization is optimal (deg C)
    low::T  # Lower optimum vernalization temperature (°C).
    high::T # Upper optimum vernalization temperature (°C).
end

struct BaseTemp{T} # min & max basetemp: base temperature
    low::T  # Base temperature before flowering (°C).
    high::T # Base temperature after flowering (°C).
end

struct nc_ratio{T} # N:C mass ratio
    root::T # Root nitrogen-to-carbon mass ratio (gN gC⁻¹).
    sto::T  # Storage-organ nitrogen-to-carbon mass ratio (gN gC⁻¹).
    pool::T # Mobile-pool nitrogen-to-carbon mass ratio (gN gC⁻¹).
end

struct ratio{T} # relative C:N ratios
    root::T # Root-to-leaf relative C:N ratio.
    sto::T  # Storage-organ-to-leaf relative C:N ratio.
    pool::T # Mobile-pool-to-leaf relative C:N ratio.
end

struct ncleaf{T} # relative C:N ratios
    low::T    # Minimum leaf nitrogen-to-carbon mass ratio (gN gC⁻¹).
    median::T # Reference leaf nitrogen-to-carbon mass ratio (gN gC⁻¹).
    high::T   # Maximum leaf nitrogen-to-carbon mass ratio (gN gC⁻¹).
end

struct K_Litter10{T}
    leaf::T # annual turnover rate at 10 °C (yr⁻¹)
    root::T # annual below-ground turnover rate at 10 °C (yr⁻¹)
end

struct NuptakeKinetics{T}
    vmax::T # Maximum uptake per unit fine-root carbon (gN kgC-1 day-1)
    kmin::T # Uptake term independent of Michaelis-Menten saturation
    Km::T   # Half-saturation concentration (gN m-3)
end

"""LPJmL-compatible sowing-date methods and daily trigger thresholds.

The method values follow LPJmL's ``NO_CALC_SDATE``, ``PREC_CALC_SDATE``,
``TEMP_WTYP_CALC_SDATE``, ``TEMP_PREC_CALC_SDATE`` and ``MULTICROP``
categories. ``NO_CALC_SDATE`` and ``MULTICROP`` retain the prescribed
initial date; the latter is reserved for future multi-crop/fallow events.
"""
struct SowingDateParameters{T, S <: Integer}
    method::S          # LPJmL sowing-date method category.
    temp_fall::T       # Cooling threshold for winter-type establishment (°C).
    temp_spring::T     # Warming threshold for spring establishment (°C).
    minimum_prec::T    # Minimum daily precipitation defining a wet sowing day (mm day⁻¹).
end

const SDATE_NO_CALC = Int32(0)
const SDATE_PRECIPITATION = Int32(1)
const SDATE_TEMPERATURE_WINTER = Int32(2)
const SDATE_TEMPERATURE = Int32(3)
const SDATE_TEMPERATURE_PRECIPITATION = Int32(4)
const SDATE_MULTICROP = Int32(5)

_convert_precision(::Type{T}, value::AbstractFloat) where {T <: AbstractFloat} = T(value)
_convert_precision(::Type{T}, value::Integer) where {T <: AbstractFloat} = value
_convert_precision(::Type{T}, value::Temp) where {T <: AbstractFloat} = Temp{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TempCO2) where {T <: AbstractFloat} = TempCO2{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TempPhotos) where {T <: AbstractFloat} = TempPhotos{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TvEff) where {T <: AbstractFloat} = TvEff{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::TvOpt) where {T <: AbstractFloat} = TvOpt{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::BaseTemp) where {T <: AbstractFloat} = BaseTemp{T}(T(value.low), T(value.high))
_convert_precision(::Type{T}, value::nc_ratio) where {T <: AbstractFloat} = nc_ratio{T}(T(value.root), T(value.sto), T(value.pool))
_convert_precision(::Type{T}, value::ratio) where {T <: AbstractFloat} = ratio{T}(T(value.root), T(value.sto), T(value.pool))
_convert_precision(::Type{T}, value::ncleaf) where {T <: AbstractFloat} = ncleaf{T}(T(value.low), T(value.median), T(value.high))
_convert_precision(::Type{T}, value::K_Litter10) where {T <: AbstractFloat} = K_Litter10{T}(T(value.leaf), T(value.root))
_convert_precision(::Type{T}, value::NuptakeKinetics) where {T <: AbstractFloat} =
    NuptakeKinetics{T}(T(value.vmax), T(value.kmin), T(value.Km))
_convert_precision(::Type{T}, value::SowingDateParameters) where {T <: AbstractFloat} =
    SowingDateParameters{T, typeof(value.method)}(
        value.method, T(value.temp_fall), T(value.temp_spring), T(value.minimum_prec),
    )

@kwdef struct CFTParameters{T <: AbstractFloat, S <: Integer}
    name::S                 # Numeric crop/CFT identifier.
    plant_type::S           # LPJmL plant-type category identifier.
    path::S                 # Photosynthetic pathway: 1 = C3, 2 = C4.
    temp::Temp{T}           # Bioclimatic cold-temperature limits (°C).
    temp_co2::TempCO2{T}    # Temperature limits of CO₂ response (°C).
    temp_photos::TempPhotos{T} # Optimum photosynthesis temperature interval (°C).
    tv_eff::TvEff{T}        # Effective vernalization temperature interval (°C).
    tv_opt::TvOpt{T}        # Optimum vernalization temperature interval (°C).
    sowing_date::SowingDateParameters{T, S} # Dynamic-sowing thresholds and daily trigger.
    psens::T                # Photoperiod sensitivity coefficient.
    pb::T                   # Lower/short-day photoperiod threshold (h).
    ps::T                   # Upper/long-day photoperiod threshold (h).
    basetemp::BaseTemp{T}   # Development base temperatures (°C).
    fphuc::T                # Heat-unit fraction ending initial LAI phase.
    flaimaxc::T             # LAI fraction at end of initial LAI phase.
    fphuk::T                # Heat-unit fraction at LAI sigmoid inflection.
    flaimaxk::T             # LAI fraction at sigmoid inflection.
    fphusen::T              # Heat-unit fraction at onset of senescence.
    flaimaxharvest::T       # Fraction of maximum LAI retained at harvest.
    laimax::T               # Maximum potential leaf-area index (m² m⁻²).
    laimin::T               # Minimum/reference LAI trait used by allocation (m² m⁻²).
    hlimit::S               # Maximum crop-cycle duration (days).
    pvd_max::S              # Maximum vernalizing-day requirement (days).
    b::T                    # Leaf maintenance respiration as fraction of `vcmax`.
    albedo_leaf::T          # Leaf shortwave albedo (0–1).
    albedo_litter::T        # Surface-litter shortwave albedo (0–1).
    alphaa::T               # Canopy allometry/allocation coefficient.
    lightextcoeff::T        # Lambert–Beer canopy light-extinction coefficient.
    longevity::T            # Characteristic leaf longevity (years).
    sla::T                  # Specific leaf area (m² leaf gC⁻¹).
    respcoeff::T            # Root/storage maintenance-respiration coefficient.
    shapesenescencenorm::T   # Normalized senescence curve shape parameter.
    fpc::T                   # Foliage projective cover scaling (0–1).
    nc_ratio::nc_ratio{T}   # Organ nitrogen-to-carbon ratios.
    ratio::ratio{T}         # Relative organ C:N ratios used in N partitioning.
    ncleaf::ncleaf{T}       # Minimum/reference/maximum leaf N:C ratios.
    k_litter10::K_Litter10{T} # Litter turnover rates at 10 °C.
    beta_root::T            # Exponential root-depth distribution parameter.
    intc::T                 # Canopy interception storage parameter.
    emax::T                 # Maximum transpiration/conductance scaling parameter.
    gmin::T                 # Minimum canopy conductance (mm s⁻¹).
    knstore::T              # Fraction of N demand allocated to storage reserve.
    no3_uptake::NuptakeKinetics{T} # Root NO₃ uptake kinetic parameters.
    nh4_uptake::NuptakeKinetics{T} # Root NH₄ uptake kinetic parameters.
    hiopt::T                # Optimal harvest index.
    himin::T                # Minimum harvest index under stress.
end

"""Return a CFT parameter set whose floating fields consistently use `T`."""
function convert_precision(::Type{T}, cft::CFTParameters{<:AbstractFloat, S}) where {T <: AbstractFloat, S <: Integer}
    names = fieldnames(typeof(cft))
    values = map(name -> _convert_precision(T, getfield(cft, name)), names)
    kwargs = NamedTuple{names}(values)
    return CFTParameters{T, S}(; kwargs...)
end

const CFT_NAMES = (
    "temperate cereals",
    "rice",
    "maize",
    "tropical cereals",
    "pulses",
    "temperate roots",
    "tropical roots",
    "oil crops sunflower",
    "oil crops soybean",
    "oil crops groundnut",
    "oil crops rapeseed",
    "sugarcane",
)

_crop_sla(longevity) = Float32(2e-4 * 10^(2.25 - 0.4 * log10(longevity * 12)) / 0.4763)

function _crop_cft(;
    id, path, temp_co2, temp_photos, tv_eff = (1000, 1000), tv_opt = (1000, 1000),
    psens = 1, pb = 0, ps = 24, basetemp, fphuc, flaimaxc, fphuk, flaimaxk = 0.95,
    fphusen, flaimaxharvest, laimax, laimin, hlimit, pvd_max = 0, beta_root,
    longevity, emax, gmin, shapesenescencenorm, storage_ratio, hiopt, himin,
    sowing_method = SDATE_NO_CALC, temp_fall = 1000, temp_spring = 1000,
)
    T = Float32
    return CFTParameters{T, Int32}(
        name = id,
        plant_type = 1,
        path = path,
        temp = Temp{T}(-1000, 1000),
        temp_co2 = TempCO2{T}(temp_co2...),
        temp_photos = TempPhotos{T}(temp_photos...),
        tv_eff = TvEff{T}(tv_eff...),
        tv_opt = TvOpt{T}(tv_opt...),
        sowing_date = SowingDateParameters{T, Int32}(
            Int32(sowing_method), T(temp_fall), T(temp_spring), T(0.1),
        ),
        psens = psens,
        pb = pb,
        ps = ps,
        basetemp = BaseTemp{T}(basetemp, basetemp),
        fphuc = fphuc,
        flaimaxc = flaimaxc,
        fphuk = fphuk,
        flaimaxk = flaimaxk,
        fphusen = fphusen,
        flaimaxharvest = flaimaxharvest,
        laimax = laimax,
        laimin = laimin,
        hlimit = hlimit,
        pvd_max = pvd_max,
        b = 0.031,
        albedo_leaf = 0.18,
        albedo_litter = 0.06,
        alphaa = 1,
        lightextcoeff = 0.5,
        longevity = longevity,
        sla = _crop_sla(longevity),
        respcoeff = 0.8,
        shapesenescencenorm = shapesenescencenorm,
        fpc = 1,
        nc_ratio = nc_ratio{T}(1 / 30, 1 / 100, 1 / 100),
        ratio = ratio{T}(1.16, storage_ratio, 3),
        ncleaf = ncleaf{T}(1 / 58.8, 1 / 25, 1 / 14.3),
        k_litter10 = K_Litter10{T}(0.97, 0.97),
        beta_root = beta_root,
        intc = 0.02,
        emax = emax,
        gmin = gmin,
        knstore = 0.1,
        no3_uptake = NuptakeKinetics{T}(1.5, 0.05, 0.70),
        nh4_uptake = NuptakeKinetics{T}(4.7, 0.05, 0.45),
        hiopt = hiopt,
        himin = himin,
    )
end

# The order is the LPJmL 6.1.1 CFT/management-band order. Rainfed and irrigated
# variants share the same biological parameter set.
const cft1 = _crop_cft(id=1, path=1, temp_co2=(0, 40), temp_photos=(12, 17),
    tv_eff=(-4, 17), tv_opt=(3, 10), pb=8, ps=20, basetemp=0,
    sowing_method=SDATE_TEMPERATURE_WINTER, temp_fall=12, temp_spring=5,
    fphuc=.05, flaimaxc=.05, fphuk=.45, fphusen=.70, flaimaxharvest=0,
    laimax=7, laimin=2, hlimit=360, pvd_max=70, beta_root=.94,
    longevity=.50, emax=8, gmin=1, shapesenescencenorm=2, storage_ratio=.99,
    hiopt=.50, himin=.20)
const cft2 = _crop_cft(id=2, path=1, temp_co2=(6, 55), temp_photos=(20, 45),
    pb=24, ps=0, basetemp=8, sowing_method=SDATE_PRECIPITATION, temp_spring=18,
    fphuc=.10, flaimaxc=.05, fphuk=.50,
    fphusen=.80, flaimaxharvest=0, laimax=7, laimin=5, hlimit=288,
    beta_root=.91, longevity=.33, emax=8, gmin=1, shapesenescencenorm=2,
    storage_ratio=1.30, hiopt=.50, himin=.25)
const cft3 = _crop_cft(id=3, path=2, temp_co2=(8, 42), temp_photos=(21, 26),
    basetemp=5, sowing_method=SDATE_TEMPERATURE_PRECIPITATION, temp_spring=14,
    fphuc=.10, flaimaxc=.05, fphuk=.50, fphusen=.75,
    flaimaxharvest=0, laimax=5, laimin=4, hlimit=334, beta_root=.94,
    longevity=.33, emax=10, gmin=1.2, shapesenescencenorm=2, storage_ratio=.83,
    hiopt=.50, himin=.30)
const cft4 = _crop_cft(id=4, path=2, temp_co2=(6, 55), temp_photos=(20, 45),
    basetemp=8, fphuc=.15, flaimaxc=.01, fphuk=.50, fphusen=.85,
    flaimaxharvest=0, laimax=7, laimin=5, hlimit=299, beta_root=.94,
    longevity=.50, emax=10, gmin=1.6, shapesenescencenorm=2, storage_ratio=.99,
    hiopt=.25, himin=.10)
const cft5 = _crop_cft(id=5, path=1, temp_co2=(-4, 45), temp_photos=(10, 30),
    basetemp=1, fphuc=.15, flaimaxc=.01, fphuk=.50, fphusen=.90,
    flaimaxharvest=0, laimax=4, laimin=4, hlimit=282, beta_root=.94,
    longevity=.50, emax=8, gmin=1, shapesenescencenorm=2, storage_ratio=.45,
    hiopt=.45, himin=.10)
const cft6 = _crop_cft(id=6, path=1, temp_co2=(-4, 45), temp_photos=(10, 30),
    basetemp=3, fphuc=.15, flaimaxc=.05, fphuk=.50, fphusen=.75,
    flaimaxharvest=.75, laimax=5, laimin=5, hlimit=299, beta_root=.94,
    longevity=.50, emax=7, gmin=1, shapesenescencenorm=.5, storage_ratio=1.74,
    hiopt=3.5, himin=1.25)
const cft7 = _crop_cft(id=7, path=1, temp_co2=(6, 55), temp_photos=(20, 45),
    basetemp=15, fphuc=.15, flaimaxc=.05, fphuk=.50, fphusen=.75,
    flaimaxharvest=.75, laimax=5, laimin=5, hlimit=360, beta_root=.94,
    longevity=.50, emax=10, gmin=1.6, shapesenescencenorm=.5, storage_ratio=3.27,
    hiopt=2, himin=1.10)
const cft8 = _crop_cft(id=8, path=1, temp_co2=(8, 42), temp_photos=(25, 32),
    basetemp=6, fphuc=.15, flaimaxc=.01, fphuk=.50, fphusen=.70,
    flaimaxharvest=0, laimax=5, laimin=5, hlimit=282, beta_root=.94,
    longevity=.33, emax=7, gmin=1, shapesenescencenorm=2, storage_ratio=1.04,
    hiopt=.40, himin=.20)
const cft9 = _crop_cft(id=9, path=1, temp_co2=(5, 45), temp_photos=(28, 32),
    basetemp=7, sowing_method=SDATE_PRECIPITATION, temp_spring=13, fphuc=.15, flaimaxc=.05,
    fphuk=.50, fphusen=.70,
    flaimaxharvest=0, laimax=5, laimin=5, hlimit=282, pvd_max=70,
    beta_root=.94, longevity=.66, emax=10, gmin=1.2, shapesenescencenorm=.5,
    storage_ratio=.42, hiopt=.40, himin=.10)
const cft10 = _crop_cft(id=10, path=1, temp_co2=(6, 55), temp_photos=(20, 45),
    basetemp=14, fphuc=.15, flaimaxc=.01, fphuk=.50, fphusen=.75,
    flaimaxharvest=0, laimax=5, laimin=5, hlimit=282, pvd_max=70,
    beta_root=.94, longevity=.50, emax=10, gmin=1.6, shapesenescencenorm=.5,
    storage_ratio=.68, hiopt=.40, himin=.30)
const cft11 = _crop_cft(id=11, path=1, temp_co2=(0, 40), temp_photos=(12, 17),
    tv_eff=(-4, 17), tv_opt=(3, 10), pb=8, ps=20, basetemp=0,
    sowing_method=SDATE_TEMPERATURE_WINTER, temp_fall=17, temp_spring=5,
    fphuc=.05, flaimaxc=.01, fphuk=.50, fphusen=.85, flaimaxharvest=0,
    laimax=7, laimin=7, hlimit=360, pvd_max=70, beta_root=.94,
    longevity=.41, emax=7, gmin=1, shapesenescencenorm=2, storage_ratio=.76,
    hiopt=.30, himin=.15)
const cft12 = _crop_cft(id=12, path=2, temp_co2=(8, 42), temp_photos=(18, 30),
    basetemp=12, fphuc=.01, flaimaxc=.01, fphuk=.40, fphusen=.95,
    flaimaxharvest=.50, laimax=6, laimin=2, hlimit=360, beta_root=.94,
    longevity=.66, emax=10, gmin=1.6, shapesenescencenorm=2, storage_ratio=4.57,
    hiopt=.80, himin=.80)

const CFTS = (cft1, cft2, cft3, cft4, cft5, cft6, cft7, cft8, cft9, cft10, cft11, cft12)

crop_cft(id::Integer) = 1 <= id <= length(CFTS) ? CFTS[id] : throw(ArgumentError("crop CFT id must be in 1:12"))
function crop_cft(name::AbstractString)
    id = findfirst(==(String(name)), CFT_NAMES)
    isnothing(id) && throw(ArgumentError("unknown crop CFT '$name'"))
    return crop_cft(id)
end
