"""
K_Soil10{T}

Temperature-scaled decomposition coefficients for fast and slow soil pools.
"""
struct K_Soil10{T} # lower and upper coldest monthly mean temperature(deg C)
    fast::T # Daily decomposition rate of fast SOC at 10 °C (day⁻¹).
    slow::T # Daily decomposition rate of slow SOC at 10 °C (day⁻¹).
end

"""
LPJmLParams{T}

Global model parameter set controlling photosynthesis, respiration,
water, nitrogen, and management process coefficients.
"""
@kwdef struct LPJmLParams{T}
    ko25::T = 3.0e4 # Rubisco Michaelis constant for O₂ at 25 °C (Pa).
    kc25::T = 30.0 # Rubisco Michaelis constant for CO₂ at 25 °C (Pa).
    theta::T = 0.9 # Curvature of the co-limitation transition in photosynthesis.
    alphac3::T = 0.08 # C3 intrinsic quantum efficiency (mol C mol⁻¹ photons).
    alphac4::T = 0.053 # C4 intrinsic quantum efficiency (mol C mol⁻¹ photons).
    k::T = 0.0548 # Extinction/scaling coefficient used by crop respiration and radiation.
    r_growth::T = 0.25 # Fraction of assimilated carbon consumed by growth respiration.
    e0::T = 308.56 # Lloyd–Taylor activation parameter for decomposition (K).
    temp_response::T = 56.02 # Lloyd–Taylor reference-temperature offset (K).
    residue_frac::T = 0.95 # fraction of residues to be submerged by tillage
    mixing_efficiency::T = 0.9 # Fraction of the full tillage bulk-density change applied per event.
    # LPJmL parameter `bioturbate = 0.5` is an annual transfer fraction.
    # fscanparam.c converts it to the daily fraction used by the process loop.
    bioturbate::T = 1 - (1 - T(0.5))^(T(1) / T(365)) # Daily surface-to-soil litter mixing fraction.
    # LPJmL reads annual turnover rates (0.04 and 0.001 yr⁻¹) and converts
    # them to daily rates in fscanparam.c before the process routines use them.
    k_soil10::K_Soil10{T} = K_Soil10{T}(T(0.04 / 365), T(0.001 / 365)) # Fast/slow SOC rates at 10 °C.
    fastfrac::T = 0.98 # Retained decomposed litter routed to fast soil pool (fraction).
    atmfrac::T = 0.5 # Decomposed litter carbon emitted directly to atmosphere (fraction).
    ALPHAM::T = 1.485 # Priestley–Taylor aerodynamic correction used by LPJmL.
    GM::T = 2.41 # Empirical canopy-conductance coefficient.
    LAMBDA_OPT::T = 0.8 # Initial/default internal-to-ambient CO₂ ratio.
    PRIESTLEY_TAYLOR::T = 1.32 # Priestley-Taylor coefficient
    MINERALDENS::T = 2700 # mineral density in kg/m3
    soildepth_evap::T = 300.0 # Effective soil depth available to evaporation (mm).
    p::T = 25 # Atmospheric-pressure/psychrometric scaling constant used by PET.
    k_temp::T = 0.0693 # factor of temperature dependence of nitrogen demand for Rubisco activity
    T_0::T = -25.0 # parameter in N uptake temperature function
    T_m::T = 15.0 # parameter in N uptake temperature function
    T_r::T = 15.0 # parameter in N uptake temperature function
    k_max::T = 0.10 # maximum fraction of soil->NH4 assumed to be nitrified
    k_2::T = 0.01 # fraction of nitrified N lost as N20 flux
    soil_cn_ratio::T = 15.0 # soil organic matter C:N ratio
    immobilization_k::T = 5.0e-3 # half-saturation coefficient for immobilization
    nitrification_a::T = 0.45 # LPJmL nitrification moisture-response parameter a.
    nitrification_b::T = 1.27 # LPJmL nitrification moisture-response parameter b.
    nitrification_c::T = 0.0012 # LPJmL nitrification moisture-response parameter c.
    nitrification_d::T = 2.84 # LPJmL nitrification moisture-response parameter d.
    CDN::T = 1.2 # shape factor for denitrification (LPJmL soil.h)
    n2o_denit_frac::T = 0.11 # fraction of denitrified N emitted as N2O
    volatil_wind::T = 1.5 # default wind speed (m/s) if no wind forcing is provided
    volatil_length::T = 1.0 # characteristic length scale (m)
    soil_infil::T = 2.0 # default soil infiltration
    soil_infil_litter::T = 2.0 # soil infiltration intensification by litter cover
    percthres::T = 1.0 # Percolation threshold/scaling coefficient.
    manure_cn::T = 14.5 # CN ration of manure gC/gN
    nfert_split_frac::T = 0.2 # fraction of fertilizer input at sowing
    nmanure_nh4_frac::T = 0.666667 # fraction of NH4 in manure input
    nfert_no3_frac::T = 0.5 # fraction of NO3 in fertilizer input
    maxsnowpack::T = 20000.0 # maximum snowpack (mm)
end
"""
lpjmlparams

Default `LPJmLParams{Float32}` singleton used across process routines.
"""
const lpjmlparams = LPJmLParams{Float32}()


"""
PhotoParams{T}

Photosynthesis-specific constants used by C3/C4 process functions.
"""
@kwdef struct PhotoParams{T}
    po2::T = 20.9e3     # Atmospheric O₂ partial pressure (Pa).
    p::T = 1.0e5        # Reference atmospheric pressure (Pa).
    q10ko::T = 1.2      # Q10 of Rubisco O₂ Michaelis constant.
    q10kc::T = 2.1      # Q10 of Rubisco CO₂ Michaelis constant.
    q10tau::T = 0.57    # Q10 of Rubisco CO₂/O₂ specificity.
    tau25::T = 2600.0   # Rubisco specificity factor at 25 °C.
    cmass::T = 12.0     # Molar mass of carbon (g mol⁻¹).
    cq::T = 4.6e-6      # PAR photon-to-mole conversion (mol photons J⁻¹).
    lambdamc4::T = 0.4  # Minimum/nominal C4 internal-to-ambient CO₂ ratio.
    lambdamc3::T = 0.8  # Initial/nominal C3 internal-to-ambient CO₂ ratio.
    tmc3::T = 45.0      # Upper thermal inhibition temperature for C3 photosynthesis (°C).
    tmc4::T = 55.0      # Upper thermal inhibition temperature for C4 photosynthesis (°C).
end
"""
photoparams

Default `PhotoParams{Float32}` singleton used in photosynthesis kernels.
"""
const photoparams = PhotoParams{Float32}()


"""
Soil texture and hydraulic default lookup vectors used to initialize `SoilParams`.
"""
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

"""
SoilParams{T}

Default soil parameter table for texture, saturation water content, and
thermal diffusivity references.
"""
@kwdef struct SoilParams{T}
    sand::Vector{T} = sand       # Sand fractions in the soil-type lookup table (0–1).
    silt::Vector{T} = silt       # Silt fractions in the soil-type lookup table (0–1).
    clay::Vector{T} = clay       # Clay fractions in the soil-type lookup table (0–1).
    w_sat::Vector{T} = w_sat     # Saturated water-content lookup values (fraction).
    tdiff_0::Vector{T} = tdiff_0 # Thermal diffusivity at zero moisture (mm² s⁻¹).
    tdiff_15::Vector{T} = tdiff_15 # Thermal diffusivity at moisture 0.15 (mm² s⁻¹).
    soildepth::Vector{T} = soildepth # Model soil-layer thicknesses (mm).
end

"""
soilparams

Default `SoilParams{Float32}` singleton for soil initialization.
"""
const soilparams = SoilParams{Float32}()


"""
SnowParams{T}

Physical constants for snow accumulation, insulation, and melt processes.
"""
@kwdef struct SnowParams{T}
    tsnow::T = 0.0 # Air-temperature threshold separating rain and snow (°C).
    snow_skin_depth::T = 40.0 # snow skin layer depth (mm water equivalent)
    th_diff_snow ::T = 0.2/6.3f5 # thermal diffusivity of snow [m2/s]
    lambda_snow::T = 0.2 # Snow thermal conductivity (W m⁻¹ K⁻¹).
    c_water2ice::T = 3.0e8 # the energy that is needed/released during water/ice conversion (J/m3)
    c_watertosnow::T = 6.70 # Conversion factor from water to snowdepth, i.e. 1 cm water equals 6.7 cm of snow
    c_roughness::T = 0.06 # height of vegetation below the canopy
end
"""
snowparams

Default `SnowParams{Float32}` singleton used by `snow!`.
"""
const snowparams = SnowParams{Float32}()

"""Numerical and bulk-thermal constants for the layered soil heat solver."""
@kwdef struct SoilThermalParams{T}
    seconds_per_day::T = 86400.0 # Seconds in one model day (s day⁻¹).
    diffusivity_conversion::T = 0.0864 # mm² s⁻¹ to m² day⁻¹
    soil_heat_capacity::T = 1.2e6 # J m⁻³ K⁻¹, LPJmL dry-soil baseline
    litter_carbon_fraction::T = 0.42 # Carbon fraction of dry surface litter (gC gDM⁻¹).
    litter_bulk_density::T = 71.1 # kg dry matter m⁻³
    litter_porosity::T = 0.952 # Pore-volume fraction of the surface-litter layer.
    litter_conductivity_dry::T = 0.05 # W m⁻¹ K⁻¹
    litter_conductivity_saturated_unfrozen::T = 0.554636 # Saturated unfrozen litter conductivity (W m⁻¹ K⁻¹).
    litter_conductivity_saturated_frozen::T = 2.106374 # Saturated frozen litter conductivity (W m⁻¹ K⁻¹).
    mineral_heat_capacity::T = 1.9259e6 # J m⁻³ K⁻¹
    water_heat_capacity::T = 4.2e6 # J m⁻³ K⁻¹
    ice_heat_capacity::T = 2.1e6 # J m⁻³ K⁻¹
    volumetric_fusion_heat::T = 3.0e8 # J m⁻³ water
    solid_conductivity::T = 8.0 # W m⁻¹ K⁻¹
    water_conductivity::T = 0.57 # W m⁻¹ K⁻¹
    ice_conductivity::T = 2.2 # W m⁻¹ K⁻¹
    phase_change_substeps::Int32 = 24 # Numerical substeps used by daily freeze–thaw integration.
end

"""Default `SoilThermalParams{Float32}` used by `soil_temperature!`."""
const soil_thermal_params = SoilThermalParams{Float32}()


"""
SoilDecompParams

Local LPJmL-style soil carbon and nitrogen decomposition response parameters used in this file.
"""
@kwdef struct SoilDecompParams{T <: AbstractFloat}
    e0::T = 308.56 # Lloyd–Taylor temperature-response activation parameter (K).
    intercept::T = 0.04021601 # Intercept of the LPJmL soil-moisture response polynomial.
    moist3::T = -5.00505434 # Cubic soil-moisture response coefficient.
    moist2::T = 4.26937932 # Quadratic soil-moisture response coefficient.
    moist1::T = 0.71890122 # Linear soil-moisture response coefficient.
    eps::T = 1e-7 # Numerical lower bound used in decomposition calculations.
end
"""
soil_decomp_params

Default `SoilDecompParams{Float32}` singleton used by soil decomposition routines.
"""
const soil_decomp_params = SoilDecompParams{Float32}()

"""
ModelParameters{T}

Precision-consistent collection of the global process parameter sets.  This
keeps model constants in the same floating-point type as the runtime state and
provides one object that can be passed through a simulation.
"""
struct ModelParameters{T <: AbstractFloat}
    lpjml::LPJmLParams{T}                  # Global LPJmL-derived process coefficients.
    photosynthesis::PhotoParams{T}         # C3/C4 biochemical photosynthesis constants.
    snow::SnowParams{T}                    # Snow accumulation, insulation, and melt constants.
    soil_thermal::SoilThermalParams{T}     # Layered soil heat and phase-change constants.
    soil_decomposition::SoilDecompParams{T} # Soil C–N decomposition response constants.
end

function ModelParameters(::Type{T}) where {T <: AbstractFloat}
    return ModelParameters{T}(
        LPJmLParams{T}(),
        PhotoParams{T}(),
        SnowParams{T}(),
        SoilThermalParams{T}(),
        SoilDecompParams{T}(),
    )
end

ModelParameters() = ModelParameters(Float32)

const FERTILIZER_MODES = (:no, :yes, :auto)

"""Normalize the LPJmL-style mineral-fertilizer mode."""
function fertilizer_mode(mode)
    normalized = mode isa Symbol ? mode : Symbol(lowercase(String(mode)))
    normalized in FERTILIZER_MODES || throw(ArgumentError(
        "fertilizer must be one of :no, :yes, or :auto",
    ))
    return normalized
end

_convert_parameter_value(::Type{T}, value::AbstractFloat) where {T <: AbstractFloat} = T(value)
_convert_parameter_value(::Type{T}, value::Integer) where {T <: AbstractFloat} = value
_convert_parameter_value(::Type{T}, value::Bool) where {T <: AbstractFloat} = value
_convert_parameter_value(::Type{T}, value::K_Soil10) where {T <: AbstractFloat} =
    K_Soil10{T}(T(value.fast), T(value.slow))

function _convert_parameter_struct(::Type{T}, ::Type{P}, params) where {T <: AbstractFloat, P}
    names = fieldnames(typeof(params))
    values = map(name -> _convert_parameter_value(T, getfield(params, name)), names)
    return P(; NamedTuple{names}(values)...)
end

convert_precision(::Type{T}, params::LPJmLParams) where {T <: AbstractFloat} =
    _convert_parameter_struct(T, LPJmLParams{T}, params)
convert_precision(::Type{T}, params::PhotoParams) where {T <: AbstractFloat} =
    _convert_parameter_struct(T, PhotoParams{T}, params)
convert_precision(::Type{T}, params::SnowParams) where {T <: AbstractFloat} =
    _convert_parameter_struct(T, SnowParams{T}, params)
convert_precision(::Type{T}, params::SoilThermalParams) where {T <: AbstractFloat} =
    _convert_parameter_struct(T, SoilThermalParams{T}, params)
convert_precision(::Type{T}, params::SoilDecompParams) where {T <: AbstractFloat} =
    _convert_parameter_struct(T, SoilDecompParams{T}, params)

"""Return a copy of a global parameter bundle converted to precision `T`."""
function convert_precision(::Type{T}, params::ModelParameters) where {T <: AbstractFloat}
    return ModelParameters{T}(
        convert_precision(T, params.lpjml),
        convert_precision(T, params.photosynthesis),
        convert_precision(T, params.snow),
        convert_precision(T, params.soil_thermal),
        convert_precision(T, params.soil_decomposition),
    )
end

"""Construct the default soil lookup table in floating-point precision `T`."""
function SoilParams(::Type{T}) where {T <: AbstractFloat}
    return SoilParams{T}(
        T.(sand), T.(silt), T.(clay), T.(w_sat), T.(tdiff_0), T.(tdiff_15),
        T.(soildepth),
    )
end

convert_precision(::Type{T}, params::SoilParams) where {T <: AbstractFloat} =
    SoilParams{T}(
        T.(params.sand), T.(params.silt), T.(params.clay), T.(params.w_sat),
        T.(params.tdiff_0), T.(params.tdiff_15), T.(params.soildepth),
    )
