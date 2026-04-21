""" LPJmL Parameters"""
struct K_Soil10{T} # lower and upper coldest monthly mean temperature(deg C)
    fast::T
    slow::T
end

@kwdef struct LPJmLParams{T}
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
    maxsnowpack::T = 20000.0 # maximum snowpack (mm)
end

""" Photosynthesis Parameters"""
@kwdef struct PhotoParams{T}
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

@kwdef struct SoilParams{T}
    sand::Vector{T} = sand
    silt::Vector{T} = silt
    clay::Vector{T} = clay
    w_sat::Vector{T} = w_sat
    tdiff_0::Vector{T} = tdiff_0
    tdiff_15::Vector{T} = tdiff_15
    soildepth::Vector{T} = soildepth
end

const lpjmlparams = LPJmLParams{Float32}()
const photoparams = PhotoParams{Float32}()
const soilparams = SoilParams{Float32}()


@kwdef struct SnowParams{T}
    tsnow::T = 0.0
    snow_skin_depth::T = 40.0 # snow skin layer depth (mm water equivalent)
    th_diff_snow ::T = 0.2/6.3f5 # thermal diffusivity of snow [m2/s]
    lambda_snow::T = 0.2
    c_water2ice::T = 0.3f9 # the energy that is needed/released during water/ice conversion (J/m3)
    c_watertosnow::T = 6.70 # Conversion factor from water to snowdepth, i.e. 1 cm water equals 6.7 cm of snow
    c_roughness::T = 0.06 # height of vegetation below the canopy
end

const snowparams = SnowParams{Float32}()