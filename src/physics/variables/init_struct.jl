function init_crop(cell_size::Int,
                   device; # GPU or CPU
                   carbon_pools = 4,
                   soil_layers = 5
)
    crop = Crop(
        device(zeros(Float32, cell_size)),      # phu
        device(zeros(Float32, cell_size)),      # vdsum
        device(zeros(Float32, cell_size)),      # husum
        device(zeros(Float32, cell_size)),      # fphu
        device(fill(false, cell_size)),         # senescence
        device(fill(false, cell_size)),         # senescence0
        device(fill(true, cell_size)),          # harvesting
        device(fill(true, cell_size)),          # harvesting0
        device(zeros(Int32, cell_size)),        # growingdays
        device(fill(false, cell_size)),         # wtype
        device(zeros(Float32, cell_size)),      # lai: 0.000415f0
        device(zeros(Float32, cell_size)),      # flaimax: 0.000083f0
        device(zeros(Float32, cell_size)),      # laimax_adjusted: 1.0f0
        device(zeros(Float32, cell_size)),      # lai_nppdeficit
        device(zeros(Float32, cell_size)),      # biomass
        device(zeros(Float32, cell_size)),      # leafc: 0.0113804f0
        device(zeros(Float32, cell_size)),      # rootc: 8.0f0
        device(zeros(Float32, cell_size)),      # poolc: 11.9886196f0
        device(zeros(Float32, cell_size)),      # stoc: 0.0f0
        device(zeros(Float32, (carbon_pools, cell_size))),  # vegc: 0.0f0
        device(zeros(Float32, cell_size)),      # carbon_sum
        device(zeros(Float32, cell_size)),      # nitrogen
        device(zeros(Float32, cell_size)),      # leafn: 0.0f0
        device(zeros(Float32, cell_size)),      # rootn: 0.0f0
        device(zeros(Float32, cell_size)),      # pooln: 0.0f0
        device(zeros(Float32, cell_size)),      # ston: 0.0f0
        device(zeros(Float32, cell_size)),      # total N demand (gN/m2)
        device(zeros(Float32, cell_size)),      # yield: 0.0f0
        device(zeros(Float32, cell_size)),      # npp
        device(zeros(Float32, cell_size)),      # resp
        device(zeros(Float32, cell_size)),      # phen
        device(zeros(Float32, cell_size)),      # albedo
        device(zeros(Float32, cell_size)),      # fpar
        device(zeros(Float32, cell_size)),      # apar
        device(zeros(Float32, cell_size)),      # gp
        device(zeros(Float32, cell_size)),      # transp
        device(zeros(Float32, cell_size)),      # canopy wet
        device(zeros(Float32, cell_size)),      # interception
        device(zeros(Float32, (soil_layers, cell_size))), # transpiration
        device(zeros(Float32, soil_layers)),    # rootdist
        device(zeros(Int32, cell_size)),        # isgrowing: 1 means crop growing, 0 means no crops
        device(zeros(Float32, cell_size)),      # nmanure
        device(zeros(Float32, cell_size)),      # nfertilizer
        device(zeros(Float32, cell_size)),      # vscal_sum
        device(zeros(Float32, cell_size)),      # vscal
        device(zeros(Float32, cell_size)),      # ndf
        device(zeros(Float32, cell_size)),      # wdf
        device(zeros(Float32, cell_size)),      # ndemand_leaf
        device(zeros(Float32, cell_size)),      # w_demandsum
        device(zeros(Float32, cell_size)),      # w_supplysum
        device(zeros(Float32, cell_size))       # wscal
    )
    
    crop_cal = Calendar(
        device(zeros(Int32, cell_size)),  # sdate
        device(zeros(Int32, cell_size)),  # hdate
        device(zeros(Int32, cell_size)),  # sowing day callback
        device(zeros(Int32, cell_size)),  # harvesting day callback
        device(zeros(Int32, cell_size))   # harvesting_year
    )

    managed_land = Managed_land(
        device(zeros(Float32, cell_size)),  # manure
        device(zeros(Float32, cell_size)),  # fertilizer
    )

    photos = Photos(
        device(zeros(Float32, cell_size)),  # agd (gpp)
        device(zeros(Float32, cell_size)),  # adt
        device(zeros(Float32, cell_size)),  # adtmm
        device(zeros(Float32, cell_size)),  # rd
        device(zeros(Float32, cell_size)),  # vm
        device(zeros(Float32, cell_size)),  # lambda
        device(zeros(Float32, cell_size))   # tstress
    )
    
    return crop, crop_cal, managed_land, photos

end

function init_pet(cell_size::Int,
                  device
)
    
    pet = PetPar(
        device(zeros(Float32, cell_size)),  # daylength
        device(zeros(Float32, cell_size)),  # par
        device(zeros(Float32, cell_size)),  # eeq
        device(zeros(Float32, cell_size))   # albedo
    )
    
    return pet

end

function init_weather(cell_size::Int,
                      device
)
    
    dailyWeather = DailyWeather(
        device(zeros(Float32, cell_size)),             # temp
        device(zeros(Float32, cell_size)),             # prec
        device(zeros(Float32, cell_size)),             # swr
        device(zeros(Float32, cell_size)),             # lwr
        device(zeros(Float32, cell_size)),             # temp_n
        device(zeros(Float32, cell_size)),             # swr_n
        device(zeros(Float32, cell_size)),             # lwr_n
        device(zeros(Float32, cell_size)),             # daily CO2 concentration
        device(zeros(Float32, 1))                      # annual CO2 concentration
    )
    
    return dailyWeather

end


function init_climbuf(cell_size::Int,
                      device;
                      NDAYS = 31,
                      NMONTH = 12,
                      NDAYS_YEAR = 365,
                      n = 5
)
    
    climBuf = ClimBuf(
        device(zeros(Float32, (NDAYS, cell_size))),       # temp
        device(zeros(Float32, (NMONTH, cell_size))),      # mtemp
        device(fill(-9999.0f0, (NMONTH, cell_size))),     # mtemp20
        device(zeros(Float32, (n, cell_size))),           # min_temp
        device(zeros(Float32, (NDAYS_YEAR, cell_size))),  # atemp
        device(zeros(Float32, cell_size)),                # atemp_mean
        device(zeros(Float32, cell_size)),                # V_req_a
        device(fill(-9999.0f0, cell_size))                # V_req
    )
    
    return climBuf

end

function init_soil(cell_size::Int,
                   soildepth::AbstractArray{T},
                   device;
                   litc_layers = 3,
                   soil_layers = 5
) where {T <: AbstractFloat}
    
    soil = Soil(
        device(zeros(Float32, (litc_layers, cell_size))),  # c_input
        device(zeros(Float32, (litc_layers, cell_size))),  # litc
        device(zeros(Float32, (litc_layers, cell_size))),  # decom_litc
        device(zeros(Float32, (soil_layers, cell_size))),  # fastc
        device(zeros(Float32, (soil_layers, cell_size))),  # slowc
        device(zeros(Float32, (soil_layers, cell_size))),  # decom_fastc
        device(zeros(Float32, (soil_layers, cell_size))),  # decom_slowc
        device(zeros(Float32, (soil_layers, cell_size))),  # swc
        device(zeros(Float32, (1, cell_size))),            # sand
        device(zeros(Float32, (1, cell_size))),            # clay
        device(zeros(Float32, (soil_layers, cell_size))),  # evaporation
        device(zeros(Float32, (soil_layers, cell_size))),  # w
        device(zeros(Float32, (soil_layers, cell_size))),  # w_fw
        device(zeros(Float32, (soil_layers, cell_size))),  # wpwp
        device(zeros(Float32, (soil_layers, cell_size))),  # wpwps
        device(zeros(Float32, (soil_layers, cell_size))),  # wfc
        device(zeros(Float32, (soil_layers, cell_size))),  # wsat
        device(zeros(Float32, (soil_layers, cell_size))),  # wsats
        device(zeros(Float32, (soil_layers, cell_size))),  # beta_soil
        device(zeros(Float32, (soil_layers, cell_size))),  # whc
        device(zeros(Float32, (soil_layers, cell_size))),  # whcs
        device(zeros(Float32, (soil_layers, cell_size))),  # Ks
        device(zeros(Float32, cell_size)),                 # agtop_cover
        device(soildepth),                                 # soil-layer depth
        device(zeros(Float32, (soil_layers, cell_size))),  # temp
        device(zeros(Float32, cell_size)),                 # tdiff_0
        device(zeros(Float32, cell_size)),                 # tdiff_15
        device(zeros(Float32, (soil_layers, cell_size))),  # NO3
        device(zeros(Float32, (soil_layers, cell_size))),  # NH4
        device(zeros(Float32, (litc_layers, cell_size))),  # litn
        device(zeros(Float32, (litc_layers, cell_size))),  # decom_litn
        device(zeros(Float32, (soil_layers, cell_size))),  # fastn
        device(zeros(Float32, (soil_layers, cell_size))),  # slown
        device(zeros(Float32, (soil_layers, cell_size))),  # decom_fastn
        device(zeros(Float32, (soil_layers, cell_size))),  # decom_slown
        device(zeros(Float32, (litc_layers, cell_size))),  # n_input
        device(zeros(Float32, cell_size)),                 # ph
        device(zeros(Float32, (soil_layers, cell_size))),  # w_influx
        device(zeros(Float32, (soil_layers, cell_size))),  # w_outflux
        device(zeros(Float32, cell_size)),                 # srunoff
        device(zeros(Float32, (soil_layers, cell_size))),  # lrunoff
        device(zeros(Float32, cell_size)),                 # outflux_f
        device(zeros(Float32, cell_size)),                 # infil
        device(zeros(Float32, (soil_layers, cell_size))),  # perc
        device(zeros(Float32, (litc_layers, litc_layers))), # tillage_frac
        device(zeros(Float32, (soil_layers, cell_size))),  # c_shift_fast
        device(zeros(Float32, (soil_layers, cell_size))),  # c_shift_slow
        device(zeros(Float32, litc_layers)),               # respose_litc
        device(zeros(Float32, soil_layers)),               # respose_fastc
        device(zeros(Float32, soil_layers)),               # respose_fastc
        device(zeros(Float32, (soil_layers, cell_size))),  # n_shift_fast
        device(zeros(Float32, (soil_layers, cell_size))),  # n_shift_slow
        device(zeros(Float32, litc_layers)),               # respose_litn
        device(zeros(Float32, soil_layers)),               # respose_fastn
        device(zeros(Float32, soil_layers)),               # respose_slown
        device(zeros(Float32, cell_size)),                 # rh
        device(zeros(Float32, cell_size)),                 # snowpack
        device(zeros(Float32, cell_size)),                 # snowheight
        device(zeros(Float32, cell_size))                  # snowfraction
    )
    
    return soil

end

function init_data_norm(cell_size::Int,
                        device;
                        carbon_pools = 4,
                        litc_layers = 3,
                        soil_layers = 5
)

    data_norm = DataNorm(
        device(zeros(Float32, cell_size)),  # npp
        device(zeros(Float32, cell_size)),  # leafc
        device(zeros(Float32, cell_size)),  # rootc
        device(zeros(Float32, cell_size)),  # poolc
        device(zeros(Float32, cell_size)),  # soc
        device(zeros(Float32, (carbon_pools, cell_size))), # vegc
        device(zeros(Float32, (litc_layers, cell_size))),  # litc
        device(zeros(Float32, (soil_layers, cell_size))),  # fastc
        device(zeros(Float32, (soil_layers, cell_size))),  # slowc
        device(zeros(Float32, (soil_layers, cell_size))),  # swc
    )

    return data_norm

end

function init_output(cell_size::Int,
                     device;
                     vegc_pools = 4,
                     litc_layers = 3,
                     soil_layers = 5
)

    output = Output(
        device(zeros(Float32, (1, cell_size))),              # gpp
        device(zeros(Float32, (1, cell_size))),              # lambda
        device(zeros(Float32, (1, cell_size))),              # vmax
        device(zeros(Float32, (1, cell_size))),              # resp
        device(zeros(Float32, (1, cell_size))),              # biomass
        device(zeros(Float32, (1, cell_size))),              # lai
        device(zeros(Float32, (1, cell_size))),              # stoc
        device(zeros(Float32, (1, cell_size))),              # yield
        device(zeros(Float32, (1, cell_size))),              # reco
        device(zeros(Float32, (1, vegc_pools*cell_size))),   # vegc
        device(zeros(Float32, (1, litc_layers*cell_size))),  # litc
        device(zeros(Float32, (1, soil_layers*cell_size))),  # fastc
        device(zeros(Float32, (1, soil_layers*cell_size))),  # slowc
        device(zeros(Float32, (1, soil_layers*cell_size))),  # swc
        device(zeros(Float32, (1, vegc_pools*cell_size))),   # vegc
        device(zeros(Float32, (1, litc_layers*cell_size))),  # litn
        device(zeros(Float32, (1, soil_layers*cell_size))),  # fastn
        device(zeros(Float32, (1, soil_layers*cell_size))),  # slown
        device(zeros(Float32, (1, cell_size))),              # rh
        device(zeros(Float32, (1, cell_size))),              # eeq (PET)
        device(zeros(Float32, (1, cell_size))),              # et (evapotranspiration)
        device(zeros(Float32, (1, cell_size))),              # prec
        device(zeros(Float32, (1, cell_size))),              # temp
        device(zeros(Float32, (1, cell_size))),              # fphu
        device(zeros(Int32, (1, cell_size))),                # growing_mask
        device(zeros(Int32, (1, cell_size))),                # harvesting mask
        device(zeros(Int32, (1, cell_size)))                 # harvesting year
    )

    return output

end

