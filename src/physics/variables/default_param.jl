lpjmlparam = LPJmLParam{Float32}()
photopar = PhotoPar{Float32}()
soilpar = SoilPar{Float32}()

function root_distribution(beta_root::AbstractFloat)

    layerbound = Float32.([200.0, 500.0, 1000.0, 2000.0, 3000.0])

    BOTTOMLAYER = length(layerbound)
    totalroots = 1 - beta_root^(layerbound[BOTTOMLAYER] / 10)
    rootdist = zeros(BOTTOMLAYER)
    rootdist[1] = (1 - beta_root^(layerbound[1] / 10)) / totalroots
    for l in 2:BOTTOMLAYER
        rootdist[l] = (beta_root^(layerbound[l-1] / 10) - beta_root^(layerbound[l] / 10)) / totalroots
    end
    
    return rootdist
end

function init_K(T, PFT::PftParameters)
    k1 = 2 * log(1 / 0.99 - 1) / (PFT.temp_co2.low - PFT.temp_photos.low)
    k2 = PFT.temp_co2.low + PFT.temp_photos.low * 0.5
    k3 = log(0.99 / 0.01) / (PFT.temp_co2.high - PFT.temp_photos.high)
    
    return K{T}(k1, k2, k3)
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
k1 = init_K(Float32, cft1)

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
k2 = init_K(Float32, cft2)

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
k3 = init_K(Float32, cft3)

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
k4 = init_K(Float32, cft4)