# using CUDA
function photosynthesis_C3!(PFT::PftParameters,
                            param::LPJmLParam,
                            photopar::PhotoPar,
                            photos::Photos,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T},
                            co2::AbstractArray{T}
) where {T <: AbstractFloat}
    
    @unpack b = PFT
    @unpack ko25, kc25, alphac3, theta = param
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p = photopar
    
    ko = ko25 * q10ko .^ ((temp .- T(25.0)) * T(0.1))
    kc = kc25 * q10kc .^ ((temp .- T(25.0)) * T(0.1))
    fac = kc .* (T(1.0) .+ po2 ./ ko)
    tau = tau25 * q10tau .^ ((temp .- T(25.0)) * T(0.1)) #reflects the abiltiy of Rubisco to discriminate between CO2 and O2
    gammastar = po2 ./ (T(2.0) * tau)

    #Choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996

    #     Intercellular CO2 partial pressure in Pa
    #     Eqn 7, Haxeltine & Prentice 1996

    pi = photos.lambda .* co2

    # calculation of C1C3, C2C3 with actual pi

    c1 = photos.tstress * alphac3 .* ((pi .- gammastar) ./ (pi .+ T(2.0) * gammastar))

    c2 = (pi .- gammastar) ./ (pi .+ fac)

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ T(1e-5))

    #   Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
    #   Eqn 5, Haxeltine & Prentice 1996

    jc = c2 .* hour2day(photos.vmax)

    #    Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #    Eqn 2, Haxeltine & Prentice 1996

    photos.agd = (je .+ jc .- sqrt.(max.(1f-7, (je .+ jc) .* (je .+ jc) .- T(4.0) * theta * je .* jc))) ./ (T(2.0) * theta) .* pet_daylength

    #    Daily dark respiration, Rd, gC/m2/day
    #    Eqn 10, Haxeltine & Prentice 1996
    # rd = b * photos.vmax

    #    Total daytime net photosynthesis, Adt, gC/m2/day
    #    Eqn 19, Haxeltine & Prentice 1996

    #    Daily dark respiration, Rd, gC/m2/day
    #    Eqn 10, Haxeltine & Prentice 1996
    # photos.rd .= ifelse.(photos.tstress .< 1e-2, zero(T), b * photos.vmax)
    gate = sigmoid.(T(50.0) * (photos.tstress .- T(1e-2)))
    photos.rd = gate * b .* photos.vmax
    photos.adt = photos.agd .- hour2day(pet_daylength) .* photos.rd

    #    Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    photos.adt = max.(photos.adt, zero(T))

    photos.adtmm = photos.adt / cmass * T(8.314) .* degCtoK(temp) / p * T(1000.0)

    # Zygote.ignore() do
    #     #    Daily dark respiration, Rd, gC/m2/day
    #     #    Eqn 10, Haxeltine & Prentice 1996
    #     photos.rd .= ifelse.(photos.tstress .< 1e-2, zero(T), b * photos.vmax)
    #     photos.adt .= photos.agd .- hour2day(pet_daylength) .* photos.rd

    #     #    Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    #     photos.adt .= max.(photos.adt, zero(T))

    #     photos.adtmm .= photos.adt / cmass * 8.314f0 .* degCtoK(temp) / p * 1000.0f0  
    #     #photos.adtmm = (photos.adt <= 0) ? 0 : photos.adt / cmass * 8.314 * degCtoK(temp)/ p * 1000.0 

    #     # idx = agd .< 0 
    #     # agd[idx] .= zero(T) #in rare occasions, agd(=GPP) can be negative, but shouldn't

    #     # idx = photos.tstress .< 1e-2
    #     # agd[idx] .= zero(T)
    #     # rd[idx] .= zero(T)
    # end

end


function photosynthesis_C4!(PFT::PftParameters,
                            param::LPJmLParam,
                            photopar::PhotoPar,
                            photos::Photos,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T}
) where {T <: AbstractFloat}
    
    @unpack b = PFT
    @unpack alphac4, theta = param
    @unpack lambdamc4, cmass, cq, p = photopar
    
    #       Parameter accounting for effect of reduced intercellular CO2
    #       concentration on photosynthesis, Phipi.
    #       Eqn 14,16, Haxeltine & Prentice 1996
    #       Fig 1b, Collatz et al 1992

    gate = sigmoid.(T(-30.0) * (photos.lambda/lambdamc4 .- one(T)))
    phipi = gate .* photos.lambda/lambdamc4 .+ (one(T) .- gate)
    # phipi = min.(one(T), photos.lambda/lambdamc4)
    # phipi = photos.lambda/lambdamc4
    c1 = photos.tstress .* phipi * alphac4
    # c2 = device(ones(T, size(c1)))

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ T(1e-5))

    # jc = c2 .* hour2day(photos.vmax)
    jc = hour2day(photos.vmax)

    #    Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #    Eqn 2, Haxeltine & Prentice 1996

    photos.agd = (je .+ jc .- sqrt.(max.(1f-7, (je .+ jc) .* (je .+ jc) .- T(4.0) * theta * je .* jc))) ./ (T(2.0) * theta) .* pet_daylength

    #    Daily dark respiration, Rd, gC/m2/day
    #    Eqn 10, Haxeltine & Prentice 1996
    # rd = b * photos.vmax

    #    Total daytime net photosynthesis, Adt, gC/m2/day
    #    Eqn 19, Haxeltine & Prentice 1996

    gate = sigmoid.(T(50.0) * (photos.tstress .- T(1e-2)))
    photos.rd = gate * b .* photos.vmax
    photos.adt = photos.agd .- hour2day(pet_daylength) .* photos.rd

    #    Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    photos.adt = max.(photos.adt, zero(T))

    photos.adtmm = photos.adt / cmass * T(8.314) .* degCtoK(temp) / p * T(1000.0)

    # Zygote.ignore() do
    #     #    Daily dark respiration, Rd, gC/m2/day
    #     #    Eqn 10, Haxeltine & Prentice 1996
    #     photos.rd .= ifelse.(photos.tstress .< 1e-2, zero(T), b * photos.vmax)
    #     photos.adt .= photos.agd .- hour2day(pet_daylength) .* photos.rd

    #     #    Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    #     photos.adt .= max.(photos.adt, zero(T))

    #     photos.adtmm .= photos.adt / cmass * 8.314f0 .* degCtoK(temp) / p * 1000.0f0  
    #     #photos.adtmm = (photos.adt <= 0) ? 0 : photos.adt / cmass * 8.314 * degCtoK(temp)/ p * 1000.0 

    #     # idx = agd .< 0 
    #     # agd[idx] .= zero(T) #in rare occasions, agd(=GPP) can be negative, but shouldn't

    #     # idx = photos.tstress .< 1e-2
    #     # agd[idx] .= zero(T)
    #     # rd[idx] .= zero(T)
    # end

end