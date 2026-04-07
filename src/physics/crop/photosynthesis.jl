# using CUDA
function photosynthesis_C3!(PFT::PftParameters,
                            param::LPJmLParam,
                            photopar::PhotoPar,
                            photos::Photos,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T},
                            co2::AbstractArray{T};
                            comp_vmax = false # compute vmax internally
) where {T <: AbstractFloat}
    
    @unpack b = PFT
    @unpack ko25, kc25, alphac3, theta = param
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p, lambdamc3 = photopar
    
    ko = ko25 * q10ko .^ ((temp .- T(25.0)) * T(0.1))
    kc = kc25 * q10kc .^ ((temp .- T(25.0)) * T(0.1))
    fac = kc .* (T(1.0) .+ po2 ./ ko)
    tau = tau25 * q10tau .^ ((temp .- T(25.0)) * T(0.1)) #reflects the abiltiy of Rubisco to discriminate between CO2 and O2
    gammastar = po2 ./ (T(2.0) * tau)

    if comp_vmax
        p_i= lambdamc3 * co2
        c1 = photos.tstress * alphac3 .* ((p_i .- gammastar) ./ (p_i .+ T(2.0) * gammastar))
        # Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996
        c2 = (p_i .- gammastar) ./ (p_i .+ fac)
        s = (24 ./ pet_daylength) * b
        sigma = 1.0f0 .- (c2 .- s) ./ (c2 .- theta * s)
        sigma = sqrt.(max.(1f-7, sigma))
        Zygote.ignore() do
            photos.lambda .= 0.8f0
        end  
        photos.vmax = (1.0f0 / b) * (c1 ./ c2) .* ((2.0f0 * theta - 1.0f0) .* s .- (2.0f0 * theta .* s .- c2) .* sigma) .* apar * cmass * cq
    end

    # calculation of C1C3, C2C3 with actual p_i (leaf internal partial pressure of CO2)
    p_i = photos.lambda .* co2

    c1 = photos.tstress * alphac3 .* ((p_i .- gammastar) ./ (p_i .+ T(2.0) * gammastar))

    c2 = (p_i .- gammastar) ./ (p_i .+ fac)

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ 1f-5)

    #   Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
    #   Eqn 5, Haxeltine & Prentice 1996

    jc = c2 .* hour2day(photos.vmax)

    #   Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #   Eqn 2, Haxeltine & Prentice 1996

    photos.agd = (je .+ jc .- sqrt.(max.(1f-7, (je .+ jc) .* (je .+ jc) .- T(4.0) * theta * je .* jc))) ./ (T(2.0) * theta) .* pet_daylength

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996

    #   Total daytime net photosynthesis, Adt, gC/m2/day
    #   Eqn 19, Haxeltine & Prentice 1996

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996
    # photos.rd .= ifelse.(photos.tstress .< 1e-2, zero(T), b * photos.vmax)
    gate = sigmoid.(T(50.0) * (photos.tstress .- T(1e-2)))
    photos.rd = gate * b .* photos.vmax
    photos.adt = photos.agd .- hour2day(pet_daylength) .* photos.rd

    #   Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    photos.adt = max.(photos.adt, zero(T))

    photos.adtmm = photos.adt / cmass * T(8.314) .* degCtoK(temp) / p * T(1000.0)

end


function photosynthesis_C4!(PFT::PftParameters,
                            param::LPJmLParam,
                            photopar::PhotoPar,
                            photos::Photos,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T};
                            comp_vmax = false # compute vmax internally
) where {T <: AbstractFloat}
    
    @unpack b = PFT
    @unpack alphac4, theta = param
    @unpack lambdamc4, cmass, cq, p = photopar
    
    #   Parameter accounting for effect of reduced intercellular CO2
    #   concentration on photosynthesis, Phipi.
    #   Eqn 14,16, Haxeltine & Prentice 1996
    #   Fig 1b, Collatz et al 1992
    if comp_vmax
        c1 = photos.tstress * alphac4
        c2 = 1.0f0
        s = (24 ./ pet_daylength) * b
        sigma = 1.0f0 .- (c2 .- s) ./ (c2 .- theta * s)
        # sigma = sqrt.(0.5f0 * (sigma .+ sqrt(sigma .* sigma .+ (1f-3)^2)))
        sigma = sqrt.(max.(1f-7, sigma))
        Zygote.ignore() do
            photos.lambda .= 0.4f0
        end  
        photos.vmax = (1.0f0 / b) * (c1 ./ c2) .* ((2.0f0 * theta - 1.0f0) .* s .- (2.0f0 * theta .* s .- c2) .* sigma) .* apar * cmass * cq
    end

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

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ 1f-5)
    
    # jc = c2 .* hour2day(photos.vmax)
    jc = hour2day(photos.vmax)

    #   Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #   Eqn 2, Haxeltine & Prentice 1996

    photos.agd = (je .+ jc .- sqrt.(max.(1f-7, (je .+ jc) .* (je .+ jc) .- T(4.0) * theta * je .* jc))) ./ (T(2.0) * theta) .* pet_daylength

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996

    #   Total daytime net photosynthesis, Adt, gC/m2/day
    #   Eqn 19, Haxeltine & Prentice 1996

    gate = sigmoid.(T(50.0) * (photos.tstress .- T(1e-2)))
    photos.rd = gate * b .* photos.vmax
    photos.adt = photos.agd .- hour2day(pet_daylength) .* photos.rd

    #   Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    photos.adt = max.(photos.adt, zero(T))

    photos.adtmm = photos.adt / cmass * T(8.314) .* degCtoK(temp) / p * T(1000.0)

end