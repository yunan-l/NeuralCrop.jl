function c3_adtmm_scalar_impl(lambda::T,
                              vcmax::T,
                              tstress::T,
                              b::T,
                              co2::T,
                              temp::T,
                              apar::T,
                              daylength::T,
                              lpjmlparams::LPJmLParams,
                              photoparams::PhotoParams) where {T <: AbstractFloat}
    if tstress < T(1e-2)
        return zero(T)
    end

    @unpack ko25, kc25, alphac3, theta = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p = photoparams

    ko = ko25 * q10ko^((temp - T(25)) * T(0.1))
    kc = kc25 * q10kc^((temp - T(25)) * T(0.1))
    fac = kc * (one(T) + po2 / ko)
    tau = tau25 * q10tau^((temp - T(25)) * T(0.1))
    gammastar = po2 / (T(2) * tau)
    internal_co2 = lambda * co2
    c1 = tstress * alphac3 *
         ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
    c2 = (internal_co2 - gammastar) / (internal_co2 + fac)

    je = c1 * apar * cmass * cq / daylength
    jc = c2 * vcmax / T(24)
    agd = (je + jc - sqrt(max(zero(T), (je + jc)^2 - T(4) * theta * je * jc))) / (T(2) * theta) * daylength
    rd = b * vcmax
    adt = agd - daylength / T(24) * rd

    return adt <= zero(T) ? zero(T) :
           adt / cmass * T(8.314) * (temp + T(273.15)) / p * T(1000)
end

function c3_adtmm_scalar(lambda::T,
                         vcmax::T,
                         tstress::T,
                         b::T,
                         co2::T,
                         temp::T,
                         apar::T,
                         daylength::T;
                         lpjmlparams::LPJmLParams = lpjmlparams,
                         photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    return c3_adtmm_scalar_impl(
        lambda, vcmax, tstress, b, co2, temp, apar, daylength,
        lpjmlparams, photoparams,
    )
end

function c4_adtmm_scalar_impl(lambda::T,
                              vcmax::T,
                              tstress::T,
                              b::T,
                              temp::T,
                              apar::T,
                              daylength::T,
                              lpjmlparams::LPJmLParams,
                              photoparams::PhotoParams) where {T <: AbstractFloat}
    if tstress < T(1e-2)
        return zero(T)
    end

    @unpack alphac4, theta = lpjmlparams
    @unpack lambdamc4, cmass, cq, p = photoparams

    phipi = min(one(T), lambda / T(lambdamc4))
    c1 = tstress * phipi * T(alphac4)
    je = c1 * apar * T(cmass) * T(cq) / daylength
    jc = vcmax / T(24)
    agd = (je + jc - sqrt(max(zero(T), (je + jc)^2 - T(4) * T(theta) * je * jc))) /
          (T(2) * T(theta)) * daylength
    rd = b * vcmax
    adt = agd - daylength / T(24) * rd

    return adt <= zero(T) ? zero(T) :
           adt / T(cmass) * T(8.314) * (temp + T(273.15)) / T(p) * T(1000)
end

function c4_adtmm_scalar(lambda::T,
                         vcmax::T,
                         tstress::T,
                         b::T,
                         temp::T,
                         apar::T,
                         daylength::T;
                         lpjmlparams::LPJmLParams = lpjmlparams,
                         photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    return c4_adtmm_scalar_impl(
        lambda, vcmax, tstress, b, temp, apar, daylength,
        lpjmlparams, photoparams,
    )
end


"""
    solve_lambda_c3_lpj(fac, vcmax, tstress, b, co2, temp, apar, daylength;
                        lower=0.02, upper=0.85, tolerance=0.001,
                        max_iterations=30)

Solve LPJmL's C3 water-stress equation
`fac * (1 - lambda) - adtmm(lambda) = 0` using the compatible CPU bisection
algorithm. `co2` is atmospheric CO₂ partial pressure in Pa, while `fac` is the
water-limited conductance term constructed by the water-balance routine.

Returns `(lambda, iterations, residual)`.
"""
function solve_lambda_c3_lpj(fac::T,
                             vcmax::T,
                             tstress::T,
                             b::T,
                             co2::T,
                             temp::T,
                             apar::T,
                             daylength::T;
                             lower::T = T(0.02),
                             upper::T = T(0.85),
                             tolerance::T = T(0.001),
                             max_iterations::Integer = 30,
                             lpjmlparams::LPJmLParams = lpjmlparams,
                             photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    objective(lambda) = fac * (one(T) - lambda) -
                        c3_adtmm_scalar_impl(
                            lambda, vcmax, tstress, b, co2, temp, apar, daylength,
                            lpjmlparams, photoparams,
                        )

    lambda, iterations = lpj_bisect(
        objective,
        lower,
        upper;
        x_accuracy = zero(T),
        y_accuracy = tolerance,
        max_iterations = max_iterations,
    )

    return lambda, iterations, objective(lambda)
end

"""
    solve_lambda_c4_lpj(fac, vcmax, tstress, b, temp, apar, daylength)

CPU reference for LPJmL's C4 water-stress lambda equation. Returns
`(lambda, iterations, residual)`.
"""
function solve_lambda_c4_lpj(fac::T,
                             vcmax::T,
                             tstress::T,
                             b::T,
                             temp::T,
                             apar::T,
                             daylength::T;
                             lower::T = T(0.02),
                             upper::T = T(0.85),
                             tolerance::T = T(0.001),
                             max_iterations::Integer = 30,
                             lpjmlparams::LPJmLParams = lpjmlparams,
                             photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    objective(lambda) = fac * (one(T) - lambda) -
                        c4_adtmm_scalar_impl(
                            lambda, vcmax, tstress, b, temp, apar, daylength,
                            lpjmlparams, photoparams,
                        )

    lambda, iterations = lpj_bisect(
        objective,
        lower,
        upper;
        x_accuracy = zero(T),
        y_accuracy = tolerance,
        max_iterations = max_iterations,
    )

    return lambda, iterations, objective(lambda)
end

"""
    solve_lambda_c3!(CFT, photos, crop, pet, temp, co2)

Solve the LPJmL water-stress equation independently for every grid cell. The
fixed 30-step loop and scalar, allocation-free objective are compatible with
both CPU and GPU backends. `co2` must be atmospheric partial pressure in Pa and
`crop_canopy_auxiliary(crop).canopy_conductance` must contain actual canopy conductance after water limitation.
"""
function solve_lambda_c3!(CFT::CFTParameters,
                          crop,
                          pet::PetPar,
                          temp::AbstractArray{T},
                          co2::AbstractArray{T};
                          lpjmlparams::LPJmLParams = lpjmlparams,
                          photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    kernel_params = (
        b = T(CFT.b),
        gmin = T(CFT.gmin),
        lpjmlparams = lpjmlparams,
        photoparams = photoparams,
    )

    launch_1D!(
        solve_lambda_c3_kernel!,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_canopy_auxiliary(crop).fpar,
        crop_canopy_auxiliary(crop).apar,
        pet.daylength,
        temp,
        co2,
        kernel_params,
    )
end


"""
    solve_lambda_c4!(CFT, photos, crop, pet, temp, co2)

GPU/CPU backend implementation of LPJmL's C4 water-stress lambda solve.
`co2` is atmospheric partial pressure in Pa.
"""
function solve_lambda_c4!(CFT::CFTParameters,
                          crop,
                          pet::PetPar,
                          temp::AbstractArray{T},
                          co2::AbstractArray{T};
                          lpjmlparams::LPJmLParams = lpjmlparams,
                          photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    kernel_params = (
        b = T(CFT.b),
        gmin = T(CFT.gmin),
        lpjmlparams = lpjmlparams,
        photoparams = photoparams,
    )

    launch_1D!(
        solve_lambda_c4_kernel!,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        crop_canopy_auxiliary(crop).canopy_conductance,
        crop_canopy_auxiliary(crop).fpar,
        crop_canopy_auxiliary(crop).apar,
        pet.daylength,
        temp,
        co2,
        kernel_params,
    )
end

"""Dispatch the water-limited lambda solve to the compile-time crop pathway."""
solve_lambda!(::Val{:C3}, CFT, crop, pet, temperature, co2; kwargs...) =
    solve_lambda_c3!(CFT, crop, pet, temperature, co2; kwargs...)
solve_lambda!(::Val{:C4}, CFT, crop, pet, temperature, co2; kwargs...) =
    solve_lambda_c4!(CFT, crop, pet, temperature, co2; kwargs...)

"""
    compute_canopy_water_supply(daylength, conductance, minimum_conductance,
                                fpar, co2)

Construct LPJmL's scalar water-supply term for the internal-CO₂ bisection.
`co2` is stored as Pa and therefore uses the existing `1e-5` Pa-to-bar
conversion at this one boundary.
"""
@inline function compute_canopy_water_supply(daylength::T,
                                             conductance::T,
                                             minimum_conductance::T,
                                             fpar::T,
                                             co2::T) where {T <: AbstractFloat}
    daily_conductance = daylength * T(3600) *
        (conductance - minimum_conductance * fpar)
    return daily_conductance, daily_conductance / T(1.6) * co2 * T(1e-5)
end

"""
    compute_lambda_c3_solution(fac, vcmax, stress, b, co2, temperature, apar,
                               daylength, params)

Fixed-iteration, allocation-free LPJmL C3 bisection. Keeping the complete
scalar solve outside the kernel makes its stopping rule testable while every
grid cell retains the exact same 30-step numerical algorithm on CPU and GPU.
"""
@inline function compute_lambda_c3_solution(fac::T,
                                             vcmax::T,
                                             stress::T,
                                             b::T,
                                             co2::T,
                                             temperature::T,
                                             apar::T,
                                             daylength::T,
                                             lpjmlparams::LPJmLParams,
                                             photoparams::PhotoParams) where {T <: AbstractFloat}
    lower = T(0.02)
    upper = T(0.85)
    lower_residual = fac * (one(T) - lower) - c3_adtmm_scalar_impl(
        lower, vcmax, stress, b, co2, temperature, apar, daylength,
        lpjmlparams, photoparams,
    )
    best = (lower + upper) * T(0.5)
    best_residual = typemax(T)
    for _ in 1:30
        midpoint = (lower + upper) * T(0.5)
        residual = fac * (one(T) - midpoint) - c3_adtmm_scalar_impl(
            midpoint, vcmax, stress, b, co2, temperature, apar, daylength,
            lpjmlparams, photoparams,
        )
        if abs(residual) < best_residual
            best_residual = abs(residual)
            best = midpoint
        end
        abs(residual) < T(0.001) && break
        if lower_residual * residual <= zero(T)
            upper = midpoint
        else
            lower = midpoint
            lower_residual = residual
        end
    end
    return best
end

"""
    compute_lambda_c4_solution(fac, vcmax, stress, b, temperature, apar,
                               daylength, params)

The C4 counterpart of `compute_lambda_c3_solution`, with the C4 assimilation
relation but the same LPJmL bracket and residual-selection rule.
"""
@inline function compute_lambda_c4_solution(fac::T,
                                             vcmax::T,
                                             stress::T,
                                             b::T,
                                             temperature::T,
                                             apar::T,
                                             daylength::T,
                                             lpjmlparams::LPJmLParams,
                                             photoparams::PhotoParams) where {T <: AbstractFloat}
    lower = T(0.02)
    upper = T(0.85)
    lower_residual = fac * (one(T) - lower) - c4_adtmm_scalar_impl(
        lower, vcmax, stress, b, temperature, apar, daylength,
        lpjmlparams, photoparams,
    )
    best = (lower + upper) * T(0.5)
    best_residual = typemax(T)
    for _ in 1:30
        midpoint = (lower + upper) * T(0.5)
        residual = fac * (one(T) - midpoint) - c4_adtmm_scalar_impl(
            midpoint, vcmax, stress, b, temperature, apar, daylength,
            lpjmlparams, photoparams,
        )
        if abs(residual) < best_residual
            best_residual = abs(residual)
            best = midpoint
        end
        abs(residual) < T(0.001) && break
        if lower_residual * residual <= zero(T)
            upper = midpoint
        else
            lower = midpoint
            lower_residual = residual
        end
    end
    return best
end

@kernel inbounds = true function solve_lambda_c4_kernel!(
    lambda::AbstractArray{T},
    vcmax::AbstractArray{T},
    tstress::AbstractArray{T},
    conductance::AbstractArray{T},
    fpar::AbstractArray{T},
    apar::AbstractArray{T},
    daylength::AbstractArray{T},
    temp::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack b, gmin, lpjmlparams, photoparams = kernel_params
    co2_cell = co2[length(co2) == 1 ? 1 : cell]

    gpd, fac = compute_canopy_water_supply(
        daylength[cell], conductance[cell], gmin, fpar[cell], co2_cell,
    )

    if gpd > T(1e-5) && tstress[cell] >= T(1e-2) &&
       daylength[cell] > zero(T) && co2_cell > zero(T)
        lambda[cell] = compute_lambda_c4_solution(
            fac, vcmax[cell], tstress[cell], b, temp[cell], apar[cell],
            daylength[cell], lpjmlparams, photoparams,
        )
    else
        # LPJmL bypasses photosynthesis and returns zero GPP here. Lambda zero
        # reproduces that result when the vector photosynthesis routine follows.
        lambda[cell] = zero(T)
    end
end

@kernel inbounds = true function solve_lambda_c3_kernel!(
    lambda::AbstractArray{T},
    vcmax::AbstractArray{T},
    tstress::AbstractArray{T},
    conductance::AbstractArray{T},
    fpar::AbstractArray{T},
    apar::AbstractArray{T},
    daylength::AbstractArray{T},
    temp::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack b, gmin, lpjmlparams, photoparams = kernel_params
    co2_cell = co2[length(co2) == 1 ? 1 : cell]

    # LPJmL receives ppm and converts it to bar. NeuralCrop stores partial
    # pressure in Pa, so the equivalent conversion is Pa * 1e-5.
    gpd, fac = compute_canopy_water_supply(
        daylength[cell], conductance[cell], gmin, fpar[cell], co2_cell,
    )

    if gpd > T(1e-5) && tstress[cell] >= T(1e-2) &&
       daylength[cell] > zero(T) && co2_cell > zero(T)
        lambda[cell] = compute_lambda_c3_solution(
            fac, vcmax[cell], tstress[cell], b, co2_cell, temp[cell], apar[cell],
            daylength[cell], lpjmlparams, photoparams,
        )
    else
        # LPJmL bypasses photosynthesis and returns zero GPP here.
        lambda[cell] = zero(T)
    end
end
