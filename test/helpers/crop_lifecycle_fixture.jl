function lifecycle_initial_data(::Type{T}, device = identity) where {T <: AbstractFloat}
    cells = 1
    layers = 5
    to_device(values) = device(values)
    return (
        latitude = to_device(T[45]),
        soilparams = (
            ph = to_device(T[6.5]),
            w_sat = to_device(fill(T(0.45), layers, cells)),
            sand = to_device(reshape(T[0.4], 1, cells)),
            clay = to_device(reshape(T[0.2], 1, cells)),
            tdiff_0 = to_device(T[0.7]),
            tdiff_15 = to_device(T[0.75]),
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        ModelState = (
            crop = (
                sdate = to_device(Int32[100]),
                phu = to_device(T[543]),
                manure = to_device(zeros(T, cells)),
                fertilizer = to_device(T[20]),
                residuefrac = to_device(T[0.67]),
            ),
            u0 = (
                swc = to_device(reshape(T[57.41, 55.32, 126.13, 274.59, 285.71], layers, cells)),
                litc = to_device(reshape(T[0.13, 187.5, 225.36], 3, cells)),
                fastc = to_device(fill(T(10), layers, cells)),
                slowc = to_device(fill(T(100), layers, cells)),
                litn = to_device(reshape(T[0.0047, 6.47, 9.47], 3, cells)),
                fastn = to_device(fill(T(1), layers, cells)),
                slown = to_device(fill(T(10), layers, cells)),
            ),
        ),
    )
end

function lifecycle_climate(::Type{T}, days::Integer, device = identity) where {T <: AbstractFloat}
    cells = 1
    to_device(values) = device(values)
    return (
        temp = to_device(fill(T(15), days, cells)),
        prec = to_device(fill(T(2), days, cells)),
        sw = to_device(fill(T(180), days, cells)),
        lw = to_device(fill(T(-40), days, cells)),
        wind = to_device(fill(T(2), days, cells)),
        co2 = to_device(fill(T(400), cld(days, 365))),
    )
end

function run_lifecycle_fixture(device = identity; T = Float32, days = 730)
    simulation = initialize_simulation(
        cft1,
        lifecycle_initial_data(T, device);
        device = device,
        T = T,
        days = days,
        diagnostics = false,
        fertilizer = :yes,
    )
    run_simulation!(simulation, lifecycle_climate(T, days, device); spinup = false)
    return simulation
end

event_days(values) = findall(!iszero, vec(Array(values)))
