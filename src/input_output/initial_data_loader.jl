"""
InitialDataLoader(data, data_index, device;
                  load_mineral_nitrogen_restart=false,
                  load_c_shift_restart=false)

Build initial model-state inputs from forcing/parameter datasets. New inputs
use the top-level `initial_state`. Mineral-N restart pools are omitted by
default because `init_states!` reconstructs NO₃ and NH₄ from slow organic N.
Set `load_mineral_nitrogen_restart=true` only when explicitly restoring a
nitrogen-limited restart state.

`c_shift` is also omitted by default. Set `load_c_shift_restart=true` only
when the supplied `initial_state` contains a calibrated or restart routing
distribution.
"""
function InitialDataLoader(data::NamedTuple,
                           data_index::Vector{Int},
                           device;
                           T::Type{<:AbstractFloat} = Float32,
                           load_mineral_nitrogen_restart::Bool = false,
                           load_c_shift_restart::Bool = false
)


    @unpack latitude, crop, soilparam = data
    hasproperty(data, :initial_state) || throw(ArgumentError(
        "initial data require the top-level `initial_state` field",
    ))
    initial_state = data.initial_state

    latitude_set = T.(latitude[data_index]) |> device

    crop = _adapt_to_device(device, (
        sdate = Int32.(crop.sdate[data_index]),
        phu = T.(crop.phu[data_index]),
        manure = T.(crop.manure[data_index]),
        fertilizer = T.(crop.fertilizer[data_index]),
        residuefrac = T.(crop.residuefrac[data_index]),
    ))

    soilparam_set = _adapt_to_device(device, (
        ph = T.(soilparam.soilph[data_index]),
        w_sat = T.(soilparam.w_sat[:, data_index]),
        sand = reshape(T.(soilparam.sand[data_index]), (1, :)),
        clay = reshape(T.(soilparam.clay[data_index]), (1, :)),
        # silt = soilparam.silt[data_index],
        tdiff_0 = T.(soilparam.tdiff_0[data_index]),
        tdiff_15 = T.(soilparam.tdiff_15[data_index]),
        soildepth = T.(soilparam.soildepth),
    ))

    u0_set = (
        swc = T.(initial_state.swc[:, data_index]),
        litc = T.(initial_state.litc[:, data_index]),
        fastc = T.(initial_state.fastc[:, data_index]),
        slowc = T.(initial_state.slowc[:, data_index]),
        litn = T.(initial_state.litn[:, data_index]),
        fastn = T.(initial_state.fastn[:, data_index]),
        slown = T.(initial_state.slown[:, data_index]),
    )
    if load_mineral_nitrogen_restart
        u0_set = merge(u0_set, (
            soil_NH4 = T.(initial_state.soil_NH4[:, data_index]),
            soil_NO3 = T.(initial_state.soil_NO3[:, data_index]),
        ))
    end
    u0_set = _adapt_to_device(device, u0_set)

    model_state = (crop = crop, u0 = u0_set)
    if load_c_shift_restart
        shift_source = if hasproperty(initial_state, :c_shift_fast) &&
                          hasproperty(initial_state, :c_shift_slow)
            initial_state
        elseif hasproperty(data, :c_shift_fast) && hasproperty(data, :c_shift_slow)
            data
        else
            throw(ArgumentError(
                "load_c_shift_restart=true requires initial_state.c_shift_fast and initial_state.c_shift_slow",
            ))
        end
        model_state = merge(model_state, (
            c_shift_fast = T.(shift_source.c_shift_fast[:, data_index]),
            c_shift_slow = T.(shift_source.c_shift_slow[:, data_index]),
        ))
    end
    model_state = _adapt_to_device(device, model_state)

    InitialData = (
        latitude = latitude_set,
        soilparams = soilparam_set,
        ModelState = model_state
    )

    return InitialData
end
