# space validation
function train_loop_rollout!(daily_crop, rollout, nn_model, ps, st, parameters, data, train_i, valid_i, year, loss_func, opt_state, η_schedule, device, save_path; batch_size=10, N_epochs=1, scheduler_offset::Int=0, save_mode::Symbol=:valid)
    
    @assert save_mode in [:valid, :train] "save_mode has to be :valid or :train"
    
    best_ps = deepcopy(ps)
    results = (i_epoch = Int[], train_loss=Float32[], learning_rate=Float32[], duration=Float32[], valid_loss=Float32[], loss_min=[Inf32], i_epoch_min=[1])
    
    progress = Progress(N_epochs, 1)
    
     # initial error 
    lowest_train_err = Inf
    lowest_valid_err = Inf
    # max_loss_threshold = 1f15
    
    for i_epoch in 1:N_epochs

        Optimisers.adjust!(opt_state, η_schedule(i_epoch + scheduler_offset)) 

        epoch_start_time = time()

        ### Training Loop ###
        loss_trian = []

        for i in 1:batch_size:length(train_i)
            batch_i = i:min(i+batch_size-1, length(train_i))
            data_index = train_i[batch_i]
            data_batch = DataLoader(data, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch

            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_trian_rollout = []
            for day in 1:rollout:365*year
                day_start = day
                day_end = min(day+rollout-1, 365*year)
                loss_p(ps) = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                l, gs = Zygote.withgradient(loss_p, ps)
                if !isnan(l) && !isinf(l)
                    push!(loss_trian_rollout, l)
                    opt_state, ps = Optimisers.update(opt_state, ps, gs[1])
                else
                    @warn "Loss too large ($l), skipping update to prevent NaN."
                end
            end
            push!(loss_trian, mean(loss_trian_rollout))
        end
        train_err = mean(loss_trian)
        epoch_time = time() - epoch_start_time
        push!(results[:i_epoch], i_epoch)
        push!(results[:train_loss], train_err)
        push!(results[:learning_rate], η_schedule(i_epoch))
        push!(results[:duration], epoch_time)
        
        ### Validation Loop ###
        loss_valid = []
        for i in 1:batch_size:length(valid_i)
            batch_i = i:min(i+batch_size-1, length(valid_i))
            data_index = valid_i[batch_i]
            data_batch = DataLoader(data, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch

            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_valid_rollout = []
            for day in 1:rollout:365*year
                day_start = day
                day_end = min(day+rollout-1, 365*year)
                l = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                push!(loss_valid_rollout, l)
            end
            push!(loss_valid, mean(loss_valid_rollout))
        end
        valid_err = mean(loss_valid)
        push!(results[:valid_loss], valid_err)
 
        next!(progress; showvalues = [(:i_epoch, i_epoch), (:train_loss, train_err), (:valid_loss, valid_err)])
        
        if i_epoch == N_epochs
            plot_loss_curve(results[:i_epoch], results[:train_loss], results[:valid_loss], save_path.fig_path)
        end
        
        if save_mode==:valid
            if valid_err < lowest_valid_err
                lowest_valid_err = valid_err 
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_valid_err
                results[:i_epoch_min] .= i_epoch
            end
        else
            if train_err < lowest_train_err
                lowest_train_err = train_err
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_train_err
                results[:i_epoch_min] .= i_epoch
            end
        end

        ps_save = adapt(Array, best_ps)
        @save save_path.ps_path ps_save

    end

    return nn_model, best_ps, st, results
    
end

# space validation
function train_loop_winter_wheat_rollout!(daily_crop, rollout, nn_model, ps, st, parameters, data, train_i, valid_i, year, loss_func, opt_state, η_schedule, device, save_path; batch_size=10, N_epochs=1, scheduler_offset::Int=0, save_mode::Symbol=:valid)
    
    @assert save_mode in [:valid, :train] "save_mode has to be :valid or :train"
    
    best_ps = deepcopy(ps)
    results = (i_epoch = Int[], train_loss=Float32[], learning_rate=Float32[], duration=Float32[], valid_loss=Float32[], loss_min=[Inf32], i_epoch_min=[1])
    
    progress = Progress(N_epochs, 1)
    
     # initial error 
    lowest_train_err = Inf
    lowest_valid_err = Inf
    # max_loss_threshold = 1f15
    

    for i_epoch in 1:N_epochs

        Optimisers.adjust!(opt_state, η_schedule(i_epoch + scheduler_offset)) 

        epoch_start_time = time()

        ### Training Loop ###
        loss_trian = []

        for i in 1:batch_size:length(train_i)
            batch_i = i:min(i+batch_size-1, length(train_i))
            data_index = train_i[batch_i]
            data_batch = DataLoader_winter_wheat(data, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch

            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_trian_rollout = []
            for day in 1:rollout:365*year
                day_start = day
                day_end = min(day+rollout-1, 365*year)
                loss_p(ps) = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                l, gs = Zygote.withgradient(loss_p, ps)
                if !isnan(l) && !isinf(l)
                    push!(loss_trian_rollout, l)
                    opt_state, ps = Optimisers.update(opt_state, ps, gs[1])
                else
                    @warn "Loss too large ($l), skipping update to prevent NaN."
                end
            end
            push!(loss_trian, mean(loss_trian_rollout))
        end
        train_err = mean(loss_trian)
        epoch_time = time() - epoch_start_time
        push!(results[:i_epoch], i_epoch)
        push!(results[:train_loss], train_err)
        push!(results[:learning_rate], η_schedule(i_epoch))
        push!(results[:duration], epoch_time)
        
        ### Validation Loop ###
        loss_valid = []
        for i in 1:batch_size:length(valid_i)
            batch_i = i:min(i+batch_size-1, length(valid_i))
            data_index = valid_i[batch_i]
            data_batch = DataLoader_winter_wheat(data, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch

            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_valid_rollout = []
            for day in 1:rollout:365*year
                day_start = day
                day_end = min(day+rollout-1, 365*year)
                l = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                push!(loss_valid_rollout, l)
            end
            push!(loss_valid, mean(loss_valid_rollout))
        end
        valid_err = mean(loss_valid)
        push!(results[:valid_loss], valid_err)
 
        next!(progress; showvalues = [(:i_epoch, i_epoch), (:train_loss, train_err), (:valid_loss, valid_err)])
        
        if i_epoch == N_epochs
            plot_loss_curve(results[:i_epoch], results[:train_loss], results[:valid_loss], save_path.fig_path)
        end
        
        if save_mode==:valid
            if valid_err < lowest_valid_err
                lowest_valid_err = valid_err 
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_valid_err
                results[:i_epoch_min] .= i_epoch
            end
        else
            if train_err < lowest_train_err
                lowest_train_err = train_err
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_train_err
                results[:i_epoch_min] .= i_epoch
            end
        end

        ps_save = adapt(Array, best_ps)
        @save save_path.ps_path ps_save

    end

    return nn_model, best_ps, st, results
    
end


# time validation
function train_loop_rollout!(daily_crop, rollout, nn_model, ps, st, parameters, data, data_valid, train_i, valid_i, year, year_valid, loss_func, opt_state, η_schedule, device, save_path; batch_size=10, batch_size_valid=10, N_epochs=1, scheduler_offset::Int=0, save_mode::Symbol=:valid)
    
    @assert save_mode in [:valid, :train] "save_mode has to be :valid or :train"
    
    best_ps = deepcopy(ps)
    results = (i_epoch = Int[], train_loss=Float32[], learning_rate=Float32[], duration=Float32[], valid_loss=Float32[], loss_min=[Inf32], i_epoch_min=[1])
    
    progress = Progress(N_epochs, 1)
    
     # initial error 
    lowest_train_err = Inf
    lowest_valid_err = Inf
    # max_loss_threshold = 1f15
    
    for i_epoch in 1:N_epochs

        Optimisers.adjust!(opt_state, η_schedule(i_epoch + scheduler_offset)) 

        epoch_start_time = time()

        ### Training Loop ###
        loss_trian = []

        for i in 1:batch_size:length(train_i)
            batch_i = i:min(i+batch_size-1, length(train_i))
            data_index = train_i[batch_i]
            data_batch = DataLoader(data, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch

            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_trian_rollout = []
            for day in 1:rollout:365*year
                day_start = day
                day_end = min(day+rollout-1, 365*year)
                loss_p(ps) = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                l, gs = Zygote.withgradient(loss_p, ps)
                if !isnan(l) && !isinf(l)
                    push!(loss_trian_rollout, l)
                    opt_state, ps = Optimisers.update(opt_state, ps, gs[1])
                else
                    @warn "Loss too large ($l), skipping update to prevent NaN."
                end
            end
            push!(loss_trian, mean(loss_trian_rollout))
        end
        train_err = mean(loss_trian)
        epoch_time = time() - epoch_start_time
        push!(results[:i_epoch], i_epoch)
        push!(results[:train_loss], train_err)
        push!(results[:learning_rate], η_schedule(i_epoch))
        push!(results[:duration], epoch_time)
        
        ### Validation Loop ###
        loss_valid = []
        for i in 1:batch_size_valid:length(valid_i)
            batch_i = i:min(i+batch_size_valid-1, length(valid_i))
            data_index = valid_i[batch_i]
            data_batch = DataLoader(data_valid, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch

            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_valid_rollout = []
            for day in 1:rollout:365*year_valid
                day_start = day
                day_end = min(day+rollout-1, 365*year_valid)
                l = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                push!(loss_valid_rollout, l)
            end
            push!(loss_valid, mean(loss_valid_rollout))
        end
        valid_err = mean(loss_valid)
        push!(results[:valid_loss], valid_err)
 
        next!(progress; showvalues = [(:i_epoch, i_epoch), (:train_loss, train_err), (:valid_loss, valid_err)])
        
        if i_epoch == N_epochs
            plot_loss_curve(results[:i_epoch], results[:train_loss], results[:valid_loss], save_path.fig_path)
        end
        
        if save_mode==:valid
            if valid_err < lowest_valid_err
                lowest_valid_err = valid_err 
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_valid_err
                results[:i_epoch_min] .= i_epoch
            end
        else
            if train_err < lowest_train_err
                lowest_train_err = train_err
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_train_err
                results[:i_epoch_min] .= i_epoch
            end
        end

        ps_save = adapt(Array, best_ps)
        @save save_path.ps_path ps_save

    end

    return nn_model, best_ps, st, results
    
end


# time validation
function train_loop_winter_wheat_rollout!(daily_crop, rollout, nn_model, ps, st, parameters, data, data_valid, train_i, valid_i, year, year_valid, loss_func, opt_state, η_schedule, device, save_path; batch_size=10, batch_size_valid=10, N_epochs=1, scheduler_offset::Int=0, save_mode::Symbol=:valid)
    
    @assert save_mode in [:valid, :train] "save_mode has to be :valid or :train"
    
    best_ps = deepcopy(ps)
    results = (i_epoch = Int[], train_loss=Float32[], learning_rate=Float32[], duration=Float32[], valid_loss=Float32[], loss_min=[Inf32], i_epoch_min=[1])
    
    progress = Progress(N_epochs, 1)
    
     # initial error 
    lowest_train_err = Inf
    lowest_valid_err = Inf
    # max_loss_threshold = 1f15
    

    for i_epoch in 1:N_epochs

        Optimisers.adjust!(opt_state, η_schedule(i_epoch + scheduler_offset)) 

        epoch_start_time = time()

        ### Training Loop ###
        loss_trian = []

        for i in 1:batch_size:length(train_i)
            batch_i = i:min(i+batch_size-1, length(train_i))
            data_index = train_i[batch_i]
            data_batch = DataLoader_winter_wheat(data, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate = data_batch
            
            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_trian_rollout = []
            for day in 1:rollout:365*year
                day_start = day
                day_end = min(day+rollout-1, 365*year)
                loss_p(ps) = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                l, gs = Zygote.withgradient(loss_p, ps)
                if !isnan(l) && !isinf(l)
                    push!(loss_trian_rollout, l)
                    opt_state, ps = Optimisers.update(opt_state, ps, gs[1])
                else
                    @warn "Loss too large ($l), skipping update to prevent NaN."
                end
            end
            push!(loss_trian, mean(loss_trian_rollout))
        end
        train_err = mean(loss_trian)
        epoch_time = time() - epoch_start_time
        push!(results[:i_epoch], i_epoch)
        push!(results[:train_loss], train_err)
        push!(results[:learning_rate], η_schedule(i_epoch))
        push!(results[:duration], epoch_time)
        
        ### Validation Loop ###
        loss_valid = []
        for i in 1:batch_size_valid:length(valid_i)
            batch_i = i:min(i+batch_size_valid-1, length(valid_i))
            data_index = valid_i[batch_i]
            data_batch = DataLoader_winter_wheat(data_valid, data_index, device)
            @unpack cft, lpjmlparam  = parameters
            @unpack climate, lpjml, soilparam = data_batch
            
            InitialData = InitilDataLoader(data_batch, data_index, device)
            climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, InitialData, length(data_index), device)
            # climbuf, crop, crop_cal, photos, pet, soil, managed_land, output = init_structs!(lpjmlparam, cft, lpjml.crop.phu, lpjml.crop.sdate, lpjml.crop.manure, lpjml.crop.fertilizer, lpjml.c_shift_fast, lpjml.c_shift_slow, lpjml.u0, soilparam, length(data_index), device)
            
            spin_up_climbuf!(cft, climate.temp_spinup, climbuf, 1, device)
            loss_valid_rollout = []
            for day in 1:rollout:365*year_valid
                day_start = day
                day_end = min(day+rollout-1, 365*year_valid)
                l = loss_func(daily_crop, day_start, day_end, nn_model, ps, st, parameters, data_batch, length(data_index),  climbuf, crop, crop_cal, photos, pet, soil, managed_land, output, device)
                push!(loss_valid_rollout, l)
            end
            push!(loss_valid, mean(loss_valid_rollout))
        end
        valid_err = mean(loss_valid)
        push!(results[:valid_loss], valid_err)
 
        next!(progress; showvalues = [(:i_epoch, i_epoch), (:train_loss, train_err), (:valid_loss, valid_err)])
        
        if i_epoch == N_epochs
            plot_loss_curve(results[:i_epoch], results[:train_loss], results[:valid_loss], save_path.fig_path)
        end
        
        if save_mode==:valid
            if valid_err < lowest_valid_err
                lowest_valid_err = valid_err 
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_valid_err
                results[:i_epoch_min] .= i_epoch
            end
        else
            if train_err < lowest_train_err
                lowest_train_err = train_err
                best_ps = deepcopy(ps)
                results[:loss_min] .= lowest_train_err
                results[:i_epoch_min] .= i_epoch
            end
        end

        ps_save = adapt(Array, best_ps)
        @save save_path.ps_path ps_save

    end

    return nn_model, best_ps, st, results
    
end