
function init_structs!(param::LPJmLParam,
                       PFT::PftParameters,
                       InitialData::NamedTuple,
                       cell_size::Int,
                       device,
)

    @unpack residue_frac, fastfrac, atmfrac, k_soil10 = param
    @unpack k_litter10, beta_root = PFT

    @unpack lpjml, soilparam = InitialData

    phu = lpjml.crop.phu
    sdate = lpjml.crop.sdate
    manure = lpjml.crop.manure
    fertilizer = lpjml.crop.fertilizer
    c_shift_fast = lpjml.c_shift_fast
    c_shift_slow = lpjml.c_shift_slow
    u0 = lpjml.u0
    
    climbuf = init_climbuf(cell_size, device)
    crop, crop_cal, managed_land, photos = init_crop(cell_size, device)
    crop.phu = copy(phu)
    Zygote.ignore() do
        rootdist = root_distribution(beta_root)
        # idx = crop.phu .< 0
        # crop.wtype[idx] .= true
        # crop.phu[idx] .= -crop.phu[idx]
        crop.wtype .= ifelse.(crop.phu .< 0, true, crop.wtype)
        crop.phu .= ifelse.(crop.phu .< 0, -crop.phu, crop.phu)
        crop.rootdist .= device(rootdist)
    end
    crop_cal.sdate = sdate
    managed_land.manure = manure
    managed_land.fertilizer = fertilizer
    pet = init_pet(cell_size, device)
    soil = init_soil(cell_size, soilparam.soildepth, device)
    soil.litc = copy(u0.litc)
    soil.fastc = copy(u0.fastc)
    soil.slowc = copy(u0.slowc)
    soil.litn = copy(u0.litn)
    soil.fastn = copy(u0.fastn)
    soil.slown = copy(u0.slown)
    soil.swc = copy(u0.swc)
    soil.NO3 = copy(u0.soil_NO3)
    soil.NH4 = copy(u0.soil_NH4)
    soil.ph = soilparam.ph
    soil.wsat = soilparam.w_sat
    soil.sand = soilparam.sand
    soil.clay = soilparam.clay
    soil.tdiff_0 = soilparam.tdiff_0
    soil.tdiff_15 = soilparam.tdiff_15

    soil.tillage_frac = device([(1 - residue_frac) 0.0f0 0.0f0; residue_frac 1.0f0 0.0f0; 0.0f0 0.0f0 1.0f0])
    soil.c_shift_fast = device(c_shift_fast * fastfrac * (1.0f0 - atmfrac))
    soil.c_shift_slow = device(c_shift_slow * (1.0f0 - fastfrac) * (1.0f0 - atmfrac))
    soil.respose_litc = device([k_litter10.leaf, k_litter10.leaf, k_litter10.root])
    soil.respose_fastc = device(fill(k_soil10.fast, 5))
    soil.respose_slowc = device(fill(k_soil10.slow, 5))

    soil.n_shift_fast = device(c_shift_fast * fastfrac * (1.0f0 - atmfrac))
    soil.n_shift_slow = device(c_shift_slow * (1.0f0 - fastfrac) * (1.0f0 - atmfrac))
    soil.respose_litn = device([k_litter10.leaf, k_litter10.leaf, k_litter10.root])
    soil.respose_fastn = device(fill(k_soil10.fast, 5))
    soil.respose_slown = device(fill(k_soil10.slow, 5))

    # var_norm = init_data_norm(cell_size, device)
    # var_norm.swc = copy(u0.swc_n)
    output = init_output(cell_size, device)

    return climbuf, crop, crop_cal, photos, pet, soil, managed_land, output
end