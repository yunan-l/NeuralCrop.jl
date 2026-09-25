using NCDatasets

@testset "HWSD C/N preprocessing" begin
    organic_carbon = fill(1.0f0, 7, 2)
    total_nitrogen = fill(1.0f0, 7, 2)
    bulk_density = fill(1.0f0, 7, 2)
    coarse_fragments = fill(0.0f0, 7, 2)
    stocks = hwsd_layer_stocks(
        organic_carbon, total_nitrogen, bulk_density, coarse_fragments,
    )
    @test stocks.carbon[1, 1] == 2000.0f0
    @test stocks.carbon[6, 1] == 5000.0f0
    @test stocks.nitrogen[1, 1] == 200.0f0
    @test stocks.nitrogen[6, 1] == 500.0f0

    component_carbon = Float32[
        10 20 30
        10 20 30
        NaN 20 NaN
        NaN 20 NaN
        NaN 20 NaN
        NaN 20 NaN
        NaN 20 NaN
    ]
    component_nitrogen = component_carbon ./ 10
    mixed = mix_hwsd_components(
        component_carbon,
        component_nitrogen,
        Float32[60, 20, 20],
        Int[3, 1, 3],
    )
    @test mixed.carbon[1:2] == Float32[16, 16]
    @test mixed.carbon[3:7] == fill(4.0f0, 5)
    @test mixed.nitrogen == mixed.carbon ./ 10
    @test mixed.uncertain == Bool[0, 0, 1, 1, 1, 1, 1]

    true_missing = copy(component_carbon)
    true_missing[2, 2] = NaN
    unresolved = mix_hwsd_components(
        true_missing,
        component_nitrogen,
        Float32[60, 20, 20],
        Int[3, 1, 3],
    )
    @test isnan(unresolved.carbon[2])

    carbon_layers = remap_hwsd_layers(stocks.carbon)
    @test carbon_layers.values[:, 1] == Float32[2000, 3000, 5000, 10000, 10000]
    @test sum(carbon_layers.values[1:4, 1]) == sum(stocks.carbon[:, 1])
    @test carbon_layers.uncertain[:, 1] == Bool[0, 0, 0, 0, 1]
    missing_deep = remap_hwsd_layers(stocks.carbon; deep_rule = :missing)
    @test isnan(missing_deep.values[5, 1])
    @test missing_deep.uncertain[5, 1]

    organic_carbon[:, 2] .= 2
    selection = CellSelection(1:1, Int32[42])
    targets = preprocess_hwsd_cn(
        organic_carbon,
        total_nitrogen,
        bulk_density,
        coarse_fragments,
        Int32[1, 1],
        Float64[1, 3],
        selection,
    )
    @test targets.soil_organic_carbon[:, 1] ==
        Float32[3500, 5250, 8750, 17500, 17500]
    @test targets.total_nitrogen[:, 1] ==
        Float32[200, 300, 500, 1000, 1000]
    @test targets.coverage == ones(Float32, 5, 1)
    @test targets.uncertain[:, 1] == Bool[0, 0, 0, 0, 1]
    @test targets.provenance.source_version == "HWSD v2.01"
    @test targets.provenance.deep_rule == :extend_deepest_density
    @test targets.provenance.conservation.carbon_g ==
        vec(sum(targets.soil_organic_carbon .* 4; dims = 2))

    grid = GridIndex(
        Float64[-0.25, 0.25], Float64[0.25, 0.75],
        Int32[0 2; 1 3], Int32[0, 1, 2, 3],
        Int32[1, 2, 1, 2], Int32[1, 1, 2, 2],
    )
    mapping = hwsd_tile_mapping(
        Float64[-0.375, -0.125, 0.125, 0.375],
        Float64[0.125, 0.375, 0.625, 0.875],
        grid,
    )
    @test reshape(mapping.target_indices, 4, 4) == Int32[
        1 1 3 3
        1 1 3 3
        2 2 4 4
        2 2 4 4
    ]
    @test all(mapping.pixel_area .> 0)
    @test isapprox(
        sum(reshape(mapping.pixel_area, 4, 4)[1:2, 1:2]),
        NeuralCropData._EARTH_RADIUS_M^2 * deg2rad(0.5) *
            (sin(deg2rad(0.5)) - sin(deg2rad(0.0)));
        rtol = 1e-14,
    )

    missing_carbon = copy(organic_carbon)
    missing_carbon[:, 1] .= NaN
    partial = preprocess_hwsd_cn(
        missing_carbon,
        total_nitrogen,
        bulk_density,
        coarse_fragments,
        Int32[1, 1],
        Float64[1, 3],
        selection;
        minimum_coverage = 0.7,
    )
    @test partial.coverage[:, 1] == fill(0.75f0, 5)
    @test partial.soil_organic_carbon[:, 1] ==
        Float32[4000, 6000, 10000, 20000, 20000]
    @test all(partial.uncertain)

    allocation = SoilPoolAllocation(
        selection,
        fill(0.25f0, 5, 1),
        fill(0.50f0, 5, 1),
        reshape(Float32[0.40, 0.25, 0.15, 0.10, 0.10], 5, 1),
        reshape(Float32[0.50, 0.20, 0.15, 0.10, 0.05], 5, 1);
        cft_id = 5,
        irrigated = true,
        provenance = (source = "test",),
    )
    mktempdir() do directory
        filled_targets = SoilCNTargets(
            targets.selection,
            targets.layer_bounds,
            targets.soil_organic_carbon,
            targets.total_nitrogen,
            targets.coverage,
            targets.uncertain,
            merge(targets.provenance, (
                fill_policy = "nearest complete 0.5-degree HWSD cell",
                donor_longitude = 10.25,
                donor_latitude = 20.25,
                fill_distance_km = 42.0,
                original_minimum_coverage = 0.0,
            )),
        )
        path = write_soil_cn_targets(joinpath(directory, "hwsd_cn.nc"), filled_targets)
        restored = read_soil_cn_targets(path)
        @test restored.selection.cell_ids == targets.selection.cell_ids
        @test restored.layer_bounds == targets.layer_bounds
        @test restored.soil_organic_carbon == targets.soil_organic_carbon
        @test restored.total_nitrogen == targets.total_nitrogen
        @test restored.coverage == targets.coverage
        @test restored.uncertain == targets.uncertain
        @test restored.provenance.source_version == targets.provenance.source_version
        @test restored.provenance.fill_policy == "nearest complete 0.5-degree HWSD cell"
        @test restored.provenance.donor_longitude == 10.25
        @test restored.provenance.fill_distance_km == 42.0

        allocation_path = write_soil_pool_allocation(
            joinpath(directory, "pool_allocation.nc"), allocation,
        )
        restored_allocation = read_soil_pool_allocation(allocation_path)
        NCDataset(allocation_path, "r") do dataset
            @test size(dataset["fast_carbon_fraction"]) == (5, 1, 1)
            @test dataset["cft_id"][:] == Int32[5]
            @test dataset["irrigated"][:] == Int8[1]
        end
        @test restored_allocation.selection.cell_ids == allocation.selection.cell_ids
        @test restored_allocation.cft_id == 5
        @test restored_allocation.irrigated
        @test restored_allocation.fast_carbon_fraction == allocation.fast_carbon_fraction
        @test restored_allocation.fast_nitrogen_fraction == allocation.fast_nitrogen_fraction
        @test restored_allocation.c_shift_fast == allocation.c_shift_fast
        @test restored_allocation.c_shift_slow == allocation.c_shift_slow

        legacy_path = joinpath(directory, "legacy_pool_allocation.nc")
        NCDataset(legacy_path, "c") do dataset
            defDim(dataset, "layer", 5)
            defDim(dataset, "cell", 1)
            defDim(dataset, "patch", 1)
            defVar(dataset, "cell_id", Int32, ("cell",))[:] = allocation.selection.cell_ids
            defVar(dataset, "pft_id", Int32, ("patch",))[:] = Int32[5]
            defVar(dataset, "irrigated", Int8, ("patch",))[:] = Int8[1]
            for name in (:fast_carbon_fraction, :fast_nitrogen_fraction, :c_shift_fast, :c_shift_slow)
                values = getproperty(allocation, name)
                defVar(dataset, String(name), Float32, ("layer", "cell", "patch"))[:, :, 1] = values
            end
            dataset.attrib["schema_version"] = string(DATA_SCHEMA_VERSION)
        end
        @test read_soil_pool_allocation(legacy_path).cft_id == 5
    end

    @test_throws ArgumentError hwsd_layer_stocks(
        fill(101.0, 7, 1), fill(1.0, 7, 1), fill(1.0, 7, 1), fill(0.0, 7, 1),
    )
    @test_throws ArgumentError remap_hwsd_layers(stocks.carbon; deep_rule = :invent)
end
