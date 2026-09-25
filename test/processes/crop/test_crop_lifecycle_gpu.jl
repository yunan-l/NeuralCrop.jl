using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

include("../../helpers/crop_lifecycle_fixture.jl")

@testset "CUDA two-year prescribed crop lifecycle" begin
    cpu = run_lifecycle_fixture(identity)
    gpu = run_lifecycle_fixture(CuArray)
    synchronize()

    for field in (:sowing_event, :harvest_event, :harvesting_mask)
        @test Array(getproperty(gpu.output.calendar, field)) ==
            Array(getproperty(cpu.output.calendar, field))
    end
    @test event_days(gpu.output.calendar.sowing_event) == [100, 465]
    @test event_days(gpu.output.calendar.harvest_event) == [137, 502]
    @test Array(gpu.output.crop.growing_mask) == Array(cpu.output.crop.growing_mask)
    @test Array(gpu.output.crop.fphu) ≈ Array(cpu.output.crop.fphu) rtol = 2.0f-5

    inactive_days = vcat(138:464, 503:730)
    for field in (
        :gpp, :npp, :lambda, :potential_vcmax, :vcmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu,
    )
        @test all(iszero, Array(getproperty(gpu.output.crop, field))[inactive_days, :])
    end
    @test all(iszero, Array(gpu.output.crop.growing_mask)[inactive_days, :])
    for field in (:harvesting_mask, :sowing_event, :harvest_event)
        @test all(
            iszero,
            Array(getproperty(gpu.output.calendar, field))[inactive_days, :],
        )
    end

    crop_state = gpu.state.prognostic.crop
    crop_fluxes = gpu.state.fluxes.crop
    crop_auxiliary = gpu.state.auxiliary.crop

    for (container, fields) in (
        (crop_state.phenology, (:vdsum, :husum, :growing_days, :is_growing)),
        (crop_auxiliary.phenology, (:fphu,)),
        (crop_state.canopy, (:lai, :lai_previous_potential, :laimax_adjusted, :lai_npp_deficit)),
        (crop_auxiliary.canopy,
         (:actual_lai, :flaimax, :fpar, :apar, :canopy_conductance, :canopy_wet)),
        (crop_state.carbon,
         (:biomass, :leaf, :root, :pool, :storage)),
        (crop_fluxes.carbon,
         (:yield, :harvest_export, :npp, :respiration, :gross_assimilation, :net_assimilation,
          :water_limited_assimilation, :leaf_respiration)),
        (crop_state.nitrogen,
         (:total, :leaf, :root, :pool, :storage, :pending_manure,
          :pending_fertilizer, :stress_sum)),
        (crop_fluxes.nitrogen,
         (:uptake, :auto_fertilizer, :seed_input, :prescribed_manure_input,
          :prescribed_fertilizer_input, :harvest_export)),
        (crop_state.water,
         (:demand_sum, :supply_sum)),
        (crop_fluxes.water,
         (:interception, :transpiration_layer)),
        (crop_auxiliary.stress,
         (:nitrogen_demand_total, :nitrogen_demand_leaf,
          :nitrogen_deficit, :water_deficit)),
        (crop_auxiliary.photosynthesis,
         (:potential_vcmax, :vcmax, :nitrogen_limitation, :lambda)),
    )
        for field in fields
            @test all(iszero, Array(getproperty(container, field)))
        end
    end
    @test Array(crop_state.phenology.harvesting) == Bool[false]
    @test Array(crop_state.phenology.harvesting_previous) == Bool[false]
    @test Array(gpu.output.annual.yield) == Float32[0]
    @test Array(crop_state.nitrogen.sufficiency) == Float32[1]
    @test Array(crop_state.water.sufficiency) == Float32[1]
end
