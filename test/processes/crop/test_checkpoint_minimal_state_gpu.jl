using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA minimal crop checkpoint classification" begin
    crop = init_crop(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)
    managed_land = init_managed_land(1, CuArray)
    pet = init_pet(1, CuArray)
    state = test_model_state(crop, soil; managed_land, pet)
    restart = crop_restart_payload(crop)

    @test propertynames(crop.state) == (:phenology, :canopy, :carbon, :nitrogen, :water)
    @test :fphu ∉ propertynames(crop.state.phenology)
    @test :phenology_fraction ∉ propertynames(crop.state.canopy)
    @test propertynames(restart) == (:state, :process_memory)
    @test restart.state === crop.state
    @test restart.process_memory.phenology.phu === crop.auxiliary.phenology.phu
    @test restart.process_memory.calendar.sowing_date === crop.auxiliary.calendar.sowing_date
    @test restart.process_memory.calendar.prescribed_sowing_date ===
          crop.auxiliary.calendar.prescribed_sowing_date
    @test :workspace ∉ propertynames(restart)

    crop.state.nitrogen.pending_fertilizer .= 10.0f0
    crop.state.phenology.husum .= 300.0f0
    crop.auxiliary.phenology.phu .= 1000.0f0
    crop.auxiliary.phenology.fphu .= 0.0f0
    fertilizer!(state, managed_land, state, 2)
    synchronize()
    @test Array(crop.fluxes.nitrogen.prescribed_fertilizer_input) == Float32[10]
    @test Array(crop.state.nitrogen.pending_fertilizer) == Float32[0]

    crop.state.canopy.lai .= 0.4f0 * cft1.laimax
    crop.state.phenology.is_growing .= 1
    albedo!(cft1, state, state, pet)
    synchronize()
    green_fraction = 1.0f0 - exp(-cft1.lightextcoeff * 0.4f0 * cft1.laimax)
    expected_canopy = green_fraction * cft1.albedo_leaf +
        (1.0f0 - green_fraction) * 0.3f0
    @test Array(crop.auxiliary.canopy.albedo) ≈ Float32[expected_canopy]
    @test Array(pet.albedo) ≈ Float32[expected_canopy]
end
