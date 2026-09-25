using Test

@testset "C3/C4 daily process-order contract" begin
    driver = "daily_crop.jl"
    source = read(joinpath(@__DIR__, "..", "..", "src", "simulations", driver), String)
    daily_start = findfirst("function _daily_crop!", source)
    daily_start === nothing && error("missing common daily crop driver")
    source = source[first(daily_start):end]

    position = call -> begin
        pattern = Regex("(?m)^[ \\t]*" * call * "[ \\t]*\\(")
        location = findfirst(pattern, source)
        location === nothing && error("missing $call call in $driver")
        first(location)
    end

    @test position("update_climbuf!") < position("cultivate!")
    @test position("litter_tillage!") < position("tillage_hydraulics!") <
          position("pedotransfer!")
    @test position("_pathway_albedo!") < position("petpar!") < position("snow!")
    @test position("pedotransfer!") <
          position("update_surface_litter_properties!") <
          position("soil_temperature!")
    @test position("soil_temperature!") < position("soil_cn_decomposition!")
    @test position("soil_cn_decomposition!") < position("nitrogen_deposition!") <
          position("phenology_crop!")
    @test position("cultivate!") < position("phenology_crop!") < position("harvest_crop!")
    @test position("harvest_crop!") < position("route_harvest_residues!") <
          position("interception!") < position("soil_infiltration!")
    @test position("soil_infiltration!") < position("_pathway_apar!") <
          position("temp_stress") < position("photosynthesis!") <
          position("transpiration!") < position("solve_lambda!")
    @test position("solve_lambda!") < position("crop_carbon!") <
          position("terminate_failed_crop!") < position("evaporation!")
    @test position("evaporation!") <
          position("soil_evapotranspiration!") < position("post_crop_nitrogen_losses!")
end
