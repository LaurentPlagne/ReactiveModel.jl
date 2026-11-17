using Test
using TOML
include(joinpath(@__DIR__, "reproduce.jl"))
const B = ReproducibleBenchmarks

@testset "Standalone benchmark environment and validation" begin
    project = TOML.parsefile(joinpath(@__DIR__, "Project.toml"))
    @test !haskey(project["deps"], "ClpJL")
    @test !haskey(project["deps"], "TinyHiGHS")
    @test project["sources"] == Dict("ReactiveModel"=>Dict("path"=>".."))
    @test !B.BENCH_TINY
    @test !(:rm_tiny in B.PATHS)
    data = B.SV.valley_data(B.SV.VALLEY_SMALL; nregions=1)
    points = B.SV.valley_points(data; nscenarios=2, order=:refinement, stocks_per_scenario=3)
    models = B.build_models(data, points)
    @test models.rm_tiny === nothing
    @test models.rm_clp.backend isa B.RM.ClpBackend
    validation = B.validate_models(models, points)
    @test validation.points == 6
    @test validation.worst_scaled_primal_residual <= B.FEASIBILITY_TOL
    point = first(points)
    value = B.solve_point!(models,:rm_highs,point; read_dual=false)
    x = B.primal_vector(models,:rm_highs,point)
    @test_throws ErrorException B.residuals(data,models.l,point,x,value+100.0)
    x[models.l.initial_volume[1]] += 100.0
    @test_throws ErrorException B.residuals(data,models.l,point,x,value)
    @test_throws ArgumentError B.campaign(mktempdir(); passes=0)
    occupied = mktempdir()
    write(joinpath(occupied,"keep"),"user data")
    @test_throws ArgumentError B.campaign(occupied; passes=1)
    @test read(joinpath(occupied,"keep"),String) == "user data"
    @test_throws ArgumentError B.scenario_chunks(1, 2)
end
