using Test
using TOML
include("parallel_comparison.jl")
const P = ParallelSolverComparison

@testset "Direct JLL backends and independent parallel workspaces" begin
    @test haskey(TOML.parsefile(joinpath(@__DIR__, "Project.toml"))["deps"], "Clp")
    @test haskey(TOML.parsefile(joinpath(@__DIR__, "Project.toml"))["deps"], "ThreadTeams")
    @test realpath(String(P.RM._CLP_LIB)) == realpath(String(P.Clp.libClp))
    data = P.SV.valley_data(P.SV.VALLEY_SMALL; nregions=1)
    l = P.SV.layout(data)
    chunks = [P.B.production_points(data, r; stocks_per_scenario=3)
              for r in P.B.scenario_chunks(2, 2)]
    reference_model = P.make_path(Val(:rm_highs), data, l, first(first(chunks)))
    reference = Dict(p=>P.solve_point!(reference_model, p) for points in chunks for p in points)
    clp = P.build_group(Val(:rm_clp), data, l, chunks, reference)
    @test all(w->w.state.problem.backend isa P.RM.ClpBackend, clp.workers)
    @test clp.workers[1].state.problem.backend !== clp.workers[2].state.problem.backend
    @test isapprox(P.run_group!(clp.execution), sum(values(reference)); atol=1e-7)
    @test clp.validation["worst_scaled_primal_residual"] <= P.B.FEASIBILITY_TOL
    @test clp.model_shape == Dict("columns"=>l.ncols, "explicit_affine_rows"=>l.nrows)
    # Every alternate frontend is checked point-by-point for objective and primal feasibility.
    for path in filter(!=(:rm_clp), P.PATHS)
        group = P.build_group(Val(path), data, l, chunks, reference)
        @test isapprox(P.run_group!(group.execution), clp.expected; atol=1e-7)
        @test group.validation["worst_scaled_primal_residual"] <= P.B.FEASIBILITY_TOL
    end
    # Team dispatch covers each chunk exactly once and reuses its partial buffer.
    old_partials = clp.execution.partials
    checksum = P.run_group!(clp.execution)
    @test old_partials === clp.execution.partials
    @test isapprox(sum(old_partials), checksum; atol=1e-7)
    @test_throws ArgumentError P.campaign(mktempdir(); workers=1)
    @test_throws ArgumentError P.campaign(mktempdir(); passes=0, workers=2)
    @test_throws ArgumentError P.campaign(mktempdir(); stocks=1, workers=2)
    occupied = mktempdir()
    write(joinpath(occupied, "keep.txt"), "preserve")
    @test_throws ArgumentError P.campaign(occupied; workers=2)
    @test read(joinpath(occupied, "keep.txt"), String) == "preserve"
    output = mktempdir()
    P.campaign(output; passes=1, workers=2, scenarios=2, stocks=3,
               presets=(P.SV.VALLEY_SMALL,))
    report = TOML.parsefile(joinpath(output, "results.toml"))
    @test report["metadata"]["timing_workers"] == [1, 2]
    @test report["metadata"]["same_clp_library"]
    @test occursin("ThreadTeams", report["metadata"]["parallel_runtime"])
    @test length(only(report["runs"])["measurements"]) == 2length(P.PATHS)
    cache = only(report["runs"])["jump_clp_cache"]
    @test !cache["incremental_interface"]
    @test cache["after_solve"] == "ATTACHED_OPTIMIZER"
    @test cache["after_update"] == "EMPTY_OPTIMIZER"
    @test cache["after_reoptimize"] == "ATTACHED_OPTIMIZER"
    @test isfile(joinpath(output, "source", "benchmarks", "parallel_comparison.jl"))
end
