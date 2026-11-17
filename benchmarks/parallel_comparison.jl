"""
Serial/parallel LP comparison with direct RM Clp_jll and HiGHS_jll backends.
JuMP comparisons use Clp.jl/MOI and HiGHS.jl/MOI, each backed by the same JLL as RM.

    julia --startup-file=no --project=benchmarks -t 8 \
        benchmarks/parallel_comparison.jl OUTPUT [passes=11] [workers=8] [scenarios=32] [stocks=9]

One persistent problem per worker; contiguous scenario chunks; same total
points for serial and parallel. ThreadTeams reuses persistent workers. Raw
samples include team dispatch/join, Julia allocations and GC time. Model construction, validation and warmup are excluded.
"""
module ParallelSolverComparison
include("reproduce.jl")
const B = ReproducibleBenchmarks
const RM = B.RM
const SV = B.SV
const JuMP = B.JuMP
using Statistics: median
using TOML
using Dates
import ThreadTeams
import Clp

const PATHS = (:rm_clp, :jump_clp_flat, :jump_clp_reset_inner, :jump_clp_named_guarded, :jump_clp_named_unconditional,
               :jump_clp_default, :rm_highs, :jump_flat, :jump_named_guarded, :jump_named_unconditional)
const CLP_WARM_OPTIONS = ("LogLevel"=>0, "SolveType"=>0, "PresolveType"=>1,
                          "Perturbation"=>50, "PrimalTolerance"=>1e-7, "DualTolerance"=>1e-7)

# SolveType=0 selects dual simplex; PresolveType=1 disables presolve in ClpSolve.
# Also retain an unmodified Clp.Optimizer path to expose the option choice.
clp_optimizer() = JuMP.optimizer_with_attributes(Clp.Optimizer, CLP_WARM_OPTIONS...)

struct RMPath{P}
    problem::P
end
struct FlatJuMPPath{S}
    state::S
end
struct ResetClpPath{S}
    state::S
end
struct NamedJuMPPath{Guarded,S}
    state::S
end
NamedJuMPPath{G}(state::S) where {G,S} = NamedJuMPPath{G,S}(state)

struct Worker{S,P}
    state::S
    points::P
end

make_path(::Val{:rm_clp}, data, l, point) =
    RMPath(RM.LinearProblem(B.reactive_const_data(data), point; backend=RM.ClpBackend()))
make_path(::Val{:rm_highs}, data, l, point) =
    RMPath(RM.LinearProblem(B.reactive_const_data(data), point;
        backend=RM.HiGHSBackend(; (Symbol(k)=>v for (k,v) in B.SOLVER_OPTIONS)...)))
make_path(::Val{:jump_flat}, data, l, point) = FlatJuMPPath(B.jump_reference(data, l, point))
make_path(::Val{:jump_named_guarded}, data, l, point) =
    NamedJuMPPath{true}(B.jump_user_model(data, point))
make_path(::Val{:jump_named_unconditional}, data, l, point) =
    NamedJuMPPath{false}(B.jump_user_model(data, point))

make_path(::Val{:jump_clp_flat}, data, l, point) =
    FlatJuMPPath(B.jump_reference(data, l, point; optimizer=clp_optimizer()))
make_path(::Val{:jump_clp_reset_inner}, data, l, point) =
    ResetClpPath(B.jump_reference(data, l, point; optimizer=clp_optimizer()))
make_path(::Val{:jump_clp_named_guarded}, data, l, point) =
    NamedJuMPPath{true}(B.jump_user_model(data, point; optimizer=clp_optimizer()))
make_path(::Val{:jump_clp_named_unconditional}, data, l, point) =
    NamedJuMPPath{false}(B.jump_user_model(data, point; optimizer=clp_optimizer()))
make_path(::Val{:jump_clp_default}, data, l, point) =
    FlatJuMPPath(B.jump_reference(data, l, point; optimizer=Clp.Optimizer))

solve_point!(path::RMPath, point) = RM.objective_value(path.problem, point)
solve_point!(path::FlatJuMPPath, point) = B.jump_solve!(path.state, point)
function solve_point!(path::ResetClpPath, point)
    # Detach before setters to avoid the non-incremental adapter's exception
    # path. The complete native LP is still reloaded on the next optimization.
    B.MOI.Utilities.reset_optimizer(JuMP.backend(path.state.model).optimizer.model)
    return B.jump_solve!(path.state, point)
end
solve_point!(path::NamedJuMPPath{G}, point) where G =
    B.jump_user_solve!(path.state, point; predicated=G)

# Function barrier: state and point types remain concrete inside the point loop.
function run_worker!(worker::Worker)
    checksum = 0.0
    for point in worker.points
        checksum += solve_point!(worker.state, point)
    end
    checksum
end

# A reusable callable and partials buffer: no per-traversal task/vector creation.
struct TeamJob{W}
    workers::W
    partials::Vector{Float64}
end
TeamJob(workers) = TeamJob(workers, zeros(length(workers)))
function (job::TeamJob)(index)
    job.partials[index] = run_worker!(job.workers[index])
    nothing
end
function run_group!(job::TeamJob)
    n = length(job.workers)
    n == 1 && return run_worker!(only(job.workers))
    ThreadTeams.foreach_chunk(job, n, n)
    sum(job.partials)
end

function primal_values(path::RMPath, l)
    solution = RM.solve(path.problem, path.problem.parameter_snapshot)
    x = zeros(l.ncols)
    for r in eachindex(l.initial_volume)
        x[l.initial_volume[r]] = RM.primal(solution, :initial_volumes, r)
    end
    for k in CartesianIndices(l.volume)
        x[l.volume[k]] = RM.primal(solution, :volumes, Tuple(k))
    end
    for k in CartesianIndices(l.debit)
        x[l.debit[k]] = RM.primal(solution, :debits, Tuple(k))
    end
    for k in CartesianIndices(l.group_debit)
        l.group_debit[k] == 0 && continue
        x[l.group_debit[k]] = RM.primal(solution, :group_debits, Tuple(k))
    end
    x[l.vB] = RM.primal(solution, :vB)
    x
end
primal_values(path::Union{FlatJuMPPath,ResetClpPath}, l) = JuMP.value.(path.state.x)
function primal_values(path::NamedJuMPPath, l)
    s = path.state
    x = zeros(l.ncols)
    for r in eachindex(l.initial_volume)
        x[l.initial_volume[r]] = JuMP.value(s.V0[r])
    end
    for k in CartesianIndices(l.volume)
        x[l.volume[k]] = JuMP.value(s.V[Tuple(k)...])
    end
    for k in CartesianIndices(l.debit)
        x[l.debit[k]] = JuMP.value(s.D[Tuple(k)...])
    end
    for k in CartesianIndices(l.group_debit)
        l.group_debit[k] == 0 && continue
        x[l.group_debit[k]] = JuMP.value(s.q[Tuple(k)...])
    end
    x[l.vB] = JuMP.value(s.vB)
    x
end

shape(path::RMPath, l) = Dict("columns"=>l.ncols, "explicit_affine_rows"=>l.nrows)
shape(path::Union{FlatJuMPPath,ResetClpPath,NamedJuMPPath}, l) = Dict(
    "columns"=>JuMP.num_variables(path.state.model),
    "explicit_affine_rows"=>JuMP.num_constraints(path.state.model; count_variable_in_set_constraints=false))

kernel_allocations(worker::Worker) = @allocated run_worker!(worker)

function build_group(path::Val{P}, data, l, chunks, reference) where P
    workers = [Worker(make_path(path, data, l, first(points)), points) for points in chunks]
    all(w -> isconcretetype(typeof(w)) && isconcretetype(eltype(w.points)), workers) ||
        error("worker state must be concrete")
    isconcretetype(eltype(workers)) || error("worker container must retain its state type")
    worst_objective, worst_residual = 0.0, 0.0
    for worker in workers, point in worker.points
        value = solve_point!(worker.state, point)
        ref = reference[point]
        isapprox(value, ref; atol=B.OBJECTIVE_ATOL, rtol=B.OBJECTIVE_RTOL) ||
            error("objective mismatch on $P")
        worst_objective = max(worst_objective, abs(value-ref))
        worst_residual = max(worst_residual, B.residuals(data, l, point,
            primal_values(worker.state, l), value))
    end
    model_shape = shape(first(workers).state, l)
    model_shape["columns"] == l.ncols || error("column count differs")
    model_shape["explicit_affine_rows"] == l.nrows || error("row count differs")
    # Warm serial and task paths before either allocation or time measurements.
    foreach(run_worker!, workers)
    execution = TeamJob(workers)
    expected = sum(reference[p] for w in workers for p in w.points)
    isapprox(run_group!(execution), expected; atol=B.OBJECTIVE_ATOL, rtol=B.OBJECTIVE_RTOL) ||
        error("parallel checksum mismatch")
    kernel_bytes = [kernel_allocations(w) for w in workers]
    (; path=P, workers, execution, expected, kernel_bytes, model_shape,
       validation=Dict("points"=>sum(length(w.points) for w in workers),
                       "worst_objective"=>worst_objective,
                       "worst_scaled_primal_residual"=>worst_residual),
       samples=Dict{String,Any}[])
end

function sample!(group)
    GC.gc()
    sample = @timed run_group!(group.execution)
    isapprox(sample.value, group.expected; atol=B.OBJECTIVE_ATOL, rtol=B.OBJECTIVE_RTOL) ||
        error("unstable timed checksum on $(group.path)")
    push!(group.samples, Dict("seconds"=>sample.time, "julia_bytes"=>sample.bytes,
          "gc_seconds"=>sample.gctime, "checksum"=>sample.value))
    nothing
end

function summarize(group, npoints)
    seconds = median(s["seconds"] for s in group.samples)
    Dict("path"=>String(group.path), "workers"=>length(group.workers),
         "points_per_worker"=>[length(w.points) for w in group.workers],
         "median_us_per_lp"=>1e6*seconds/npoints,
         "median_julia_bytes_per_lp"=>median(s["julia_bytes"] for s in group.samples)/npoints,
         "median_gc_seconds"=>median(s["gc_seconds"] for s in group.samples),
         "kernel_bytes_per_worker_sweep"=>group.kernel_bytes,
         "model_shape"=>group.model_shape, "validation"=>group.validation,
         "samples"=>group.samples)
end

function archive_source(output)
    source = joinpath(output, "source")
    mkpath(source)
    root = dirname(@__DIR__)
    for name in ("src", "ext", "examples")
        cp(joinpath(root, name), joinpath(source, name))
    end
    for name in ("Project.toml", "LICENSE")
        cp(joinpath(root, name), joinpath(source, name))
    end
    bench = joinpath(source, "benchmarks")
    mkpath(bench)
    for file in readdir(@__DIR__)
        (endswith(file, ".jl") || endswith(file, ".py") || file in ("Project.toml", "Manifest.toml", "PROVENANCE.md")) || continue
        cp(joinpath(@__DIR__, file), joinpath(bench, file))
    end
end

# Clp.jl 1.4 reports a non-incremental MOI interface. Record the observed
# nested caching transition instead of treating a retained JuMP model as native reuse.
function clp_cache_observation(data, l, points)
    path = make_path(Val(:jump_clp_flat), data, l, first(points))
    solve_point!(path, first(points))
    cache = JuMP.backend(path.state.model).optimizer.model
    cache isa B.MOI.Utilities.CachingOptimizer || error("unexpected Clp cache layout")
    before = string(cache.state)
    B.jump_update!(path.state, points[2])
    after_update = string(cache.state)
    JuMP.optimize!(path.state.model)
    Dict("incremental_interface"=>B.MOI.supports_incremental_interface(Clp.Optimizer()),
         "after_solve"=>before, "after_update"=>after_update,
         "after_reoptimize"=>string(cache.state),
         "native_copy_policy"=>"Clp.jl MOI.copy_to calls Clp_loadProblem when reattaching")
end

function campaign(output; passes=11, workers=8, scenarios=32, stocks=9,
                  presets=(SV.VALLEY_SMALL, SV.VALLEY_MEDIUM, SV.VALLEY_LARGE))
    passes > 0 || throw(ArgumentError("passes must be positive"))
    2 <= workers <= min(ThreadTeams.nworkers(), scenarios) ||
        throw(ArgumentError("need 2 <= workers <= min(ThreadTeams workers, scenarios)"))
    stocks >= 2 || throw(ArgumentError("need at least two stocks to measure updates"))
    scenarios > 0 || throw(ArgumentError("scenarios must be positive"))
    ispath(output) && !isempty(readdir(output)) && throw(ArgumentError("output directory must be new or empty"))
    mkpath(output)
    metadata = B.metadata()
    metadata["timing_workers"] = [1, workers]
    metadata["paths"] = collect(String.(PATHS))
    metadata["timing_scope"] = "update and solve; parallel includes persistent ThreadTeams dispatch/join, no per-traversal task creation; construction, validation, warmup and pre-sample GC excluded"
    metadata["parallel_runtime"] = "ThreadTeams.foreach_chunk with reusable callable/partials"
    metadata["threadteams_version"] = string(pkgversion(ThreadTeams))
    metadata["threadteams_workers_available"] = ThreadTeams.nworkers()
    metadata["threadteams_spin_seconds"] = ThreadTeams.spin_time()
    metadata["jump_clp_warm_options"] = Dict(CLP_WARM_OPTIONS)
    metadata["jump_clp_default_options"] = "Clp.Optimizer defaults, silent; no warm option overrides"
    same_clp_library = realpath(String(RM._CLP_LIB)) == realpath(String(Clp.libClp))
    same_clp_library || error("RM and JuMP/Clp must use the same native library")
    metadata["same_clp_library"] = same_clp_library
    metadata["clp_library_sha256"] = bytes2hex(B.sha256(read(String(Clp.libClp))))
    metadata["order"] = "contiguous scenarios per worker; bisection stock order within scenario; rotate configuration order per trial"
    metadata["clp_interface"] = "RM direct Clp_jll; JuMP via Clp.jl/MOI using that exact Clp_jll binary"
    metadata["highs_interfaces"] = "RM direct HiGHS_jll; JuMP via existing HiGHS.jl/MOI using HiGHS_jll"
    metadata["kernel_allocation_scope"] = "separate warmed serial sweep per worker, excludes task harness; aggregate samples include team dispatch"
    metadata["gc_scope"] = "process-wide GC seconds observed by @timed; not a per-task or causal attribution"
    metadata["input_provenance"] = "independent-synthetic-2026-11-13; see benchmarks/PROVENANCE.md"
    report = Dict{String,Any}("metadata"=>metadata, "runs"=>Any[])
    for spec in presets
        data = SV.valley_data(spec; nregions=1)
        l = SV.layout(data)
        all_points = B.production_points(data, 1:scenarios; stocks_per_scenario=stocks)
        reference_model = make_path(Val(:rm_highs), data, l, first(all_points))
        reference = Dict(p=>solve_point!(reference_model, p) for p in all_points)
        clp_cache = clp_cache_observation(data, l, all_points)
        configurations = Any[]
        for count in (1, workers)
            chunks = [B.production_points(data, r; stocks_per_scenario=stocks)
                      for r in B.scenario_chunks(scenarios, count)]
            for path in PATHS
                push!(configurations, build_group(Val(path), data, l, chunks, reference))
            end
        end
        for trial in 1:passes
            for offset in eachindex(configurations)
                index = mod1(offset + trial - 1, length(configurations))
                sample!(configurations[index])
            end
            println(spec.name, ": sample ", trial, "/", passes)
            flush(stdout)
        end
        measurements = [summarize(c, length(all_points)) for c in configurations]
        run = Dict("preset"=>spec.name, "ncols"=>l.ncols, "nrows"=>l.nrows,
                   "scenarios"=>scenarios, "stocks_per_scenario"=>stocks, "npoints"=>length(all_points),
                   "nregions"=>1, "seed"=>spec.seed, "passes"=>passes,
                   "jump_clp_cache"=>clp_cache, "measurements"=>measurements)
        push!(report["runs"], run)
        # Save completed presets even if a later preset is interrupted.
        open(joinpath(output, "results.toml"), "w") do io
            TOML.print(io, report; sorted=true)
        end
        for m in measurements
            println(rpad(m["path"], 26), " workers=", m["workers"], " ",
                    round(m["median_us_per_lp"]; digits=2), " us/LP; ",
                    round(m["median_julia_bytes_per_lp"]; digits=2), " Julia B/LP")
        end
        flush(stdout)
    end
    archive_source(output)
    println("Results: ", abspath(output))
    report
end

function main(args=ARGS)
    isempty(args) && throw(ArgumentError("provide a new output directory"))
    passes = length(args) > 1 ? parse(Int, args[2]) : 11
    workers = length(args) > 2 ? parse(Int, args[3]) : 8
    scenarios = length(args) > 3 ? parse(Int, args[4]) : 32
    stocks = length(args) > 4 ? parse(Int, args[5]) : 9
    campaign(args[1]; passes, workers, scenarios, stocks)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ParallelSolverComparison.main()
end
