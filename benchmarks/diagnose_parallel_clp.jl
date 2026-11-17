"""
Controlled diagnostics for scheduler placement and JuMP/Clp model reuse.

    julia --startup-file=no --project=benchmarks -t 8 \
        benchmarks/diagnose_parallel_clp.jl OUTPUT [passes=11]

Uses the same models, chunks and JLL binaries for every scheduler. Results are
diagnostics on one machine, not estimates of universal runtime overhead.
"""
module ParallelClpDiagnostics
include("parallel_comparison.jl")
const P = ParallelSolverComparison
const B = P.B
const TT = P.ThreadTeams
const MOI = B.MOI
using Statistics: median
using TOML

const MODES = (:spawn, :team_spin100, :team_spin0, :range_spin0)

function dispatch!(job, ::Val{:spawn})
    tasks = [Threads.@spawn job(i) for i in eachindex(job.workers)]
    foreach(fetch, tasks)
    sum(job.partials)
end
dispatch!(job, ::Union{Val{:team_spin100},Val{:team_spin0}}) = P.run_group!(job)
function dispatch!(job, ::Val{:range_spin0})
    TT.foreach_range(length(job.workers); threads=length(job.workers), schedule=:static) do lo, hi
        for i in lo:hi
            job(i)
        end
    end
    sum(job.partials)
end

# Stop the parked pool before @spawn, so its workers cannot compete with tasks.
# For team modes, startup and one complete warm traversal are outside the timer.
function prepare!(job, mode)
    TT.shutdown!()
    TT.set_spin_time!(mode === :team_spin100 ? 1e-4 : 0.0)
    dispatch!(job, Val(mode))
end

struct PlacementProbe{J}
    job::J
    tids::Vector{Int}
    seconds::Vector{Float64}
end
function (probe::PlacementProbe)(i)
    probe.tids[i] = Threads.threadid()
    start = time_ns()
    probe.job(i)
    probe.seconds[i] = (time_ns()-start)*1e-9
    nothing
end

function placement(job, mode)
    probe = PlacementProbe(job, zeros(Int,length(job.workers)), zeros(length(job.workers)))
    if mode === :spawn
        foreach(fetch, [Threads.@spawn probe(i) for i in eachindex(job.workers)])
    elseif mode === :range_spin0
        TT.foreach_range(length(job.workers); threads=length(job.workers), schedule=:static) do lo, hi
            foreach(probe, lo:hi)
        end
    else
        TT.foreach_chunk(probe, length(job.workers), length(job.workers))
    end
    Dict("thread_ids"=>probe.tids, "chunk_seconds"=>probe.seconds)
end

function paired_schedulers(data, l, points, reference; passes, workers=8)
    chunks = [B.production_points(data, r; stocks_per_scenario=9)
              for r in B.scenario_chunks(32, workers)]
    results = Any[]
    for path in (:rm_clp, :rm_highs, :jump_flat)
        group = P.build_group(Val(path), data, l, chunks, reference)
        samples = Dict(mode=>Any[] for mode in MODES)
        for trial in 1:passes
            for offset in eachindex(MODES)
                mode = MODES[mod1(offset+trial-1,length(MODES))]
                prepare!(group.execution, mode)
                GC.gc()
                sample = @timed dispatch!(group.execution, Val(mode))
                @assert isapprox(sample.value, group.expected; atol=B.OBJECTIVE_ATOL,rtol=B.OBJECTIVE_RTOL)
                push!(samples[mode], Dict("seconds"=>sample.time, "bytes"=>sample.bytes,
                      "gc_seconds"=>sample.gctime, "checksum"=>sample.value, "trial"=>trial))
            end
        end
        for mode in MODES
            prepare!(group.execution,mode)
            probe = placement(group.execution,mode)
            values = samples[mode]
            result = Dict("path"=>String(path), "mode"=>String(mode), "samples"=>values,
                "median_us_per_lp"=>1e6*median(s["seconds"] for s in values)/length(points),
                "median_bytes_per_lp"=>median(s["bytes"] for s in values)/length(points),
                "placement_probe"=>probe, "validation"=>group.validation,
                "kernel_bytes_per_worker_sweep"=>group.kernel_bytes)
            push!(results,result)
            println(path," / ",mode,": ",round(result["median_us_per_lp"];digits=2)," us/LP")
        end
        TT.shutdown!()
        flush(stdout)
    end
    results
end

# Explicit reattachment isolates the native reload and index-map construction.
# This uses MOI's attach_optimizer, not a patched wrapper. The ordinary optimize!
# control below checks that decomposing this operation does not change totals.
function clp_phases!(path, points, reference)
    state = path.state
    cache = P.JuMP.backend(state.model).optimizer.model
    seconds, bytes = zeros(4), zeros(Int,4)
    iterations = Int[]
    detached = 0
    for point in points
        update = @timed B.jump_update!(state,point)
        detached += cache.state == MOI.Utilities.EMPTY_OPTIMIZER
        attach = @timed if cache.state == MOI.Utilities.EMPTY_OPTIMIZER
            MOI.Utilities.attach_optimizer(cache)
        end
        solve = @timed P.JuMP.optimize!(state.model)
        query = @timed P.JuMP.objective_value(state.model)
        @assert P.JuMP.termination_status(state.model) == MOI.OPTIMAL
        @assert isapprox(query.value,reference[point];atol=B.OBJECTIVE_ATOL,rtol=B.OBJECTIVE_RTOL)
        push!(iterations,Int(P.Clp.Clp_numberIterations(cache.optimizer)))
        for (i,sample) in enumerate((update,attach,solve,query))
            seconds[i] += sample.time
            bytes[i] += sample.bytes
        end
    end
    Dict("phase_seconds"=>seconds,"phase_bytes"=>bytes,"iterations"=>iterations,
         "native_reloads"=>detached)
end

function clp_decomposition(data,l,points,reference;passes)
    path = P.make_path(Val(:jump_clp_flat),data,l,first(points))
    worker = P.Worker(path,points)
    P.run_worker!(worker)
    clp_phases!(path,points,reference)
    ordinary, decomposed = Any[], Any[]
    for trial in 1:passes
        # Alternate ordinary and instrumented traversals in each pair.
        for mode in (isodd(trial) ? (:ordinary,:phases) : (:phases,:ordinary))
            GC.gc()
            if mode === :ordinary
                sample = @timed P.run_worker!(worker)
                @assert isapprox(sample.value,sum(values(reference));atol=B.OBJECTIVE_ATOL,rtol=B.OBJECTIVE_RTOL)
                push!(ordinary,Dict("seconds"=>sample.time,"bytes"=>sample.bytes,"gc_seconds"=>sample.gctime))
            else
                push!(decomposed,clp_phases!(path,points,reference))
            end
        end
    end
    # No update/reload: retained JuMP/Clp state, reoptimized at the last point.
    # Deliberately not compared with RM's repeated-point cache, which skips solve.
    retained = Any[]
    for _ in 1:passes
        sample = @timed for _ in eachindex(points)
            P.JuMP.optimize!(path.state.model)
        end
        push!(retained,Dict("seconds"=>sample.time,"bytes"=>sample.bytes))
    end
    rm = P.make_path(Val(:rm_clp),data,l,first(points))
    rm_iterations = Int[]
    for point in points
        value = P.solve_point!(rm,point)
        @assert isapprox(value,reference[point];atol=B.OBJECTIVE_ATOL,rtol=B.OBJECTIVE_RTOL)
        push!(rm_iterations,Int(ccall((:Clp_numberIterations,P.RM._CLP_LIB),Cint,(Ptr{Cvoid},),rm.problem.backend.inner)))
    end
    result = Dict("phases"=>["update_and_detach","copy_reattach","optimize_attached","objective_query"],
        "ordinary_samples"=>ordinary,"instrumented_samples"=>decomposed,
        "retained_no_update_samples"=>retained,"rm_iterations_first_traversal"=>rm_iterations,
        "ordinary_median_us_per_lp"=>1e6*median(s["seconds"] for s in ordinary)/length(points),
        "instrumented_median_us_per_lp"=>1e6*median(sum(s["phase_seconds"]) for s in decomposed)/length(points),
        "median_phase_us_per_lp"=>[1e6*median(s["phase_seconds"][i] for s in decomposed)/length(points) for i in 1:4],
        "median_phase_bytes_per_lp"=>[median(s["phase_bytes"][i] for s in decomposed)/length(points) for i in 1:4],
        "retained_median_us_per_lp"=>1e6*median(s["seconds"] for s in retained)/length(points))
    println("Clp ordinary / instrumented / retained us/LP: ",result["ordinary_median_us_per_lp"]," / ",result["instrumented_median_us_per_lp"]," / ",result["retained_median_us_per_lp"])
    println("Clp phase us/LP: ",result["median_phase_us_per_lp"])
    result
end

function campaign(output;passes=11)
    passes > 0 || throw(ArgumentError("passes must be positive"))
    TT.nworkers() >= 8 || error("need at least 8 default-pool workers")
    ispath(output) && !isempty(readdir(output)) && error("output must be new or empty")
    mkpath(output)
    metadata = B.metadata()
    metadata["timing_scope"] = "same-process paired scheduler diagnostics; same 8 models and 288 points; startup, warm traversal and GC outside each timer"
    metadata["scheduler_modes"] = collect(String.(MODES))
    metadata["threadteams_version"] = string(pkgversion(TT))
    metadata["default_pool_tids"] = collect(Threads.threadpooltids(:default))
    metadata["interactive_pool_tids"] = collect(Threads.threadpooltids(:interactive))
    metadata["caller_tid"] = Threads.threadid()
    metadata["same_clp_library"] = realpath(String(P.RM._CLP_LIB)) == realpath(String(P.Clp.libClp))
    @assert metadata["same_clp_library"]
    metadata["clp_library_sha256"] = bytes2hex(B.sha256(read(String(P.Clp.libClp))))
    metadata["clp_phase_scope"] = "instrumented serial update/detach, explicit MOI attach_optimizer, JuMP optimize! on attached cache, objective query; point-wise timers; ordinary sweep control"
    report = Dict("metadata"=>metadata,"runs"=>Any[])
    oldspin = TT.spin_time()
    try
        for spec in (P.SV.VALLEY_SMALL,P.SV.VALLEY_MEDIUM,P.SV.VALLEY_LARGE)
            data = P.SV.valley_data(spec;nregions=1)
            l = P.SV.layout(data)
            points = B.production_points(data,1:32;stocks_per_scenario=9)
            model = P.make_path(Val(:rm_highs),data,l,first(points))
            reference = Dict(p=>P.solve_point!(model,p) for p in points)
            println("Preset: ",spec.name); flush(stdout)
            schedulers = paired_schedulers(data,l,points,reference;passes)
            TT.shutdown!()
            clp = clp_decomposition(data,l,points,reference;passes)
            push!(report["runs"],Dict("preset"=>spec.name,"ncols"=>l.ncols,"nrows"=>l.nrows,
                "points"=>length(points),"passes"=>passes,"schedulers"=>schedulers,"clp"=>clp))
            open(joinpath(output,"results.toml"),"w") do io
                TOML.print(io,report;sorted=true)
            end
            flush(stdout)
        end
        P.archive_source(output)
    finally
        TT.shutdown!()
        TT.set_spin_time!(oldspin)
    end
    report
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ParallelClpDiagnostics.campaign(ARGS[1];passes=length(ARGS)>1 ? parse(Int,ARGS[2]) : 11)
end
