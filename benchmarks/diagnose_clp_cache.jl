"""
Test explicit cache detachment for the non-incremental Clp.jl adapter.

    julia --startup-file=no --project=benchmarks \
        benchmarks/diagnose_clp_cache.jl OUTPUT [passes=11]

Ordinary JuMP updates are the control. The other paths explicitly reset an MOI
cache before updating. They still reload the full native LP, without patching
JuMP, MOI or Clp.jl. The inner-cache path depends on the recorded adapter layout.
"""
module ClpCacheDiagnostics
include("parallel_comparison.jl")
const P = ParallelSolverComparison
const B = P.B
const MOI = B.MOI
using Statistics: median
using TOML
using Profile
const MODES = (:ordinary, :reset_inner, :reset_outer)

reset_before_update!(state, ::Val{:ordinary}) = nothing
function reset_before_update!(state, ::Val{:reset_inner})
    MOI.Utilities.reset_optimizer(P.JuMP.backend(state.model).optimizer.model)
end
function reset_before_update!(state, ::Val{:reset_outer})
    MOI.Utilities.reset_optimizer(P.JuMP.backend(state.model))
end
function sweep!(path, points, mode; reference=nothing, data=nothing, l=nothing)
    checksum = 0.0
    for point in points
        reset_before_update!(path.state, mode)
        value = P.solve_point!(path,point)
        if reference !== nothing
            @assert isapprox(value,reference[point];atol=B.OBJECTIVE_ATOL,rtol=B.OBJECTIVE_RTOL)
            B.residuals(data,l,point,P.primal_values(path,l),value)
        end
        checksum += value
    end
    checksum
end

function campaign(output;passes=11)
    passes > 0 || throw(ArgumentError("passes must be positive"))
    ispath(output) && !isempty(readdir(output)) && error("output must be new or empty")
    mkpath(output)
    P.ThreadTeams.shutdown!()
    metadata = B.metadata()
    metadata["timing_scope"] = "serial paired ordinary updates versus explicit MOI cache reset before every update; reset, update, reload, solve and query included; warmup, validation and GC excluded"
    metadata["modes"] = collect(String.(MODES))
    metadata["same_clp_library"] = realpath(String(P.RM._CLP_LIB)) == realpath(String(P.Clp.libClp))
    @assert metadata["same_clp_library"]
    report = Dict("metadata"=>metadata,"runs"=>Any[])
    for spec in (P.SV.VALLEY_SMALL,P.SV.VALLEY_MEDIUM,P.SV.VALLEY_LARGE)
        data = P.SV.valley_data(spec;nregions=1)
        l = P.SV.layout(data)
        points = B.production_points(data,1:32;stocks_per_scenario=9)
        model = P.make_path(Val(:rm_highs),data,l,first(points))
        reference = Dict(p=>P.solve_point!(model,p) for p in points)
        expected = sum(values(reference))
        paths = Dict(mode=>P.make_path(Val(:jump_clp_flat),data,l,first(points)) for mode in MODES)
        for mode in MODES
            sweep!(paths[mode],points,Val(mode);reference,data,l)
            sweep!(paths[mode],points,Val(mode))
        end
        samples = Dict(mode=>Any[] for mode in MODES)
        for trial in 1:passes, offset in eachindex(MODES)
            mode = MODES[mod1(offset+trial-1,length(MODES))]
            GC.gc()
            sample = @timed sweep!(paths[mode],points,Val(mode))
            @assert isapprox(sample.value,expected;atol=B.OBJECTIVE_ATOL,rtol=B.OBJECTIVE_RTOL)
            push!(samples[mode],Dict("seconds"=>sample.time,"bytes"=>sample.bytes,
                 "gc_seconds"=>sample.gctime,"checksum"=>sample.value,"trial"=>trial))
        end
        results = Any[]
        for mode in MODES
            observations = samples[mode]
            result = Dict("mode"=>String(mode),"samples"=>observations,
                "median_us_per_lp"=>1e6*median(s["seconds"] for s in observations)/length(points),
                "median_bytes_per_lp"=>median(s["bytes"] for s in observations)/length(points))
            push!(results,result)
            println(spec.name," ",mode,": ",round(result["median_us_per_lp"];digits=2)," us/LP; ",round(result["median_bytes_per_lp"];digits=2)," B/LP")
        end
        # Qualitative profile: keep caller thread only, since parked Julia threads
        # otherwise dominate the report. Profiling is outside all timing samples.
        Profile.clear()
        Profile.@profile sweep!(paths[:ordinary],points,Val(:ordinary))
        open(joinpath(output,spec.name*"-profile.txt"),"w") do io
            Profile.print(io;format=:flat,sortedby=:count,C=true,mincount=2,
                          threads=Threads.threadid(),maxdepth=40)
        end
        push!(report["runs"],Dict("preset"=>spec.name,"points"=>length(points),
            "passes"=>passes,"ncols"=>l.ncols,"nrows"=>l.nrows,"measurements"=>results))
        open(joinpath(output,"results.toml"),"w") do io
            TOML.print(io,report;sorted=true)
        end
        flush(stdout)
    end
    P.archive_source(output)
    report
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ClpCacheDiagnostics.campaign(ARGS[1];passes=length(ARGS)>1 ? parse(Int,ARGS[2]) : 11)
end
