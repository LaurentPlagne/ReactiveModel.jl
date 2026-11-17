"""
Small, serial, independently checked benchmark campaign.

    julia --startup-file=no --project=benchmarks benchmarks/reproduce.jl [output_dir] [passes]

Writes raw samples, validation residuals, versions, solver options, and the
source revision. Model construction, numerical verification and warmup are
outside timings. Run a larger campaign separately before drawing performance
conclusions; this command also serves as an installation smoke test.
"""
module ReproducibleBenchmarks
include(joinpath(@__DIR__, "valley_comparison.jl"))
using Dates
using Pkg
using SHA
using TOML

const BASELINE_PATHS = (:rm_highs, :rm_clp, :jump, :moi, :poi)
const OBJECTIVE_ATOL = 1e-6
const OBJECTIVE_RTOL = 1e-8
const FEASIBILITY_TOL = 1e-7

function git_output(args...; fallback="unavailable")
    root = dirname(@__DIR__)
    try
        read(pipeline(`git -C $root $args`; stderr=devnull), String)
    catch error
        error isa ProcessFailedException || rethrow()
        fallback
    end
end

function primal_vector(models, path, point)
    l = models.l
    if path in (:rm_highs, :rm_clp)
        sol = RM.solve(getproperty(models, path), point)
        x = zeros(l.ncols)
        for r in eachindex(l.initial_volume)
            x[l.initial_volume[r]] = RM.primal(sol, :initial_volumes, r)
        end
        for k in CartesianIndices(l.volume)
            x[l.volume[k]] = RM.primal(sol, :volumes, Tuple(k))
        end
        for k in CartesianIndices(l.debit)
            x[l.debit[k]] = RM.primal(sol, :debits, Tuple(k))
        end
        for k in CartesianIndices(l.group_debit)
            l.group_debit[k] == 0 && continue
            x[l.group_debit[k]] = RM.primal(sol, :group_debits, Tuple(k))
        end
        x[l.vB] = RM.primal(sol, :vB)
        return x
    end
    state = getproperty(models, path)
    path === :moi ? MOI.get(state.optimizer, MOI.VariablePrimal(), state.x) : JuMP.value.(state.x)
end

function residuals(data, l, point, x, value)
    all(isfinite, x) || error("non-finite primal solution")
    lower, upper = zeros(l.ncols), zeros(l.ncols)
    rhs = zeros(l.nrows)
    costs = zeros(length(l.objective_columns))
    SV.column_bounds!(lower, upper, data, l, point)
    SV.row_rhs!(rhs, data, l, point)
    SV.objective_costs!(costs, data, l, point)
    violation = 0.0
    for j in eachindex(x)
        scale = max(1.0, abs(x[j]), isfinite(lower[j]) ? abs(lower[j]) : 0.0,
                    isfinite(upper[j]) ? abs(upper[j]) : 0.0)
        violation = max(violation, (lower[j]-x[j])/scale, (x[j]-upper[j])/scale)
    end
    for (i, (indices, coefficients, sense)) in enumerate(l.rows)
        isfinite(rhs[i]) || continue
        activity = sum(c*x[j] for (j,c) in zip(indices, coefficients))
        scale = max(1.0, abs(rhs[i]), sum(abs(c*x[j]) for (j,c) in zip(indices, coefficients)))
        delta = activity-rhs[i]
        violation = max(violation, (sense === :eq ? abs(delta) : sense === :le ? delta : -delta)/scale)
    end
    recomputed = x[l.vB] + sum(costs[i]*x[l.objective_columns[i]] for i in eachindex(costs))
    isapprox(value, recomputed; atol=OBJECTIVE_ATOL, rtol=OBJECTIVE_RTOL) ||
        error("objective does not match primal values: $value versus $recomputed")
    violation <= FEASIBILITY_TOL || error("primal feasibility residual $violation exceeds $FEASIBILITY_TOL")
    violation
end

function validate_models(models, points; paths=BASELINE_PATHS)
    worst_residual, worst_objective = 0.0, 0.0
    for point in points
        reference = RM.objective_value(models.rm_highs, point)
        for path in paths
            value = solve_point!(models, path, point; read_dual=false)
            isapprox(reference, value; atol=OBJECTIVE_ATOL, rtol=OBJECTIVE_RTOL) ||
                error("objective mismatch on $path at $point: $value versus $reference")
            worst_objective = max(worst_objective, abs(value-reference))
            x = primal_vector(models, path, point)
            worst_residual = max(worst_residual, residuals(models.data, models.l, point, x, value))
        end
    end
    (; points=length(points), worst_objective, worst_scaled_primal_residual=worst_residual)
end

function samples!(models, path, points; passes)
    checksum = run_path!(models, path, points; read_dual=false)
    samples = map(1:passes) do _
        GC.gc()
        m = @timed run_path!(models, path, points; read_dual=false)
        isapprox(m.value, checksum; atol=OBJECTIVE_ATOL, rtol=OBJECTIVE_RTOL) || error("unstable timing checksum for $path")
        Dict("seconds"=>m.time, "julia_bytes"=>m.bytes, "gc_seconds"=>m.gctime, "checksum"=>m.value)
    end
    Dict("path"=>String(path), "median_us_per_solve"=>1e6*median(s["seconds"] for s in samples)/length(points),
         "median_julia_bytes_per_solve"=>median(s["julia_bytes"] for s in samples)/length(points),
         "samples"=>samples)
end

function metadata()
    root = dirname(@__DIR__)
    versions = Dict(info.name=>string(info.version) for info in values(Pkg.dependencies()) if info.version !== nothing)
    Dict("utc"=>string(Dates.now(Dates.UTC)), "julia"=>string(VERSION),
         "cpu"=>Sys.cpu_info()[1].model, "architecture"=>string(Sys.ARCH), "os"=>string(Sys.KERNEL),
         "julia_threads"=>Threads.nthreads(), "timing_workers"=>1,
         "revision"=>strip(git_output("rev-parse", "HEAD")),
         "dirty"=>!isempty(git_output("status", "--porcelain")),
         "project_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__, "Project.toml")))),
         "manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml")))),
         "versions"=>versions, "highs_options"=>Dict(SOLVER_OPTIONS),
         "clp_options"=>"ReactiveModel.ClpBackend defaults; see recorded package source",
         "objective_atol"=>OBJECTIVE_ATOL, "objective_rtol"=>OBJECTIVE_RTOL,
         "scaled_primal_feasibility_tolerance"=>FEASIBILITY_TOL,
         "dual_policy"=>"No equality assertion between different valid dual optima; primal feasibility and objective checked",
         "allocation_scope"=>"Julia allocations; native solver allocations excluded",
         "timing_scope"=>"serial update and solve; construction, validation, JIT and first sweep excluded")
end

function campaign(output_dir; passes=5)
    passes > 0 || throw(ArgumentError("passes must be positive"))
    ispath(output_dir) && !isempty(readdir(output_dir)) && throw(ArgumentError("output directory must be new or empty"))
    mkpath(output_dir)
    report = Dict{String,Any}("metadata"=>metadata(), "runs"=>Any[])
    for spec in (SV.VALLEY_SMALL, SV.VALLEY_MEDIUM, SV.VALLEY_LARGE)
        data = SV.valley_data(spec; nregions=1)
        points = SV.valley_points(data; nscenarios=2, order=:refinement, stocks_per_scenario=3)
        models = build_models(data, points)
        validation = validate_models(models, points)
        shape = SV.describe(data)
        measurements = [samples!(models, path, points; passes) for path in BASELINE_PATHS]
        push!(report["runs"], Dict("preset"=>spec.name, "ncols"=>shape.ncols, "nrows"=>shape.nrows,
            "nnz"=>shape.nnz, "seed"=>spec.seed, "nregions"=>1, "nscenarios"=>2,
            "stocks_per_scenario"=>3, "order"=>"refinement", "passes"=>passes,
            "validation"=>Dict(string(k)=>v for (k,v) in pairs(validation)), "measurements"=>measurements))
        println(spec.name, ": ", validation)
        for m in measurements
            println("  ", rpad(m["path"], 12), round(m["median_us_per_solve"]; digits=2), " us/solve; ",
                    round(m["median_julia_bytes_per_solve"]; digits=1), " Julia B/solve")
        end
    end
    open(joinpath(output_dir,"results.toml"), "w") do io
        TOML.print(io, report; sorted=true)
    end
    for file in ("Project.toml", "Manifest.toml")
        cp(joinpath(@__DIR__, file), joinpath(output_dir,file))
    end
    # Keep exact source even for an uncommitted checkout or a source archive.
    source = joinpath(output_dir, "source")
    mkpath(source)
    for folder in ("src", "ext", "examples")
        cp(joinpath(dirname(@__DIR__),folder), joinpath(source,folder))
    end
    cp(joinpath(dirname(@__DIR__),"Project.toml"), joinpath(source,"Project.toml"))
    bench_source = joinpath(source,"benchmarks")
    mkpath(bench_source)
    for file in readdir(@__DIR__)
        (endswith(file,".jl") || file in ("Project.toml","Manifest.toml")) || continue
        cp(joinpath(@__DIR__,file),joinpath(bench_source,file))
    end
    println("Results: ", abspath(output_dir))
    report
end

function main(args=ARGS)
    output = isempty(args) ? joinpath(@__DIR__,"results","reproduce-"*Dates.format(Dates.now(Dates.UTC), "yyyymmdd-HHMMSS")) : args[1]
    passes = length(args) >= 2 ? parse(Int,args[2]) : 5
    campaign(output; passes)
end
end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ReproducibleBenchmarks.main()
end
