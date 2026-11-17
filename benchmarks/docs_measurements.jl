"""
The measurements quoted in the manual (`docs/src/why/performance.md` and
`docs/src/backends.md`).

    julia --startup-file=no --project=benchmarks -t 8 \\
        benchmarks/docs_measurements.jl [nworkers]

For each preset of `synthetic_valley.jl`, the production traversal (three
concave regions, nine stock points per scenario, all scenarios) is timed
serially and with `nworkers` threads (default 8) for:

* ReactiveModel on `HiGHSBackend` and `ClpBackend`, plus `TinyHiGHSBackend`
  when RM_BENCH_TINYHIGHS=1 in a separate environment;
* JuMP written as a user writes it, with naive updates (`jump_user`) and with
  hand-written change predicates (`jump_user_pred`);
* direct MOI and POI on HiGHS.

Before any timing, enabled ReactiveModel backends, JuMP, MOI and POI are
checked against each other on objective and on a dual at every point of a
parity sequence.
"""

include(joinpath(@__DIR__, "valley_comparison.jl"))

const DOC_PATHS = benchmark_paths((:rm_highs, :rm_tiny, :rm_clp, :jump_user, :jump_user_pred,
                   :moi, :poi))

function doc_run(spec; nscenarios=spec.nscenarios, nregions=3, nworkers=8)
    pdata = SV.valley_data(spec)
    ppoints = SV.valley_points(pdata; nscenarios=8, order=:refinement,
                               stocks_per_scenario=3)
    parity = check_parity!(build_models(pdata, ppoints), ppoints)

    data = SV.valley_data(spec; nregions)
    points = production_points(data, 1:nscenarios)
    models = build_region_models(data, points, nregions)
    schedule = production_schedule(models, data, points; nregions)
    serial = [timed_production!(models, path, schedule) for path in DOC_PATHS]

    chunks = scenario_chunks(nscenarios, nworkers)
    worker_points = [production_points(data, chunk) for chunk in chunks]
    worker_models = [build_region_models(data, worker_points[w], nregions)
                     for w in 1:nworkers]
    worker_schedules = [production_schedule(worker_models[w], data,
        worker_points[w]; nregions) for w in 1:nworkers]
    parallel = [timed_production_parallel!(worker_models, worker_schedules, path)
                for path in DOC_PATHS]
    (; shape=SV.describe(data), nsolves=solves(schedule), parity, serial, parallel)
end

function docs_main(args=ARGS)
    nworkers = isempty(args) ? min(Threads.nthreads(), 8) : parse(Int, args[1])
    for name in ("small", "medium", "large")
        spec = PRESETS[name]
        result = doc_run(spec; nworkers)
        println("### ", spec.name, "  ", result.shape, "  solves=", result.nsolves)
        println("  parity ", result.parity)
        for (s, p) in zip(result.serial, result.parallel)
            println(rpad(string(s.path), 16),
                " serial ", round(s.per_solve_us; digits=1), " us ",
                round(Int, s.per_solve_bytes), " B | ", nworkers, "w ",
                round(p.per_solve_us; digits=1), " us ",
                round(Int, p.per_solve_bytes), " B gc ",
                round(100p.gc_share; digits=1), "%")
        end
        flush(stdout)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    docs_main()
end
