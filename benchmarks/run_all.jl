"""
    benchmarks/run_all.jl

Master benchmark runner for ReactiveModel.jl.
Runs all comparisons in clean isolated processes and outputs formatted tables reproducible by any user.

Usage:
    julia --startup-file=no --project=benchmarks -t 8 benchmarks/run_all.jl [nworkers]
"""

using Printf
using Dates

println("="^70)
println("ReactiveModel.jl Benchmark Suite")
println("Date: ", Dates.today())
println("Julia: ", VERSION)
println("CPU: ", Sys.CPU_NAME)
println("Threads: ", Threads.nthreads())
println("="^70)

nworkers = isempty(ARGS) ? min(Threads.nthreads(), 8) : parse(Int, ARGS[1])
cmd_base = Base.julia_cmd()
bench_dir = @__DIR__

println("\n>>> Step 1/3: Checking correctness and parity across all front-ends (test_frontends.jl)...")
run(`$cmd_base --startup-file=no --project=$bench_dir $(joinpath(bench_dir, "test_frontends.jl"))`)

println("\n>>> Step 2/3: Running economic dispatch comparison (Tutorial 1 LP, dispatch_comparison.jl)...")
run(`$cmd_base --startup-file=no --project=$bench_dir $(joinpath(bench_dir, "dispatch_comparison.jl"))`)

println("\n>>> Step 3/3: Running synthetic cascade benchmarks (Docs measurements, docs_measurements.jl)...")
run(`$cmd_base --startup-file=no --project=$bench_dir -t $nworkers $(joinpath(bench_dir, "docs_measurements.jl")) $nworkers`)

println("\n" * "="^70)
println("All benchmarks completed successfully.")
println("="^70)
