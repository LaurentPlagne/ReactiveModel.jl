# Synthetic benchmarks and diagnostics

The presets are constructed examples with seeds 11, 22 and 33. Their graph
construction, dimensions and coefficient formulas are specified in
[PROVENANCE.md](PROVENANCE.md). All study arrays are generated in code.

## Repeated solves

The [serial/parallel protocol](../docs/src/parallel-results.md) gives the commands,
numerical checks, timing boundaries and output format. It compares RM's direct
Clp_jll/HiGHS_jll adapters with explicit JuMP/Clp and JuMP/HiGHS implementations,
using persistent models and ThreadTeams dispatch. Review the
[comparison code](../docs/src/benchmarks-jump.md) when choosing a baseline.

Start with the consistency checks from the repository root:

```sh
julia --startup-file=no --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=benchmarks benchmarks/test_reproduce.jl
julia --startup-file=no --project=benchmarks benchmarks/test_frontends.jl small medium large
julia --startup-file=no --project=benchmarks -t 2 benchmarks/test_parallel_comparison.jl
```

`reproduce.jl OUTPUT 5` provides a small serial comparison against JuMP, MOI and
POI on six points per preset. `diagnose_parallel_clp.jl` and
`diagnose_clp_cache.jl` inspect scheduler placement and adapter phases;
`show_jump_case.jl` displays a concrete formulation and update path.

For each timing run, keep raw output with the checkout, dependency manifest,
backend options and machine configuration. Run experiments sequentially on an
otherwise available machine. Compare validated numerical outputs before timing,
and distinguish point-loop allocations from complete team dispatch and native
solver memory. Synthetic fixed-cut LPs do not describe complete SDDP training
or establish a general solver/frontend ranking.

## First use and compilation

The [compilation protocol](../docs/src/compilation-results.md) covers
`compilation_scaling.jl` and `reduction_compilation.jl`. These generate algebraic
models and numerical helper inputs, separate cold calls from warmed probes,
and retain their numerical assertions. Run each first-use configuration in a
fresh process.

Benchmark dependencies belong to this environment, whose relative source is
ReactiveModel (`..`). Optional ClpJL and TinyHiGHS are not required. Keep the
manifest used for a run; resolving another environment can change the experiment.
A populated depot does not establish empty-depot installation.
