# What does the JuMP benchmark solve?

The flat JuMP/HiGHS case in [the synthetic benchmark workflow](parallel-results.md) is a
persistent, **flat-buffer** formulation of a synthetic reservoir LP. It is built
once, uses preallocated numeric buffers, and checks parameter changes and
individual numerical values before calling JuMP setters. It is not a rebuild
at every point and it is not a claim about a typical JuMP application.

The [ThreadTeams workflow](parallel-results.md) compares that same approach
and two named-variable versions with both HiGHS and Clp, plus a default-option
Clp control and an explicit inner-cache-reset variant. This page shows their source.
JuMP uses HiGHS.jl/MOI or Clp.jl/MOI; RM calls Clp_jll or HiGHS_jll directly.
Clp.jl is a dependency only of the benchmark project. Both Clp paths use the
same native Clp_jll binary in this workflow.

The [adapter diagnostics](benchmark-diagnostics.md) explain the Clp paths.
The driver times ordinary updates and explicit cache reset separately;
native reload is included in both cases.

## Mathematical problem

The presets operate constructed cascades over **8, 12 or 16 periods**. Variables are
initial and final-period reservoir stocks, unit discharges, group discharges
for generating units, and a scalar future value. Constraints enforce water
balance, discharge aggregation, final-stock conditions and piecewise-linear
future-value cuts. The maximization objective is generation revenue plus
future value.

Between points, a scenario selects inflows and prices, while initial stock
changes the fixed initial-volume bounds and, through the supplied functions,
some revenue coefficients. The tested region stays fixed. Data are generated
from a recorded seed; no external study inputs are needed.

| Preset | Reservoirs | Units | Variables | Explicit affine rows |
| :--- | ---: | ---: | ---: | ---: |
| Small | 1 | 2 | 50 | 44 |
| Medium | 3 | 5 | 172 | 102 |
| Large | 5 | 9 | 390 | 236 |

The counts of columns and explicit affine rows are checked for each
RM and JuMP path. Variable-bound sets are not counted as affine rows.

## Run the JuMP cases yourself

From the repository root, after instantiating `--project=benchmarks`:

```sh
julia --startup-file=no --project=benchmarks benchmarks/show_jump_case.jl
```

With no argument, this prints objectives for both solvers. An optional `clp` or
`highs` argument selects one. Flat-buffer and named/block-guarded models are
checked at six points against the corresponding direct RM backend. The full measurement command is
on the [serial/parallel results page](parallel-results.md).

The source entry points are:

| Case name | Construction | Update function |
| :--- | :--- | :--- |
| Flat buffers | `jump_reference` | `jump_update!` |
| Clp, explicit reset | `jump_reference` | Inner-cache reset, then `jump_update!` |
| Named, guarded | `jump_user_model` | `jump_user_update_predicated!` |
| Named, unconditional | `jump_user_model` | `jump_user_update!` |

The first pair lives in `benchmarks/valley_comparison.jl`; the named versions
live in `benchmarks/jump_user_model.jl`. The code excerpts below are read from
those files during the documentation build. They need the benchmark's data,
layout and helper functions; use the command above to run the complete example.

## Flat-buffer construction

`x[1:ncols]` stores variable references; `rowrefs` stores constraint references.
Affine expressions are assembled once from the shared layout. These construction
allocations are outside the warmed update-and-solve timing.

```@eval
using Markdown, ReactiveModel
source = read(joinpath(pkgdir(ReactiveModel), "benchmarks", "valley_comparison.jl"), String)
fragment = split(split(source, "function jump_reference("; limit=2)[2], "\nfunction jump_update!"; limit=2)[1]
Markdown.parse("```julia\nfunction jump_reference(" * fragment * "\n```")
```

## Flat-buffer updates

The application retains copies of current bounds, RHS values and costs. It
computes affected components into scratch buffers and suppresses setter calls
when individual values compare equal. A new point is optimized using the same
JuMP model, followed by a termination-status and objective query.

```@eval
using Markdown, ReactiveModel
source = read(joinpath(pkgdir(ReactiveModel), "benchmarks", "valley_comparison.jl"), String)
fragment = split(split(source, "function jump_update!("; limit=2)[2], "\nfunction moi_reference"; limit=2)[1]
Markdown.parse("```julia\nfunction jump_update!(" * fragment * "\n```")
```

## Named variables and equations

The alternative model keeps reservoir, unit and time indices in the JuMP
containers. Here is its declaration and water-balance excerpt; the complete
function also adds group sums, final-stock conditions, cuts and the objective.

```@eval
using Markdown, ReactiveModel
source = read(joinpath(pkgdir(ReactiveModel), "benchmarks", "jump_user_model.jl"), String)
fragment = split(split(source, "function jump_user_model("; limit=2)[2], "\n    outside_cuts"; limit=2)[1]
fragment = join(filter(line -> !startswith(strip(line), "#"), split(fragment, '\n')), '\n')
Markdown.parse("```julia\nfunction jump_user_model(" * fragment * "\n    # ... remaining constraints and objective in the full source ...\nend\n```")
```

The named/guarded update checks which parameter blocks changed before calling
setters through references such as `state.V0[r]` and
`state.water_balance[r, t]`. It guards blocks, while the flat version also
compares individual numbers. The unconditional version reapplies the
parameter-dependent quantities even when a block's dependencies stayed equal.
These are distinct update policies, rather than interchangeable baselines.

```@eval
using Markdown, ReactiveModel
source = read(joinpath(pkgdir(ReactiveModel), "benchmarks", "jump_user_model.jl"), String)
fragment = split(split(source, "function jump_user_update_predicated!("; limit=2)[2], "\nfunction jump_user_solve!"; limit=2)[1]
Markdown.parse("```julia\nfunction jump_user_update_predicated!(" * fragment * "\n```")
```

## Clp: the wrapper also matters

The same construction functions accept `optimizer=Clp.Optimizer` or a configured
optimizer factory. The dual/no-presolve case uses
`JuMP.optimizer_with_attributes(Clp.Optimizer, "SolveType"=>0, "PresolveType"=>1,
"Perturbation"=>50, "PrimalTolerance"=>1e-7, "DualTolerance"=>1e-7)`.
A separate flat-buffer case uses Clp's default solver options.

Clp.jl 1.4 reports no incremental MOI interface. Its inner cache is
attached after a solve, empty after a numerical update, and attached after the
next optimization. Its reattachment loads the complete native problem through
`Clp_loadProblem`. Thus, keeping the algebraic JuMP model does not guarantee
native model/basis reuse for this adapter. RM/Clp sends numerical C API updates
and reoptimizes with `Clp_dual`. The [benchmark workflow](parallel-results.md)
keeps this adapter behavior separate from algebraic syntax when comparing cases.

The additional flat Clp variant calls
`MOI.Utilities.reset_optimizer(JuMP.backend(model).optimizer.model)` before
`jump_update!`. Its reset is part of the timed path. Access to that inner cache
depends on the Clp.jl/MOI layout and still leads to a full native reload.
The implementation is `ResetClpPath` in `benchmarks/parallel_comparison.jl`.

## Interpretation

The paths deliberately expose formulation and update-policy choices. Measurements
on one do not characterize every implementation using JuMP. RM/HiGHS versus
JuMP/HiGHS holds the native solver fixed; RM/Clp versus JuMP/HiGHS additionally
changes the solver. Construction, compilation, result requirements and any
rebuilds should be included when evaluating a complete application.
