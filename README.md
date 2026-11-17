# ReactiveModel.jl

[![CI](https://github.com/LaurentPlagne/ReactiveModel.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/LaurentPlagne/ReactiveModel.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-manual-blue.svg)](https://laurentplagne.github.io/ReactiveModel.jl/)

**[Read the documentation](https://laurentplagne.github.io/ReactiveModel.jl/)** ·
[Getting started](https://laurentplagne.github.io/ReactiveModel.jl/tutorials/first-model/) ·
[API reference](https://laurentplagne.github.io/ReactiveModel.jl/api/)

ReactiveModel is an experimental Julia library for repeatedly solving a family
of linear programs. You describe how bounds, matrix coefficients, right-hand
sides and objective coefficients depend on study data and a parameter point.
ReactiveModel keeps a solver model and updates its numerical values between
solves.

The intended use is a sequence of related problems with the same variable and
constraint families, for example in simulation, rolling-horizon optimization
or dynamic programming. The formulation and its parameter dependencies stay
together in one definition.

## Why ReactiveModel?

- **One definition to maintain.** Bounds, equations and their parameter-dependent
  formulas serve both construction and subsequent updates, including business
  functions of `(study_data, study_parameters)`.
- **Checked updates.** Start with `safe_mode=true`: every numerical formula is
  reevaluated and actual changes are applied even if a condition misses them.
  The report identifies missed changes and reevaluations without numerical change.
- **Reusable numerical workspaces.** Persistent model and scratch storage serve
  successive solves. Measure callbacks, result handling and solver state reuse
  on the execution path your application needs.
- **Applications compiled with trimming.** A supplied workload, including safe
  updates, is compiled and executed at runtime-selected sizes with JuliaC and
  HiGHS in the documented configuration.

## Try it from a checkout

Julia 1.12 or later is required by the package. See the backend and compilation
guides for validation scope. From the repository root:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. examples/getting_started.jl
```

The example below is also available in `examples/getting_started.jl`.
HiGHS and Clp are installed with the package; no adjacent solver checkout is
needed.

```julia
using ReactiveModel

@model Dispatch(data, point) begin
    @variable production[g in 1:2] in [0.0, data.capacity[g]]
    @constraint balance sum(production[g] for g in 1:2) == data.demand[point.t]
    @objective Min sum(data.cost[g] * production[g] for g in 1:2)
end

data = (capacity=[80.0, 60.0], cost=[40.0, 70.0], demand=[100.0, 110.0])
problem = LinearProblem(ConstData(Dispatch(data), data), (t=1,);
                        backend=HiGHSBackend(), safe_mode=true)

first_cost = objective_value(problem, (t=1,))
solution = solve(problem, (t=2,))
second_cost = objective_value(solution)
production = [primal(solution, :production, g) for g in 1:2]

@assert isapprox(first_cost, 4600.0)
@assert isapprox(second_cost, 5300.0)
@assert production ≈ [80.0, 30.0]
(first_cost, second_cost, production)
```

For realistic functions and a branching control loop, see
[business functions and safe updates](docs/src/tutorials/business-functions.md),
or run `examples/business_functions.jl`. `update_report(problem)` exposes the
checks; observed unchanged values are optimization leads, not proof that a
dependency can be removed. The constructor default is `safe_mode=false`.

## Usage contract

- Study data can contain arrays, dictionaries and mutable objects. Keep their
  contents constant between rebuilds. After mutation, call `rebuild!` on every
  problem sharing the data before solving again.
- Parameter points must be `isbits` values of the same concrete type for one
  problem, such as `(t=1,)`. Dependency inference follows their use in the DSL;
  helper functions must depend on explicit inputs.
- Changing family keys or dimensions requires a new definition and problem.
- A `Solution` is a reusable live view: a later solve updates it. Copy values
  you need to retain. A failed solve or rebuild invalidates earlier results.
- Use an independent `LinearProblem` for each concurrent task.

## Scope and related tools

The primary target is continuous LPs. HiGHS also supports integer and binary
variables through `domain=:integer` and `domain=:binary`. Backend capabilities
and validation limits are described in the [backend guide](docs/src/backends.md).
The HiGHS backend uses `HiGHS_jll` directly, without HiGHS.jl or MathOptInterface.

JuMP, MathOptInterface and ParametricOptInterface provide other ways to formulate
and update optimization models. ReactiveModel explores a narrower contract for
repeated solves. See [modeling approaches](docs/src/why/ergonomics.md) for their
roles and references.

Performance depends on the model, update pattern, solver and requested results.
The [serial/parallel protocol](docs/src/parallel-results.md) describes comparisons
with explicit JuMP/Clp and JuMP/HiGHS implementations on generated inputs. The
[JuMP comparison code](docs/src/benchmarks-jump.md) explains the baseline, and the
[input specification](docs/src/benchmark-inputs.md) defines the synthetic presets.
Use the [benchmark method](docs/src/benchmarks.md) to collect validated numerical
outputs, timings and allocations for your checkout and machine. Julia allocation
counts exclude native solver memory.

The trimming workload has been compiled with Julia 1.13.2, JuliaC 0.3.9 and
HiGHS on macOS aarch64, then executed at several runtime sizes. This evidence
is specific to the tested application. See [compilation](docs/src/guide/trimming.md).
The [design guide](docs/src/guide/design.md) explains family specialization,
execution plans, numerical storage and backend interfaces. The
[compilation diagnostic](docs/src/compilation-results.md) separates first use
from warmed execution and examines reductions in business functions.

## Documentation and verification

The [online manual](https://laurentplagne.github.io/ReactiveModel.jl/) is built and
published from `main` after the CI tests and executable documentation pass.
Its [Markdown sources](docs/src/index.md) are also available in this repository.

Start with [installation](docs/src/installation.md) and the
[first model](docs/src/tutorials/first-model.md). The manual also covers the
[functional API](docs/src/tutorials/first-model-functional.md),
[business functions and checked updates](docs/src/tutorials/business-functions.md),
[indexed families](docs/src/tutorials/indexed-families.md),
[a reservoir stage problem](docs/src/tutorials/bellman-recursion.md),
[a small SDDP training loop](docs/src/tutorials/sddp.md),
[data and dependencies](docs/src/guide/concepts.md), and
[result lifetime and errors](docs/src/guide/solve-cycle.md).

From the repository root:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```

The documentation build executes tutorial examples and writes `docs/build/`.
Open `docs/build/index.html` locally. See [contribution instructions](CONTRIBUTING.md).

## Author and license

Laurent Plagne — <laurent.plagne@gmail.com>. Released under the [MIT license](LICENSE).
