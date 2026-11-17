# Compiling applications with trimming

Trimming produces an application containing the methods reachable from its
entry point. The workflow here validates a supplied numerical application;
it is separate from installing or using the package in a Julia session.
Complete [installation](../installation.md) first.

ReactiveModel retains concrete model families, callbacks and update plans inside
a fixed execution boundary. Their types no longer propagate into `LinearProblem`
through its study definition or updater parameters. The generated adapter keeps
their methods reachable for the supplied application. Study arrays can still have
sizes determined at runtime: their lengths and
numerical values are not encoded in the model's type. A model with one generator
and a model with five generators can therefore use the same compiled family.
Specialization is per family and callback. Adding many different families can
increase compilation cost; adding keys to the same family changes its runtime
dimensions rather than creating a type per variable.

The supported contract remains unchanged: study data stays read-only between
rebuilds, and successive parameter points are `isbits` values of the same type.

## What is validated

The application in `test/trim/app.jl` is compiled with JuliaC's `--trim=safe` and
then executed with different problem sizes supplied on the command line. Its
workload checks:

- the DSL and functional constructors;
- changing bounds, matrix coefficients, right-hand sides, and objective costs;
- explicit field lists and two-argument reactivity predicates, alongside inference;
- transitions between signed zeros in inferred and explicit dependencies;
- integer, symbol, tuple, and scalar family keys, and typed empty families;
- primal values, duals, activities, and lazy result collections;
- rebuilding after study-data mutation and creating an independent problem;
- an infeasible point followed by a feasible point;
- a MILP containing both integer and binary variables;
- transitive parameter aliases and affine constants on both constraint sides
  and in objectives;
- filtered dependent families, nested sums and typed empty families;
- distinct custom keys with identical hashes and models with the same name
  declared in different modules.
- safe-mode correction of deliberately incomplete conditions, including bounds,
  matrix, RHS and objective changes, same-point checking and report reset.
- traversal across a plan group boundary, with nine distinct variable and
  constraint families and runtime-selected key counts;
- callback exceptions and recovery, GC rooting, independent clones, and replacing
  study data through the execution adapter;
- five distinct result key types and retained family snapshots.

This validation uses **HiGHSBackend**, Julia **1.13.2**, JuliaC **0.3.9**, on
**macOS aarch64**. It establishes compatibility for these application paths and
versions. Other backends, platforms, Julia versions, and application callbacks
need their own compilation and execution checks. The unit tests also exercise
the continuous workload with Clp in normal Julia execution.
Read [design and compilation](design.md) and the
[client compilation diagnostic](../compilation-results.md) for specialization costs.

## Reproduce the compilation

From the repository root, use Julia 1.13.2 and JuliaC 0.3.9 for the recorded
configuration. Follow [JuliaC installation](https://github.com/JuliaLang/JuliaC.jl)
to make its `juliac` command available, then run:

```sh
rm_bundle=$(mktemp -d)
juliac --project . --output-exe rm-typed-checks \
  --bundle "$rm_bundle" --trim=safe --experimental \
  test/trim/app.jl

"$rm_bundle/bin/rm-typed-checks" 1
"$rm_bundle/bin/rm-typed-checks" 5
"$rm_bundle/bin/rm-typed-checks" 9
# Invalid size: expect exit code 1.
"$rm_bundle/bin/rm-typed-checks" 0
```

The positive-size runs return exit code zero only when the numerical checks
pass. The same binary handles all sizes. Size zero returns exit code 1. Compilation alone is not sufficient: always execute
the produced application and check its results.
Use a fresh bundle directory for each build: JuliaC 0.3.9's artifact bundler
does not overwrite artifacts already present in the output directory.

Keep the bundle together when moving the application: the executable relies on
the Julia runtime and native solver libraries shipped alongside it. Trimming
does not imply an application without allocations or native dependencies.
See [JuliaC](https://github.com/JuliaLang/JuliaC.jl) for installation and bundling.

## Writing an inferable application

Define the model and callbacks in normal Julia source before compiling the
application. Keep the model family types known on each reachable construction
path. Dynamically sized `Vector` and `Matrix` study data are allowed; fields such
as `Any`, abstract callback fields, runtime `eval`, and loading new model code
after compilation can introduce unresolved calls.

The public constructors prepare a concrete update plan and its execution adapter,
while leaving numerical assembly and solver initialization until the first solve. `rebuild!` prepares a
new numerical generation with the same family types. Changing the schema of
families can require a different compiled construction path.

Scalar `primal`, `dual`, and `activity` queries read current numerical buffers
without materializing collections. Accessing `solution.variables` or
`solution.constraints` builds a typed snapshot collection on demand. The
`Solution` object itself remains a live view, reused by later solves.

## Integer and binary variables

Domains are part of the fixed definition and are passed through the DSL:

```@example trim_domains
using ReactiveModel

@model CompilationDomainExample(data, point) begin
    @variable units[i in 1:data.n] in [0.0, data.limit] domain=:integer
    @variable enabled in [0.0, 1.0] domain=:binary
    @constraint floor[i in 1:data.n] units[i] >= point.minimum
    @constraint enabled_floor enabled >= 0.5
    @objective Min sum(units[i] for i in 1:data.n) + enabled
end

data = (n=2, limit=10.0)
problem = LinearProblem(ConstData(CompilationDomainExample(data), data),
                        (minimum=1.5,); backend=HiGHSBackend())
solution = solve(problem, (minimum=1.5,))
@assert isapprox(objective_value(solution), 5.0)
@assert isapprox(primal(solution, :units, 1), 2.0)
@assert isapprox(primal(solution, :enabled), 1.0)
@assert isnan(dual(solution, :enabled_floor))
(objective_value(solution), primal(solution, :units, 1), primal(solution, :enabled))
```

The objective is `5.0`: two integer units at `2.0` each and the binary variable
at `1.0`. With continuous domains, that formulation would have a fractional
solution. The domain must stay fixed across solves; parameter-dependent bounds
can still change. Clp rejects integer and binary domains explicitly.

MILP solutions do not have LP row duals; HiGHS stores `NaN` in those result
buffers. Claims about reusing a simplex basis or interpreting LP dual values do
not extend to the full integer solve.

## Diagnose inference before compiling

JET's optimization analysis helps locate unresolved calls in construction and
update paths. The optional script `test/jet_checks.jl` runs these analyses on the
same numerical workloads, including explicit checks of concrete execution
adapters whose payload types are hidden from client code. Run it from an
environment containing a compatible
JET version and ReactiveModel; JET is not a runtime dependency of ReactiveModel.

The snippet below assumes that separate environment is already active and the
working directory is the repository root:

```julia
include("test/jet_checks.jl")
```

A clean inference report is useful evidence, but the safe compiler and the
executed application remain the acceptance test. See
[JET's optimization analysis](https://aviatesk.github.io/JET.jl/stable/optanalysis/).
