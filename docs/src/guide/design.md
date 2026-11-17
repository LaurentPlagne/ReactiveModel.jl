# Design and compilation costs

ReactiveModel targets repeated numerical changes within a fixed LP family.
Initial construction and later updates use the same formulas. A driver supplies
the next point without reproducing the formulation as a sequence of setters.
This page explains the implementation choices and their costs; the
[usage contract](concepts.md) remains the reference for application code.

## Layers and source layout

| Layer | Responsibility | Source |
| :--- | :--- | :--- |
| Definition | Families, keys, domains and numerical formulas | `src/core/definitions.jl` |
| DSL | Validate syntax, normalize affine expressions, emit typed callbacks | `src/DSL/ast_utils.jl`, `src/DSL/codegen.jl` |
| Execution | Fixed client boundary, typed callbacks and shared numeric kernels | `src/core/execution.jl`, `src/core/callbacks.jl`, `src/core/updates.jl` |
| Workspace | Generations, initial assembly and independent clones | `src/core/problem.jl`, `src/core/assembly.jl` |
| Results | Live views, scalar queries, solve/recovery policy | `src/core/results.jl`, `src/core/solve.jl` |
| Native interface | Common arrays, solver-specific operations and result accessors | `src/Backends/` |
| Safe checks | Missed changes and unchanged reevaluation counters | `src/UpdateSafety.jl` |

`src/ReactiveModel.jl` declares the module/API and includes these layers in type
dependency order. Files share one module. Optional solvers remain extensions
under `ext/`; file separation does not add runtime dispatch.

## What enters a type

The mathematical definition uses tuples of concrete blocks and a named view of
the same objects. Adding a family or changing its callback type produces a new
definition type.

`LinearProblem` retains the original `ConstData` object through a field typed by
the study-data type, without exposing its mathematical schema in the problem's
type parameters. A fixed execution handle roots the concrete definition, update
plan and workspace. Its adapter is generated with a concrete `@cfunction`
signature for that payload. Numerical assembly, updates, lazy result layout,
data replacement and same-backend cloning pass through this boundary. Callback
methods remain statically reachable for the tested trimming paths.

Different formulas can share a `LinearProblem` type when study-data, parameter,
backend and distinct result key types match. Result families are grouped by
key type; their names and counts are runtime metadata. Adding another family
with an existing key type does not enlarge the public problem, workspace or
solution type. A new distinct key type can still change those types. This
reduces actual public type parameters, separately from abbreviated printing.

The handle owns a GC root for its payload and a reusable argument slot. Calls
preserve the handle, and the slot is cleared on success and exceptions. Each
problem and clone owns its own mutable execution state. The adapter is called
directly from Julia, with no external library executing between caller and
callee; callback exceptions keep their type and backtrace in the checked paths.

Keys, numerical values and study-vector lengths remain runtime values. A
thousand columns sharing one rule should normally be one indexed family,
rather than a thousand distinct scalar declarations. The former reuses one
callback specialization; the latter requires more compilation work.

```julia
@model IndexedPlant(data, point) begin
    @variable output[g in 1:data.n] in [0.0, available_capacity(data, point)]
    @constraint demand sum(output[g] for g in 1:data.n) >= required_output(data, point)
    @objective Min unit_cost(data, point) * sum(output[g] for g in 1:data.n)
end
```

This schematic fragment assumes the three business functions. `data.n` changes
loop lengths, not the number of family types. The executed
[business-function example](../tutorials/business-functions.md) supplies complete
functions of `(study_data, study_parameters)`.

## Bounded specialization

Concrete tuples are useful here, but unbounded inlining is expensive. The
execution plan stores groups of at most eight families, constructed once, with
call boundaries between groups. Shared numerical kernels are separate from
callbacks, limiting the machinery compiled into every family's code. Eight is
an internal engineering choice, not an optimum established for all Julia versions.

The visitor explicitly specializes its forwarded callable as `f::F where F`.
Julia can avoid specializing a function argument that is only forwarded; see
the official [specialization guidance](https://docs.julialang.org/en/v1/manual/performance-tips/#Be-aware-of-when-Julia-avoids-specializing).
Retaining that specialization avoids dynamic calls and boxing at these boundaries.

Mathematical tuples remain typed inside the execution payload. Compilation
still grows with distinct rules; hiding their types from client signatures does
not eliminate callback compilation.
The [compilation protocol](../compilation-results.md) describes how to investigate this tradeoff.

A coefficient emitter needs column positions, not the complete model type.
Its dictionaries are grouped by key type, with `(family_name, key)` as the
lookup key. Another integer-keyed family does not change this index/emitter
type. Key equality, rather than a hash alone, resolves collisions.
Updaters and assembly visitors retain study payloads and compact numerical
storage, excluding unrelated callback/schema types from their kernels.

The small generated index factory inspects schema field types to choose the
distinct key-type maps before traversal. It never evaluates study values or
callbacks. Starting with a fixed group type keeps mixed-key construction
inferable; growing a tuple of new key types recursively can lose that inference.

## Readable client diagnostics

The full parametric type can still be long. RM abbreviates its own aggregate
types in REPL output, contextual nested signatures, stacktraces and `MethodError`
messages. A problem type shows its parameter type name and backend, for example
`LinearProblem{P=NamedTuple, B=ClpBackend, …}`. The model name is retained as a
value and appears when displaying the problem instance. Models with different
names can share a type. Internal updater labels retain the family name when the
DSL provides it.

These abbreviations are a display policy, separate from the execution boundary.
Different hidden parameters can still distinguish types. For full investigation:

```julia
ReactiveModel.full_type_string(typeof(problem))
ReactiveModel.full_type_string(typeof(problem.const_data.definition))
showerror(IOContext(stderr, :reactivemodel_full_types => true), error)
```

These schematic inspection calls assume a problem and caught exception. The
IOContext flag also works on nested types. Compact labels are descriptions,
not Julia type expressions to paste into declarations. The compilation
diagnostic script records full spelling size with abbreviation disabled.

Formatting is scoped to `IOContext` streams. Bare `IO` streams retain Base's
full type printing; use `show(IOContext(io, :compact => true), typeof(problem))`
for compact stream output. This boundary preserves statically resolvable error
paths in the tested safe-trimmed application.

## Numerical storage and ownership

| Representation | Meaning and owner |
| :--- | :--- |
| Definition | Complete mathematical rules and fixed key spaces |
| Workspace arrays and row dictionaries | Complete current bounds, RHS, costs and rows for comparisons and activities |
| Setter packets and scratch | Changed entries for one update; shared scratch is overwritten by the next component |
| Backend/native model | Solver ABI arrays, imported matrix, solution and solver-specific state |

`_AssemblyBuffers` holds references to existing arrays, not another model copy.
Sparse row and objective comparisons share a kernel. Removed nonzero entries
produce explicit zero setters. Retained row packets own their delta arrays so
the next row cannot overwrite them through shared scratch.

A retained result collection owns its numerical snapshot. A retained family
view keeps only that family's block alive. Names, typed keys and global offsets
come from a lazy result layout, independent of callback types. Named access
(`solution.variables.output[key]`) and declaration-order iteration use this
metadata.

Native CSC arrays may become stale after direct coefficient updates. Clp
refreshes its native matrix before a reload; importing the original arrays
would undo changes that RM already considers applied. Distinct storage
lifetimes are therefore intentional and must not be mutated independently.

## Reductions and business functions

`sum` can incur compilation work for its generator, callable and input types.
A new closure expression or callback type can require another specialization.
Increasing a `Vector{Float64}`'s length normally preserves those types.

Supported DSL sums involving decision variables are lowered to emission loops
over runtime keys. Many distinct family terms still create a larger syntax
tree. Replacing an indexed family with many scalar declarations preserves that cost.

Keep numerical reductions in named helpers with explicit inputs and stable
types. A loop accumulator can reduce first-use cost for a particular helper;
replacing every `sum` is not an established optimization. The included
[sum/loop diagnostic protocol](../compilation-results.md) separates cold calls
from warmed allocation probes. Shared native assembly uses reductions over fixed
concrete buffer types, reusing their specialization across client schemas.

## Safe and selective evaluation

With `safe_mode=true`, every formula is evaluated on every solve, including an
identical point. Conditions are observed but cannot suppress evaluation. Only
numerical differences reach the solver; an unchanged LP can reuse its outcome
after checking. With `safe_mode=false`, the point cache and component conditions
can skip evaluation. Both paths share formulas and numerical comparison kernels.

The safe report counts complete component comparisons, not scalar setters.
A dictionary maps `(family, component)` to its counter slot, avoiding repeated
scans over all counters. Equal values after a positive condition identify an
observed reevaluation cost, not proof of global invariance. Conditions are not
rewritten automatically.

Safe mode does not validate business mathematics. A custom constraint must
explicitly supply its RHS, including an intended zero. Study data stays
read-only between rebuilds; helpers must be deterministic from explicit inputs.

## Backend interfaces and cloning

The six adapters share mathematical assembly: bounds/cost arrays and a Julia
CSC matrix. Pure Julia adapters use its one-based indices; native adapters
convert them to zero-based buffers. Each adapter keeps its own coefficient
import policy, solve, options and recovery behavior. Sharing assembly does not
equate native state histories. Numerical values are `Float64`; native integer
buffers follow the solver ABI. See the [policy table](../backends.md#Shared-assembly-and-solver-policies).

Result reads use `backend_status`, `backend_objective_value`,
`backend_primal_value` and `backend_dual_value`. Default methods use conventional
fields; a backend with other storage overrides them. This protocol does not
depend on MOI. See [backend implementation](../backends.md).

`similar(problem)` creates independent native state. HiGHS options configured
through RM are recorded as concrete string pairs and copied. The model, basis,
solution and safe report are not shared. Native option calls outside RM are not
tracked. Other configurable backends must supply a clone factory; the generic
fallback uses a zero-argument constructor.

## Failures and result lifetime

A nonoptimal solve can follow a complete numerical update. RM retains that
model, invalidates prior results and raises `SolveError`. A callback/backend
exception can leave a partial update: results are invalidated and the next call
rebuilds. These paths have different recovery costs.

`Solution` is a reusable live view guarded by a generation token. Scalar queries
use the result protocol. Collections are materialized lazily and can be retained
as snapshots. Its result layout retains family names, typed keys and positions,
without retaining mathematical formulas or model identity. The layout is cached
on the first result or structural query. Models with the same result layout,
backend, parameter and key types can share `Solution` and `Workspace` types even
when their formulas differ. The typed definition and update plan remain in the
execution payload, outside client problem and result type parameters.
Applications copy the values they want to keep. See
[the solve cycle](solve-cycle.md).

For representation changes, compare incremental values with fresh models,
including removed coefficients, objective constants and failure/recovery.
Check warmed allocations and first-use latency separately. Execute the trimmed
application across a group boundary and at several runtime key counts; an
inference report alone is insufficient. The supplied workload includes nine
distinct DSL families, safe mode and several key types. Compilation and trimming
evidence remains specific to the tested application, versions and platform.
