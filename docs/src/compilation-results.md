# Compilation diagnostic protocol

The supplied diagnostics separate model definition, construction, first use
and warmed execution. They use generated algebraic LPs and ordinary numerical
helpers. The [design guide](guide/design.md) explains the typed schemas,
execution boundaries and runtime key spaces being examined.

## Family count and key count

`benchmarks/compilation_scaling.jl` creates N variable families and N constraint
families with K keys per family: N×K columns and N×K rows. Its objective contains
an indexed term for each variable family. Changing N changes the generated
client syntax and distinct rules; changing K changes runtime dimensions within
the same family structure.

Run each configuration in a fresh Julia process. The script records package
loading, DSL definition, construction, first solve, first update, first result
and collection materialization separately. A warmed sweep changes RHS values
and checks objectives and Julia allocations. Dependency caches may already be
populated; state their condition when reporting a run.

From the repository root, after instantiating the package environment:

```sh
rm_compile=$(mktemp -d)
julia --startup-file=no --project=. benchmarks/compilation_scaling.jl "$rm_compile/families-4.toml" 4 1
julia --startup-file=no --project=. benchmarks/compilation_scaling.jl "$rm_compile/families-32.toml" 32 1
julia --startup-file=no --project=. benchmarks/compilation_scaling.jl "$rm_compile/keys-128.toml" 4 128
```

Repeat configurations into distinct output files. For revision comparisons,
use separate frozen checkouts with matching Julia and dependency versions,
identical diagnostic scripts and alternating execution order. Run processes
sequentially and record competing machine load. A single sample cannot
characterize compilation variability.

## Type spelling and diagnostics

The script records complete problem, workspace, solution and client type
spellings after timing, with compact formatting disabled. Text length describes
a diagnostic representation, not compiler memory or type identity. Compact
printing does not change the compiler's types.

Use `ReactiveModel.full_type_string(problem)` or the
`:reactivemodel_full_types => true` IOContext flag for complete inspection.
The [design guide](guide/design.md#Readable-client-diagnostics) explains the
display policy. Distinct formulas still need compilation even when their
public problem and result types match.

## Numerical reductions

`benchmarks/reduction_compilation.jl` compares a `sum` generator and an explicit
loop computing the same thresholded weighted sum. Inputs are a `Vector{Float64}`
and a fixed named-tuple point type. Each helper is first invoked through
`invokelatest` so its compilation occurs inside the first-call timer.

```sh
rm_reduction=$(mktemp -d)
julia --startup-file=no --project=. benchmarks/reduction_compilation.jl "$rm_reduction/sum.toml" sum 4096
julia --startup-file=no --project=. benchmarks/reduction_compilation.jl "$rm_reduction/loop.toml" loop 4096
```

The scripts check numerical agreement and record first-call and warmed costs
separately. Vary vector length and repeat fresh processes to distinguish effects
of runtime dimensions from new callable or input types. Results for one helper
do not justify replacing every reduction with a loop. Ordinary numerical
reductions and DSL sums involving decision variables follow different paths.

## Reporting scope

Retain raw TOML, the source and manifest used, Julia version, machine details,
thread count and cache conditions. State which operation reaches compilation
inside its timer. Package loading, first model construction, result extraction
and warmed updates are different workloads; Julia allocation counts exclude
native solver memory.

These diagnostics explore first use and specialization. For repeated solves
and parallel throughput, use the [serial/parallel protocol](parallel-results.md).
