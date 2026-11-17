# Optional developer check. Requires JET in the active diagnostic environment.
using JET
using ReactiveModel
include("trim/workload.jl")

println("Julia ", VERSION, "; JET ", pkgversion(JET))
for (label, f, types) in (
    ("LP construction, updates, results and rebuild", ReactiveModelTrimWorkload.verify_lp, Tuple{Int}),
    ("Functional API", ReactiveModelTrimWorkload.verify_functional, Tuple{Int}),
    ("Affine aliases, colliding keys and module identities", ReactiveModelTrimWorkload.verify_math, Tuple{Int}),
    ("Filtered dependent families and nested sums", ReactiveModelTrimWorkload.verify_grammar, Tuple{Int}),
    ("MILP DSL", ReactiveModelTrimWorkload.verify_milp, Tuple{}),
    ("Safe updates and report", ReactiveModelTrimWorkload.verify_safe_updates, Tuple{Int}),
    ("Bounded family groups", ReactiveModelTrimWorkload.verify_grouped_rules, Tuple{Int}),
    ("Compact DSL execution, callback errors, clones and data replacement", ReactiveModelTrimWorkload.verify_execution, Tuple{Int}),
    ("Interrupted rebuild recovery and stale scalar properties", ReactiveModelTrimWorkload.verify_rebuild_recovery, Tuple{Int}),
    ("Numeric bound pairs, compact callbacks and mixed precisions", ReactiveModelTrimWorkload.verify_numeric_bounds, Tuple{Int}),
    ("Five result key types and retained family snapshots", ReactiveModelTrimWorkload.verify_many_result_keys, Tuple{Int}),
    ("Composed objectives, constants, rebuild and independent data", ReactiveModelTrimWorkload.verify_composed_objective, Tuple{Int}),
)
    println("\n", label)
    result = JET.report_opt(f, types; target_modules=(ReactiveModel, ReactiveModelTrimWorkload))
    show(stdout, MIME"text/plain"(), result)
    println()
end

# The execution boundary deliberately hides payload types from client methods.
# Inspect each concrete adapter too; a clean client report alone would no longer
# demonstrate inference inside the callback engine.
compact_data=(scale=1.0, capacity=10.0)
grouped_data=(n=1,)
for (label, problem) in (
    ("Typed adapter: mixed keys and rebuild", ReactiveModelTrimWorkload.lp_problem(1)),
    ("Typed adapter: functional callbacks", ReactiveModelTrimWorkload.functional_problem(1)),
    ("Typed adapter: composed objectives", ReactiveModelTrimWorkload.composed_problem(1)),
    ("Typed adapter: five key types", ReactiveModelTrimWorkload.many_key_problem(1)),
    ("Typed adapter: grouped DSL rules", LinearProblem(
        ConstData(ReactiveModelTrimWorkload.GroupedRules(grouped_data), grouped_data),
        (limit=10.0, rhs=1.0),
    )),
    ("Typed adapter: DSL exceptions and data replacement", LinearProblem(
        ConstData(ReactiveModelTrimWorkload.CompactExecution(compact_data), compact_data),
        (rhs=1.0, fail=false); safe_mode=true,
    )),
)
    println("\n", label)
    payload = problem.numerical_updates.payload[]
    types = Tuple{typeof(payload),UInt8,ReactiveModel._ExecutionArgument}
    result = JET.report_opt(ReactiveModel._ExecutionThunk(), types;
        target_modules=(ReactiveModel, ReactiveModelTrimWorkload))
    show(stdout, MIME"text/plain"(), result)
    println()
end
