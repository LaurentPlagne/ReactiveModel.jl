module ReactiveModel

"""
    ReactiveModel

A small reactive linear-modeling layer.  A model is assembled once from named
variable blocks, constraint blocks and an objective.  The only operation an
application normally performs afterwards is `solve(problem, parameters)`.
Parameter snapshots and the `has_changed` predicates decide which numerical
parts are regenerated; the backend is modified only when a value really
changed.
With `safe_mode=true`, all numerical formulas are reevaluated independently of
the predicates, and a report records missed changes and unchanged evaluations.

The model description is independent from its numerical backend. HiGHS and
Clp implement the direct backend contract without an intermediate modeling
interface. Pure-Julia solver integrations are optional extensions.
"""

using Clp_jll: Clp_jll
using HiGHS_jll: HiGHS_jll
using SparseArrays: SparseArrays

public backend_status, backend_objective_value, backend_primal_value, backend_dual_value
public full_type_string

export AbstractBackend,
    HiGHSBackend,
    TinyHiGHSBackend,
    ClpBackend,
    ClpJLBackend,
    XpressBackend,
    CPLEXBackend,
    initialize!,
    backend_set_verbose_flag!,
    backend_change_col_bounds!,
    backend_change_cols_bounds!,
    update_variable_block_bounds!,
    update_constraint_rhs!,
    update_constraint_variable_coefficients!,
    backend_update_objective_coefficients!,
    backend_set_objective_sense!,
    backend_solve!,
    backend_has_duals,
    backend_write_to_file,
    backend_change_cols_bounds_scratch,
    backend_flush_row_bounds!,
    update_constraint_block_rhs!,
    backend_reset_basis!,
    AbstractVariable,
    AbstractConstraint,
    AbstractObjective,
    VariableTag,
    ConstraintTag,
    ObjectiveTag,
    lower_bound,
    upper_bound,
    variable_bounds,
    has_changed,
    lhs_has_changed,
    rhs_has_changed,
    lhs_terms!,
    rhs_value,
    objective_terms!,
    VariableBlock,
    ConstraintBlock,
    Objective,
    compose_objectives,
    ProblemDefinition,
    ConstData,
    Parameters,
    Workspace,
    LinearProblem,
    SolveError,
    UpdateCheck,
    update_report,
    reset_update_report!,
    @coeff,
    set_coeff!,
    solve,
    rebuild!,
    objective_value,
    termination_status,
    is_optimal,
    primal,
    dual,
    activity,
    iteration_space,
    constraint_sense,
    @model,
    @variable,
    @external_variables,
    @constraint,
    @objective,
    variable_block_offset,
    constraint_block_offset,
    variable_block_keys,
    constraint_block_keys

# Includes follow type dependencies; backend implementations remain optional
# where their external solver packages are optional. See the design manual.
include("Backends/types.jl")
include("Backends/results.jl")
include("core/definitions.jl")
include("UpdateSafety.jl")
include("core/callbacks.jl")
include("core/results.jl")
include("core/problem.jl")
include("Backends/interface.jl")
include("Backends/assembly.jl")
include("Backends/HiGHSBackend.jl")
include("Backends/TinyHiGHSBackend.jl")
include("Backends/ClpBackend.jl")
include("Backends/ClpJLBackend.jl")
include("Backends/XpressBackend.jl")
include("Backends/CPLEXBackend.jl")
include("core/assembly.jl")
include("core/updates.jl")
include("core/execution.jl")
include("core/solve.jl")
include("DSL/ast_utils.jl")
include("DSL/codegen.jl")
include("core/display.jl")

end # module ReactiveModel
