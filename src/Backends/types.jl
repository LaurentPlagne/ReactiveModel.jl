"""
    AbstractBackend

Abstract supertype of persistent linear-programming backends.

A concrete backend implements the update/solve operations below and the result
protocol [`backend_status`](@ref), [`backend_objective_value`](@ref),
[`backend_primal_value`](@ref) and [`backend_dual_value`](@ref). The default result
methods use the conventional native-backend fields; another representation can
override them. `HiGHSBackend` is the reference implementation.

| Function | Called | Frequency |
|---|---|---|
| `backend_set_verbose_flag!` | before `initialize!` | once per model build |
| `initialize!` | first solve, rebuild, or recovery from a backend exception | once per model build |
| `backend_change_col_bounds!` | 1-3 columns of one block move | per changed block |
| `update_variable_block_bounds!` | more than 3 columns of one block move | per changed block |
| `update_constraint_rhs!` | 1-3 rows of one block move | per changed block |
| `update_constraint_block_rhs!` | more than 3 rows of one block move | per changed block |
| `update_constraint_variable_coefficients!` | one row's coefficients move | per changed row |
| `backend_update_objective_coefficients!` | the objective moves | per solve where it does |
| `backend_set_objective_sense!` | the objective sense moves | never, in practice |
| `backend_solve!` | every solve where anything changed | once per solve |
| `backend_write_to_file` | on demand, for debugging | — |

`backend_has_duals` defaults to `true` and needs overriding only when it can
be `false` (a MIP, for `HiGHSBackend`). `backend_change_cols_bounds!`,
`backend_change_cols_bounds_scratch` and `backend_flush_row_bounds!` are part
of the exported surface but are not called by this package's own update path
today; a backend may leave them at their defaults.
"""
abstract type AbstractBackend end

# These small mutable records are the RM-side equivalent of the numerical
# records consumed by the LP backends.  The backend contract below deliberately
# operates on them rather than on the mathematical definition objects.
mutable struct _BackendVariableBlock
    domain::Symbol
    offset::Int
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    values::Vector{Float64}
end

mutable struct _BackendConstraint
    variable_indices::Vector{Int}
    variable_coefficients::Vector{Float64}
    rhs::Float64
    op::Symbol
end

mutable struct _BackendConstraintBlock
    offset::Int
    constraints::Vector{_BackendConstraint}
    values::Vector{Float64}
end

mutable struct _BackendObjectiveCoefficients
    variable_indices::Vector{Int}
    variable_coefficients::Vector{Float64}
end

mutable struct _BackendObjective
    coefficients::_BackendObjectiveCoefficients
    op::Symbol
end
