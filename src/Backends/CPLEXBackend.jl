# ---------------------------------------------------------------------------
# CPLEXBackend: IBM ILOG CPLEX backend for MILP and LP
# ---------------------------------------------------------------------------

"""
    CPLEXBackend(; kwargs...)

Persistent IBM ILOG CPLEX backend for Mixed-Integer Linear Programs (MILP) and LPs,
wrapping `CPLEX.jl` native C API.

Requires `using CPLEX` to load the `ReactiveModelCPLEXExt` extension.

This optional integration has not been qualified against a licensed solver.
Incremental numerical updates and export are unsupported and raise
`ArgumentError`; use `rebuild!` to load changed data.
"""
mutable struct CPLEXBackend <: AbstractBackend
    env::Any
    lp::Any
    col_value::Vector{Float64}
    col_dual::Vector{Float64}
    row_value::Vector{Float64}
    row_dual::Vector{Float64}
    is_mip::Bool
    objective_value::Float64
    status::Symbol
    options::Dict{Symbol, Any}

    function CPLEXBackend(env, lp, col_value, col_dual, row_value, row_dual,
                          is_mip, objective_value, status, options)
        new(env, lp, col_value, col_dual, row_value, row_dual,
            is_mip, objective_value, status, options)
    end
end

function CPLEXBackend(args...; kwargs...)
    error("`CPLEXBackend` requires `CPLEX.jl` and IBM ILOG CPLEX Studio to be installed and loaded.\nPlease run `using CPLEX` or `import CPLEX` first.")
end

# Status codes from CPLEX.jl's generated C API (CPX_STAT_* and CPXMIP_*).
# Unknown outcomes remain errors rather than being reported as infeasible.
function _cplex_status(is_mip::Bool, status::Integer)
    if is_mip
        status in (101, 102) && return :optimal
        status == 103 && return :infeasible
        status == 118 && return :unbounded
        status in (119, 133) && return :infeasible_or_unbounded
        status in (107, 108, 131, 132) && return :time_limit
        status in (105, 106) && return :node_limit
        status in (104, 128) && return :solution_limit
        status in (111, 112, 117) && return :memory_limit
        status in (113, 114) && return :interrupted
        status == 127 && return :feasible
    else
        status == 1 && return :optimal
        status == 3 && return :infeasible
        status == 2 && return :unbounded
        status == 4 && return :infeasible_or_unbounded
        status in (11, 25) && return :time_limit
        status == 10 && return :iteration_limit
        status in (12, 21, 22) && return :objective_limit
        status == 13 && return :interrupted
        status in (6, 41) && return :numerical_error
        status == 23 && return :feasible
    end
    return :other_error
end

_cplex_unsupported(operation::Symbol) = throw(ArgumentError(
    "CPLEXBackend does not implement $operation; use rebuild! to load changed data " *
    "or choose a backend supporting this operation",
))

backend_change_col_bounds!(::CPLEXBackend, ::Int, ::Float64, ::Float64) =
    _cplex_unsupported(:backend_change_col_bounds!)
backend_change_cols_bounds!(::CPLEXBackend, ::Vector{Int}, ::Vector{Float64}, ::Vector{Float64}) =
    _cplex_unsupported(:backend_change_cols_bounds!)
update_variable_block_bounds!(::CPLEXBackend, ::_BackendVariableBlock, ::Int) =
    _cplex_unsupported(:update_variable_block_bounds!)
update_constraint_rhs!(::CPLEXBackend, ::_BackendConstraint, ::Int) =
    _cplex_unsupported(:update_constraint_rhs!)
update_constraint_block_rhs!(::CPLEXBackend, ::_BackendConstraintBlock, ::Int) =
    _cplex_unsupported(:update_constraint_block_rhs!)
update_constraint_variable_coefficients!(::CPLEXBackend, ::Int, ::_BackendConstraint) =
    _cplex_unsupported(:update_constraint_variable_coefficients!)
backend_update_objective_coefficients!(::CPLEXBackend, ::_BackendObjectiveCoefficients) =
    _cplex_unsupported(:backend_update_objective_coefficients!)
backend_write_to_file(::CPLEXBackend, ::String) = _cplex_unsupported(:backend_write_to_file)
