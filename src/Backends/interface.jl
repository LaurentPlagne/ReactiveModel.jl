function initialize! end
function backend_set_verbose_flag! end
function backend_change_col_bounds! end
function backend_change_cols_bounds! end
function backend_change_cols_bounds_scratch end
function update_variable_block_bounds! end
function update_constraint_rhs! end
function update_constraint_block_rhs! end
function update_constraint_variable_coefficients! end
function backend_update_objective_coefficients! end
function backend_set_objective_sense! end
function backend_solve! end
function backend_has_duals end
function backend_write_to_file end
function backend_flush_row_bounds! end

backend_has_duals(::Any) = true
backend_change_cols_bounds_scratch(::Any) = nothing
function backend_change_cols_bounds!(backend, indices, lowers, uppers, ::Any)
    return backend_change_cols_bounds!(backend, indices, lowers, uppers)
end
backend_flush_row_bounds!(::Any) = nothing
function update_constraint_block_rhs!(
    backend::AbstractBackend, block::_BackendConstraintBlock, offset::Int
)
    for (i, constraint) in enumerate(block.constraints)
        update_constraint_rhs!(backend, constraint, offset + i)
    end
    return nothing
end
_backend_has_duals(backend::AbstractBackend) = backend_has_duals(backend)

function _backend_constraint_bounds(constraint::_BackendConstraint)
    constraint.op === :eq && return constraint.rhs, constraint.rhs
    constraint.op === :le && return -Inf, constraint.rhs
    return constraint.rhs, Inf
end

