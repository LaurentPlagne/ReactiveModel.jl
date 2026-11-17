"""
    backend_status(backend::AbstractBackend) -> Symbol

Return the current solve status (`:optimal`, `:infeasible`, or another documented
nonoptimal status). The default reads `backend.status`. A backend with a different
representation overrides this method; RM does not require that field to exist.
"""
backend_status(backend::AbstractBackend) = getfield(backend, :status)

"""
    backend_objective_value(backend::AbstractBackend) -> Float64

Return the solver objective without RM's separately stored affine constant.
The default reads `backend.objective_value`. Values are meaningful only after
an optimal solve; public solution queries enforce their generation/status checks.
"""
backend_objective_value(backend::AbstractBackend) = getfield(backend, :objective_value)

"""
    backend_primal_value(backend::AbstractBackend, column::Int) -> Float64

Read one primal value by one-based RM column index without materializing a
collection. The default reads `backend.col_value[column]`. Override this method
to use another storage representation. A backend must finish extracting its
current solution before `backend_solve!` returns.
"""
backend_primal_value(backend::AbstractBackend, column::Int) = getfield(backend, :col_value)[column]

"""
    backend_dual_value(backend::AbstractBackend, row::Int) -> Float64

Read one dual of a normalized constraint row, using one-based RM indices. The
default reads `backend.row_dual[row]`. Backends without row duals return `NaN`
and report `backend_has_duals(backend) == false`.
"""
backend_dual_value(backend::AbstractBackend, row::Int) = getfield(backend, :row_dual)[row]
