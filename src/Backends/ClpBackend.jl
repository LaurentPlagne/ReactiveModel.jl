# ---------------------------------------------------------------------------
# ClpBackend: Direct C-API backend for COIN-OR Clp solver
# ---------------------------------------------------------------------------

"""
Direct C-API backend for the COIN-OR Clp LP solver, implementing ReactiveModel's
`AbstractBackend` interface.

The package calls the native C API through `Clp_jll` directly. `Clp.jl` is
used only in the separate benchmark environment for comparisons through JuMP;
it is not a dependency of the ReactiveModel package.

Column and row bounds are updated locally in buffer arrays and flushed lazily to Clp
upon `solve!`, amortizing the `ccall` cost.
"""

const _CLP_LIB = get(ENV, "REACTIVEMODEL_CLP_LIB", Clp_jll.libClp)

mutable struct ClpBackend <: AbstractBackend
    inner::Ptr{Cvoid}
    col_value::Vector{Cdouble}
    col_dual::Vector{Cdouble}
    row_value::Vector{Cdouble}
    row_dual::Vector{Cdouble}
    a_start::Vector{Cint}
    a_index::Vector{Cint}
    a_value::Vector{Cdouble}
    cc::Vector{Cdouble}
    cl::Vector{Cdouble}
    cu::Vector{Cdouble}
    rl::Vector{Cdouble}
    ru::Vector{Cdouble}
    cols_dirty::Bool
    rows_dirty::Bool
    solved::Bool
    status::Symbol
    objective_value::Float64
    saved_basis::Vector{Cuchar}
    has_saved_basis::Bool
    # Coefficient setters update the native model. Refresh cached CSC arrays
    # from Clp before a reload when they no longer describe that matrix.
    matrix_modified::Bool
end

"""
    ClpBackend()

A persistent Clp model, empty until the first `initialize!`.
"""
function ClpBackend()
    backend = ClpBackend(
        C_NULL,
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cint[],
        Cint[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        false,
        false,
        false,
        :not_solved,
        NaN,
        Cuchar[],
        false,
        false,
    )
    _clp_empty!(backend)
    finalizer(_clp_finalize!, backend)
    return backend
end

Base.unsafe_convert(::Type{Ptr{Cvoid}}, backend::ClpBackend) = backend.inner

function _clp_finalize!(backend::ClpBackend)
    if backend.inner != C_NULL
        ccall((:Clp_deleteModel, _CLP_LIB), Cvoid, (Ptr{Cvoid},), backend.inner)
        backend.inner = C_NULL
    end
    return nothing
end

function _clp_empty!(backend::ClpBackend)
    _clp_finalize!(backend)
    backend.inner = ccall((:Clp_newModel, _CLP_LIB), Ptr{Cvoid}, ())
    backend.inner == C_NULL && error("Clp_newModel returned a null model")
    ccall((:Clp_setLogLevel, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Cint), backend.inner, Cint(0))
    # Retain the existing perturbation policy and native iteration/time limits.
    ccall((:Clp_setPerturbation, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Cint), backend.inner, Cint(50))
    backend.cols_dirty = false
    backend.rows_dirty = false
    backend.solved = false
    return nothing
end

function backend_set_verbose_flag!(backend::ClpBackend; verbose::Bool)
    ccall(
        (:Clp_setLogLevel, _CLP_LIB),
        Cvoid,
        (Ptr{Cvoid}, Cint),
        backend.inner,
        verbose ? Cint(1) : Cint(0),
    )
    return nothing
end

backend_has_duals(::ClpBackend) = true

function _clp_variable_domain_error(variable_blocks)
    all(block -> block.domain === :continuous, variable_blocks) ||
        throw(ArgumentError("ClpBackend supports continuous variables only"))
    return nothing
end

function initialize!(
    backend::ClpBackend,
    variable_blocks::Vector{_BackendVariableBlock},
    constraint_blocks::Vector{_BackendConstraintBlock},
    objective::_BackendObjective,
)
    _clp_variable_domain_error(variable_blocks)

    arrays=_backend_arrays(variable_blocks, constraint_blocks, objective)
    num_col, num_row=arrays.ncols, arrays.nrows
    backend.a_start=arrays.starts
    backend.a_index=arrays.indices
    backend.a_value=arrays.values
    backend.cc=arrays.cost
    backend.cl=arrays.lower
    backend.cu=arrays.upper
    backend.rl=arrays.row_lower
    backend.ru=arrays.row_upper

    _clp_empty!(backend)
    a_start, a_index, a_value = backend.a_start, backend.a_index, backend.a_value
    cl, cu, cc, rl, ru = backend.cl, backend.cu, backend.cc, backend.rl, backend.ru
    GC.@preserve a_start a_index a_value cl cu cc rl ru begin
        ccall(
            (:Clp_loadProblem, _CLP_LIB),
            Cvoid,
            (
                Ptr{Cvoid},
                Cint,
                Cint,
                Ptr{Cint},
                Ptr{Cint},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
            ),
            backend.inner,
            Cint(num_col),
            Cint(num_row),
            a_start,
            a_index,
            a_value,
            cl,
            cu,
            cc,
            rl,
            ru,
        )
    end
    sense = objective.op === :max ? -1.0 : 1.0
    ccall((:Clp_setObjSense, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Cdouble), backend.inner, sense)

    backend.col_value = zeros(Cdouble, num_col)
    backend.col_dual = zeros(Cdouble, num_col)
    backend.row_value = zeros(Cdouble, num_row)
    backend.row_dual = zeros(Cdouble, num_row)
    backend.cols_dirty = false
    backend.rows_dirty = false
    backend.solved = false
    backend.status = :not_solved
    backend.saved_basis = zeros(Cuchar, num_col + num_row)
    backend.has_saved_basis = false
    backend.matrix_modified = false
    return nothing
end

function backend_change_col_bounds!(backend::ClpBackend, index::Int, lower::Float64, upper::Float64)
    backend.cl[index] = lower
    backend.cu[index] = upper
    backend.cols_dirty = true
    return nothing
end

function backend_change_cols_bounds!(
    backend::ClpBackend, indices::Vector{Int}, lowers::Vector{Float64}, uppers::Vector{Float64}
)
    for i in eachindex(indices)
        backend.cl[indices[i]] = lowers[i]
        backend.cu[indices[i]] = uppers[i]
    end
    backend.cols_dirty = true
    return nothing
end

function update_variable_block_bounds!(
    backend::ClpBackend, block::_BackendVariableBlock, offset::Int
)
    copyto!(backend.cl, offset + 1, block.lower_bounds, 1, length(block.lower_bounds))
    copyto!(backend.cu, offset + 1, block.upper_bounds, 1, length(block.upper_bounds))
    backend.cols_dirty = true
    return nothing
end

function _clp_flush_col_bounds!(backend::ClpBackend)
    backend.cols_dirty || return nothing
    cl, cu = backend.cl, backend.cu
    GC.@preserve cl begin
        ccall((:Clp_chgColumnLower, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), backend.inner, cl)
    end
    GC.@preserve cu begin
        ccall((:Clp_chgColumnUpper, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), backend.inner, cu)
    end
    backend.cols_dirty = false
    return nothing
end

function update_constraint_rhs!(backend::ClpBackend, constraint::_BackendConstraint, row::Int)
    backend.rl[row], backend.ru[row] = _backend_constraint_bounds(constraint)
    backend.rows_dirty = true
    return nothing
end

function _clp_flush_row_bounds!(backend::ClpBackend)
    backend.rows_dirty || return nothing
    rl, ru = backend.rl, backend.ru
    GC.@preserve rl begin
        ccall((:Clp_chgRowLower, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), backend.inner, rl)
    end
    GC.@preserve ru begin
        ccall((:Clp_chgRowUpper, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), backend.inner, ru)
    end
    backend.rows_dirty = false
    return nothing
end

function update_constraint_variable_coefficients!(
    backend::ClpBackend, row::Int, constraint::_BackendConstraint
)
    for (column, coefficient) in zip(constraint.variable_indices, constraint.variable_coefficients)
        ccall(
            (:Clp_modifyCoefficient, _CLP_LIB),
            Cvoid,
            (Ptr{Cvoid}, Cint, Cint, Cdouble, Cuchar),
            backend.inner,
            Cint(row - 1),
            Cint(column - 1),
            coefficient,
            Cuchar(0),
        )
    end
    isempty(constraint.variable_indices) || (backend.matrix_modified = true)
    return nothing
end

function backend_update_objective_coefficients!(
    backend::ClpBackend, coefficients::_BackendObjectiveCoefficients
)
    for i in eachindex(coefficients.variable_indices)
        backend.cc[coefficients.variable_indices[i]] = coefficients.variable_coefficients[i]
    end
    cc = backend.cc
    GC.@preserve cc begin
        ccall(
            (:Clp_chgObjCoefficients, _CLP_LIB),
            Cvoid,
            (Ptr{Cvoid}, Ptr{Cdouble}),
            backend.inner,
            cc,
        )
    end
    return nothing
end

function backend_set_objective_sense!(backend::ClpBackend, sense::Symbol)
    csense = sense === :max ? -1.0 : 1.0
    ccall((:Clp_setObjSense, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Cdouble), backend.inner, csense)
    return nothing
end

# Secondary statuses 2, 3 and 4 indicate that optimality applies only to the
# scaled model while the original model retains infeasibilities. Check this
# status as well as Clp_isProvenOptimal before accepting the solution.
function _clp_proven_optimal(backend::ClpBackend)
    ccall((:Clp_isProvenOptimal, _CLP_LIB), Cint, (Ptr{Cvoid},), backend.inner) != 0 || return false
    secondary = ccall((:Clp_secondaryStatus, _CLP_LIB), Cint, (Ptr{Cvoid},), backend.inner)
    return !(2 <= secondary <= 4)
end

function backend_solve!(
    backend::ClpBackend,
    variable_blocks::Vector{_BackendVariableBlock},
    variable_offsets::Vector{Int},
    constraint_blocks::Vector{_BackendConstraintBlock},
    constraint_offsets::Vector{Int},
)
    _clp_flush_col_bounds!(backend)
    _clp_flush_row_bounds!(backend)

    optimal = false
    if backend.solved
        # Reoptimize the existing native LP with its current basis.
        ccall((:Clp_dual, _CLP_LIB), Cint, (Ptr{Cvoid}, Cint), backend.inner, Cint(0))
        optimal = _clp_proven_optimal(backend)
        if !optimal
            # Retry the initial solve on the current matrix after a failed warm solve.
            ccall((:Clp_initialSolve, _CLP_LIB), Cint, (Ptr{Cvoid},), backend.inner)
            optimal =
                ccall((:Clp_isProvenOptimal, _CLP_LIB), Cint, (Ptr{Cvoid},), backend.inner) != 0
        end
    else
        ccall((:Clp_initialSolve, _CLP_LIB), Cint, (Ptr{Cvoid},), backend.inner)
        optimal = ccall((:Clp_isProvenOptimal, _CLP_LIB), Cint, (Ptr{Cvoid},), backend.inner) != 0
    end

    if !optimal
        backend.solved = false
        backend.has_saved_basis = false
        backend.status = :infeasible
        backend.objective_value = NaN
        return false, NaN
    end

    backend.solved = true

    col_ptr = ccall((:Clp_getColSolution, _CLP_LIB), Ptr{Cdouble}, (Ptr{Cvoid},), backend.inner)
    row_ptr = ccall((:Clp_getRowActivity, _CLP_LIB), Ptr{Cdouble}, (Ptr{Cvoid},), backend.inner)
    dual_ptr = ccall((:Clp_getRowPrice, _CLP_LIB), Ptr{Cdouble}, (Ptr{Cvoid},), backend.inner)
    (col_ptr == C_NULL || row_ptr == C_NULL || dual_ptr == C_NULL) &&
        error("Clp returned a null solution pointer")
    unsafe_copyto!(pointer(backend.col_value), col_ptr, length(backend.col_value))
    unsafe_copyto!(pointer(backend.row_value), row_ptr, length(backend.row_value))
    unsafe_copyto!(pointer(backend.row_dual), dual_ptr, length(backend.row_dual))
    fill!(backend.col_dual, 0.0)

    for block in variable_blocks
        unsafe_copyto!(
            pointer(block.values),
            pointer(backend.col_value, block.offset + 1),
            length(block.values),
        )
    end
    for block in constraint_blocks
        unsafe_copyto!(
            pointer(block.values), pointer(backend.row_dual, block.offset + 1), length(block.values)
        )
    end

    backend.status = :optimal
    backend.objective_value = ccall(
        (:Clp_getObjValue, _CLP_LIB), Cdouble, (Ptr{Cvoid},), backend.inner
    )
    return true, backend.objective_value
end

function backend_write_to_file(backend::ClpBackend, filename)
    status = ccall(
        (:Clp_writeMps, _CLP_LIB),
        Cint,
        (Ptr{Cvoid}, Cstring, Cint, Cint, Cdouble),
        backend.inner,
        filename,
        Cint(0),
        Cint(2),
        Cdouble(1.0),
    )
    status == 0 || @warn("Clp_writeMps failed with status $status")
    return nothing
end

"""
    _clp_refresh_matrix!(backend::ClpBackend)

Copy the current native matrix into compact CSC buffers before resetting the
solver. Incremental coefficient changes update Clp directly, while cached
`a_start`, `a_index` and `a_value` may still describe an earlier matrix. Loading
those stale arrays would undo changes that RM's current row dictionaries already
consider applied. This refresh keeps reload/recovery consistent with the latest LP.
"""
function _clp_refresh_matrix!(backend::ClpBackend)
    m = backend.inner
    num_col = length(backend.cl)
    num_col == 0 && return nothing
    starts = unsafe_wrap(
        Array, ccall((:Clp_getVectorStarts, _CLP_LIB), Ptr{Cint}, (Ptr{Cvoid},), m), num_col
    )
    lengths = unsafe_wrap(
        Array, ccall((:Clp_getVectorLengths, _CLP_LIB), Ptr{Cint}, (Ptr{Cvoid},), m), num_col
    )
    storage_size = maximum(Int(starts[j]) + Int(lengths[j]) for j in 1:num_col)
    indices = unsafe_wrap(
        Array, ccall((:Clp_getIndices, _CLP_LIB), Ptr{Cint}, (Ptr{Cvoid},), m), storage_size
    )
    elements = unsafe_wrap(
        Array, ccall((:Clp_getElements, _CLP_LIB), Ptr{Cdouble}, (Ptr{Cvoid},), m), storage_size
    )
    nnz = sum(Int, lengths)
    a_start = Vector{Cint}(undef, num_col + 1)
    a_index = Vector{Cint}(undef, nnz)
    a_value = Vector{Cdouble}(undef, nnz)
    k = 0
    for j in 1:num_col
        a_start[j] = Cint(k)
        for p in (Int(starts[j]) + 1):(Int(starts[j]) + Int(lengths[j]))
            k += 1
            a_index[k] = indices[p]
            a_value[k] = elements[p]
        end
    end
    a_start[num_col + 1] = Cint(k)
    backend.a_start, backend.a_index, backend.a_value = a_start, a_index, a_value
    backend.matrix_modified = false
    return nothing
end

function _clp_reload!(backend::ClpBackend)
    # Preserve the current objective sense before replacing the native model.
    objective_sense = ccall((:Clp_getObjSense, _CLP_LIB), Cdouble, (Ptr{Cvoid},), backend.inner)
    # Refresh native coefficient changes before loading the cached matrix.
    backend.matrix_modified && _clp_refresh_matrix!(backend)
    _clp_empty!(backend)
    num_col = length(backend.cl)
    num_row = length(backend.rl)
    a_start, a_index, a_value = backend.a_start, backend.a_index, backend.a_value
    cl, cu, cc, rl, ru = backend.cl, backend.cu, backend.cc, backend.rl, backend.ru
    GC.@preserve a_start a_index a_value cl cu cc rl ru begin
        ccall(
            (:Clp_loadProblem, _CLP_LIB),
            Cvoid,
            (
                Ptr{Cvoid},
                Cint,
                Cint,
                Ptr{Cint},
                Ptr{Cint},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
                Ptr{Cdouble},
            ),
            backend.inner,
            Cint(num_col),
            Cint(num_row),
            a_start,
            a_index,
            a_value,
            cl,
            cu,
            cc,
            rl,
            ru,
        )
    end
    ccall(
        (:Clp_setObjSense, _CLP_LIB), Cvoid, (Ptr{Cvoid}, Cdouble), backend.inner, objective_sense
    )
    backend.cols_dirty = false
    backend.rows_dirty = false
    backend.solved = false
    backend.has_saved_basis = false
    return nothing
end

function backend_reset_basis!(backend::ClpBackend)
    if backend.inner != C_NULL && !isempty(backend.a_start)
        _clp_reload!(backend)
    else
        backend.solved = false
        backend.has_saved_basis = false
    end
    return nothing
end
