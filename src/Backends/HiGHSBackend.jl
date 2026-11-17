include("HiGHSCAPI.jl")

# ---------------------------------------------------------------------------
# HiGHSBackend: Direct C-API backend for HiGHS
# ---------------------------------------------------------------------------

"""
    HiGHSBackend(; options...)

Persistent LP/MILP backend using the HiGHS C API from `HiGHS_jll` directly.
Neither `HiGHS.jl` nor MathOptInterface is loaded by this backend. Keywords
are native HiGHS options, for example `time_limit=10.0` or `presolve="off"`.
RM records successfully applied option values in a concrete string-pair vector
so an independent problem created by `similar` retains that configuration.
RM's initialization policy still selects sequential simplex and logging options;
native calls made outside RM are not tracked in this configuration record.
"""
mutable struct HiGHSBackend <: AbstractBackend
    inner::Ptr{Cvoid}
    col_value::Vector{Float64}
    col_dual::Vector{Float64}
    row_value::Vector{Float64}
    row_dual::Vector{Float64}
    a_start::Vector{Cint}
    a_index::Vector{Cint}
    a_value::Vector{Cdouble}
    cc::Vector{Cdouble}
    cl::Vector{Cdouble}
    cu::Vector{Cdouble}
    rl::Vector{Cdouble}
    ru::Vector{Cdouble}
    integrality::Vector{Cint}
    cint_indices::Vector{Cint}
    is_mip::Bool
    objective_value::Float64
    status::Symbol
    options::Vector{Pair{String,String}}
end

function HiGHSBackend(; kwargs...)
    inner = _HiGHSCAPI.Highs_create()
    inner == C_NULL && error("HiGHS model allocation failed")
    if _HiGHSCAPI.Highs_getSizeofHighsInt(inner) != sizeof(Cint)
        _HiGHSCAPI.Highs_destroy(inner)
        error("HiGHSBackend requires a HiGHS library with 32-bit HighsInt")
    end
    backend = HiGHSBackend(
        inner,
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Cint[],
        Cint[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cdouble[],
        Cint[],
        Cint[],
        false,
        NaN,
        :not_solved,
        Pair{String,String}[],
    )
    finalizer(_highs_finalize!, backend)
    # HiGHS shares its scheduler between instances, including JuMP models.
    # Enforcing a different thread count here causes Highs_run to fail after a JuMP MILP.
    # The simplex remains sequential via `parallel=off` below.
    try
        _highs_set_option!(backend, "output_flag", false)
        for (name, value) in pairs(kwargs)
            _highs_set_option!(backend, String(name), value)
        end
    catch
        _highs_finalize!(backend)
        rethrow()
    end
    return backend
end

function _highs_finalize!(backend::HiGHSBackend)
    if backend.inner != C_NULL
        _HiGHSCAPI.Highs_destroy(backend.inner)
        backend.inner = C_NULL
    end
    return nothing
end

backend_has_duals(backend::HiGHSBackend) = !backend.is_mip

Base.unsafe_convert(::Type{Ptr{Cvoid}}, backend::HiGHSBackend) = backend.inner

function _highs_check_status(status, context::AbstractString)
    status == _HiGHSCAPI.kHighsStatusError && error("$context failed with status $status")
    return nothing
end

function _highs_set_option!(backend::HiGHSBackend, name::AbstractString, value)
    # Let HiGHS parse the value according to the option's native type. This
    # preserves, for example, time_limit=10 as well as time_limit=10.0.
    status = _HiGHSCAPI.Highs_setStringOptionValue(backend, name, string(value))
    _highs_check_status(status, "HiGHS option $name=$value")
    setting=String(name)=>string(value)
    index=findfirst(pair->pair.first==name, backend.options)
    if index===nothing
        push!(backend.options, setting)
    else
        backend.options[index]=setting
    end
    return nothing
end

"""Clone native option values set through RM, without copying a model or basis."""
function _similar_backend(backend::HiGHSBackend)
    cloned=HiGHSBackend()
    try
        for (name, value) in backend.options
            _highs_set_option!(cloned, name, value)
        end
    catch
        _highs_finalize!(cloned)
        rethrow()
    end
    return cloned
end

function backend_set_verbose_flag!(backend::HiGHSBackend; verbose::Bool)
    _highs_set_option!(backend, "parallel", "off")
    _highs_set_option!(backend, "solver", "simplex")
    _highs_set_option!(backend, "log_to_console", false)
    _highs_set_option!(backend, "simplex_max_concurrency", 1)
    _highs_set_option!(backend, "output_flag", verbose)
    return nothing
end

function _backend_empty!(backend::HiGHSBackend)
    backend.inner == C_NULL && (backend.inner = _HiGHSCAPI.Highs_create())
    status = _HiGHSCAPI.Highs_clearModel(backend)
    _highs_check_status(status, "HiGHS model clear")
    return nothing
end

function initialize!(
    backend::HiGHSBackend,
    variable_blocks::Vector{_BackendVariableBlock},
    constraint_blocks::Vector{_BackendConstraintBlock},
    objective::_BackendObjective,
)
    _backend_empty!(backend)
    arrays=_backend_arrays(variable_blocks, constraint_blocks, objective; protect_tiny=true)
    num_col, num_row=arrays.ncols, arrays.nrows
    backend.a_start=arrays.starts
    backend.a_index=arrays.indices
    backend.a_value=arrays.values
    backend.cc=arrays.cost
    backend.cl=arrays.lower
    backend.cu=arrays.upper
    backend.rl=arrays.row_lower
    backend.ru=arrays.row_upper
    backend.integrality=arrays.integrality
    backend.is_mip=arrays.is_mip
    is_mip=arrays.is_mip
    sense = if objective.op === :max
        _HiGHSCAPI.kHighsObjSenseMaximize
    else
        _HiGHSCAPI.kHighsObjSenseMinimize
    end
    status = if is_mip
        _HiGHSCAPI.Highs_passMip(
            backend,
            Cint(num_col),
            Cint(num_row),
            Cint(length(backend.a_value)),
            Cint(_HiGHSCAPI.kHighsMatrixFormatColwise),
            Cint(sense),
            Cdouble(0),
            backend.cc,
            backend.cl,
            backend.cu,
            backend.rl,
            backend.ru,
            backend.a_start,
            backend.a_index,
            backend.a_value,
            backend.integrality,
        )
    else
        _HiGHSCAPI.Highs_passLp(
            backend,
            Cint(num_col),
            Cint(num_row),
            Cint(length(backend.a_value)),
            Cint(_HiGHSCAPI.kHighsMatrixFormatColwise),
            Cint(sense),
            Cdouble(0),
            backend.cc,
            backend.cl,
            backend.cu,
            backend.rl,
            backend.ru,
            backend.a_start,
            backend.a_index,
            backend.a_value,
        )
    end
    _highs_check_status(status, "HiGHS model import")
    for (row, column, coefficient) in arrays.tiny
        status = _HiGHSCAPI.Highs_changeCoeff(backend, Cint(row - 1), Cint(column - 1), coefficient)
        _highs_check_status(status, "HiGHS tiny coefficient restoration")
    end
    backend.col_value = zeros(Cdouble, num_col)
    backend.col_dual = zeros(Cdouble, num_col)
    backend.row_value = zeros(Cdouble, num_row)
    backend.row_dual = zeros(Cdouble, num_row)
    backend.cint_indices = Cint[]
    backend.status = :not_solved
    return nothing
end

function backend_change_col_bounds!(
    backend::HiGHSBackend, index::Int, lower::Float64, upper::Float64
)
    status = _HiGHSCAPI.Highs_changeColBounds(backend, Cint(index - 1), lower, upper)
    status == _HiGHSCAPI.kHighsStatusError && error(
        "HiGHS column bounds update failed: column=$index, lower=$lower, upper=$upper, status=$status",
    )
    return nothing
end

function backend_change_cols_bounds!(
    backend::HiGHSBackend, indices::Vector{Int}, lowers::Vector{Float64}, uppers::Vector{Float64}
)
    n = length(indices)
    length(backend.cint_indices) < n && resize!(backend.cint_indices, n)
    for i in 1:n
        backend.cint_indices[i] = Cint(indices[i] - 1)
    end
    status = _HiGHSCAPI.Highs_changeColsBoundsBySet(
        backend, Cint(n), backend.cint_indices, lowers, uppers
    )
    if status == _HiGHSCAPI.kHighsStatusError
        first_updates = [(indices[i], lowers[i], uppers[i]) for i in 1:min(n, 5)]
        error(
            "HiGHS column bounds update failed: count=$n, first (column, lower, upper)=$first_updates, status=$status",
        )
    end
    return nothing
end

function update_variable_block_bounds!(
    backend::HiGHSBackend, block::_BackendVariableBlock, offset::Int
)
    from_col = Cint(offset)
    to_col = Cint(offset + length(block.lower_bounds) - 1)
    status = _HiGHSCAPI.Highs_changeColsBoundsByRange(
        backend, from_col, to_col, block.lower_bounds, block.upper_bounds
    )
    _highs_check_status(status, "HiGHS variable block bounds update")
    return nothing
end

function update_constraint_rhs!(backend::HiGHSBackend, constraint::_BackendConstraint, row::Int)
    lower, upper = _backend_constraint_bounds(constraint)
    status = _HiGHSCAPI.Highs_changeRowBounds(backend, Cint(row - 1), lower, upper)
    _highs_check_status(status, "HiGHS row bounds update")
    return nothing
end

function update_constraint_block_rhs!(
    backend::HiGHSBackend, block::_BackendConstraintBlock, offset::Int
)
    n = length(block.constraints)
    length(backend.rl) < n && (resize!(backend.rl, n); resize!(backend.ru, n))
    for i in 1:n
        backend.rl[i], backend.ru[i] = _backend_constraint_bounds(block.constraints[i])
    end
    status = _HiGHSCAPI.Highs_changeRowsBoundsByRange(
        backend, Cint(offset), Cint(offset + n - 1), backend.rl, backend.ru
    )
    _highs_check_status(status, "HiGHS constraint block bounds update")
    return nothing
end

function update_constraint_variable_coefficients!(
    backend::HiGHSBackend, row::Int, constraint::_BackendConstraint
)
    for (column, coefficient) in zip(constraint.variable_indices, constraint.variable_coefficients)
        status = _HiGHSCAPI.Highs_changeCoeff(backend, Cint(row - 1), Cint(column - 1), coefficient)
        _highs_check_status(status, "HiGHS coefficient update")
    end
    return nothing
end

function backend_update_objective_coefficients!(
    backend::HiGHSBackend, coefficients::_BackendObjectiveCoefficients
)
    n = length(coefficients.variable_indices)
    length(backend.cint_indices) < n && resize!(backend.cint_indices, n)
    for i in 1:n
        backend.cint_indices[i] = Cint(coefficients.variable_indices[i] - 1)
    end
    status = _HiGHSCAPI.Highs_changeColsCostBySet(
        backend, Cint(n), backend.cint_indices, coefficients.variable_coefficients
    )
    _highs_check_status(status, "HiGHS objective update")
    return nothing
end

function backend_set_objective_sense!(backend::HiGHSBackend, sense::Symbol)
    hsense = sense === :max ? _HiGHSCAPI.kHighsObjSenseMaximize : _HiGHSCAPI.kHighsObjSenseMinimize
    status = _HiGHSCAPI.Highs_changeObjectiveSense(backend, hsense)
    _highs_check_status(status, "HiGHS objective sense update")
    return nothing
end

function backend_solve!(
    backend::HiGHSBackend,
    variable_blocks::Vector{_BackendVariableBlock},
    variable_offsets::Vector{Int},
    constraint_blocks::Vector{_BackendConstraintBlock},
    constraint_offsets::Vector{Int},
)
    run_status = _HiGHSCAPI.Highs_run(backend)
    model_status = _HiGHSCAPI.Highs_getModelStatus(backend)
    if run_status != _HiGHSCAPI.kHighsStatusOk &&
        model_status == _HiGHSCAPI.kHighsModelStatusUnknown
        _HiGHSCAPI.Highs_clearSolver(backend)
        run_status = _HiGHSCAPI.Highs_run(backend)
        model_status = _HiGHSCAPI.Highs_getModelStatus(backend)
    end
    if run_status != _HiGHSCAPI.kHighsStatusOk ||
        model_status != _HiGHSCAPI.kHighsModelStatusOptimal
        if get(task_local_storage(), :reactivemodel_diag_highs, false)
            lp_path = joinpath(tempdir(), "reactivemodel-highs-$(getpid())-$(time_ns()).lp")
            write_status = _HiGHSCAPI.Highs_writeModel(backend, lp_path)
            @warn "HiGHS nonoptimal" run_status model_status write_status lp_path columns=length(
                backend.col_value
            ) rows=length(backend.row_value)
        end
        backend.status =
            model_status == _HiGHSCAPI.kHighsModelStatusInfeasible ? :infeasible : :other_error
        backend.objective_value = NaN
        return false, NaN
    end
    solution_status = _HiGHSCAPI.Highs_getSolution(
        backend, backend.col_value, backend.col_dual, backend.row_value, backend.row_dual
    )
    _highs_check_status(solution_status, "HiGHS solution extraction")
    backend.is_mip && (fill!(backend.col_dual, NaN); fill!(backend.row_dual, NaN))
    for block in variable_blocks
        for i in eachindex(block.values)
            block.values[i] = backend.col_value[block.offset + i]
        end
    end
    for block in constraint_blocks
        for i in eachindex(block.values)
            block.values[i] = backend.row_dual[block.offset + i]
        end
    end
    backend.status = :optimal
    backend.objective_value = _HiGHSCAPI.Highs_getObjectiveValue(backend)
    return true, backend.objective_value
end

function backend_write_to_file(backend::HiGHSBackend, filename)
    return _HiGHSCAPI.Highs_writeModel(backend, filename)
end
