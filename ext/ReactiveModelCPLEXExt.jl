module ReactiveModelCPLEXExt

using ReactiveModel
import ReactiveModel: CPLEXBackend,
    _BackendVariableBlock, _BackendConstraintBlock, _BackendObjective,
    _backend_arrays, _backend_row_data, _backend_column_types, _cplex_status,
    backend_has_duals, backend_set_verbose_flag!, backend_set_objective_sense!,
    initialize!, backend_solve!
using CPLEX

function _check(code, operation::Symbol)
    code == 0 || error("CPLEX $operation failed with code $code")
    return nothing
end

function ReactiveModel.CPLEXBackend(; silent::Bool=true, time_limit=nothing, kwargs...)
    isempty(kwargs) || throw(ArgumentError(
        "CPLEXBackend supports only silent and time_limit constructor options; " *
        "unsupported options: $(join(keys(kwargs), ", "))",
    ))
    opts = Dict{Symbol,Any}(:silent => silent, :time_limit => time_limit)
    env_stat = Ref{Cint}(0)
    env = CPLEX.CPXopenCPLEX(env_stat)
    _check(env_stat[], :CPXopenCPLEX)
    lp = C_NULL
    try
        lp_stat = Ref{Cint}(0)
        lp = CPLEX.CPXcreateprob(env, lp_stat, "ReactiveModelCPLEX")
        _check(lp_stat[], :CPXcreateprob)
        _check(CPLEX.CPXsetintparam(env, CPLEX.CPXPARAM_ScreenOutput,
            silent ? CPLEX.CPX_OFF : CPLEX.CPX_ON), :CPXsetintparam)
        if time_limit !== nothing
            _check(CPLEX.CPXsetdblparam(env, CPLEX.CPXPARAM_TimeLimit,
                Float64(time_limit)), :CPXsetdblparam)
        end
    catch
        lp == C_NULL || CPLEX.CPXfreeprob(env, Ref(lp))
        CPLEX.CPXcloseCPLEX(Ref(env))
        rethrow()
    end
    backend = CPLEXBackend(
        env, lp, Float64[], Float64[], Float64[], Float64[],
        false, NaN, :not_solved, opts,
    )
    finalizer(backend) do b
        if b.lp != C_NULL
            CPLEX.CPXfreeprob(b.env, Ref(b.lp))
            b.lp = C_NULL
        end
        if b.env != C_NULL
            CPLEX.CPXcloseCPLEX(Ref(b.env))
            b.env = C_NULL
        end
    end
    return backend
end

ReactiveModel._similar_backend(backend::CPLEXBackend) =
    ReactiveModel.CPLEXBackend(; backend.options...)

backend_has_duals(backend::CPLEXBackend) = !backend.is_mip

function backend_set_verbose_flag!(backend::CPLEXBackend; verbose::Bool)
    _check(CPLEX.CPXsetintparam(backend.env, CPLEX.CPXPARAM_ScreenOutput,
        verbose ? CPLEX.CPX_ON : CPLEX.CPX_OFF), :CPXsetintparam)
    return nothing
end

function backend_set_objective_sense!(backend::CPLEXBackend, sense::Symbol)
    objsen = sense === :max ? CPLEX.CPX_MAX : CPLEX.CPX_MIN
    _check(CPLEX.CPXchgobjsen(backend.env, backend.lp, objsen), :CPXchgobjsen)
    return nothing
end

function initialize!(backend::CPLEXBackend,
                     variable_blocks::Vector{_BackendVariableBlock},
                     constraint_blocks::Vector{_BackendConstraintBlock},
                     objective::_BackendObjective)
    arrays = _backend_arrays(variable_blocks, constraint_blocks, objective;
                             drop_tolerance=1e-12)
    row_sense, rhs = _backend_row_data(constraint_blocks, Cchar)
    coltype = _backend_column_types(variable_blocks, Cchar)
    counts = Cint[arrays.starts[k + 1] - arrays.starts[k] for k in 1:arrays.ncols]
    objsen = objective.op === :max ? CPLEX.CPX_MAX : CPLEX.CPX_MIN
    _check(CPLEX.CPXcopylp(
        backend.env, backend.lp, Cint(arrays.ncols), Cint(arrays.nrows), objsen,
        arrays.cost, rhs, row_sense, arrays.starts, counts, arrays.indices,
        arrays.values, arrays.lower, arrays.upper, C_NULL,
    ), :CPXcopylp)
    if arrays.is_mip
        _check(CPLEX.CPXcopyctype(backend.env, backend.lp, coltype), :CPXcopyctype)
    end
    backend.col_value = zeros(arrays.ncols)
    backend.col_dual = zeros(arrays.ncols)
    backend.row_value = zeros(arrays.nrows)
    backend.row_dual = zeros(arrays.nrows)
    backend.is_mip = arrays.is_mip
    backend.objective_value = NaN
    backend.status = :not_solved
    return nothing
end

function backend_solve!(backend::CPLEXBackend,
                        variable_blocks::Vector{_BackendVariableBlock},
                        variable_offsets::Vector{Int},
                        constraint_blocks::Vector{_BackendConstraintBlock},
                        constraint_offsets::Vector{Int})
    backend.objective_value = NaN
    backend.status = :not_solved
    if backend.is_mip
        _check(CPLEX.CPXmipopt(backend.env, backend.lp), :CPXmipopt)
    else
        _check(CPLEX.CPXlpopt(backend.env, backend.lp), :CPXlpopt)
    end
    status = _cplex_status(backend.is_mip, CPLEX.CPXgetstat(backend.env, backend.lp))
    if status !== :optimal
        backend.status = status
        return false, NaN
    end

    obj = Ref{Cdouble}(NaN)
    _check(CPLEX.CPXgetobjval(backend.env, backend.lp, obj), :CPXgetobjval)
    if !isempty(backend.col_value)
        _check(CPLEX.CPXgetx(backend.env, backend.lp, backend.col_value,
            Cint(0), Cint(length(backend.col_value) - 1)), :CPXgetx)
    end
    if backend.is_mip
        fill!(backend.col_dual, NaN)
        fill!(backend.row_dual, NaN)
    elseif !isempty(backend.row_dual)
        _check(CPLEX.CPXgetpi(backend.env, backend.lp, backend.row_dual,
            Cint(0), Cint(length(backend.row_dual) - 1)), :CPXgetpi)
    end
    backend.objective_value = obj[]
    backend.status = :optimal
    for block in variable_blocks, i in eachindex(block.values)
        block.values[i] = backend.col_value[block.offset + i]
    end
    for block in constraint_blocks, i in eachindex(block.values)
        block.values[i] = backend.row_dual[block.offset + i]
    end
    return true, backend.objective_value
end

end # module ReactiveModelCPLEXExt
