module ReactiveModelXpressExt

using ReactiveModel
import ReactiveModel: XpressBackend,
    _BackendVariableBlock, _BackendConstraintBlock, _BackendConstraint,
    _BackendObjective, _BackendObjectiveCoefficients, _backend_constraint_bounds,
    _backend_arrays, _backend_row_data, _backend_column_types, _xpress_status,
    backend_has_duals, backend_set_verbose_flag!, backend_set_objective_sense!,
    backend_write_to_file, initialize!, backend_change_col_bounds!,
    backend_change_cols_bounds!, update_variable_block_bounds!,
    update_constraint_rhs!, update_constraint_block_rhs!,
    update_constraint_variable_coefficients!, backend_update_objective_coefficients!,
    backend_solve!
using Xpress

function _check(code, operation::Symbol)
    code == 0 || error("Xpress $operation failed with code $code")
    return nothing
end

function ReactiveModel.XpressBackend(; silent::Bool=true, time_limit=nothing, mip_rel_gap=nothing, kwargs...)
    for (key, value) in kwargs
        value isa Real || throw(ArgumentError("XpressBackend option $key must be numeric"))
    end
    if time_limit !== nothing
        time_limit isa Real && isfinite(time_limit) && time_limit >= 0 ||
            throw(ArgumentError("XpressBackend time_limit must be finite and nonnegative"))
    end
    if mip_rel_gap !== nothing
        mip_rel_gap isa Real && isfinite(mip_rel_gap) && mip_rel_gap >= 0 ||
            throw(ArgumentError("XpressBackend mip_rel_gap must be finite and nonnegative"))
    end
    opts = Dict{Symbol,Any}(kwargs...)
    opts[:silent] = silent
    opts[:time_limit] = time_limit
    opts[:mip_rel_gap] = mip_rel_gap
    prob_ref = Ref{Xpress.XPRSprob}()
    _check(Xpress.XPRScreateprob(prob_ref), :XPRScreateprob)
    prob = prob_ref[]
    try
        if silent
            _check(Xpress.XPRSsetlogfile(prob, ""), :XPRSsetlogfile)
        end
        _check(Xpress.XPRSsetintcontrol(prob, Xpress.XPRS_OUTPUTLOG,
            silent ? 0 : 1), :XPRSsetintcontrol)
        if time_limit !== nothing
            # MAXTIME < 0 stops LP or MIP search unconditionally; 0 is unlimited.
            seconds = time_limit == 0 ? 0 : -max(1, ceil(Int, time_limit))
            _check(Xpress.XPRSsetintcontrol(prob, Xpress.XPRS_MAXTIME,
                seconds), :XPRSsetintcontrol)
        end
        if mip_rel_gap !== nothing
            _check(Xpress.XPRSsetdblcontrol(prob, Xpress.XPRS_MIPRELSTOP,
                Float64(mip_rel_gap)), :XPRSsetdblcontrol)
        end
        for (key, value) in kwargs
            control = getproperty(Xpress, Symbol("XPRS_", uppercase(string(key))))
            if value isa Integer
                _check(Xpress.XPRSsetintcontrol(prob, control, Int(value)), :XPRSsetintcontrol)
            else
                _check(Xpress.XPRSsetdblcontrol(prob, control, Float64(value)), :XPRSsetdblcontrol)
            end
        end
    catch
        Xpress.XPRSdestroyprob(prob)
        rethrow()
    end
    backend = XpressBackend(
        prob, Float64[], Float64[], Float64[], Float64[],
        Cint[], Cint[], Cdouble[], Cdouble[], Cdouble[], Cdouble[],
        Cdouble[], Cdouble[], UInt8[], Cdouble[], UInt8[], Cint[],
        false, NaN, :not_solved, opts,
    )
    finalizer(backend) do b
        if b.prob != C_NULL
            Xpress.XPRSdestroyprob(b.prob)
            b.prob = C_NULL
        end
    end
    return backend
end

ReactiveModel._similar_backend(backend::XpressBackend) =
    ReactiveModel.XpressBackend(; backend.options...)

backend_has_duals(backend::XpressBackend) = !backend.is_mip

function backend_set_verbose_flag!(backend::XpressBackend; verbose::Bool)
    _check(Xpress.XPRSsetintcontrol(backend.prob, Xpress.XPRS_OUTPUTLOG,
        verbose ? 1 : 0), :XPRSsetintcontrol)
    return nothing
end

function backend_set_objective_sense!(backend::XpressBackend, sense::Symbol)
    obj_sense = sense === :max ? Xpress.XPRS_OBJ_MAXIMIZE : Xpress.XPRS_OBJ_MINIMIZE
    _check(Xpress.XPRSchgobjsense(backend.prob, obj_sense), :XPRSchgobjsense)
    return nothing
end

function backend_write_to_file(backend::XpressBackend, filename::String)
    _check(Xpress.XPRSwriteprob(backend.prob, filename, ""), :XPRSwriteprob)
    return nothing
end

function initialize!(backend::XpressBackend,
                     variable_blocks::Vector{_BackendVariableBlock},
                     constraint_blocks::Vector{_BackendConstraintBlock},
                     objective::_BackendObjective)
    prob = backend.prob::Xpress.XPRSprob
    arrays = _backend_arrays(variable_blocks, constraint_blocks, objective;
                             drop_tolerance=1e-12)
    rowtype, rhs = _backend_row_data(constraint_blocks, UInt8)
    coltype = _backend_column_types(variable_blocks, UInt8)
    int_indices = Cint[i - 1 for i in eachindex(coltype) if coltype[i] != UInt8('C')]
    if arrays.is_mip
        ent_types = UInt8[coltype[i + 1] for i in int_indices]
        _check(Xpress.XPRSloadmip(
            prob, "ReactiveModelXpress", Cint(arrays.ncols), Cint(arrays.nrows),
            rowtype, rhs, C_NULL, arrays.cost, arrays.starts, C_NULL,
            arrays.indices, arrays.values, arrays.lower, arrays.upper,
            Cint(length(int_indices)), Cint(0), ent_types, int_indices, C_NULL,
            C_NULL, C_NULL, C_NULL, C_NULL,
        ), :XPRSloadmip)
    else
        _check(Xpress.XPRSloadlp(
            prob, "ReactiveModelXpress", Cint(arrays.ncols), Cint(arrays.nrows),
            rowtype, rhs, C_NULL, arrays.cost, arrays.starts, C_NULL,
            arrays.indices, arrays.values, arrays.lower, arrays.upper,
        ), :XPRSloadlp)
    end
    backend_set_objective_sense!(backend, objective.op)
    backend.col_value = zeros(arrays.ncols)
    backend.col_dual = zeros(arrays.ncols)
    backend.row_value = zeros(arrays.nrows)
    backend.row_dual = zeros(arrays.nrows)
    backend.a_start = arrays.starts
    backend.a_index = arrays.indices
    backend.a_value = arrays.values
    backend.col_cost = arrays.cost
    backend.col_lower = arrays.lower
    backend.col_upper = arrays.upper
    backend.row_lower = arrays.row_lower
    backend.row_upper = arrays.row_upper
    backend.rowtype = rowtype
    backend.rhs = rhs
    backend.integrality = coltype
    backend.int_indices = int_indices
    backend.is_mip = arrays.is_mip
    backend.objective_value = NaN
    backend.status = :not_solved
    return nothing
end

function backend_change_col_bounds!(backend::XpressBackend, index::Int,
                                    lower::Float64, upper::Float64)
    col = Cint[index - 1]
    _check(Xpress.XPRSchgbounds(backend.prob, Cint(1), col,
        UInt8['L'], Cdouble[lower]), :XPRSchgbounds)
    _check(Xpress.XPRSchgbounds(backend.prob, Cint(1), col,
        UInt8['U'], Cdouble[upper]), :XPRSchgbounds)
    backend.col_lower[index] = lower
    backend.col_upper[index] = upper
    return nothing
end

function backend_change_cols_bounds!(backend::XpressBackend, indices::Vector{Int},
                                     lowers::Vector{Float64}, uppers::Vector{Float64})
    columns = Cint[index - 1 for index in indices]
    n = length(columns)
    _check(Xpress.XPRSchgbounds(backend.prob, Cint(n), columns,
        fill(UInt8('L'), n), lowers), :XPRSchgbounds)
    _check(Xpress.XPRSchgbounds(backend.prob, Cint(n), columns,
        fill(UInt8('U'), n), uppers), :XPRSchgbounds)
    for k in eachindex(indices)
        backend.col_lower[indices[k]] = lowers[k]
        backend.col_upper[indices[k]] = uppers[k]
    end
    return nothing
end

function update_variable_block_bounds!(backend::XpressBackend,
                                      block::_BackendVariableBlock, offset::Int)
    indices = collect(offset + 1:offset + length(block.lower_bounds))
    backend_change_cols_bounds!(backend, indices, block.lower_bounds, block.upper_bounds)
    return nothing
end

function update_constraint_rhs!(backend::XpressBackend,
                               constraint::_BackendConstraint, row::Int)
    lower, upper = _backend_constraint_bounds(constraint)
    value = backend.rowtype[row] == UInt8('G') ? lower : upper
    _check(Xpress.XPRSchgrhs(backend.prob, Cint(1), Cint[row - 1],
        Cdouble[value]), :XPRSchgrhs)
    backend.rhs[row] = value
    backend.row_lower[row] = lower
    backend.row_upper[row] = upper
    return nothing
end

function update_constraint_block_rhs!(backend::XpressBackend,
                                     block::_BackendConstraintBlock, offset::Int)
    n = length(block.constraints)
    rows = Cint[offset + i - 1 for i in 1:n]
    rhs = Float64[c.rhs for c in block.constraints]
    _check(Xpress.XPRSchgrhs(backend.prob, Cint(n), rows, rhs), :XPRSchgrhs)
    for (i, constraint) in enumerate(block.constraints)
        row = offset + i
        backend.rhs[row] = rhs[i]
        backend.row_lower[row], backend.row_upper[row] = _backend_constraint_bounds(constraint)
    end
    return nothing
end

function update_constraint_variable_coefficients!(backend::XpressBackend,
                                                  row::Int,
                                                  constraint::_BackendConstraint)
    for (column, coefficient) in zip(constraint.variable_indices, constraint.variable_coefficients)
        _check(Xpress.XPRSchgcoef(backend.prob, Cint(row - 1), Cint(column - 1),
            coefficient), :XPRSchgcoef)
    end
    return nothing
end

function backend_update_objective_coefficients!(backend::XpressBackend,
                                                coefficients::_BackendObjectiveCoefficients)
    columns = Cint[column - 1 for column in coefficients.variable_indices]
    values = coefficients.variable_coefficients
    _check(Xpress.XPRSchgobj(backend.prob, Cint(length(columns)), columns,
        values), :XPRSchgobj)
    for (column, value) in zip(coefficients.variable_indices, values)
        backend.col_cost[column] = value
    end
    return nothing
end

function backend_solve!(backend::XpressBackend,
                        variable_blocks::Vector{_BackendVariableBlock},
                        variable_offsets::Vector{Int},
                        constraint_blocks::Vector{_BackendConstraintBlock},
                        constraint_offsets::Vector{Int})
    prob = backend.prob::Xpress.XPRSprob
    backend.status = :not_solved
    backend.objective_value = NaN
    if backend.is_mip
        _check(Xpress.XPRSmipoptimize(prob, ""), :XPRSmipoptimize)
    else
        _check(Xpress.XPRSlpoptimize(prob, ""), :XPRSlpoptimize)
    end
    raw_status = Ref{Cint}(0)
    stop = Ref{Cint}(0)
    _check(Xpress.XPRSgetintattrib(prob, backend.is_mip ? Xpress.XPRS_MIPSTATUS : Xpress.XPRS_LPSTATUS,
        raw_status), :XPRSgetintattrib)
    _check(Xpress.XPRSgetintattrib(prob, Xpress.XPRS_STOPSTATUS, stop), :XPRSgetintattrib)
    status = _xpress_status(backend.is_mip, raw_status[], stop[])
    if status !== :optimal
        backend.status = status
        return false, NaN
    end
    obj = Ref{Cdouble}(NaN)
    if backend.is_mip
        _check(Xpress.XPRSgetdblattrib(prob, Xpress.XPRS_MIPOBJVAL, obj), :XPRSgetdblattrib)
        _check(Xpress.XPRSgetmipsol(prob, backend.col_value, backend.row_value), :XPRSgetmipsol)
        fill!(backend.col_dual, NaN)
        fill!(backend.row_dual, NaN)
    else
        _check(Xpress.XPRSgetdblattrib(prob, Xpress.XPRS_LPOBJVAL, obj), :XPRSgetdblattrib)
        _check(Xpress.XPRSgetlpsol(prob, backend.col_value, backend.row_value,
            backend.row_dual, backend.col_dual), :XPRSgetlpsol)
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

end # module ReactiveModelXpressExt
