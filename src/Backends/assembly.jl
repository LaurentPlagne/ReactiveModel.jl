"""
    _backend_data(variable_blocks, constraint_blocks, objective;
                  drop_tolerance=nothing, protect_tiny=false)

Build the common Float64 column/row arrays and a one-based Julia CSC matrix.
Allocation happens during model construction,
not during repeated solves. Column order follows block offsets; row senses are
normalized to lower/upper bounds. Binary bounds have already been restricted
by the mathematical assembly.

With `drop_tolerance=nothing`, retain all entries, including explicit zeros.
A numerical tolerance discards entries with `abs(value) <= drop_tolerance`
before sparse assembly, matching the optional solvers' import policies.
With `protect_tiny=true`, return the original tiny coefficients separately and
import placeholders in the matrix. HiGHS restores them through its coefficient
API after import. This solver-specific policy is opt-in and does not alter Clp.
"""
function _backend_data(variable_blocks, constraint_blocks, objective;
                       drop_tolerance=nothing, protect_tiny=false)
    ncols=sum(b->length(b.lower_bounds), variable_blocks; init=0)
    nrows=sum(b->length(b.constraints), constraint_blocks; init=0)
    cost=zeros(Float64, ncols)
    for (column, value) in
        zip(objective.coefficients.variable_indices, objective.coefficients.variable_coefficients)
        cost[column]=value
    end
    lower=zeros(Float64, ncols)
    upper=zeros(Float64, ncols)
    for block in variable_blocks, i in eachindex(block.lower_bounds)
        column=block.offset+i
        lower[column]=block.lower_bounds[i]
        upper[column]=block.upper_bounds[i]
    end
    row_lower=zeros(Float64, nrows)
    row_upper=zeros(Float64, nrows)
    rows, columns, values=Int[], Int[], Float64[]
    tiny=Tuple{Int,Int,Float64}[]
    for block in constraint_blocks, (i, constraint) in enumerate(block.constraints)
        row=block.offset+i
        row_lower[row], row_upper[row]=_backend_constraint_bounds(constraint)
        for (column, value) in zip(constraint.variable_indices, constraint.variable_coefficients)
            drop_tolerance !== nothing && abs(value) <= drop_tolerance && continue
            if protect_tiny && abs(value)<1e-9
                push!(tiny, (row, column, value))
                value=1.0
            end
            push!(rows, row)
            push!(columns, column)
            push!(values, value)
        end
    end
    matrix=SparseArrays.sparse(rows, columns, values, nrows, ncols)
    return (;
        ncols,
        nrows,
        cost,
        lower,
        upper,
        row_lower,
        row_upper,
        is_mip=any(b -> b.domain !== :continuous && !isempty(b.lower_bounds), variable_blocks),
        matrix,
        tiny,
    )
end

# Native solvers consume the same mathematical data with zero-based C indices.
function _backend_arrays(variable_blocks, constraint_blocks, objective; kwargs...)
    data = _backend_data(variable_blocks, constraint_blocks, objective; kwargs...)
    integrality = zeros(Cint, data.ncols)
    for block in variable_blocks, i in eachindex(block.lower_bounds)
        integrality[block.offset + i] = block.domain === :continuous ? Cint(0) : Cint(1)
    end
    return (; data...,
        integrality,
        starts=Cint.(data.matrix.colptr .- 1),
        indices=Cint.(data.matrix.rowval .- 1),
        values=data.matrix.nzval,
    )
end

# Commercial APIs encode row senses and variable domains as character arrays.
# Use the declared sense directly; large finite right-hand sides are not infinity.
function _backend_row_data(constraint_blocks, ::Type{T}) where {T<:Integer}
    nrows = sum(b -> length(b.constraints), constraint_blocks; init=0)
    sense = Vector{T}(undef, nrows)
    rhs = Vector{Float64}(undef, nrows)
    for block in constraint_blocks, (i, constraint) in enumerate(block.constraints)
        row = block.offset + i
        sense[row] = T(constraint.op === :eq ? 'E' : constraint.op === :le ? 'L' : 'G')
        rhs[row] = constraint.rhs
    end
    return sense, rhs
end

function _backend_column_types(variable_blocks, ::Type{T}) where {T<:Integer}
    ncols = sum(b -> length(b.lower_bounds), variable_blocks; init=0)
    types = Vector{T}(undef, ncols)
    for block in variable_blocks
        code = T(block.domain === :binary ? 'B' : block.domain === :integer ? 'I' : 'C')
        for i in eachindex(block.lower_bounds)
            types[block.offset + i] = code
        end
    end
    return types
end
