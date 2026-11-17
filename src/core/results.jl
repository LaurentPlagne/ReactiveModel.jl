struct SolutionVariableBlock{K}
    name::Symbol
    keys::Vector{K}
    positions::Dict{K,Int}
    offset::Int
    values::Vector{Float64}
end

struct SolutionConstraintBlock{K}
    name::Symbol
    keys::Vector{K}
    positions::Dict{K,Int}
    offset::Int
    duals::Vector{Float64}
    activities::Vector{Float64}
end

struct _SolutionGroup{B}
    blocks::Vector{B}
end
_SolutionGroup(::Type{B}) where {B} = _SolutionGroup(B[])

"""Result families grouped by key type; names and family counts are runtime data."""
struct SolutionCollection{G<:Tuple}
    groups::G
    names::Vector{Symbol}
    positions::Dict{Symbol,Tuple{Int,Int}}
end

"""Named family view into an owned snapshot collection."""
struct _ResultFamily{B}
    block::B
    _ResultFamily{B}(block) where {B} = new{B}(block::B)
end

@generated function _result_family_type(::Type{G}) where {G}
    B = Union{(group.parameters[1] for group in G.parameters)...}
    return :(_ResultFamily{$B})
end

struct _ResultBlocks{C,T} <: AbstractVector{T}
    collection::C
end
_ResultBlocks(collection::SolutionCollection{G}) where {G} =
    _ResultBlocks{typeof(collection),_result_family_type(G)}(collection)
Base.IndexStyle(::Type{<:_ResultBlocks}) = IndexLinear()
Base.size(view::_ResultBlocks) = (length(getfield(getfield(view, :collection), :names)),)
function Base.getindex(view::_ResultBlocks, index::Int)
    collection = getfield(view, :collection)
    return collection[getfield(collection, :names)[index]]
end

_result_block(::Tuple{}, group, index) = throw(BoundsError())
@inline function _result_block(groups::Tuple, group::Int, index::Int)
    group == 1 && return first(groups).blocks[index]
    return _result_block(Base.tail(groups), group-1, index)
end

# Inline the lookup and property bridge together so scalar reads can eliminate
# temporary union-backed family views on Julia 1.12 as well as newer compilers.
@inline function Base.getindex(collection::SolutionCollection{G}, name::Symbol) where {G}
    group, index = getfield(collection, :positions)[name]
    block = _result_block(getfield(collection, :groups), group, index)
    return _result_family_type(G)(block)
end
Base.keys(collection::SolutionCollection) = copy(getfield(collection, :names))
Base.length(collection::SolutionCollection) = length(getfield(collection, :names))
Base.haskey(collection::SolutionCollection, name::Symbol) = haskey(getfield(collection, :positions), name)
Base.getindex(collection::SolutionCollection, index::Int) = collection[getfield(collection, :names)[index]]
function Base.iterate(collection::SolutionCollection, state::Int=1)
    state > length(collection) && return nothing
    return collection[state], state+1
end
Base.values(collection::SolutionCollection) = _ResultBlocks(collection)
Base.pairs(collection::SolutionCollection) = (name=>collection[name] for name in getfield(collection, :names))
@inline function Base.getproperty(collection::SolutionCollection, name::Symbol)
    name === :blocks && return _ResultBlocks(collection)
    name === :by_name && return collection
    return collection[name]
end
Base.propertynames(collection::SolutionCollection, private::Bool=false) =
    (Tuple(getfield(collection, :names))..., :blocks, :by_name)
function Base.show(io::IO, collection::SolutionCollection)
    print(io, "SolutionCollection(")
    show(io, getfield(collection, :names))
    print(io, ")")
end

@generated function Base.propertynames(::_ResultFamily{B}, private::Bool=false) where {B}
    return QuoteNode(fieldnames(first(Base.uniontypes(B))))
end

# A retained family owns only its block, rather than keeping every other
# family's numerical snapshot alive. Explicit branches preserve scalar inference
# even for applications with more key types than Julia's union-splitting limit.
@generated function _result_field(family::_ResultFamily{B}, ::Val{F}) where {B,F}
    cases = [:(block isa $T && return getfield(block::$T, $(QuoteNode(F)))) for T in Base.uniontypes(B)]
    return quote
        block = getfield(family, :block)
        $(cases...)
        throw(ArgumentError("invalid result-family storage"))
    end
end
@generated function Base.show(io::IO, family::_ResultFamily{B}) where {B}
    cases = [:(block isa $T && return show(io, block::$T)) for T in Base.uniontypes(B)]
    return quote
        block = getfield(family, :block)
        $(cases...)
        throw(ArgumentError("invalid result-family storage"))
    end
end

@inline function Base.getproperty(family::_ResultFamily, name::Symbol)
    name === :name && return _result_field(family, Val(:name))
    name === :keys && return _result_field(family, Val(:keys))
    name === :positions && return _result_field(family, Val(:positions))
    name === :offset && return _result_field(family, Val(:offset))
    name === :values && return _result_field(family, Val(:values))
    name === :duals && return _result_field(family, Val(:duals))
    name === :activities && return _result_field(family, Val(:activities))
    throw(ArgumentError("unknown result-family property: $name"))
end

_result_block_value(block::SolutionVariableBlock, key) = block.values[block.positions[key]]
_result_block_value(block::SolutionConstraintBlock, key) = block.duals[block.positions[key]]
@generated function Base.getindex(family::_ResultFamily{B}, key) where {B}
    cases = [:(block isa $T && return _result_block_value(block::$T, key)) for T in Base.uniontypes(B)]
    return quote
        block = getfield(family, :block)
        $(cases...)
        throw(ArgumentError("invalid result-family storage"))
    end
end

# These factories inspect only field types, never callback/study values. Their
# output has one group per key type, independent of family names and counts.
function _solution_group_type(F, kind)
    families = fieldtypes(fieldtype(fieldtype(F, kind), :blocks))
    key_types = Type[]
    for family in families
        K = eltype(fieldtype(family, :keys))
        K in key_types || push!(key_types, K)
    end
    block_type = kind === :variables ? SolutionVariableBlock : SolutionConstraintBlock
    return Tuple{(_SolutionGroup{block_type{K}} for K in key_types)...}
end
@generated function _solution_collection_type(::Type{F}, ::Val{K}) where {F,K}
    return :(SolutionCollection{$(_solution_group_type(F, K))})
end
@generated function _solution_groups(::Type{F}, ::Val{K}) where {F,K}
    G = _solution_group_type(F, K)
    return Expr(:tuple, [:(_SolutionGroup($(group.parameters[1]))) for group in G.parameters]...)
end

_empty_result_block(block::VariableBlock, offset) =
    SolutionVariableBlock(block.name, block.keys, block.positions, offset, Float64[])
_empty_result_block(block::ConstraintBlock, offset) =
    SolutionConstraintBlock(block.name, block.keys, block.positions, offset, Float64[], Float64[])

@inline function _append_result_block!(
    groups::Tuple{_SolutionGroup{B},Vararg}, collection, block::B, group::Int
) where {B}
    push!(first(groups).blocks, block)
    push!(getfield(collection, :names), block.name)
    getfield(collection, :positions)[block.name] = (group, length(first(groups).blocks))
    return nothing
end
@inline function _append_result_block!(groups::Tuple, collection, block, group::Int)
    return _append_result_block!(Base.tail(groups), collection, block, group+1)
end

struct _ResultLayoutBuilder{C}
    collection::C
    offset::Base.RefValue{Int}
end
function (builder::_ResultLayoutBuilder)(family)
    block = _empty_result_block(family, builder.offset[])
    _append_result_block!(getfield(builder.collection, :groups), builder.collection, block, 1)
    builder.offset[] += length(family.keys)
    return nothing
end
function _empty_solution_collection(definition::F, ::Val{K}) where {F,K}
    collection = SolutionCollection(_solution_groups(F, Val(K)), Symbol[], Dict{Symbol,Tuple{Int,Int}}())
    _foreach_family(_ResultLayoutBuilder(collection, Ref(0)), getfield(definition, K).blocks)
    return collection
end
_empty_solution_variables(definition::ProblemDefinition) = _empty_solution_collection(definition, Val(:variables))
_empty_solution_constraints(definition::ProblemDefinition) = _empty_solution_collection(definition, Val(:constraints))

"""
Result metadata independent of mathematical formulas and model identity.
The empty collections retain names, typed keys and positions; numerical result
arrays are materialized separately so retained collections remain snapshots.
"""
struct _ResultLayout{V,C}
    variables::V
    constraints::C
end

_ResultLayout(definition::ProblemDefinition) = _ResultLayout(
    _empty_solution_variables(definition), _empty_solution_constraints(definition)
)

"""
Reusable view of the current solve. The workspace updates this same object on
each numerical solve and clears its lazy variable and constraint collections.
A collection retained separately by the caller remains a snapshot, but a
retained `Solution` follows the model's latest successful solve.
Only result metadata is retained here, not the mathematical definition or its
callbacks. Models with the same result layout can share this solution type.
"""
mutable struct Solution{B,V,C,I}
    layout::_ResultLayout{V,C}
    backend::B
    variable_indices::_VariableIndices{I}
    row_coefficients::Vector{Dict{Int,Float64}}
    generation_token::Base.RefValue{UInt}
    generation::UInt
    objective_value::Float64
    termination_status::Symbol
    variables_cache::Union{Nothing,V}
    constraints_cache::Union{Nothing,C}
end

function Base.getproperty(solution::Solution, name::Symbol)
    name === :variables && return _variables(solution)
    name === :constraints && return _constraints(solution)
    name === :objective_value && return objective_value(solution)
    name === :termination_status && return termination_status(solution)
    return getfield(solution, name)
end

function Base.propertynames(::Solution, private::Bool=false)
    return (:variables, :constraints, :objective_value, :termination_status)
end

function _check_current(solution::Solution)
    getfield(solution, :generation) == getfield(solution, :generation_token)[] ||
        throw(ErrorException("solution is no longer current"))
    return nothing
end

function _variables(solution::Solution)
    _check_current(solution)
    cache = getfield(solution, :variables_cache)
    cache === nothing || return cache
    cache = _materialize_variables(solution)
    setfield!(solution, :variables_cache, cache)
    return cache
end

function _constraints(solution::Solution)
    _check_current(solution)
    cache = getfield(solution, :constraints_cache)
    cache === nothing || return cache
    cache = _materialize_constraints(solution)
    setfield!(solution, :constraints_cache, cache)
    return cache
end

"""Return the objective value of the current solution."""
objective_value(s::Solution) = (_check_current(s); getfield(s, :objective_value))
"""Return the termination status of the current solution."""
termination_status(s::Solution) = (_check_current(s); getfield(s, :termination_status))
"""Return `true` when the solution was produced by an optimal backend solve."""
is_optimal(s::Solution) = termination_status(s) === :optimal

_result_position(::Tuple{}, group, index, key) = throw(BoundsError())
@inline function _result_position(groups::Tuple, group::Int, index::Int, key)
    if group == 1
        block = first(groups).blocks[index]
        actual_key = key === nothing && length(block.keys) == 1 ? first(block.keys) : key
        return block.offset + block.positions[actual_key]
    end
    return _result_position(Base.tail(groups), group-1, index, key)
end
function _family_position(collection::SolutionCollection, name::Symbol, key)
    group, index = getfield(collection, :positions)[name]
    return _result_position(getfield(collection, :groups), group, index, key)
end

"""Return one primal value directly, without materializing result collections."""
function primal(s::Solution, block::Symbol, key=nothing)
    _check_current(s)
    column = _column_position(getfield(s, :variable_indices), block, key; normalize=true)
    return backend_primal_value(getfield(s, :backend), column)
end

function _constraint_position(s::Solution, block::Symbol, key)
    return _family_position(getfield(s, :layout).constraints, block, key)
end
_dual_current(s::Solution, row::Int) = backend_dual_value(getfield(s, :backend), row)

"""Return one dual value without materializing result collections."""
function dual(s::Solution, block::Symbol, key=nothing)
    _check_current(s)
    return _dual_current(s, _constraint_position(s, block, key))
end

"""
Return the normalized variable expression of one constraint without materializing
collections. Constants moved into the RHS are excluded: for `x + 1 == demand`,
the activity is `x`, consistent with the normalized RHS `demand - 1`.
"""
function activity(s::Solution, block::Symbol, key=nothing)
    _check_current(s)
    row = _constraint_position(s, block, key)
    return sum(
        coefficient * backend_primal_value(getfield(s, :backend), column) for
        (column, coefficient) in getfield(s, :row_coefficients)[row];
        init=0.0,
    )
end
