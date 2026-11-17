"""
    SolveError(status, message)

Raised when a backend does not return an optimal solution. `status` is the
backend's compact symbolic status. The numerical model is kept: the next point
is updated from the failed one, as after a successful solve. No previous lazy
solution remains valid, and asking again for the failed point raises again.
"""
struct SolveError <: Exception
    status::Symbol
    message::String
end

function Base.showerror(io::IO, error::SolveError)
    return print(io, error.message, " (termination status: ", error.status, ")")
end

"""
    VarKey(block::Symbol, key)

A variable family name and its complete key. Equality resolves hash collisions
using the key itself.
"""
struct VarKey{K}
    block::Symbol
    key::K
end
Base.isequal(a::VarKey, b::VarKey) = a.block === b.block && isequal(a.key, b.key)
Base.:(==)(a::VarKey, b::VarKey) = isequal(a, b)
Base.hash(k::VarKey, seed::UInt) = hash(k.key, hash(k.block, seed))

"""Column lookup for all families sharing key type `K`, independent of callbacks."""
struct _KeyIndex{K}
    positions::Dict{VarKey{K},Int}
    scalar_keys::Dict{Symbol,K}
end
_KeyIndex(::Type{K}) where {K} = _KeyIndex(Dict{VarKey{K},Int}(), Dict{Symbol,K}())

function _index_family!(index::_KeyIndex{K}, name::Symbol, keys::Vector{K}, offset) where {K}
    for (i, key) in enumerate(keys)
        index.positions[VarKey(name, key)] = offset+i
    end
    length(keys)==1 && (index.scalar_keys[name]=only(keys))
    return index
end

@inline function _populate_index!(
    groups::Tuple{_KeyIndex{K},Vararg}, name, keys::Vector{K}, offset
) where {K}
    _index_family!(first(groups), name, keys, offset)
    return nothing
end
@inline function _populate_index!(groups::Tuple, name, keys, offset)
    _populate_index!(Base.tail(groups), name, keys, offset)
    return nothing
end

struct _PopulateIndices{G<:Tuple}
    groups::G
    offset::Base.RefValue{Int}
end
function (builder::_PopulateIndices)(family)
    _populate_index!(builder.groups, family.name, family.keys, builder.offset[])
    builder.offset[] += length(family.keys)
    return nothing
end

"""
Compact typed column maps. The tuple has one entry per distinct key type,
not per mathematical family. Its type excludes callbacks and model identity,
so every coefficient emitter does not specialize on the complete definition.
"""
struct _VariableIndices{G<:Tuple}
    groups::G
    count::Int
end
@generated function _VariableIndices(definition::F) where {F}
    # Inspect field types only. Family names, key values and dimensions stay
    # runtime values. Fix groups before traversal rather than growing an
    # inferred tuple recursively: the latter loses inference for mixed key types.
    families = fieldtypes(fieldtype(fieldtype(F, :variables), :blocks))
    key_types = Type[]
    for family in families
        key_type = eltype(fieldtype(family, :keys))
        key_type in key_types || push!(key_types, key_type)
    end
    constructors = Expr(:tuple, [:(_KeyIndex($K)) for K in key_types]...)
    return quote
        groups = $constructors
        offset = Ref(0)
        _foreach_family(_PopulateIndices(groups, offset), definition.variables.blocks)
        _VariableIndices(groups, offset[])
    end
end

_lookup_index(::Tuple{}, block, key, normalize) = 0
@noinline function _lookup_index(groups::Tuple, block, key, normalize)
    group=first(groups)
    column=if normalize && key===nothing && haskey(group.scalar_keys, block)
        get(group.positions, VarKey(block, group.scalar_keys[block]), 0)
    else
        get(group.positions, VarKey(block, key), 0)
    end
    return column==0 ? _lookup_index(Base.tail(groups), block, key, normalize) : column
end
function _column_position(indices::_VariableIndices, block, key; normalize=false)
    column=_lookup_index(indices.groups, block, key, normalize)
    column==0 && throw(KeyError((block, key)))
    return column
end
@inline Base.getindex(vi::_VariableIndices, k::VarKey) = _column_position(vi, k.block, k.key)
@inline Base.getindex(vi::_VariableIndices, k::Tuple{Symbol,Any}) = vi[VarKey(k...)]
@inline Base.haskey(vi::_VariableIndices, k::VarKey) =
    _lookup_index(vi.groups, k.block, k.key, false)!=0
@inline Base.haskey(vi::_VariableIndices, k::Tuple{Symbol,Any}) = haskey(vi, VarKey(k...))
Base.length(vi::_VariableIndices) = vi.count

mutable struct _Emitter{F}
    variable_indices::_VariableIndices{F}
    coefficients::Dict{Int,Float64}
    constant::Float64
end
_Emitter(indices::_VariableIndices, coefficients) = _Emitter(indices, coefficients, 0.0)
function _add_constant!(emit::_Emitter, value)
    value = Float64(value)
    isfinite(value) || throw(ArgumentError("non-finite affine constant"))
    emit.constant += value
    return nothing
end

_always_changed(old, new) = true
const _UNNAMED_BLOCK = Symbol("#ReactiveModel unnamed block#")

"""
    AbstractVariable

Supertype for the mathematical definition of one variable family. Extend
both [`lower_bound`](@ref) and [`upper_bound`](@ref) for a concrete subtype.
The methods receive the application payloads stored inside `ConstData` and
`Parameters`. Neither method has a numerical fallback.
"""
abstract type AbstractVariable end

"""
    VariableTag{ModelIdentity, VarName} <: AbstractVariable

Singleton tag type representing a variable family declared in a `@model` definition,
using an identity type emitted in the caller module for each model declaration.
"""
struct VariableTag{M,V} <: AbstractVariable end

"""
    iteration_space(mathematical_object, data)

Return an iterable collection of keys defining the index space of a variable family
or constraint family. The element type `K` of the space must be concrete. Keys need stable equality
and hashing for the lifetime of the definition; they need not be `isbits`.
"""
function iteration_space end

"""`lower_bound(variable, study, point, key)` returns one variable lower bound."""
function lower_bound end

"""`upper_bound(variable, study, point, key)` returns one variable upper bound."""
function upper_bound end

"""`variable_bounds(variable, study, point, key)` returns `(lower, upper)` bounds for one variable key."""
@inline function variable_bounds(variable::AbstractVariable, const_data, parameters, key)
    if applicable(variable_bounds, variable, key)
        return variable_bounds(variable, key)
    elseif applicable(variable_bounds, variable)
        return variable_bounds(variable)
    else
        return (
            Float64(lower_bound(variable, const_data, parameters, key)),
            Float64(upper_bound(variable, const_data, parameters, key)),
        )
    end
end
"""
    has_changed(mathematical_object, previous, current)

Return whether the numerical values defined by a variable family, constraint
family, or objective may differ between two parameter points. The conservative
default is `true`: ReactiveModel then reevaluates every number in that family,
compares it with the current numerical model, and forwards only actual changes
to the backend. Return `false` only when the complete family is known to be
independent of the parameter change.
"""
has_changed(::AbstractVariable, old_parameters, new_parameters) = true

_formula_value(value, const_data, parameters, key) = value
function _formula_value(formula::F, const_data, parameters, key) where {F<:Function}
    return formula(const_data, parameters, key)
end

struct _VariableFormula{L,U,H} <: AbstractVariable
    lower::L
    upper::U
    predicate::H
end

function lower_bound(variable::_VariableFormula, const_data, parameters, key)
    return _formula_value(variable.lower, const_data, parameters, key)
end
function upper_bound(variable::_VariableFormula, const_data, parameters, key)
    return _formula_value(variable.upper, const_data, parameters, key)
end
@inline variable_bounds(variable::_VariableFormula, const_data, parameters, key) = (
    Float64(lower_bound(variable, const_data, parameters, key)),
    Float64(upper_bound(variable, const_data, parameters, key)),
)
function has_changed(variable::_VariableFormula, old_parameters, new_parameters)
    return variable.predicate(old_parameters, new_parameters)
end

"""
    VariableBlock(variable, index_space=(nothing,); domain=:continuous)
    VariableBlock(name, index_space, variable; domain=:continuous)
    VariableBlock(name, index_space; lower=-Inf, upper=Inf,
                  domain=:continuous, has_changed=(previous, current) -> true)

An indexed family of variables. `variable` is either a concrete
[`AbstractVariable`](@ref) value or its singleton type. In the compact form,
the surrounding [`ProblemDefinition`](@ref) takes the family name from its
`NamedTuple` key. The explicit `name` form remains available and is checked
against that key. All keys in `index_space` must have one concrete type. A
typed empty index space is allowed and represents a structural hole.
Supported domains are `:continuous`, `:integer`, and `:binary`.
In the function-based form, a bound is a number or a function of
`(data, point, key)`. The default predicate reevaluates at distinct points;
there is no syntax-based dependency inference for these callbacks.
"""
mutable struct VariableBlock{K,V<:AbstractVariable}
    name::Symbol
    keys::Vector{K}
    positions::Dict{K,Int}
    variable::V
    domain::Symbol
end

function VariableBlock(
    name::Symbol, index_space, variable::V; domain::Symbol=:continuous
) where {V<:AbstractVariable}
    keys = vec([x for x in index_space])
    K = eltype(keys)
    isconcretetype(K) ||
        throw(ArgumentError("keys of variable block $name must have one concrete type"))
    length(unique(keys)) == length(keys) ||
        throw(ArgumentError("duplicate key in variable block $name"))
    domain ∈ (:continuous, :integer, :binary) ||
        throw(ArgumentError("unknown variable domain $domain"))
    return VariableBlock{K,V}(
        name, keys, Dict{K,Int}(k => i for (i, k) in enumerate(keys)), variable, domain
    )
end

function VariableBlock(name::Symbol, index_space, ::Type{V}; kwargs...) where {V<:AbstractVariable}
    return VariableBlock(name, index_space, V(); kwargs...)
end

function VariableBlock(variable::V, index_space=(nothing,); kwargs...) where {V<:AbstractVariable}
    return VariableBlock(_UNNAMED_BLOCK, index_space, variable; kwargs...)
end

function VariableBlock(::Type{V}, index_space=(nothing,); kwargs...) where {V<:AbstractVariable}
    return VariableBlock(V(), index_space; kwargs...)
end

function VariableBlock(
    name::Symbol,
    index_space;
    lower=(-Inf),
    upper=Inf,
    domain::Symbol=:continuous,
    has_changed=_always_changed,
)
    return VariableBlock(name, index_space, _VariableFormula(lower, upper, has_changed); domain)
end

function VariableBlock(name::Symbol, variable::V, data; kwargs...) where {V<:AbstractVariable}
    return VariableBlock(name, iteration_space(variable, data), variable; kwargs...)
end

function VariableBlock(name::Symbol, variable::AbstractVariable; kwargs...)
    return VariableBlock(name, (nothing,), variable; kwargs...)
end
function VariableBlock(name::Symbol, ::Type{V}; kwargs...) where {V<:AbstractVariable}
    return VariableBlock(name, (nothing,), V(); kwargs...)
end
VariableBlock(name::Symbol; kwargs...) = VariableBlock(name, (nothing,); kwargs...)

"""
    AbstractConstraint

Supertype for one mathematical family of linear constraints. Extend
[`lhs_terms!`](@ref) and [`rhs_value`](@ref) for a concrete subtype.
The methods receive the application payloads, and neither has a numerical
fallback.
"""
abstract type AbstractConstraint end

"""
    ConstraintTag{ModelIdentity, CtrName} <: AbstractConstraint

Singleton tag type representing a constraint family declared in a `@model` definition,
using the unique identity type of its model declaration.
"""
struct ConstraintTag{M,C} <: AbstractConstraint end

"""Emit the variable terms of one indexed constraint into `equation`."""
function lhs_terms! end
"""
    rhs_value(constraint, study_data, point, key) -> Float64

Return the complete normalized RHS for one row. A custom constraint must define
this method, or a compact `rhs_value(constraint, key)` / `rhs_value(constraint)`
method for values that do not depend on the supplied study/point. An omitted
method raises `MethodError`; zero must be declared explicitly. Safe mode checks
updates of supplied formulas, not whether an absent formula was intended as zero.
"""
@inline function rhs_value(constraint::AbstractConstraint, const_data, parameters, key)
    if applicable(rhs_value, constraint, key)
        return Float64(rhs_value(constraint, key))
    elseif applicable(rhs_value, constraint)
        return Float64(rhs_value(constraint))
    else
        throw(MethodError(rhs_value, (constraint, const_data, parameters, key)))
    end
end
has_changed(::AbstractConstraint, old_parameters, new_parameters) = true

"""Return whether a constraint family's left-hand side may have changed."""
function lhs_has_changed(constraint::AbstractConstraint, old_parameters, new_parameters)
    return has_changed(constraint, old_parameters, new_parameters)
end

"""Return whether a constraint family's right-hand side may have changed."""
function rhs_has_changed(constraint::AbstractConstraint, old_parameters, new_parameters)
    return has_changed(constraint, old_parameters, new_parameters)
end

"""
    constraint_sense(constraint) -> Symbol

Return the relational operator (`:le`, `:eq`, `:ge`) associated with a constraint family.
Defaults to `:le`.
"""
constraint_sense(::AbstractConstraint) = :le

struct _ConstraintFormula{T,R,H} <: AbstractConstraint
    terms::T
    rhs::R
    predicate::H
end

function lhs_terms!(constraint::_ConstraintFormula, emit, const_data, parameters, key)
    return constraint.terms(emit, const_data, parameters, key)
end
function rhs_value(constraint::_ConstraintFormula, const_data, parameters, key)
    return _formula_value(constraint.rhs, const_data, parameters, key)
end
function has_changed(constraint::_ConstraintFormula, old_parameters, new_parameters)
    return constraint.predicate(old_parameters, new_parameters)
end

"""
    ConstraintBlock(index_space, sense, constraint)
    ConstraintBlock(sense, constraint)
    ConstraintBlock(name, index_space, sense, constraint)
    ConstraintBlock(name, index_space, sense, terms, rhs; has_changed=(previous, current) -> true)

An indexed family of constraints. `sense` is `:eq`, `:le`, or `:ge`, and
`constraint` is either a concrete [`AbstractConstraint`](@ref) value or its
singleton type. In the compact form, the surrounding [`ProblemDefinition`](@ref)
takes the family name from its `NamedTuple` key. Omit `index_space` for one
scalar constraint. A typed empty index space is allowed and represents a
structural hole.
In the function-based form, `terms(emit, data, point, key)` emits the complete
current variable expression with `set_coeff!`; `rhs` is a number or a function
of `(data, point, key)`. The predicate guards both components and defaults to
reevaluation at distinct points.
"""
mutable struct ConstraintBlock{K,C<:AbstractConstraint}
    name::Symbol
    keys::Vector{K}
    positions::Dict{K,Int}
    sense::Symbol
    constraint::C
end

function ConstraintBlock(
    name::Symbol, index_space, sense::Symbol, constraint::C
) where {C<:AbstractConstraint}
    sense ∈ (:eq, :le, :ge) || throw(ArgumentError("constraint sense must be :eq, :le or :ge"))
    keys = vec([x for x in index_space])
    K = eltype(keys)
    isconcretetype(K) ||
        throw(ArgumentError("keys of constraint block $name must have one concrete type"))
    length(unique(keys)) == length(keys) ||
        throw(ArgumentError("duplicate key in constraint block $name"))
    return ConstraintBlock{K,C}(
        name, keys, Dict{K,Int}(k => i for (i, k) in enumerate(keys)), sense, constraint
    )
end

function ConstraintBlock(
    name::Symbol, index_space, sense::Symbol, ::Type{C}
) where {C<:AbstractConstraint}
    return ConstraintBlock(name, index_space, sense, C())
end

function ConstraintBlock(index_space, sense::Symbol, constraint::C) where {C<:AbstractConstraint}
    return ConstraintBlock(_UNNAMED_BLOCK, index_space, sense, constraint)
end

function ConstraintBlock(index_space, sense::Symbol, ::Type{C}) where {C<:AbstractConstraint}
    return ConstraintBlock(index_space, sense, C())
end

function ConstraintBlock(sense::Symbol, constraint::C) where {C<:AbstractConstraint}
    return ConstraintBlock((nothing,), sense, constraint)
end

function ConstraintBlock(sense::Symbol, ::Type{C}) where {C<:AbstractConstraint}
    return ConstraintBlock(sense, C())
end

function ConstraintBlock(
    name::Symbol, constraint::C, data; sense::Symbol=constraint_sense(constraint)
) where {C<:AbstractConstraint}
    return ConstraintBlock(name, iteration_space(constraint, data), sense, constraint)
end

function ConstraintBlock(
    name::Symbol, index_space, sense::Symbol, terms::Function, rhs; has_changed=_always_changed
)
    return ConstraintBlock(name, index_space, sense, _ConstraintFormula(terms, rhs, has_changed))
end

"""
    AbstractObjective

Supertype for a linear objective definition. Extend
[`objective_terms!`](@ref) for a concrete subtype.
The method receives the application payloads and has no numerical fallback.
"""
abstract type AbstractObjective end

"""
    ObjectiveTag{ModelIdentity} <: AbstractObjective

Singleton tag type representing an objective declared in a `@model` definition,
using the unique identity type of its model declaration.
"""
struct ObjectiveTag{M} <: AbstractObjective end

"""Emit the linear objective terms into `expression`."""
function objective_terms! end
has_changed(::AbstractObjective, old_parameters, new_parameters) = true

struct _ObjectiveFormula{T,H} <: AbstractObjective
    terms::T
    predicate::H
end

function objective_terms!(objective::_ObjectiveFormula, emit, const_data, parameters)
    return objective.terms(emit, const_data, parameters)
end
function has_changed(objective::_ObjectiveFormula, old_parameters, new_parameters)
    return objective.predicate(old_parameters, new_parameters)
end

"""
    Objective(sense, objective)
    Objective(sense, terms; has_changed=(previous, current) -> true)

A minimization (`:min`) or maximization (`:max`) objective whose coefficients
are defined by a concrete [`AbstractObjective`](@ref) value or its singleton
type. In the function-based form, `terms(emit, data, point)` emits the
complete current expression. Its predicate defaults to reevaluation at distinct
points; it must cover every dependency when supplied explicitly.
"""
struct Objective{O<:AbstractObjective}
    sense::Symbol
    objective::O
    function Objective{O}(sense::Symbol, objective::O) where {O<:AbstractObjective}
        sense ∈ (:min, :max) || throw(ArgumentError("objective sense must be :min or :max"))
        return new{O}(sense, objective)
    end
end

Objective(sense::Symbol, objective::O) where {O<:AbstractObjective} = Objective{O}(sense, objective)

Objective(sense::Symbol, ::Type{O}) where {O<:AbstractObjective} = Objective(sense, O())

function Objective(sense::Symbol, terms::Function; has_changed=_always_changed)
    return Objective(sense, _ObjectiveFormula(terms, has_changed))
end

struct _CompositeObjective{T<:Tuple} <: AbstractObjective
    fragments::T
end

"""
    compose_objectives(a::Objective, b::Objective, rest::Objective...)
    compose_objectives(a::Objective)

Combine objective fragments with the same sense (`:min` or `:max`). Different
senses throw `ArgumentError`; the unary form returns `a` unchanged. At least one
objective is required. This combines objectives only, without merging variable
or constraint definitions.

Fragments emit in argument order into one expression. Repeated variable terms
and objective constants are added, including cancellations. Each fragment must
describe its complete current contribution. Change detection is a short-circuit
OR of the child predicates: when any child changes, every fragment is reevaluated.
There is no per-fragment numerical cache. With `safe_mode=true`, checks and
diagnostics concern the aggregate objective, not individual fragments.

DSL fragments can be objective-only `@model` definitions declaring shared names
with `@external_variables`; compose their `.objective` fields with other DSL or
function-based objectives. All referenced variables must exist in the final
`ProblemDefinition`. Keep family keys fixed. After changing study data, invalidate
any application memoization and call `rebuild!`; changed keys require a new
definition. For independent clones use `similar(problem; data=independent_data)`;
data and any private caches must be independent for concurrent solves.
"""
compose_objectives(a::Objective) = a

function compose_objectives(a::Objective, b::Objective, rest::Objective...)
    objectives = (a, b, rest...)
    all(objective -> objective.sense === a.sense, objectives) ||
        throw(ArgumentError("composed objectives must have the same sense"))
    return Objective(a.sense, _CompositeObjective(map(objective -> objective.objective, objectives)))
end

"""
Typed variable-family schema. `blocks` and `by_name` refer to the same block
objects in traversal and named-access forms. Keys and cardinalities remain
runtime data; the set of families and callback types is part of this schema's
type. Numerical lookup/emitter types intentionally do not retain this schema.
"""
struct VariableCollection{B<:Tuple,N<:NamedTuple}
    blocks::B
    by_name::N
end

"""Typed constraint-family schema, with traversal and named views of the same blocks."""
struct ConstraintCollection{B<:Tuple,N<:NamedTuple}
    blocks::B
    by_name::N
end

function _assign_name!(block::Union{VariableBlock,ConstraintBlock}, name::Symbol)
    if block.name === _UNNAMED_BLOCK
        block.name = name
    elseif block.name !== name
        throw(ArgumentError("block name $(block.name) does not match NamedTuple key $name"))
    end
    return block
end

_collection(::Type{VariableCollection}, value::VariableCollection) = value
_collection(::Type{ConstraintCollection}, value::ConstraintCollection) = value

function _merge_family_names(a::NamedTuple, b::NamedTuple)
    for name in keys(a)
        haskey(b, name) && throw(ArgumentError("duplicate family name $name"))
    end
    return merge(a, b)
end

function Base.merge(a::VariableCollection, b::VariableCollection)
    return VariableCollection((a.blocks..., b.blocks...), _merge_family_names(a.by_name, b.by_name))
end
function Base.merge(
    a::VariableCollection, b::VariableCollection, c::VariableCollection, rest::VariableCollection...
)
    return foldl(merge, (c, rest...); init=merge(a, b))
end

function Base.merge(a::ConstraintCollection, b::ConstraintCollection)
    return ConstraintCollection(
        (a.blocks..., b.blocks...), _merge_family_names(a.by_name, b.by_name)
    )
end
function Base.merge(
    a::ConstraintCollection,
    b::ConstraintCollection,
    c::ConstraintCollection,
    rest::ConstraintCollection...,
)
    return foldl(merge, (c, rest...); init=merge(a, b))
end

function Base.filter(f::Function, c::ConstraintCollection)
    blocks = filter(b -> f(b.name), c.blocks)
    nt = NamedTuple{map(b -> b.name, blocks)}(blocks)
    return ConstraintCollection(blocks, nt)
end

function _collection(::Type{VariableCollection}, nt::NamedTuple)
    blocks = map((block, name) -> _assign_name!(block, name), values(nt), keys(nt))
    return VariableCollection(blocks, nt)
end

function _collection(::Type{ConstraintCollection}, nt::NamedTuple)
    blocks = map((block, name) -> _assign_name!(block, name), values(nt), keys(nt))
    return ConstraintCollection(blocks, nt)
end

function _collection(::Type{VariableCollection}, value)
    throw(ArgumentError("variables must be a NamedTuple or VariableCollection"))
end
function _collection(::Type{ConstraintCollection}, value)
    throw(ArgumentError("constraints must be a NamedTuple or ConstraintCollection"))
end

"""
    ProblemDefinition(; variables, constraints, objective)

The fixed mathematical structure of an LP. `variables` and `constraints` are
named tuples of blocks; their keys become the public mathematical family names.
`objective` is an [`Objective`](@ref).

Construct a family for each distinct mathematical rule, and use its index space
for repeated rows/columns with that rule. A thousand keys in one family share
one callback specialization. A thousand distinct family/callback types can
require substantially more client compilation. This representation preserves
static callback reachability for trimming; it does not promise compilation cost
independent of the number of distinct rules. Once a problem is built, family
names, keys, domains and callback definitions must remain fixed.
"""
struct ProblemDefinition{V<:VariableCollection,C<:ConstraintCollection,O<:Objective}
    variables::V
    constraints::C
    objective::O
end

function ProblemDefinition(; variables, constraints, objective)
    return ProblemDefinition(
        _collection(VariableCollection, variables),
        _collection(ConstraintCollection, constraints),
        objective,
    )
end

"""
    ConstData(definition, data=nothing)

Study data held by reference and treated as read-only between rebuilds.
`definition` contains the fixed variable and constraint families; `data` can
contain dynamically sized arrays, dictionaries, and mutable structs. It does
not need to be `isbits`.

ReactiveModel does not copy, freeze, or monitor these objects. The application
must keep their contents constant while solving. After modifying study data,
call [`rebuild!`](@ref) on every problem sharing those data before solving again,
even if the parameter point also changes. Otherwise cached formulas and solver
state may be stale. Changing family keys or dimensions requires a new definition.
"""
struct ConstData{D,F<:ProblemDefinition}
    definition::F
    data::D
end

ConstData(definition::ProblemDefinition) = ConstData(definition, nothing)

function Base.getproperty(data::ConstData, name::Symbol)
    return if name === :definition || name === :data
        getfield(data, name)
    else
        getproperty(getfield(data, :data), name)
    end
end

"""
    Parameters(value)

Immutable per-solve parameter point. The wrapped value describes what changes
between solves: selectors such as time-step or scenario indices, or numerical
state such as an initial stock. The constructor requires an `isbits` value,
which ReactiveModel retains by value as the previous parameter snapshot.
All points passed to one problem must have the same concrete type. Dynamically
sized arrays belong in [`ConstData`](@ref), whose contents stay constant between
calls to [`rebuild!`](@ref).
"""
struct Parameters{P}
    value::P
    function Parameters(value::P) where {P}
        isbitstype(P) || throw(
            ArgumentError(
                "Parameters must wrap a small isbits selector; put study arrays in ConstData"
            ),
        )
        return new{P}(value)
    end
end

function Base.getproperty(parameters::Parameters, name::Symbol)
    return if name === :value
        getfield(parameters, name)
    else
        getproperty(getfield(parameters, :value), name)
    end
end

Base.isequal(a::Parameters, b::Parameters) = isequal(a.value, b.value)

_parameters(parameters::Parameters) = parameters
_parameters(parameters) = Parameters(parameters)
