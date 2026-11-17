Base.keys(c::VariableCollection) = keys(c.by_name)
Base.keys(c::ConstraintCollection) = keys(c.by_name)

Base.getindex(c::VariableCollection, name::Symbol) = c.by_name[name]
Base.getindex(c::ConstraintCollection, name::Symbol) = c.by_name[name]

Base.haskey(c::VariableCollection, name::Symbol) = haskey(c.by_name, name)
Base.haskey(c::ConstraintCollection, name::Symbol) = haskey(c.by_name, name)

function Base.getproperty(c::VariableCollection, name::Symbol)
    return name === :blocks || name === :by_name ? getfield(c, name) : c.by_name[name]
end

function Base.getproperty(c::ConstraintCollection, name::Symbol)
    return name === :blocks || name === :by_name ? getfield(c, name) : c.by_name[name]
end

@inline function _variable_bounds(
    block::VariableBlock{K,V}, const_data::C, parameters::P, key::K
) where {K,V,C,P}
    return variable_bounds(block.variable, const_data.data, parameters.value, key)
end

@inline _has_changed(block::VariableBlock{K,V}, old::P, new::P) where {K,V,P} =
    has_changed(block.variable, old.value, new.value)
@inline _has_changed(block::ConstraintBlock{K,C}, old::P, new::P) where {K,C,P} =
    has_changed(block.constraint, old.value, new.value)
@inline _has_changed(objective::Objective{O}, old::P, new::P) where {O,P} =
    has_changed(objective.objective, old.value, new.value)

@inline _constraint_rhs(
    block::ConstraintBlock{K,C}, const_data::D, parameters::P, key::K
) where {K,C,D,P} = Float64(rhs_value(block.constraint, const_data.data, parameters.value, key))

@inline function _emit_objective!(
    objective::Objective{O}, emit::_Emitter, const_data::C, parameters::P
) where {O,C,P}
    return _emit_objective_fragment!(objective.objective, emit, const_data.data, parameters.value)
end

# Assembly, updates and composition share this boundary. Explicit O/E/D/P
# specialization retains concrete callbacks and lets matching fragment, emitter,
# data and point types reuse the same native kernel across all three callers.
@noinline function _emit_objective_fragment!(
    fragment::O, emit::E, data::D, point::P
) where {O<:AbstractObjective,E,D,P}
    objective_terms!(fragment, emit, data, point)
    return nothing
end

"""
    set_coeff!(equation, variable_block::Symbol, key, coefficient)
    set_coeff!(equation, variable_block::Symbol, coefficient)

Set or accumulate the linear coefficient of `variable_block[key]` in `equation`.
"""
@inline function set_coeff!(eq::_Emitter, block::Symbol, key, coefficient)
    idx = eq.variable_indices[VarKey(block, key)]
    c = Float64(coefficient)
    isfinite(c) || throw(ArgumentError("non-finite coefficient for $block at column $idx"))
    eq.coefficients[idx] = get(eq.coefficients, idx, 0.0) + c
    return nothing
end

@inline function set_coeff!(eq::_Emitter, block::Symbol, coefficient)
    if haskey(eq.variable_indices, VarKey(block, 0))
        set_coeff!(eq, block, 0, coefficient)
    else
        set_coeff!(eq, block, nothing, coefficient)
    end
end

mutable struct _UpdateScratch{F}
    emitter::_Emitter{F}
    columns::Vector{Int}
    lowers::Vector{Float64}
    uppers::Vector{Float64}
    rows::Vector{Int}
    row_lowers::Vector{Float64}
    row_uppers::Vector{Float64}
    changed_columns::Vector{Int}
    changed_costs::Vector{Float64}
    objective_coefficients::_BackendObjectiveCoefficients
    changed_coefficients::Vector{Int}
    changed_coefficient_values::Vector{Float64}
end

function _UpdateScratch(variable_indices)
    changed_columns = Int[]
    changed_costs = Float64[]
    changed_coefficients = Int[]
    changed_coefficient_values = Float64[]
    return _UpdateScratch(
        _Emitter(variable_indices, Dict{Int,Float64}()),
        Int[],
        Float64[],
        Float64[],
        Int[],
        Float64[],
        Float64[],
        changed_columns,
        changed_costs,
        _BackendObjectiveCoefficients(changed_columns, changed_costs),
        changed_coefficients,
        changed_coefficient_values,
    )
end

# Mutable parameter state is independent of the concrete update plan. Updaters
# share this small record, not a container which recursively contains themselves.
abstract type AbstractBlockUpdater end

struct _UpdateState{P}
    previous::Base.RefValue{Parameters{P}}
    current::Base.RefValue{Parameters{P}}
    changed::Base.RefValue{Bool}
    safety::_UpdateSafety
end

mutable struct _ExecutionArgument
    value::Any
end

"""
Fixed execution boundary. `payload` roots the concrete implementation behind
the C-callable adapter; its callback/schema types do not enter client types.
The adapter is generated statically for that payload during construction.
"""
struct _ExecutionHandle{P}
    function_pointer::Ptr{Cvoid}
    object_pointer::Ptr{Cvoid}
    payload::Any
    argument::_ExecutionArgument
    state::_UpdateState{P}
    variable_count::Int
    constraint_count::Int
    model_name::Symbol
end

abstract type _AbstractFamilyPlan end
struct _EmptyFamilyPlan <: _AbstractFamilyPlan end

"""Fixed execution node: at most eight typed families and an already-built tail."""
struct _FamilyPlan{H<:Tuple,T<:_AbstractFamilyPlan} <: _AbstractFamilyPlan
    head::H
    tail::T
end

_family_plan(blocks::Tuple) = _family_plan(identity, blocks)
_family_plan(f::F, ::Tuple{}) where {F} = _EmptyFamilyPlan()

# Map only one bounded group behind a call boundary. In particular, never
# instantiate Base.map on the full heterogeneous model tuple during construction.
# Explicit F retains concrete factories and static callback reachability.
_map_family_group(::typeof(identity), blocks::Tuple) = blocks

@noinline function _map_family_group(f::F, blocks::Tuple) where {F}
    return map(f, blocks)
end

@noinline function _family_plan(f::F, blocks::Tuple) where {F}
    if length(blocks)<=8
        return _FamilyPlan(_map_family_group(f, blocks), _EmptyFamilyPlan())
    end
    head = _map_family_group(f, ntuple(i->getfield(blocks, i), Val(8)))
    # Factories share their index/offset counters across groups. Finish the head
    # before recursing so block order, including empty families, stays unchanged.
    tail = _family_plan(f, ntuple(i->getfield(blocks, i+8), Val(length(blocks)-8)))
    return _FamilyPlan(head, tail)
end
Base.length(::_EmptyFamilyPlan) = 0
Base.length(plan::_FamilyPlan) = length(plan.head)+length(plan.tail)

"""
Reusable execution plan. Model formulas stay concrete, but traversal is bounded
by preconstructed family groups. Updaters retain study payloads and compact
numeric/index storage, never the complete ProblemDefinition or LinearProblem.
"""
struct _NumericalUpdates{P,V<:_AbstractFamilyPlan,C<:_AbstractFamilyPlan,O<:AbstractBlockUpdater}
    state::_UpdateState{P}
    variables::V
    constraints::C
    objective::O
end

# Inline only a bounded group. Inlining an entire heterogeneous model duplicates
# the visitor and its callees in the client's build/solve methods.
@inline _foreach_small(f, ::Tuple{}) = nothing
@inline function _foreach_small(f::F, blocks::Tuple) where {F}
    f(first(blocks))
    _foreach_small(f, Base.tail(blocks))
    return nothing
end

"""
Visit concrete family types in groups of at most eight. The call boundary between
groups limits the size of client IR while retaining statically reachable callback
methods for safe trimming. Family cardinalities remain runtime values.
"""
_foreach_family(f, ::Tuple{}) = nothing
@noinline function _foreach_family(f::F, blocks::Tuple) where {F}
    if length(blocks)<=8
        _foreach_small(f, blocks)
    else
        _foreach_small(f, ntuple(i->getfield(blocks, i), Val(8)))
        _foreach_family(f, ntuple(i->getfield(blocks, i+8), Val(length(blocks)-8)))
    end
    return nothing
end

_foreach_family(f, ::_EmptyFamilyPlan) = nothing
# F is explicit because this visitor forwards the callable to another kernel.
# Julia otherwise may avoid specializing a forwarded Function argument.
@noinline function _foreach_family(f::F, plan::_FamilyPlan) where {F}
    _foreach_small(f, plan.head)
    _foreach_family(f, plan.tail)
    return nothing
end

struct _EmitObjectiveFragment{E,D,P}
    emit::E
    data::D
    point::P
end

@inline function (visitor::_EmitObjectiveFragment)(fragment::O) where {O<:AbstractObjective}
    _emit_objective_fragment!(fragment, visitor.emit, visitor.data, visitor.point)
    return nothing
end

function objective_terms!(objective::_CompositeObjective, emit, data, point)
    _foreach_family(_EmitObjectiveFragment(emit, data, point), objective.fragments)
    return nothing
end

@inline _any_objective_changed_small(::Tuple{}, old::P, new::P) where {P} = false
@inline function _any_objective_changed_small(fragments::Tuple, old::P, new::P) where {P}
    return has_changed(first(fragments), old, new) ||
           _any_objective_changed_small(Base.tail(fragments), old, new)
end

# Pure tuple recursion: no mutable query state or numerical fragment caches.
# Bound inlining to the same eight-element groups as emission.
@noinline function _any_objective_changed(fragments::Tuple, old::P, new::P) where {P}
    if length(fragments) <= 8
        return _any_objective_changed_small(fragments, old, new)
    end
    return _any_objective_changed_small(ntuple(i -> getfield(fragments, i), Val(8)), old, new) ||
           _any_objective_changed(ntuple(i -> getfield(fragments, i + 8), Val(length(fragments) - 8)), old, new)
end

function has_changed(objective::_CompositeObjective, old::P, new::P) where {P}
    return _any_objective_changed(objective.fragments, old, new)
end

"""
    @emit equation variable[index] = coefficient

Emit one term from a constraint or objective definition. The block name is a
mathematical name, and its key is either the sole index or a tuple of indices.
The syntax is not an assignment: the right-hand side is the coefficient of the
variable in `equation`. Repeated emissions for the same variable are added.

The shorter `@emit variable[index] = coefficient` form uses a local emitter
named exactly `emit`; the explicit form above is canonical and is robust to
renaming that argument.
"""
function _parse_coeff_expr(expr, target_sym)
    expr isa Expr || throw(ArgumentError("@coeff expects an expression, got $expr"))

    if expr.head === :call && expr.args[1] === :←
        lhs, coeff = expr.args[2], expr.args[3]
    elseif expr.head === :call && expr.args[1] === :(=>)
        lhs, coeff = expr.args[2], expr.args[3]
    elseif expr.head === :(=)
        lhs, coeff = expr.args[1], expr.args[2]
    elseif expr.head === :call && expr.args[1] === :*
        coeff, lhs = expr.args[2], expr.args[3]
    else
        throw(
            ArgumentError(
                "@coeff: unsupported syntax `$expr`. Use `var[k] ← coeff` or `var[k] = coeff`."
            ),
        )
    end

    name_symbol = if lhs isa Symbol
        lhs
    elseif lhs isa QuoteNode && lhs.value isa Symbol
        lhs.value
    elseif lhs isa Expr && lhs.head === :ref
        root = lhs.args[1]
        s = if root isa Symbol
            root
        else
            (root isa QuoteNode && root.value isa Symbol ? root.value : nothing)
        end
        s === nothing &&
            throw(ArgumentError("left side of @coeff must be a variable name, got $root"))
        indices = lhs.args[2:end]
        s, (length(indices) == 1 ? indices[1] : Expr(:tuple, indices...))
    else
        throw(ArgumentError("left side of @coeff must be a variable name, got $lhs"))
    end

    name, key = name_symbol isa Tuple ? name_symbol : (name_symbol, :nothing)
    return :($(@__MODULE__).set_coeff!($target_sym, $(QuoteNode(name)), $key, $coeff))
end

"""
    @coeff [equation] variable[index] ← coefficient
    @coeff [equation] variable ← coefficient
    @coeff [equation] begin ... end

Set linear coefficients in `lhs_terms!` or `objective_terms!`.
Lowers to `set_coeff!(equation, :variable, index, coefficient)`.

Supports the arrow `←` (`\\leftarrow<TAB>`), `=`, `=>`, and `begin ... end` blocks.
The default equation target is automatically detected as `eq` or `emit`.
"""
macro coeff(args...)
    target = if length(args) == 2 && args[1] isa Symbol
        args[1]
    else
        :(@isdefined(eq) ? eq : emit)
    end
    expr = length(args) == 2 ? args[2] : args[1]

    if expr isa Expr && expr.head === :block
        lines = [line for line in expr.args if !(line isa LineNumberNode)]
        return esc(Expr(:block, [_parse_coeff_expr(l, target) for l in lines]...))
    else
        return esc(_parse_coeff_expr(expr, target))
    end
end
