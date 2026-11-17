# Compact diagnostics do not change type identity, dispatch or compilation.
# Scope the formatting to RM-owned types; retain an IOContext escape hatch for
# inspecting the complete parameter graph in compiler/debugging investigations.
const _CompactRMTypes = Union{
    LinearProblem,
    Workspace,
    Solution,
    ConstData,
    ProblemDefinition,
    VariableCollection,
    ConstraintCollection,
    _NumericalUpdates,
    _ExecutionHandle,
    _ExecutionPayload,
    _FamilyPlan,
    _VariableIndices,
    _UpdateScratch,
    VariableUpdater,
    ConstraintUpdater,
    ObjectiveUpdater,
}

_model_display_name(::Type) = :functional
_formula_display_name(::Type) = :functional
_formula_display_name(::Type{<:ObjectiveTag{M}}) where {M} = _model_display_name(M)
_family_display_name(::Type{T}) where {T} = nameof(T)
_family_display_name(::Type{<:VariableTag{M,N}}) where {M,N} = N
_family_display_name(::Type{<:ConstraintTag{M,N}}) where {M,N} = N

function _definition_display_name(::Type{C}) where {C<:ConstData}
    definition = fieldtype(C, :definition)
    objective = fieldtype(fieldtype(definition, :objective), :objective)
    return _formula_display_name(objective)
end

function _show_compact_type(io, ::Type{T}) where {T}
    return print(io, "ReactiveModel.", nameof(T), "{…}")
end
function _show_compact_type(io, ::Type{LinearProblem{C,P,B,W}}) where {C,P,B,W}
    return print(
        io,
        "ReactiveModel.LinearProblem{",
        "P=",
        nameof(P),
        ", B=",
        nameof(B),
        ", …}",
    )
end
function _show_compact_type(io, ::Type{VariableUpdater{P,K,V,D,B,S}}) where {P,K,V,D,B,S}
    return print(
        io, "ReactiveModel.VariableUpdater{", _family_display_name(V), "; B=", nameof(B), ", …}"
    )
end
function _show_compact_type(io, ::Type{ConstraintUpdater{P,K,C,D,B,S}}) where {P,K,C,D,B,S}
    return print(
        io, "ReactiveModel.ConstraintUpdater{", _family_display_name(C), "; B=", nameof(B), ", …}"
    )
end

# Keep bare-stream Type printing on Base's generic method. Specializing that
# method for every IO also makes some SparseArrays error paths dynamically
# dispatch on an abstract Type during safe trimming. Debug/REPL/type-signature
# output uses IOContext, where the compact representation is explicitly scoped.
function Base.show(io::IOContext, type::Type{<:_CompactRMTypes})
    # Compiler diagnostics can print a type containing free TypeVars. Such a
    # value need not bind a static parameter in a `::Type{T}` method; use the
    # actual argument so printing an unresolved call does not itself throw.
    if get(io, :reactivemodel_full_types, false) || !isconcretetype(type)
        invoke(show, Tuple{IO,Type}, io, type)
    else
        _show_compact_type(io, type)
    end
end

"""
    ReactiveModel.full_type_string(type_or_problem) -> String

Return the complete type spelling, disabling RM's compact type display. REPL
output, contextual nested signatures and `MethodError` messages abbreviate RM's
large internal parameter graphs. Abbreviation does not alter the underlying
type or make compilation cheaper. Bare streams retain Base's full formatting;
use `show(IOContext(io, :compact => true), type)` for compact stream output.
The returned text can be large and allocates memory.

For complete exception diagnostics, use
`showerror(IOContext(stderr, :reactivemodel_full_types => true), error)`.
The same IOContext key works with `show` on nested types. Compact type labels
contain ellipses/labels and are descriptions, not expressions to paste as Julia
type declarations. Use this function when investigating hidden type differences.
"""
function full_type_string(type::Type)
    io=IOBuffer()
    show(IOContext(io, :reactivemodel_full_types=>true, :limit=>false), type)
    return String(take!(io))
end
full_type_string(problem::LinearProblem) = full_type_string(typeof(problem))
