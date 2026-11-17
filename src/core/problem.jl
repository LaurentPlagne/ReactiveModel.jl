"""
    Workspace

Mutable numerical state owned by a [`LinearProblem`](@ref): persistent backend,
generated rows and columns, previous parameter point, reusable update storage,
and the current lazy solution. Applications normally do not construct or
mutate a workspace directly.

Complete row dictionaries and bound/RHS/cost arrays are the reference for
incremental comparisons. Backend packets may contain only changed entries;
solver-owned arrays and basis state belong to the backend. These representations
serve different lifetimes and must not be mutated independently by client code.
Result metadata contains names and keys without retaining mathematical callbacks.
"""
mutable struct Workspace{P,B,V,C,I}
    backend::B
    variable_indices::_VariableIndices{I}
    lower_values::Vector{Float64}
    upper_values::Vector{Float64}
    row_blocks::Vector{Tuple{Symbol,Any}}
    row_coefficients::Vector{Dict{Int,Float64}}
    row_rhs::Vector{Float64}
    objective_coefficients::Dict{Int,Float64}
    objective_offset::Base.RefValue{Float64}
    backend_variable_blocks::Vector{_BackendVariableBlock}
    variable_block_offsets::Vector{Int}
    backend_constraint_blocks::Vector{_BackendConstraintBlock}
    constraint_block_offsets::Vector{Int}
    backend_objective::Union{Nothing,_BackendObjective}
    parameter_snapshot::Union{Nothing,Parameters{P}}
    solution::Union{Nothing,Solution{B,V,C,I}}
    result_layout::Union{Nothing,_ResultLayout{V,C}}
    solution_generation::Base.RefValue{UInt}
    update_scratch::_UpdateScratch{I}
    update_safety::_UpdateSafety
    built::Bool
    # Status of the solve at `parameter_snapshot`: `:optimal`, or the backend's
    # status when that point was infeasible (the model is then kept, see
    # `_prepare_backend!`), or `:not_solved` before the first solve.
    last_status::Symbol
end

function Workspace{P}(backend::B, definition::ProblemDefinition; safe_mode::Bool=false) where {P,B}
    V = _solution_collection_type(typeof(definition), Val(:variables))
    C = _solution_collection_type(typeof(definition), Val(:constraints))
    variable_indices = _VariableIndices(definition)
    return Workspace{P,B,V,C,typeof(variable_indices.groups)}(
        backend,
        variable_indices,
        Float64[],
        Float64[],
        Tuple{Symbol,Any}[],
        Dict{Int,Float64}[],
        Float64[],
        Dict{Int,Float64}(),
        Ref(0.0),
        _BackendVariableBlock[],
        Int[],
        _BackendConstraintBlock[],
        Int[],
        nothing,
        nothing,
        nothing,
        nothing,
        Ref(UInt(0)),
        _UpdateScratch(variable_indices),
        _UpdateSafety(safe_mode),
        false,
        :not_solved,
    )
end

"""
    LinearProblem(const_data, initial_parameters; backend=HiGHSBackend(), safe_mode=false)

A fixed study paired with its reusable numerical workspace. `initial_parameters`
sets the parameter type and is retained for cloning before the first solve.
The concrete update plan is prepared at construction; numerical coefficients
and the backend model are built on the first solve. A fixed execution handle
keeps the mathematical schema and callback types out of the client problem type.
Its statically generated adapter retains the concrete implementation internally,
visiting preconstructed groups of at most eight families. Dimensions and
numerical values remain runtime data. Compilation still grows
with the number of distinct family rules; reuse indexed families for repeated
structure and construct the workspace outside the application's solve loop.

With `safe_mode=true`, every solve reevaluates all numerical formulas, even at
the same point and even when a change predicate returns false. Only actual
differences are sent to the backend. Inspect [`update_report`](@ref) for missed
updates that were corrected and reevaluations without numerical changes.
This mode is recommended while developing and validating an application.
Callbacks must be deterministic functions of the supplied study data and point;
study data and family keys retain their usual read-only/rebuild contract.
Safe mode does not validate the mathematical formulation or the backend itself.
The default retains the predicate-based execution path used by existing code.
"""
mutable struct LinearProblem{C,P,B,W<:Workspace{P,B}}
    const_data::C
    workspace::W
    numerical_updates::_ExecutionHandle{P}
end

"""Show a short summary; `typeof(problem)` still exposes the full compiler type."""
function Base.show(io::IO, problem::LinearProblem)
    return print(
        io,
        "LinearProblem(",
        problem.numerical_updates.model_name,
        "; ",
        problem.numerical_updates.variable_count,
        " variable families, ",
        problem.numerical_updates.constraint_count,
        " constraint families; backend=",
        nameof(typeof(problem.backend)),
        ", safe_mode=",
        problem.update_safety.enabled,
        ", built=",
        problem.built,
        ")",
    )
end

function Base.getproperty(lp::LinearProblem, name::Symbol)
    return if name === :const_data || name === :workspace || name === :numerical_updates
        getfield(lp, name)
    elseif name === :definition
        getfield(lp, :const_data).definition
    else
        getproperty(getfield(lp, :workspace), name)
    end
end

function Base.setproperty!(lp::LinearProblem, name::Symbol, value)
    return if name === :const_data
        _execute(lp.numerical_updates, _EXEC_REBIND, value)
        setfield!(lp, name, value)
    elseif name === :workspace || name === :numerical_updates
        setfield!(lp, name, value)
    else
        setproperty!(getfield(lp, :workspace), name, value)
    end
end

function LinearProblem(
    const_data::ConstData{D},
    parameters::Parameters{P};
    backend::B=HiGHSBackend(),
    safe_mode::Bool=false,
) where {D,P,B}
    workspace = Workspace{P}(backend, const_data.definition; safe_mode)
    updates = _direct_numerical_updates(const_data, workspace, parameters)
    execution = _execution_handle(const_data, workspace, updates)
    return LinearProblem{ConstData{D},P,B,typeof(workspace)}(
        const_data, workspace, execution
    )
end

function LinearProblem(
    const_data::ConstData, parameters; backend::B=HiGHSBackend(), safe_mode::Bool=false
) where {B}
    return LinearProblem(const_data, _parameters(parameters); backend, safe_mode)
end

function LinearProblem(definition::ProblemDefinition, parameters; kwargs...)
    return LinearProblem(ConstData(definition), parameters; kwargs...)
end

"""
    rebuild!(lp::LinearProblem, [new_const_data=lp.const_data]; parameters=lp.parameter_snapshot)

Reconstruct the numerical model from the study data and invalidate cached
solutions. Call this after any mutation of `ConstData.data`, before the next
solve, even if the parameter point changes. Rebuild every problem sharing the
mutated data, with no concurrent solves or result queries during the change.

The existing definition and its family keys are retained: changing numerical
tables of the same shape is supported; changing topology, dimensions, or index
spaces requires constructing a new definition and problem. The wrapper
`LinearProblem` is retained, but numerical buffers may be reallocated.

If assembly or backend initialization throws, the current solution is invalid
and the next solve must assemble again. A schema rejected before assembly
leaves the previous model and solution intact.
"""
function rebuild!(
    lp::LinearProblem{C,P,B};
    parameters::Union{Nothing,Parameters{P},P}=lp.parameter_snapshot,
) where {C,P,B}
    return _rebuild_const_data!(lp, lp.const_data, parameters)
end

function rebuild!(
    lp::LinearProblem{C,P,B}, new_const_data::ConstData;
    parameters::Union{Nothing,Parameters{P},P}=lp.parameter_snapshot,
) where {C,P,B}
    return _rebuild_const_data!(lp, new_const_data, parameters)
end

function _rebuild_const_data!(
    lp::LinearProblem{C,P,B},
    @nospecialize(new_const_data::ConstData),
    parameters::Union{Nothing,Parameters{P},P},
) where {C,P,B}
    _execute(lp.numerical_updates, _EXEC_REBIND, new_const_data)
    setfield!(lp, :const_data, new_const_data::C)
    if parameters !== nothing
        p = _checked_parameters(lp, parameters)
        _build!(lp, p)
    else
        _invalidate_model!(lp)
    end
    return lp
end

function rebuild!(
    lp::LinearProblem{C,P,B}, new_data::D;
    parameters::Union{Nothing,Parameters{P},P}=lp.parameter_snapshot,
) where {D,C<:ConstData{D},P,B}
    const_data = _execute(lp.numerical_updates, _EXEC_DATA, new_data)::C
    return _rebuild_const_data!(lp, const_data, parameters)
end

_similar_backend(backend::AbstractBackend) = typeof(backend)()

"""
    similar(lp::LinearProblem, [initial_parameters=lp.parameter_snapshot]; backend=_similar_backend(lp.backend), safe_mode=lp.update_safety.enabled, data=lp.const_data.data)

Create a fresh, independent `LinearProblem`. By default it shares the exact
same read-only `ConstData`, safe for use in another thread or task.
Pass `data=...` to provide independent study data of the same type, retaining
the mathematical definition and its keys. Dimensions and index spaces must
remain unchanged, as for `rebuild!`. Callbacks receive the supplied data from
the first assembly; the source study is not evaluated while building the clone.
A built source produces one assembly of the clone, while an unbuilt source
leaves the clone unbuilt until its first solve.

Unlike `deepcopy` (which is prohibited because solver backends contain non-shareable native
C pointers or stateful workspaces), `similar` provisions a new, isolated solver backend
and clean numerical workspace while sharing study data whose contents remain read-only between rebuilds.
The clone inherits safe mode by default but starts a separate empty update report.
HiGHS options configured through RM are retained, but the native model, basis
and solution are not copied. Passing `backend=...` supplies a new configuration.
Other configurable backends can specialize RM's `_similar_backend` factory;
the generic fallback calls their zero-argument constructor.
Use `similar(lp; safe_mode=false)` to construct an unchecked workspace after
validating the model and its change conditions.
"""
function Base.similar(
    lp::LinearProblem{C,P,B},
    initial_parameters::Union{Parameters{P},P}=lp.numerical_updates.state.current[];
    backend::AbstractBackend=_similar_backend(lp.backend),
    safe_mode::Bool=lp.update_safety.enabled,
    data=lp.const_data.data,
) where {C,P,B}
    p = initial_parameters === nothing ? nothing : _checked_parameters(lp, initial_parameters)
    data = _checked_clone_data(lp.const_data, data)
    if typeof(backend) === B
        return _execute(lp.numerical_updates, _EXEC_CLONE, (backend, p, safe_mode, data))::typeof(lp)
    end
    const_data = data === lp.const_data.data ? lp.const_data :
        _execute(lp.numerical_updates, _EXEC_DATA, data)::C
    new_lp = LinearProblem(const_data, p; backend, safe_mode)
    if p !== nothing && lp.built
        _build!(new_lp, p)
    end
    return new_lp
end

function _checked_clone_data(::ConstData{D}, data) where {D}
    data isa D || throw(ArgumentError("clone study data must retain type $D"))
    return data::D
end

Base.copy(lp::LinearProblem) = similar(lp)

function Base.deepcopy_internal(::LinearProblem, ::IdDict)
    throw(
        ArgumentError(
            "Cannot deepcopy a LinearProblem: solver backends contain stateful or native C resources. " *
            "To create an independent problem for another thread, instantiate a new LinearProblem " *
            "sharing the same ConstData: `LinearProblem(lp.const_data, initial_parameters; backend=...)` " *
            "or call `similar(lp)`.",
        ),
    )
end

function Base.deepcopy_internal(backend::AbstractBackend, ::IdDict)
    throw(
        ArgumentError(
            "Cannot deepcopy an $(typeof(backend)): backends contain stateful or native C resources. " *
            "Instantiate a fresh backend instead.",
        ),
    )
end

function Base.deepcopy_internal(::Workspace, ::IdDict)
    throw(
        ArgumentError(
            "Cannot deepcopy a Workspace: solver workspaces contain non-thread-safe solver state. " *
            "Instantiate a fresh Workspace instead.",
        ),
    )
end
