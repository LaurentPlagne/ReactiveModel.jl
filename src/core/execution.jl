# The C-callable adapter is a statically generated Julia entry point, not an
# external solver callback. No foreign code runs between its Julia caller and
# callee, so Julia exceptions retain their original type and backtrace.
const _EXEC_UPDATE = UInt8(0)
const _EXEC_BUILD = UInt8(1)
const _EXEC_REBIND = UInt8(2)
const _EXEC_DATA = UInt8(3)
const _EXEC_LAYOUT = UInt8(4)
const _EXEC_CLONE = UInt8(5)

mutable struct _ExecutionPayload{C,U,W}
    const_data::C
    updates::U
    workspace::W
end

struct _ExecutionThunk end

function (::_ExecutionThunk)(payload::_ExecutionPayload{C}, command::UInt8, slot::_ExecutionArgument) where {C}
    argument = slot.value
    if command == _EXEC_UPDATE
        updates = payload.updates
        _foreach_family(update -> update(), updates.variables)
        _foreach_family(update -> update(), updates.constraints)
        updates.objective()
    elseif command == _EXEC_BUILD
        parameters = payload.updates.state.current[]
        _build_direct_data!(payload.const_data, payload.workspace, parameters)
        payload.updates = _direct_numerical_updates(
            payload.const_data, payload.workspace, parameters; state=payload.updates.state
        )
    elseif command == _EXEC_REBIND
        argument isa C || throw(ArgumentError("replacement study must retain the mathematical schema and data type"))
        payload.const_data = argument::C
    elseif command == _EXEC_DATA
        return _replace_execution_data(payload, argument, payload.const_data)
    elseif command == _EXEC_LAYOUT
        return _ResultLayout(payload.const_data.definition)
    elseif command == _EXEC_CLONE
        return _clone_execution(payload, argument, payload.workspace.backend, payload.updates.state, payload.const_data)
    else
        throw(ArgumentError("invalid execution operation"))
    end
    return nothing
end

function _replace_execution_data(payload, argument, const_data::ConstData{D}) where {D}
    return ConstData(const_data.definition, argument::D)
end

function _clone_execution(payload, argument, ::B, ::_UpdateState{P}, source::ConstData{D}) where {B,P,D}
    backend, parameters, safe_mode, data = argument::Tuple{B,Parameters{P},Bool,D}
    const_data = data === source.data ? source : ConstData(source.definition, data)
    cloned = LinearProblem(const_data, parameters; backend, safe_mode)
    payload.workspace.built && _build!(cloned, parameters)
    return cloned
end

@generated function _execution_pointer(::Type{T}) where {T}
    return :(@cfunction($(_ExecutionThunk()), Any, (Ref{$T}, UInt8, Ref{_ExecutionArgument})))
end

function _execution_handle(const_data, workspace, updates)
    payload = _ExecutionPayload(const_data, updates, workspace)
    T = typeof(payload)
    root = Base.cconvert(Ref{T}, payload)
    return _ExecutionHandle(
        _execution_pointer(T), Ptr{Cvoid}(Base.unsafe_convert(Ref{T}, root)), root,
        _ExecutionArgument(nothing), updates.state,
        length(updates.variables), length(updates.constraints), _definition_display_name(typeof(const_data)),
    )
end

function _execute(handle::_ExecutionHandle, command::UInt8, argument=nothing)
    @nospecialize argument
    # A reusable, concrete Ref argument keeps the adapter ABI statically known.
    # Clear it even on exceptions so a clone's backend or a rejected replacement
    # study cannot be retained by the source problem after the call.
    slot = handle.argument
    slot.value = argument
    GC.@preserve handle begin
        try
            return ccall(
                handle.function_pointer, Any, (Ptr{Cvoid}, UInt8, Ref{_ExecutionArgument}),
                handle.object_pointer, command, slot,
            )
        finally
            slot.value = nothing
        end
    end
end

function _result_layout(lp::LinearProblem{CD,P,B,Workspace{P,B,V,C,I}}) where {CD,P,B,V,C,I}
    layout = lp.workspace.result_layout
    layout === nothing || return layout
    layout = _execute(lp.numerical_updates, _EXEC_LAYOUT)::_ResultLayout{V,C}
    lp.workspace.result_layout = layout
    return layout
end
