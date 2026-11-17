function _parse_when_args(args; allow_domain=false, context=:constraint)
    when_cond = nothing
    lhs_cond = nothing
    rhs_cond = nothing
    seen = Set{Symbol}()
    idx = 1
    while idx <= length(args)
        tok = args[idx]
        if allow_domain && tok isa Expr && tok.head in (:(=), :kw) && tok.args[1] === :domain
            idx += 1
            continue
        end
        target = if tok in (:update_when, :when, :changed_when)
            :all
        elseif tok in (:update_lhs_when, :lhs_when)
            :lhs
        elseif tok in (:update_rhs_when, :rhs_when)
            :rhs
        else
            throw(ArgumentError("unknown update clause: $tok"))
        end
        idx < length(args) || throw(ArgumentError("missing condition after $tok"))
        context === :constraint || target === :all ||
            throw(ArgumentError("$tok is only valid on constraints; use update_when on $context declarations"))
        target in seen && throw(ArgumentError("duplicate update clause for $target"))
        push!(seen, target)
        cond = args[idx+1]
        cond === nothing && throw(ArgumentError("missing condition after $tok"))
        if target === :all
            when_cond = cond
        elseif target === :lhs
            lhs_cond = cond
        else
            rhs_cond = cond
        end
        idx += 2
    end
    # Component-specific conditions override the general condition in either order.
    :lhs in seen || (lhs_cond = when_cond)
    :rhs in seen || (rhs_cond = when_cond)
    return when_cond, lhs_cond, rhs_cond
end

function _parse_variable_domain(args)
    domain = QuoteNode(:continuous)
    seen = false
    for arg in args
        if arg isa Expr && arg.head in (:(=), :kw) && arg.args[1] === :domain
            seen && throw(ArgumentError("duplicate variable domain"))
            domain = arg.args[2]
            seen = true
            if domain isa QuoteNode && domain.value ∉ (:continuous, :integer, :binary)
                throw(ArgumentError("variable domain must be :continuous, :integer or :binary"))
            end
        end
    end
    domain
end

"""
    _extract_param_deps(expr, param_sym::Symbol)

Inspect an AST expression to identify which fields of `param_sym` are accessed (e.g. `param.t`).
Returns a `Set{Symbol}` of field names, `Set{Symbol}()` if `param_sym` does not appear,
or `:opaque` if `param_sym` is passed to an arbitrary function or indexed dynamically.
"""
function _extract_param_deps(expr, param_sym::Symbol)
    deps = Set{Symbol}()
    is_opaque = false
    recognized_syms = Set([param_sym])

    function walk(e)
        is_opaque && return
        if e in recognized_syms
            is_opaque = true
            return
        elseif e isa Expr
            if e.head === :. && e.args[1] in recognized_syms
                field = e.args[2] isa QuoteNode ? e.args[2].value : e.args[2]
                if field isa Symbol
                    push!(deps, field)
                    return
                else
                    is_opaque = true
                    return
                end
            else
                for a in e.args
                    walk(a)
                end
            end
        end
    end

    walk(expr)
    return is_opaque ? :opaque : deps
end

function _formula_param_deps(expr, setup_stmts, param_sym)
    setup = _filter_setup_stmts(setup_stmts, expr)
    _extract_param_deps(Expr(:block, setup..., expr), param_sym)
end

"""
    _deps_to_condition(deps)

Compare dependencies with `isequal`, consistently with parameter snapshots.
Compare the whole point if opaque, or return `false` if invariant.
"""
function _deps_to_condition(deps)
    if deps === :opaque
        return :(!$(GlobalRef(Base, :isequal))(o, n))
    elseif isempty(deps)
        return false
    else
        fields = sort(collect(deps))
        expr = :(!$(GlobalRef(Base, :isequal))(o.$(fields[1]), n.$(fields[1])))
        for f in fields[2:end]
            expr = Expr(:||, expr, :(!$(GlobalRef(Base, :isequal))(o.$f, n.$f)))
        end
        return expr
    end
end

"""
    _normalize_when_cond(cond, param_sym::Symbol, setup_stmts=[])

Normalize user-supplied `when` conditions supporting:
- Symbol `:t` or nonempty tuple `(:t, :scenario)` -> compare fields with `isequal`
- Lambda `(old, new) -> expr`
- `changed(param.t)` or `δ(param.t)`
- Standard boolean literals or legacy `(o, n)` expressions
"""
function _normalize_when_cond(cond, param_sym::Symbol, setup_stmts=[])
    cond === false && return false
    cond === true && return true

    if cond isa QuoteNode && cond.value isa Symbol
        return _deps_to_condition(Set([cond.value]))
    end
    if cond isa Symbol && cond !== :false && cond !== :true
        return _deps_to_condition(Set([cond]))
    end
    if cond isa Expr && cond.head === :tuple
        syms = Symbol[]
        for a in cond.args
            if a isa QuoteNode && a.value isa Symbol
                push!(syms, a.value)
            elseif a isa Symbol
                push!(syms, a)
            else
                throw(ArgumentError("update field tuples must contain only field names, for example (:t, :scenario)"))
            end
        end
        isempty(syms) && throw(ArgumentError("update field tuples must not be empty; use false for an invariant component"))
        return _deps_to_condition(Set(syms))
    end
    if cond isa Expr && cond.head === :->
        args_node = cond.args[1]
        args_node isa Expr && args_node.head === :tuple && length(args_node.args) == 2 &&
            all(a -> a isa Symbol, args_node.args) && args_node.args[1] !== args_node.args[2] ||
            throw(ArgumentError("update predicates must declare two distinct arguments: (previous, current) -> condition"))
        return :(($cond)(o, n))
    end
    if cond isa Expr && cond.head === :call && cond.args[1] in (:changed, :δ, :diff)
        length(cond.args) == 2 || throw(ArgumentError("$(cond.args[1]) expects one dependency expression"))
        deps = _formula_param_deps(cond.args[2], setup_stmts, param_sym)
        return _deps_to_condition(deps)
    end
    cond isa Expr || throw(ArgumentError("expected a field name, field tuple, Boolean or two-argument update predicate"))
    # Keep legacy expressions using the implicit o/n bindings compatible.
    return cond
end

function _build_nested_comprehension(key_tuple, iters, filter_cond)
    if isempty(iters)
        return :([$key_tuple])
    end
    inner = if filter_cond !== nothing
        Expr(:generator, key_tuple, Expr(:filter, filter_cond, Expr(:(=), iters[end][1], iters[end][2])))
    else
        Expr(:generator, key_tuple, Expr(:(=), iters[end][1], iters[end][2]))
    end
    for i in reverse(iters[1:end-1])
        inner = Expr(:flatten, Expr(:generator, inner, Expr(:(=), i[1], i[2])))
    end
    return Expr(:comprehension, inner)
end

function _extract_symbols(expr, syms=Set{Symbol}())
    if expr isa Symbol
        push!(syms, expr)
    elseif expr isa Expr
        if expr.head === :. && length(expr.args) == 2 && expr.args[2] isa QuoteNode
            _extract_symbols(expr.args[1], syms)
        else
            for arg in expr.args
                _extract_symbols(arg, syms)
            end
        end
    end
    return syms
end

function _extract_lhs_symbols!(lhs, syms=Set{Symbol}())
    if lhs isa Symbol
        push!(syms, lhs)
    elseif lhs isa Expr
        if lhs.head in (:tuple, :vect, :parameters)
            for a in lhs.args
                _extract_lhs_symbols!(a, syms)
            end
        elseif lhs.head === :(::)
            _extract_lhs_symbols!(lhs.args[1], syms)
        end
    end
    return syms
end

function _assigned_symbols(stmt)
    syms = Set{Symbol}()
    if stmt isa Expr && stmt.head === :(=)
        _extract_lhs_symbols!(stmt.args[1], syms)
    end
    return syms
end

function _filter_setup_stmts(setup_stmts, expr)
    isempty(setup_stmts) && return setup_stmts
    used_in_expr = _extract_symbols(expr)
    stmt_info = []
    for s in setup_stmts
        defs = _assigned_symbols(s)
        uses = _extract_symbols(s isa Expr && s.head === :(=) ? s.args[2] : s)
        push!(stmt_info, (s, defs, uses))
    end
    needed_syms = copy(used_in_expr)
    needed_indices = Set{Int}()
    changed = true
    while changed
        changed = false
        for (idx, (s, defs, uses)) in reverse(collect(enumerate(stmt_info)))
            if idx ∉ needed_indices
                if isempty(defs) || !isempty(intersect(defs, needed_syms))
                    push!(needed_indices, idx)
                    union!(needed_syms, uses)
                    changed = true
                end
            end
        end
    end
    return [stmt_info[idx][1] for idx in 1:length(stmt_info) if idx in needed_indices]
end

# ---------------------------------------------------------------------------
# Code generation (Codegen) for ReactiveModel declarative macros
# ---------------------------------------------------------------------------

function _validate_setup_target(lhs)
    lhs isa Symbol && return [lhs]
    if lhs isa Expr && lhs.head === :tuple && !isempty(lhs.args)
        return reduce(vcat, (_validate_setup_target(a) for a in lhs.args))
    end
    throw(ArgumentError("setup assignments must bind names or tuples of names; define helpers outside @model"))
end

function _validate_dsl_value(expr)
    expr isa Expr || return nothing
    expr.head in (:(=), :+=, :-=, :*=, :/=, :block, :let, :for, :while, :try,
                  :function, :->, :return, :global, :local, :macrocall, :quote) &&
        throw(ArgumentError("unsupported statement in a DSL value: $expr; use a pure helper defined outside @model"))
    for arg in expr.args
        # Generator bindings are syntax, not assignments to model setup names.
        if expr.head in (:generator, :filter) && arg isa Expr && arg.head === :(=)
            length(arg.args) == 2 || throw(ArgumentError("invalid generator binding"))
            _validate_dsl_value(arg.args[2])
        else
            _validate_dsl_value(arg)
        end
    end
    nothing
end

function _validate_model_structure(variables, constraints, setup, known_vars, data_sym, param_sym)
    for families in (variables, constraints)
        names = first.(families)
        length(unique(names)) == length(names) || throw(ArgumentError("duplicate family name in @model"))
    end
    bindings = Set([data_sym, param_sym])
    definitions = [_validate_setup_target(s.args[1]) for s in setup]
    pending = Set(Iterators.flatten(definitions))
    for (stmt, names) in zip(setup, definitions)
        length(unique(names)) == length(names) || throw(ArgumentError("duplicate name in setup assignment"))
        any(n -> n in bindings || n in known_vars, names) &&
            throw(ArgumentError("setup names cannot be reassigned or shadow model arguments or variables"))
        rhs = stmt.args[2]
        isempty(intersect(_extract_symbols(rhs), pending)) ||
            throw(ArgumentError("setup assignments must use earlier bindings; forward or recursive aliases are unsupported"))
        _contains_variable(rhs, known_vars) && throw(ArgumentError("setup aliases must be scalar data or parameters, not decision-variable expressions"))
        _validate_dsl_value(rhs)
        union!(bindings, names)
        setdiff!(pending, names)
    end
    for variable in variables
        _, _, _, lower, upper, _, domain = variable
        for value in (lower, upper, domain)
            _contains_variable(value, known_vars) && throw(ArgumentError("variable bounds and domains cannot depend on decision variables"))
            _validate_dsl_value(value)
        end
        deps = _formula_param_deps(domain, setup, param_sym)
        (deps !== :opaque && isempty(deps)) ||
            throw(ArgumentError("variable domains must be independent of parameters"))
    end
    for family in (variables..., constraints...)
        name, iters, filter = family[1:3]
        name in bindings && throw(ArgumentError("family names cannot shadow setup names or model arguments: $name"))
        for (index, space) in iters
            (index in bindings || index in known_vars) && throw(ArgumentError("iterator names cannot shadow setup, model arguments or variables: $index"))
            _validate_dsl_value(space)
        end
        _validate_dsl_value(filter)
        structure = Expr(:tuple, (i[2] for i in iters)..., filter)
        _contains_variable(structure, known_vars) && throw(ArgumentError("family indices and filters cannot depend on decision variables"))
        deps = _formula_param_deps(structure, setup, param_sym)
        (deps !== :opaque && isempty(deps)) ||
            throw(ArgumentError("family indices and filters must be independent of parameters; create a new definition for new keys"))
    end
end

function _iteration_space_method(struct_expr, data_sym, iters, filter_cond, setup_stmts)
    if isempty(iters)
        return :(ReactiveModel.iteration_space(::$struct_expr, $data_sym) = [0])
    end

    iter_vars = [i[1] for i in iters]
    key_tuple = length(iter_vars) == 1 ? iter_vars[1] : Expr(:tuple, iter_vars...)
    comp_expr = _build_nested_comprehension(key_tuple, iters, filter_cond)
    setup_iter = _filter_setup_stmts(setup_stmts, comp_expr)
    return :(function ReactiveModel.iteration_space(::$struct_expr, $data_sym)
        $(setup_iter...)
        return $comp_expr
    end)
end

function _build_model_expr(model_name::Symbol, data_sym::Symbol, param_sym::Symbol, body::Expr)
    body.head === :block || throw(ArgumentError("@model requires a begin/end block"))
    data_sym !== param_sym || throw(ArgumentError("data and parameter argument names must differ"))
    setup_stmts = []
    variables = []
    constraints = []
    objective = nothing
    known_vars = Set{Symbol}()

    # 1. First pass: identify all variables
    for line in body.args
        line isa LineNumberNode && continue
        if line isa Expr && line.head === :macrocall
            mname = line.args[1]
            if mname === Symbol("@variable")
                length(line.args) >= 3 || throw(ArgumentError("@variable requires a declaration"))
                decl = line.args[3]
                ref_expr, l_expr, u_expr = parse_variable_bounds(decl)
                vname, iters, filter_cond = parse_domain_ref(ref_expr)
                when_cond, _, _ = _parse_when_args(line.args[4:end]; allow_domain=true, context=:variable)
                domain = _parse_variable_domain(line.args[4:end])
                push!(known_vars, vname)
                push!(variables, (vname, iters, filter_cond, l_expr, u_expr, when_cond, domain))
            elseif mname in (Symbol("@external_variables"), Symbol("@known_variables"))
                length(line.args) >= 3 || throw(ArgumentError("external variables require at least one name"))
                for v in line.args[3:end]
                    if v isa Symbol
                        push!(known_vars, v)
                    elseif v isa Expr && v.head === :tuple
                        isempty(v.args) && throw(ArgumentError("external variables require at least one name"))
                        for sub_v in v.args
                            sub_v isa Symbol || throw(ArgumentError("external variables must be names"))
                            push!(known_vars, sub_v)
                        end
                    else
                        throw(ArgumentError("external variables must be names"))
                    end
                end
            end
        end
    end

    # 2. Second pass: constraints, objective, and setup statements
    for line in body.args
        line isa LineNumberNode && continue
        if line isa Expr && line.head === :macrocall
            mname = line.args[1]
            if mname in (Symbol("@variable"), Symbol("@external_variables"), Symbol("@known_variables"))
                # Already collected
                continue
            elseif mname === Symbol("@constraint")
                length(line.args) >= 4 || throw(ArgumentError("@constraint requires a name and one relation"))
                # @constraint Name[domain] equation
                # or @constraint Name begin equation end
                ref_expr = line.args[3]
                cname, iters, filter_cond = parse_domain_ref(ref_expr)
                eq_expr = line.args[4]
                if eq_expr isa Expr && eq_expr.head === :block
                    eq_lines = [l for l in eq_expr.args if !(l isa LineNumberNode)]
                    length(eq_lines) == 1 || throw(ArgumentError("@constraint expects exactly one relation"))
                    eq_expr = only(eq_lines)
                end
                sense, lhs, rhs = split_relation(eq_expr)
                _, lhs_cond, rhs_cond = _parse_when_args(line.args[5:end])
                push!(constraints, (cname, iters, filter_cond, sense, lhs, rhs, lhs_cond, rhs_cond))
            elseif mname === Symbol("@objective")
                length(line.args) >= 4 || throw(ArgumentError("@objective requires a sense and one expression"))
                # @objective Min / Max expr
                sense = line.args[3]
                sense in (:Min, :min, :minimize, :Max, :max, :maximize) ||
                    throw(ArgumentError("objective sense must be Min or Max, got: $sense"))
                objective === nothing || throw(ArgumentError("duplicate @objective declaration"))
                obj_sense = (sense in (:Min, :min, :minimize)) ? :min : :max
                obj_expr = line.args[4]
                if obj_expr isa Expr && obj_expr.head === :block
                    obj_lines = [l for l in obj_expr.args if !(l isa LineNumberNode)]
                else
                    obj_lines = [obj_expr]
                end
                when_cond, _, _ = _parse_when_args(line.args[5:end]; context=:objective)
                objective = (obj_sense, obj_lines, when_cond)
            else
                throw(ArgumentError("unsupported macro inside @model: $mname"))
            end
        else
            line isa Expr && line.head === :(=) && length(line.args) == 2 ||
                throw(ArgumentError("@model accepts declarations and setup assignments; loops, statements and conditional declarations are unsupported"))
            push!(setup_stmts, line)
        end
    end

    _validate_model_structure(variables, constraints, setup_stmts, known_vars, data_sym, param_sym)
    for (_, _, _, _, lhs, rhs, _, _) in constraints
        _validate_dsl_value(lhs)
        _validate_dsl_value(rhs)
    end
    if objective !== nothing
        foreach(_validate_dsl_value, objective[2])
    end
    emitter = gensym(:emit)
    scalar_key = gensym(:key)

    model_identity = gensym(Symbol(model_name, "_Identity"))
    emitted_code = [:(struct $model_identity end)]
    push!(emitted_code, :(ReactiveModel._model_display_name(::Type{$model_identity}) = $(QuoteNode(model_name))))

    # 3. Variable emission
    var_block_fields = []
    for (vname, iters, filter_cond, l_expr, u_expr, when_cond, domain) in variables
        struct_expr = :(ReactiveModel.VariableTag{$model_identity, $(QuoteNode(vname))})

        # iteration_space
        push!(
            emitted_code,
            _iteration_space_method(struct_expr, data_sym, iters, filter_cond, setup_stmts),
        )

        # variable_bounds
        is_const_bounds = (l_expr isa Number) && (u_expr isa Number)
        return_bounds_expr = (u_expr === nothing) ? :($l_expr) : :(Float64($l_expr), Float64($u_expr))
        if is_const_bounds && isempty(iters)
            push!(emitted_code, :(ReactiveModel.variable_bounds(::$struct_expr) = (Float64($l_expr), Float64($u_expr))))
        elseif isempty(iters)
            setup_bounds = _filter_setup_stmts(setup_stmts, return_bounds_expr)
            push!(emitted_code, :(function ReactiveModel.variable_bounds(::$struct_expr, $data_sym, $param_sym, $scalar_key)
                $(setup_bounds...)
                return $return_bounds_expr
            end))
        else
            iter_vars = [i[1] for i in iters]
            key_arg = length(iter_vars) == 1 ? iter_vars[1] : Expr(:tuple, iter_vars...)
            setup_bounds = _filter_setup_stmts(setup_stmts, return_bounds_expr)
            push!(emitted_code, :(function ReactiveModel.variable_bounds(::$struct_expr, $data_sym, $param_sym, $key_arg)
                $(setup_bounds...)
                return $return_bounds_expr
            end))
        end

        final_when_cond = if when_cond !== nothing
            _normalize_when_cond(when_cond, param_sym, setup_stmts)
        else
            _deps_to_condition(_formula_param_deps(Expr(:tuple, l_expr, u_expr), setup_stmts, param_sym))
        end
        push!(emitted_code, :(ReactiveModel.has_changed(::$struct_expr, o, n) = $final_when_cond))

        push!(var_block_fields, (vname, :(ReactiveModel.VariableBlock(
            $(QuoteNode(vname)), $struct_expr(), $data_sym; domain=$domain))))
    end

    # 4. Constraint emission
    ctr_block_fields = []
    for (cname, iters, filter_cond, sense, lhs, rhs, lhs_cond, rhs_cond) in constraints
        struct_expr = :(ReactiveModel.ConstraintTag{$model_identity, $(QuoteNode(cname))})

        # iteration_space
        push!(
            emitted_code,
            _iteration_space_method(struct_expr, data_sym, iters, filter_cond, setup_stmts),
        )

        # constraint_sense
        push!(emitted_code, :(ReactiveModel.constraint_sense(::$struct_expr) = $(QuoteNode(sense))))

        # lhs_terms!
        lhs_stmts, lhs_constant = lower_affine_expr(lhs, known_vars, emitter)
        rhs_stmts, rhs_constant = lower_affine_expr(rhs, known_vars, emitter, -1.0)
        append!(lhs_stmts, rhs_stmts)
        # rhs_constant carries a minus sign from coefficient normalization.
        normalized_rhs = :(-($rhs_constant) - ($lhs_constant))
        iter_vars = [i[1] for i in iters]
        key_arg = isempty(iters) ? scalar_key : (length(iter_vars) == 1 ? iter_vars[1] : Expr(:tuple, iter_vars...))
        setup_lhs = _filter_setup_stmts(setup_stmts, Expr(:block, lhs_stmts...))
        push!(emitted_code, :(function ReactiveModel.lhs_terms!(::$struct_expr, $emitter, $data_sym, $param_sym, $key_arg)
            $(setup_lhs...)
            $(lhs_stmts...)
        end))

        # rhs_value
        if lhs_constant == 0.0 && rhs_constant == 0.0
            push!(emitted_code, :(ReactiveModel.rhs_value(::$struct_expr) = 0.0))
        else
            setup_rhs = _filter_setup_stmts(setup_stmts, normalized_rhs)
            push!(emitted_code, :(function ReactiveModel.rhs_value(::$struct_expr, $data_sym, $param_sym, $key_arg)
                $(setup_rhs...)
                return Float64($normalized_rhs)
            end))
        end

        final_lhs_cond = if lhs_cond !== nothing
            _normalize_when_cond(lhs_cond, param_sym, setup_stmts)
        else
            _deps_to_condition(_formula_param_deps(Expr(:block, lhs_stmts...), setup_stmts, param_sym))
        end

        final_rhs_cond = if rhs_cond !== nothing
            _normalize_when_cond(rhs_cond, param_sym, setup_stmts)
        else
            _deps_to_condition(_formula_param_deps(normalized_rhs, setup_stmts, param_sym))
        end

        push!(emitted_code, :(ReactiveModel.lhs_has_changed(::$struct_expr, o, n) = $final_lhs_cond))
        push!(emitted_code, :(ReactiveModel.rhs_has_changed(::$struct_expr, o, n) = $final_rhs_cond))

        push!(ctr_block_fields, (cname, :(ReactiveModel.ConstraintBlock($(QuoteNode(cname)), $struct_expr(), $data_sym))))
    end

    # 5. Objective emission
    obj_expr_val = :(ReactiveModel.Objective(:min, (eq, d, s) -> nothing))
    if objective !== nothing
        obj_sense, obj_lines, obj_when_cond = objective
        obj_struct_expr = :(ReactiveModel.ObjectiveTag{$model_identity})
        obj_stmts = []
        length(obj_lines) == 1 || throw(ArgumentError("@objective expects exactly one affine expression"))
        terms, constant = lower_affine_expr(only(obj_lines), known_vars, emitter)
        append!(obj_stmts, terms)
        constant == 0.0 || push!(obj_stmts, :(ReactiveModel._add_constant!($emitter, $constant)))
        setup_obj = _filter_setup_stmts(setup_stmts, Expr(:block, obj_stmts...))
        push!(emitted_code, :(function ReactiveModel.objective_terms!(::$obj_struct_expr, $emitter, $data_sym, $param_sym)
            $(setup_obj...)
            $(obj_stmts...)
        end))

        final_obj_when = if obj_when_cond !== nothing
            _normalize_when_cond(obj_when_cond, param_sym, setup_stmts)
        else
            _deps_to_condition(_formula_param_deps(Expr(:block, obj_stmts...), setup_stmts, param_sym))
        end
        push!(emitted_code, :(ReactiveModel.has_changed(::$obj_struct_expr, o, n) = $final_obj_when))

        obj_expr_val = :(ReactiveModel.Objective($(QuoteNode(obj_sense)), $obj_struct_expr()))
    end

    # 6. Model constructor function
    var_symbols = [QuoteNode(v[1]) for v in var_block_fields]
    var_calls = [v[2] for v in var_block_fields]
    ctr_symbols = [QuoteNode(c[1]) for c in ctr_block_fields]
    ctr_calls = [c[2] for c in ctr_block_fields]

    setup_ctor = _filter_setup_stmts(setup_stmts, Expr(:block, var_calls..., ctr_calls...))
    push!(emitted_code, :(function $model_name($data_sym)
        $(setup_ctor...)
        variables_nt = NamedTuple{($(var_symbols...),)}(($(var_calls...),))
        constraints_nt = NamedTuple{($(ctr_symbols...),)}(($(ctr_calls...),))
        return ReactiveModel.ProblemDefinition(
            variables = variables_nt,
            constraints = constraints_nt,
            objective = $obj_expr_val
        )
    end))

    return Expr(:block, emitted_code...)
end

"""
    @model ModelName(data, param) begin ... end
    @model ModelName begin ... end

Define a reactive linear optimization model with a JuMP-like declarative syntax
and compile-time lowering to typed Julia callbacks.

Parameters default to `(data, param)`.
Variables accept `domain=:continuous` (default), `domain=:integer` or
`domain=:binary`, for example `@variable x in [0, 1] domain=:binary`.
The domain is fixed when the definition is constructed. The selected backend
must support it; HiGHS supports integer and binary variables, Clp does not.
Reactivity conditions (`update_when`, `update_lhs_when`, `update_rhs_when`) are
**automatically inferred** by AST inspection with conservative whole-point fallback.
For explicit dependencies, use `update_when :t` or `update_rhs_when (:t, :scenario)`.
Inferred dependencies and field lists compare values with `isequal`.
For a custom rule, use `update_when ((previous, current) -> condition)` returning
a Boolean. Explicit conditions replace inference for their component and must
cover all its dependencies, including normalized affine constants.
`update_lhs_when` and `update_rhs_when` apply only to constraints; a specific
condition overrides `update_when` regardless of clause order. Use `false` for
an invariant component. Legacy expressions using implicit `o`/`n` remain accepted.
"""
macro model(sig, body)
    model_name, data_sym, param_sym = if sig isa Symbol
        sig, :data, :param
    elseif sig isa Expr && sig.head === :call
        1 <= length(sig.args) <= 3 && all(a -> a isa Symbol, sig.args) ||
            throw(ArgumentError("@model signature must name the model and at most two distinct argument names"))
        name = sig.args[1]
        d = length(sig.args) >= 2 ? sig.args[2] : :data
        p = length(sig.args) >= 3 ? sig.args[3] : :param
        name, d, p
    else
        throw(ArgumentError("Expected `@model ModelName(data, param) begin ... end`, got: $sig"))
    end

    return esc(_build_model_expr(model_name, data_sym, param_sym, body))
end
