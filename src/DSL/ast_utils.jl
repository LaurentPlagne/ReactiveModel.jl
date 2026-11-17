# ---------------------------------------------------------------------------
# AST analysis and lowering utilities for the ReactiveModel DSL
# ---------------------------------------------------------------------------

"""
    split_relation(expr) -> (sense::Symbol, lhs, rhs)

Split an equation or inequality at its relational operator.
Return the relational sense (`:le`, `:eq`, `:ge`), the left-hand side, and the right-hand side.
"""
function split_relation(expr)
    if expr isa Expr
        if expr.head === :call && length(expr.args) == 3
            op = expr.args[1]
            if op in (:<=, :le, :≤)
                return :le, expr.args[2], expr.args[3]
            elseif op in (:>=, :ge, :≥)
                return :ge, expr.args[2], expr.args[3]
            elseif op in (:(==), :eq)
                return :eq, expr.args[2], expr.args[3]
            end
        end
    end
    throw(ArgumentError("Expected relation (<=, ==, >=), got: $expr"))
end

# True if an AST refers to a declared variable. This also catches variables
# buried in unsupported calls, so no mathematical fragment is silently dropped.
_contains_variable(e::Symbol, vars) = e in vars
_contains_variable(e::Expr, vars) = any(a -> _contains_variable(a, vars), e.args)
_contains_variable(e, vars) = false

_affine_scale(a, b) = a == 0.0 || b == 0.0 ? 0.0 :
    a == 1.0 ? b : b == 1.0 ? a : :($a * $b)
_affine_add(a, b) = a == 0.0 ? b : b == 0.0 ? a : :($a + $b)

"""
    lower_affine_expr(expr, known_vars, target=:eq, scale=1.0)

Return coefficient emission statements and the scalar constant of an affine
expression. Products and quotients may scale an affine expression; nonlinear
or unsupported operations involving declared variables fail at expansion.
"""
function lower_affine_expr(expr, known_vars::Set{Symbol}, target=:eq, scale=1.0)
    if expr isa Symbol && expr in known_vars
        return [:(ReactiveModel.set_coeff!($target, $(QuoteNode(expr)), $scale))], 0.0
    elseif expr isa Expr && expr.head === :ref && expr.args[1] in known_vars
        length(expr.args) > 1 || throw(ArgumentError("use a bare scalar variable name instead of empty indexing"))
        any(a -> _contains_variable(a, known_vars), expr.args[2:end]) &&
            throw(ArgumentError("variable-dependent index in affine expression: $expr"))
        indices = expr.args[2:end]
        key = length(indices) == 1 ? indices[1] : Expr(:tuple, indices...)
        return [:(ReactiveModel.set_coeff!($target, $(QuoteNode(expr.args[1])), $key, $scale))], 0.0
    elseif !_contains_variable(expr, known_vars)
        return [], _affine_scale(scale, expr)
    elseif expr isa Expr && expr.head === :call && expr.args[1] in (:+, :-)
        op = expr.args[1]
        op === :- && !(2 <= length(expr.args) <= 3) &&
            throw(ArgumentError("subtraction expects one or two operands"))
        stmts = []
        constant = 0.0
        for (i, arg) in enumerate(expr.args[2:end])
            factor = op === :- && (i > 1 || length(expr.args) == 2) ? :(-$scale) : scale
            terms, offset = lower_affine_expr(arg, known_vars, target, factor)
            append!(stmts, terms)
            constant = _affine_add(constant, offset)
        end
        return stmts, constant
    elseif expr isa Expr && expr.head === :call && expr.args[1] === :*
        factors = expr.args[2:end]
        variable_factors = findall(a -> _contains_variable(a, known_vars), factors)
        length(variable_factors) == 1 ||
            throw(ArgumentError("nonlinear product in affine expression: $expr"))
        index = only(variable_factors)
        factor = scale
        for i in eachindex(factors)
            i == index || (factor = _affine_scale(factor, factors[i]))
        end
        return lower_affine_expr(factors[index], known_vars, target, factor)
    elseif expr isa Expr && expr.head === :call && expr.args[1] === :/ && length(expr.args) == 3
        numerator, denominator = expr.args[2:end]
        _contains_variable(denominator, known_vars) &&
            throw(ArgumentError("variable denominator in affine expression: $expr"))
        return lower_affine_expr(numerator, known_vars, target, :($scale / $denominator))
    elseif expr isa Expr && expr.head === :if && length(expr.args) == 3
        cond, a, b = expr.args
        _contains_variable(cond, known_vars) &&
            throw(ArgumentError("variable-dependent condition in affine expression: $expr"))
        terms_a, constant_a = lower_affine_expr(a, known_vars, target, scale)
        terms_b, constant_b = lower_affine_expr(b, known_vars, target, scale)
        terms = [Expr(:if, cond, Expr(:block, terms_a...), Expr(:block, terms_b...))]
        constant = constant_a == 0.0 && constant_b == 0.0 ? 0.0 : Expr(:if, cond, constant_a, constant_b)
        return terms, constant
    elseif expr isa Expr && expr.head === :call && expr.args[1] === :sum &&
           length(expr.args) == 2 && expr.args[2] isa Expr && expr.args[2].head === :generator
        gen = expr.args[2]
        raw_iters = gen.args[2:end]
        any(a -> _contains_variable(a, known_vars), raw_iters) &&
            throw(ArgumentError("variable-dependent iterator in affine expression: $expr"))
        terms, constant = lower_affine_expr(gen.args[1], known_vars, target, scale)
        if length(raw_iters) == 1 && raw_iters[1] isa Expr && raw_iters[1].head === :filter
            filter = raw_iters[1]
            body = Expr(:if, filter.args[1], Expr(:block, terms...))
            iters = filter.args[2:end]
        else
            body = Expr(:block, terms...)
            iters = raw_iters
        end
        foreach(_parse_iterator, iters)
        iter_spec = length(iters) == 1 ? iters[1] : Expr(:block, iters...)
        statements = isempty(terms) ? [] : [Expr(:for, iter_spec, body)]
        # The explicit Float64 init also handles empty filtered generators.
        offset = constant == 0.0 ? 0.0 :
            Expr(:call, :sum, Expr(:parameters, Expr(:kw, :init, 0.0)),
                 Expr(:generator, :(Float64($constant)), raw_iters...))
        return statements, offset
    end
    throw(ArgumentError("unsupported affine expression: $expr"))
end

lower_linear_expr(expr, known_vars::Set{Symbol}, target=:eq, sign=1.0) =
    first(lower_affine_expr(expr, known_vars, target, sign))

function _parse_iterator(arg)
    if arg isa Expr && arg.head === :call && length(arg.args) == 3 && arg.args[1] in (:in, :∈)
        name, space = arg.args[2:3]
    elseif arg isa Expr && arg.head in (:(=), :kw) && length(arg.args) == 2
        name, space = arg.args
    else
        throw(ArgumentError("expected an explicit iterator such as i in data.indices, got: $arg"))
    end
    name isa Symbol || throw(ArgumentError("iterator names must be single symbols, got: $name"))
    return name, space
end

"""
    parse_domain_ref(ref_expr) -> (name::Symbol, iterators, filter_cond)

Parse une expression d'indexation telle que `volumes[r in R, t in T; r > 0]` ou `vB`.
"""
function parse_domain_ref(expr)
    expr isa Symbol && return expr, [], nothing
    expr isa Expr && expr.head in (:ref, :typed_vcat) && length(expr.args) > 1 ||
        throw(ArgumentError("expected a scalar name or indexed family with explicit iterators, got: $expr"))
    name = expr.args[1]
    name isa Symbol || throw(ArgumentError("family name must be a symbol, got: $name"))
    args = expr.args[2:end]
    filter_cond = nothing
    if expr.head === :typed_vcat
        length(args) == 2 || throw(ArgumentError("use one filter expression after the semicolon"))
        filter_cond = pop!(args)
    else
        filters = filter(a -> a isa Expr && a.head === :parameters, args)
        length(filters) <= 1 || throw(ArgumentError("use one filter expression after the semicolon"))
        if !isempty(filters)
            length(only(filters).args) == 1 || throw(ArgumentError("combine family filters with &&"))
            filter_cond = only(only(filters).args)
            filter!(a -> !(a isa Expr && a.head === :parameters), args)
        end
    end
    isempty(args) && throw(ArgumentError("indexed families need at least one iterator"))
    iterators = [_parse_iterator(arg) for arg in args]
    names = first.(iterators)
    length(unique(names)) == length(names) || throw(ArgumentError("duplicate iterator name in $expr"))
    return name, iterators, filter_cond
end

"""
    parse_variable_bounds(decl_expr) -> (ref_expr, lower_expr, upper_expr)

Parse variable declaration and its bounds:
- `x[indices] in [lower, upper]`
- `lower <= x[indices] <= upper`
- `x[indices] >= lower`
- `x[indices]` -> `(-Inf, Inf)`
"""
function parse_variable_bounds(expr)
    # Form 1: var[...] in [lower, upper] or var[...] in bounds_tuple_expr
    if expr isa Expr && expr.head === :call && length(expr.args) == 3 && expr.args[1] in (:in, :∈)
        ref_expr = expr.args[2]
        bounds = expr.args[3]
        if bounds isa Expr && bounds.head in (:vect, :tuple)
            length(bounds.args) == 2 || throw(ArgumentError("variable bounds need exactly two endpoints"))
            return ref_expr, bounds.args[1], bounds.args[2]
        else
            return ref_expr, bounds, nothing
        end
    # Form 2: lower <= var[...] <= upper
    elseif expr isa Expr && expr.head === :comparison && length(expr.args) == 5
        l, op1, ref_expr, op2, u = expr.args
        if op1 in (:<=, :≤) && op2 in (:<=, :≤)
            return ref_expr, l, u
        end
    # Form 3: var[...] >= lower
    elseif expr isa Expr && expr.head === :call && length(expr.args) == 3 && expr.args[1] in (:>=, :≥)
        ref_expr = expr.args[2]
        return ref_expr, expr.args[3], :Inf
    # Form 4: var[...] <= upper
    elseif expr isa Expr && expr.head === :call && length(expr.args) == 3 && expr.args[1] in (:<=, :≤)
        ref_expr = expr.args[2]
        return ref_expr, -Inf, expr.args[3]
    # Form 5: no bounds specified
    elseif expr isa Symbol || (expr isa Expr && expr.head === :ref)
        return expr, -Inf, Inf
    end
    throw(ArgumentError("Unsupported @variable declaration syntax: $expr"))
end
