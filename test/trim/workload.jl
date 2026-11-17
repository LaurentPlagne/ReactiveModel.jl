module ReactiveModelTrimWorkload

using ReactiveModel

@model TypedLP(data, point) begin
    @variable x[i in 1:data.n] in [0.0, data.capacity[i] + point.extra] update_when :extra
    @variable anchor in [0.0, 0.0]
    @variable label[s in data.labels] in [0.0, 0.0]
    @variable pair[k in data.pairs] in [0.0, 0.0]
    @variable hole[i in data.empty] in [0.0, 1.0]
    @constraint floor[i in 1:data.n] point.a * x[i] >= point.rhs update_lhs_when :a update_rhs_when (
        :rhs,
    )
    @constraint zero anchor >= 0.0
    @constraint named[s in data.labels] label[s] == 0.0
    @constraint paired[k in data.pairs] pair[k] == 0.0
    @constraint no_rows[i in data.empty] x[i] <= 0.0
    @objective Min sum(data.cost[i] * point.cost * x[i] for i in 1:data.n) update_when (
        (previous, current) -> !isequal(previous.cost, current.cost)
    )
end

@model TypedMILP(data, point) begin
    @variable units in [0.0, point.limit] domain=:integer
    @variable active in [0.0, 1.0] domain=:binary
    @constraint floor units >= point.minimum
    @constraint activation units - data.bigM * active <= 0.0
    @objective Min units + 0.25 * active
end

function lp_problem(n::Int; backend=HiGHSBackend())
    data = (
        n=n,
        capacity=fill(10.0, n),
        cost=ones(n),
        labels=[:first, :second],
        pairs=[(:first, 1), (:second, 2)],
        empty=Int[],
    )
    point = (a=1.0, rhs=1.0, cost=2.0, extra=0.0)
    return LinearProblem(ConstData(TypedLP(data), data), point; backend)
end

function verify_lp(n::Int; backend=HiGHSBackend())
    lp = lp_problem(n; backend)
    p1 = (a=1.0, rhs=1.0, cost=2.0, extra=0.0)
    p2 = (a=2.0, rhs=3.0, cost=4.0, extra=1.0)
    objective_value(lp, p1) == 2.0n || return false
    sol = solve(lp, p2)
    objective_value(sol) == 6.0n || return false
    primal(sol, :x, 1) == 1.5 || return false
    primal(sol, :label, :first) == 0.0 || return false
    primal(sol, :pair, (:first, 1)) == 0.0 || return false
    activity(sol, :floor, 1) == 3.0 || return false
    dual(sol, :floor, 1) == 2.0 || return false
    saved = sol.variables.x
    sol.variables.label[:first] == 0.0 || return false
    sol.variables.pair[(:second, 2)] == 0.0 || return false
    sol.constraints.floor[1] == 2.0 || return false
    isempty(sol.variables.hole.values) || return false
    isempty(sol.constraints.no_rows.activities) || return false
    dual(sol, :floor, 1) == 2.0 || return false

    # Rebuild at the same point after mutating a dynamic study array.
    lp.const_data.data.cost .= 2.0
    rebuild!(lp)
    objective_value(lp, p2) == 12.0n || return false
    p3 = (a=2.0, rhs=4.0, cost=4.0, extra=1.0)
    objective_value(lp, p3) == 16.0n || return false
    saved[1] == 1.5 || return false
    independent = similar(lp)
    objective_value(independent, p3) == 16.0n || return false

    bad = (a=2.0, rhs=1000.0, cost=4.0, extra=1.0)
    failed = false
    try
        objective_value(lp, bad)
    catch error
        error isa SolveError || rethrow()
        failed = error.status === :infeasible
    end
    failed || return false
    return objective_value(lp, p2) == 12.0n
end

function verify_milp()
    data = (bigM=10.0,)
    lp = LinearProblem(ConstData(TypedMILP(data), data), (minimum=1.5, limit=10.0))
    objective_value(lp, (minimum=1.5, limit=10.0)) == 2.25 || return false
    sol = solve(lp, (minimum=2.5, limit=9.0))
    objective_value(sol) == 3.25 || return false
    primal(sol, :units) == 3.0 || return false
    primal(sol, :active) == 1.0 || return false
    return !backend_has_duals(lp.backend) && isnan(dual(sol, :floor))
end

function functional_problem(n::Int)
    definition = ProblemDefinition(;
        variables=(x=VariableBlock(:x, 1:n; lower=0.0, upper=10.0),),
        constraints=(
            floor=ConstraintBlock(
                :floor, 1:n, :ge, (eq, d, p, i) -> set_coeff!(eq, :x, i, p.a), (d, p, i) -> p.rhs
            ),
        ),
        objective=Objective(:min, (eq, d, p) -> begin
            for i in 1:n
                set_coeff!(eq, :x, i, p.cost)
            end
        end),
    )
    return LinearProblem(definition, (a=1.0, rhs=1.0, cost=2.0))
end

function verify_functional(n::Int)
    lp = functional_problem(n)
    objective_value(lp, (a=1.0, rhs=1.0, cost=2.0)) == 2.0n || return false
    return objective_value(lp, (a=2.0, rhs=3.0, cost=4.0)) == 6.0n
end

# The publication regressions are part of the compiled acceptance workload.
struct CollidingKey
    id::Int
end
Base.hash(::CollidingKey, seed::UInt) = hash(0, seed)

@model AffineAliases(data, point) begin
    selected = data.values[point.t]
    factor = selected / 2.0
    @variable x[k in data.keys] in [0.0, selected + 4.0 * length(data.keys) + 5.0]
    @constraint balance[k in data.keys] factor * (x[k] + selected) == x[k] / 2.0 + selected * k.id
    @objective Min sum(x[k] + selected for k in data.keys)
end

module OwnerA
    using ReactiveModel
    @model Shared(data, point) begin
        @variable x in [0.0, 1.0]
        @constraint cap x <= 1.0
        @objective Max x + 10.0
    end
end
module OwnerB
    using ReactiveModel
    @model Shared(data, point) begin
        @variable x in [0.0, 9.0]
        @constraint cap x <= 9.0
        @objective Min -x + 20.0
    end
end

@model SignedZeroConditions(data, point) begin
    @variable inferred in [0.0, signbit(point.zero) ? 2.0 : 1.0]
    @variable explicit in [0.0, signbit(point.zero) ? 2.0 : 1.0] update_when :zero
    @objective Max inferred + explicit
end

function verify_math(n::Int)
    data = (values=[2.0, 4.0], keys=[CollidingKey(i + 2) for i in 1:n])
    lp = LinearProblem(ConstData(AffineAliases(data), data), (t=1,))
    for t in (1, 2, 1)
        selected = data.values[t]
        sol = solve(lp, (t=t,))
        expected = 0.0
        for k in data.keys
            x = (selected * k.id - selected^2 / 2.0) / (selected / 2.0 - 0.5)
            expected += x + selected
            isapprox(primal(sol, :x, k), x) || return false
        end
        isapprox(objective_value(sol), expected) || return false
        isapprox(objective_value(lp, (t=t,)), expected) || return false
    end
    signed = LinearProblem(ConstData(SignedZeroConditions(nothing)), (zero=0.0,))
    objective_value(signed, (zero=0.0,)) == 2.0 || return false
    objective_value(signed, (zero=-0.0,)) == 4.0 || return false
    objective_value(signed, (zero=0.0,)) == 2.0 || return false
    a = LinearProblem(ConstData(OwnerA.Shared(nothing)), nothing)
    b = LinearProblem(ConstData(OwnerB.Shared(nothing)), nothing)
    objective_value(a, nothing) == 11.0 || return false
    return objective_value(b, nothing) == 11.0
end

@model StructuredDSL(key, eq) begin
    @variable anchor in [0.0, 0.0]
    @variable x[i = 1:key.n, j in 1:i; i > 0 && j > 0] in [0.0, eq.limit]
    @constraint fixed[i in 1:key.n, j in 1:i] x[i, j] == eq.value
    @objective Min anchor + sum(sum(x[i, j] + eq.shift for j in 1:i) for i in 1:key.n)
end

function verify_grammar(n::Int)
    data = (n=n,)
    p = (limit=10.0, value=2.0, shift=1.0)
    lp = LinearProblem(ConstData(StructuredDSL(data), data), p)
    objective_value(lp, p) == 3.0*n*(n+1)/2 || return false
    q = (limit=9.0, value=3.0, shift=2.0)
    objective_value(lp, q) == 5.0*n*(n+1)/2 || return false
    empty = LinearProblem(ConstData(StructuredDSL((n=0,)), (n=0,)), p)
    return objective_value(empty, p) == 0.0
end

@model CheckedCallbacks(data, point) begin
    @variable x[i in 1:data.n] in [0.0, point.limit] update_when false
    @constraint floor[i in 1:data.n] point.factor*x[i] >= point.rhs update_when false
    @objective Min sum(point.cost*x[i] for i in 1:data.n) + point.offset update_when false
end

function verify_safe_updates(n::Int)
    data=(n=n,)
    p=(limit=10.0, factor=1.0, rhs=1.0, cost=2.0, offset=0.0)
    lp=LinearProblem(ConstData(CheckedCallbacks(data), data), p; safe_mode=true)
    objective_value(lp, p)==2.0n || return false
    q=(limit=9.0, factor=2.0, rhs=3.0, cost=4.0, offset=5.0)
    objective_value(lp, q)==6.0n+5.0 || return false
    checks=update_report(lp)
    length(checks)==4 || return false
    all(c->c.checks==1 && c.missed==1, checks) || return false
    reset_update_report!(lp)
    isempty(update_report(lp)) || return false
    objective_value(lp, q)==6.0n+5.0 || return false
    return all(c->c.checks==1 && c.missed==0, update_report(lp))
end

@model GroupedRules(data, point) begin
    @variable a[i in 1:data.n] in [0.0, point.limit]
    @variable b[i in 1:data.n] in [0.0, point.limit]
    @variable c[i in 1:data.n] in [0.0, point.limit]
    @variable d[i in 1:data.n] in [0.0, point.limit]
    @variable e[i in 1:data.n] in [0.0, point.limit]
    @variable f[i in 1:data.n] in [0.0, point.limit]
    @variable g[i in 1:data.n] in [0.0, point.limit]
    @variable h[i in 1:data.n] in [0.0, point.limit]
    @variable j[i in 1:data.n] in [0.0, point.limit]
    @constraint ca[i in 1:data.n] a[i] >= point.rhs
    @constraint cb[i in 1:data.n] b[i] >= point.rhs
    @constraint cc[i in 1:data.n] c[i] >= point.rhs
    @constraint cd[i in 1:data.n] d[i] >= point.rhs
    @constraint ce[i in 1:data.n] e[i] >= point.rhs
    @constraint cf[i in 1:data.n] f[i] >= point.rhs
    @constraint cg[i in 1:data.n] g[i] >= point.rhs
    @constraint ch[i in 1:data.n] h[i] >= point.rhs
    @constraint cj[i in 1:data.n] j[i] >= point.rhs
    @objective Min sum(a[i]+b[i]+c[i]+d[i]+e[i]+f[i]+g[i]+h[i]+j[i] for i in 1:data.n)
end

function verify_grouped_rules(n::Int)
    data=(n=n,)
    p=(limit=10.0, rhs=1.0)
    lp=LinearProblem(ConstData(GroupedRules(data), data), p)
    objective_value(lp, p)==9.0n || return false
    sol=solve(lp, (limit=9.0, rhs=2.0))
    objective_value(sol)==18.0n || return false
    return primal(sol, :j, n)==2.0 && lp.numerical_updates.variable_count==9
end

@noinline function checked_execution_rhs(data, point)
    point.fail && throw(DomainError(point.rhs, "execution callback probe"))
    return data.scale*point.rhs
end

@model CompactExecution(data, point) begin
    @variable x in [0.0, data.capacity]
    @constraint floor x >= checked_execution_rhs(data, point)
    @objective Min x
end

function stale_scalar_properties(solution)
    objective_failed = false
    status_failed = false
    try
        solution.objective_value
    catch error
        error isa ErrorException || rethrow()
        objective_failed = true
    end
    try
        solution.termination_status
    catch error
        error isa ErrorException || rethrow()
        status_failed = true
    end
    return objective_failed && status_failed
end

function verify_rebuild_recovery(n::Int)
    data = (scale=Float64(n), capacity=10.0*n)
    p = (rhs=1.0, fail=false)
    bad = (rhs=1.0, fail=true)
    lp = LinearProblem(ConstData(CompactExecution(data), data), p)
    solution = solve(lp, p)
    saved = solution.variables
    independent = similar(lp)
    failed = false
    try
        rebuild!(lp, (scale=2.0*n, capacity=10.0*n); parameters=bad)
    catch error
        error isa DomainError || rethrow()
        failed = true
    end
    failed && !lp.built && lp.parameter_snapshot === nothing || return false
    lp.workspace.last_status === :not_solved || return false
    stale_scalar_properties(solution) || return false
    objective_value(lp, p) == 2n || return false
    solve(lp, p) === solution || return false
    solution.objective_value == 2n || return false
    solution.termination_status === :optimal || return false
    saved.x[0] == n || return false
    return objective_value(independent, p) == n
end

function verify_execution(n::Int)
    data=(scale=Float64(n), capacity=10.0*n)
    p=(rhs=1.0, fail=false)
    lp=LinearProblem(ConstData(CompactExecution(data), data), p; safe_mode=true)
    lp.workspace.result_layout===nothing || return false
    GC.gc()
    sol=solve(lp, p)
    objective_value(sol)==n || return false
    caught=false
    try
        objective_value(lp, (rhs=1.0, fail=true))
    catch error
        error isa DomainError || rethrow()
        caught=true
    end
    caught && !lp.built || return false
    stale_scalar_properties(sol) || return false
    lp.numerical_updates.argument.value===nothing || return false
    objective_value(lp, p)==n || return false
    independent=similar(lp)
    replacement=similar(lp;data=(scale=3.0*n,capacity=10.0*n))
    objective_value(replacement,p)==3n || return false
    objective_value(lp,p)==n || return false
    rebuild!(lp, (scale=2.0*n, capacity=10.0*n))
    objective_value(lp, p)==2n || return false
    objective_value(independent, p)==n || return false
    return primal(sol, :x)==2n && activity(sol, :floor)==2n && dual(sol, :floor)==1.0
end

@model ManyResultKeys(data, point) begin
    @variable i[k in data.ints] in [point.minimum, 10.0]
    @variable s[k in data.symbols] in [point.minimum, 10.0]
    @variable t[k in data.tuples] in [point.minimum, 10.0]
    @variable w[k in data.words] in [point.minimum, 10.0]
    @variable c[k in data.chars] in [point.minimum, 10.0]
    @objective Min sum(i[k] for k in data.ints) + sum(s[k] for k in data.symbols) +
        sum(t[k] for k in data.tuples) + sum(w[k] for k in data.words) + sum(c[k] for k in data.chars)
end

function many_key_problem(n::Int)
    data=(ints=[n], symbols=[:first], tuples=[(n,1)], words=["word"], chars=['a'])
    return LinearProblem(ConstData(ManyResultKeys(data), data), (minimum=1.0,))
end
function verify_many_result_keys(n::Int)
    problem=many_key_problem(n)
    solution=solve(problem, (minimum=1.0,))
    total=solution.variables.i[n] + solution.variables.s[:first] + solution.variables.t[(n,1)] +
        solution.variables.w["word"] + solution.variables.c['a']
    total==5.0 || return false
    saved=solution.variables.w
    objective_value(problem, (minimum=2.0,))==10.0 || return false
    return saved["word"]==1.0 && solution.variables.t[(n,1)]==2.0
end

@model ComposedVariables(data, point) begin
    @variable x in [1.0, 10.0]
end

@model ComposedPrimary(data, point) begin
    @external_variables x
    @objective Min data.scale[1] * point.a * x + point.offset
end

@model ComposedSecondary(data, point) begin
    @external_variables x
    @objective Min point.b * x + 2.0 * point.offset
end

function composed_problem(n::Int)
    data=(scale=[Float64(n)],)
    base=ComposedVariables(data)
    objective=compose_objectives(ComposedPrimary(data).objective,
        ComposedSecondary(data).objective,
        Objective(:min,(emit,d,p)->set_coeff!(emit,:x,1.0);
            has_changed=(old,new)->false))
    definition=ProblemDefinition(variables=base.variables,constraints=base.constraints,
        objective=objective)
    return LinearProblem(ConstData(definition,data),(a=2.0,b=3.0,offset=1.0))
end

function verify_composed_objective(n::Int)
    lp=composed_problem(n)
    p=(a=2.0,b=3.0,offset=1.0)
    q=(a=3.0,b=4.0,offset=2.0)
    objective_value(lp,p)==2n+7 || return false
    objective_value(lp,q)==3n+11 || return false
    lp.objective_offset[]==6.0 || return false
    independent=similar(lp;data=(scale=copy(lp.const_data.data.scale),))
    lp.const_data.data.scale[1]=2.0*n
    rebuild!(lp)
    objective_value(lp,q)==6n+11 || return false
    objective_value(independent,q)==3n+11 || return false
    replacement=similar(lp;data=(scale=[3.0*n],))
    objective_value(replacement,q)==9n+11 || return false
    rebuild!(lp,(scale=[4.0*n],))
    return objective_value(lp,p)==8n+7
end

numeric_bounds_pair(data,point)=(point.lo,point.hi)

@model NumericBounds(data,point) begin
    @variable x[i in 1:data.n] in numeric_bounds_pair(data,point)
    @constraint floor[i in 1:data.n] x[i]>=-100.0
    @objective Max sum(point.cost*x[i] for i in 1:data.n)
end

struct CompactNumericBounds{T} <: AbstractVariable
    limit::T
end
ReactiveModel.variable_bounds(v::CompactNumericBounds)=(zero(v.limit),v.limit)

function verify_numeric_bounds(n::Int)
    data=(n=n,)
    p=(lo=0,hi=4,cost=1.0)
    lp=LinearProblem(ConstData(NumericBounds(data),data),p)
    objective_value(lp,p)==4.0n || return false
    q=(lo=1,hi=3,cost=-1.0)
    objective_value(lp,q)==-1.0n || return false
    objective_value(lp,p)==4.0n || return false
    clone=similar(lp)
    rebuild!(lp;parameters=q)
    objective_value(lp,q)==-1.0n || return false
    objective_value(clone,p)==4.0n || return false

    f=(lo=-Inf32,hi=2.5f0,cost=1.0)
    floats=LinearProblem(ConstData(NumericBounds(data),data),f)
    objective_value(floats,f)==2.5n || return false
    objective_value(floats,(lo=0.0f0,hi=Inf32,cost=-1.0))==0.0 || return false

    definition=ProblemDefinition(
        variables=(x=VariableBlock(:x,1:n,CompactNumericBounds(Float32(n))),),
        constraints=NamedTuple(),
        objective=Objective(:max,(e,d,p)->foreach(i->set_coeff!(e,:x,i,1.0),1:d.n)),
    )
    compact=LinearProblem(ConstData(definition,data),nothing)
    return objective_value(compact,nothing)==Float64(n*n)
end

function verify(n::Int)
    return verify_lp(n) &&
           verify_milp() &&
           verify_functional(n) &&
           verify_math(n) &&
           verify_grammar(n) &&
           verify_safe_updates(n) &&
           verify_grouped_rules(n) &&
           verify_execution(n) &&
           verify_rebuild_recovery(n) &&
           verify_numeric_bounds(n) &&
           verify_composed_objective(n) &&
           verify_many_result_keys(n)
end

end # module
