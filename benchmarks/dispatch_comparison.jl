"""
Reproducible front-end comparison for the dispatch example.

Run from the repository root with:

    julia --startup-file=no --project=benchmarks \
        benchmarks/dispatch_comparison.jl

The benchmark deliberately separates cold construction, steady-state cached
solve, and a warm update. Solver time is included in all three timings; the
`@allocated` report is the Julia-side allocation count and is not a claim about
native allocations inside a solver.
"""

using BenchmarkTools
using JuMP
using HiGHS
import MathOptInterface
import ParametricOptInterface as POI
const MOI = MathOptInterface

import ReactiveModel
include(joinpath(@__DIR__, "..", "examples", "examples.jl"))
using .ReactiveModelExamples

const DISPATCH_DATA = dispatch_const_data(
    [100.0, 110.0],
    [80.0 60.0; 75.0 65.0],
    [40.0 70.0; 45.0 68.0],
)
const P = DispatchParameters(1)
const P2 = DispatchParameters(2)

@inline dispatch_demand(p::DispatchParameters) = DISPATCH_DATA.demand[p.time_step]
@inline dispatch_capacity(p::DispatchParameters, g) = DISPATCH_DATA.capacity[p.time_step, g]
@inline dispatch_cost(p::DispatchParameters, g) = DISPATCH_DATA.cost[p.time_step, g]

mutable struct ReactiveState{C,P,B}
    lp::ReactiveModel.LinearProblem{C,P,B}
    use_second::Bool
end

function reactive_state(backend=ReactiveModel.ClpBackend())
    lp = ReactiveModel.LinearProblem(
        DISPATCH_DATA, ReactiveModel.Parameters(P);
        backend,
    )
    ReactiveState(lp, false)
end

function reactive_cold(backend=ReactiveModel.ClpBackend())
    state = reactive_state(backend)
    ReactiveModel.objective_value(ReactiveModel.solve(state.lp, P))
end

function reactive_cached!(state::ReactiveState)
    ReactiveModel.objective_value(ReactiveModel.solve(state.lp, P))
end

function reactive_warm_update!(state::ReactiveState)
    state.use_second = !state.use_second
    ReactiveModel.objective_value(state.lp, state.use_second ? P2 : P)
end

function reactive_lazy_solution_warm_update!(state::ReactiveState)
    state.use_second = !state.use_second
    solution = ReactiveModel.solve(state.lp, state.use_second ? P2 : P)
    ReactiveModel.objective_value(solution)
end

function reactive_materialized_warm_update!(state::ReactiveState)
    state.use_second = !state.use_second
    solution = ReactiveModel.solve(state.lp, state.use_second ? P2 : P)
    solution.variables.production[1]
    solution.constraints.balance[nothing]
    solution
end

mutable struct JuMPState
    model::Model
    production::Vector{VariableRef}
    balance::ConstraintRef
    use_second::Bool
end

function jump_state()
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    for (name, value) in (("parallel", "off"),
                          ("solver", "simplex"), ("simplex_max_concurrency", 1))
        set_optimizer_attribute(model, name, value)
    end
    @variable(model, production[1:2] >= 0)
    for g in 1:2
        set_upper_bound(production[g], dispatch_capacity(P, g))
    end
    balance = @constraint(model, sum(production) == dispatch_demand(P))
    @objective(model, Min, sum(dispatch_cost(P, g) * production[g] for g in 1:2))
    JuMPState(model, collect(production), balance, false)
end

function jump_apply!(state::JuMPState, p::DispatchParameters)
    for g in 1:2
        set_upper_bound(state.production[g], dispatch_capacity(p, g))
        set_objective_coefficient(state.model, state.production[g], dispatch_cost(p, g))
    end
    set_normalized_rhs(state.balance, dispatch_demand(p))
    optimize!(state.model)
    objective_value(state.model)
end

function jump_cold()
    state = jump_state()
    optimize!(state.model)
    state
end

function jump_cached!(state::JuMPState)
    optimize!(state.model)
    state
end

function jump_warm_update!(state::JuMPState)
    state.use_second = !state.use_second
    jump_apply!(state, state.use_second ? P2 : P)
end

mutable struct MOIDispatchState
    optimizer::Any
    variables::Vector{MOI.VariableIndex}
    lower_bounds::Vector{MOI.ConstraintIndex}
    upper_bounds::Vector{MOI.ConstraintIndex}
    balance::MOI.ConstraintIndex
    objective_coefficients::Vector{Float64}
    use_second::Bool
end

function _affine(variables, coefficients)
    MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(c, v) for (c, v) in zip(coefficients, variables)], 0.0,
    )
end

function moi_state()
    optimizer = HiGHS.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("output_flag"), false)
    for (name, value) in (("parallel", "off"),
                          ("solver", "simplex"), ("simplex_max_concurrency", 1))
        MOI.set(optimizer, MOI.RawOptimizerAttribute(name), value)
    end
    variables = MOI.VariableIndex[MOI.add_variable(optimizer) for _ in 1:2]
    lower = [MOI.add_constraint(optimizer, v, MOI.GreaterThan(0.0)) for v in variables]
    upper = [MOI.add_constraint(optimizer, v, MOI.LessThan(dispatch_capacity(P, g))) for (g, v) in enumerate(variables)]
    balance = MOI.add_constraint(optimizer, _affine(variables, [1.0, 1.0]), MOI.EqualTo(dispatch_demand(P)))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    coefficients = [dispatch_cost(P, g) for g in 1:2]
    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), _affine(variables, coefficients))
    MOIDispatchState(optimizer, variables, lower, upper, balance, coefficients, false)
end

function moi_apply!(state::MOIDispatchState, p::DispatchParameters)
    for g in 1:2
        MOI.set(state.optimizer, MOI.ConstraintSet(), state.upper_bounds[g], MOI.LessThan(dispatch_capacity(p, g)))
    end
    MOI.set(state.optimizer, MOI.ConstraintSet(), state.balance, MOI.EqualTo(dispatch_demand(p)))
    for i in eachindex(state.objective_coefficients)
        value = dispatch_cost(p, i)
        if !isequal(state.objective_coefficients[i], value)
            MOI.modify(state.optimizer,
                MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                MOI.ScalarCoefficientChange(state.variables[i], value))
            state.objective_coefficients[i] = value
        end
    end
    MOI.optimize!(state.optimizer)
    MOI.get(state.optimizer, MOI.ObjectiveValue())
end

function moi_cold()
    state = moi_state()
    MOI.optimize!(state.optimizer)
    state
end

moi_cached!(state::MOIDispatchState) = (MOI.optimize!(state.optimizer); state)

function moi_warm_update!(state::MOIDispatchState)
    state.use_second = !state.use_second
    moi_apply!(state, state.use_second ? P2 : P)
end

mutable struct POIDispatchState
    model::Model
    production::Vector{VariableRef}
    capacity::Vector{VariableRef}
    demand::VariableRef
    cost::Vector{VariableRef}
    use_second::Bool
end

function poi_state()
    model = Model(() -> POI.Optimizer(HiGHS.Optimizer))
    set_silent(model)
    for (name, value) in (("parallel", "off"),
                          ("solver", "simplex"), ("simplex_max_concurrency", 1))
        set_optimizer_attribute(model, name, value)
    end
    @variable(model, production[1:2] >= 0)
    @variable(model, capacity[g=1:2] in MOI.Parameter(dispatch_capacity(P, g)))
    @variable(model, demand in MOI.Parameter(dispatch_demand(P)))
    @variable(model, cost[g=1:2] in MOI.Parameter(dispatch_cost(P, g)))
    @constraint(model, capacity_constraint[g=1:2], production[g] <= capacity[g])
    @constraint(model, balance, sum(production) == demand)
    @objective(model, Min, sum(cost[g] * production[g] for g in 1:2))
    POIDispatchState(model, collect(production), collect(capacity), demand, collect(cost), false)
end

function poi_apply!(state::POIDispatchState, p::DispatchParameters)
    for g in 1:2
        set_parameter_value(state.capacity[g], dispatch_capacity(p, g))
        set_parameter_value(state.cost[g], dispatch_cost(p, g))
    end
    set_parameter_value(state.demand, dispatch_demand(p))
    optimize!(state.model)
    objective_value(state.model)
end

function poi_cold()
    state = poi_state()
    optimize!(state.model)
    state
end

poi_cached!(state::POIDispatchState) = (optimize!(state.model); state)

function poi_warm_update!(state::POIDispatchState)
    state.use_second = !state.use_second
    poi_apply!(state, state.use_second ? P2 : P)
end

function _trial(f, state; samples=100)
    # One evaluation per sample: each call contains one complete model action.
    trial = @benchmark $f($state) samples=samples evals=1
    (median_ns=BenchmarkTools.median(trial).time,
     alloc_bytes=BenchmarkTools.median(trial).allocs == 0 ? 0 : BenchmarkTools.median(trial).memory,
     allocations=BenchmarkTools.median(trial).allocs)
end

function _cold_trial(f; samples=30)
    trial = @benchmark $f() samples=samples evals=1
    estimate = BenchmarkTools.median(trial)
    (median_ns=estimate.time, alloc_bytes=estimate.memory, allocations=estimate.allocs)
end

function type_report(reactive, jump, moi, poi)
    reactive_lazy = reactive_state()
    ReactiveModel.objective_value(ReactiveModel.solve(reactive_lazy.lp, P))
    (
        reactive_problem=string(typeof(reactive.lp)),
        reactive_backend=string(typeof(reactive.lp.backend)),
        reactive_const_data=string(typeof(reactive.lp.const_data)),
        reactive_workspace=string(typeof(reactive.lp.workspace)),
        reactive_variable_block=string(typeof(reactive.lp.const_data.definition.variables.blocks[1])),
        jump_model=string(typeof(jump.model)),
        jump_variables=string(typeof(jump.production)),
        moi_optimizer=string(typeof(moi.optimizer)),
        moi_variables=string(typeof(moi.variables)),
        reactive_lazy_type_size=Base.summarysize(reactive_lazy.lp),
        reactive_materialized_type_size=Base.summarysize(reactive.lp),
        jump_type_size=Base.summarysize(jump.model),
        moi_type_size=Base.summarysize(moi.optimizer),
        poi_model=string(typeof(poi.model)),
        poi_parameter=string(typeof(poi.demand)),
        poi_type_size=Base.summarysize(poi.model),
    )
end

function check_dispatch_parity!(reactive, reactive_highs, jump, moi, poi)
    for expected in (5755.0, 4600.0)
        values = (
            reactive_warm_update!(reactive),
            reactive_warm_update!(reactive_highs),
            jump_warm_update!(jump),
            moi_warm_update!(moi),
            poi_warm_update!(poi),
        )
        all(isapprox(value, expected; atol=1e-8, rtol=0) for value in values) ||
            error("Dispatch comparator mismatch: expected $expected, got $values")
    end
    nothing
end

function run_dispatch_comparison(; samples=100)
    # Warm up compilation before measuring steady-state paths.
    reactive = reactive_state(); ReactiveModel.solve(reactive.lp, P)
    reactive_highs = reactive_state(
        ReactiveModel.HiGHSBackend(output_flag=false),
    ); ReactiveModel.solve(reactive_highs.lp, P)
    jump = jump_state(); optimize!(jump.model)
    moi = moi_state(); MOI.optimize!(moi.optimizer)
    poi = poi_state(); optimize!(poi.model)
    check_dispatch_parity!(reactive, reactive_highs, jump, moi, poi)

    report = (
        cold=(reactive=_cold_trial(reactive_cold; samples=min(samples, 30)),
              reactive_highs=_cold_trial(
                  () -> reactive_cold(ReactiveModel.HiGHSBackend(output_flag=false));
                  samples=min(samples, 30),
              ),
              jump=_cold_trial(jump_cold; samples=min(samples, 30)),
              moi=_cold_trial(moi_cold; samples=min(samples, 30)),
              poi=_cold_trial(poi_cold; samples=min(samples, 30))),
        reactive=(cached=_trial(reactive_cached!, reactive; samples),
                  warm_update=_trial(reactive_warm_update!, reactive; samples),
                  lazy_solution_warm_update=_trial(
                      reactive_lazy_solution_warm_update!, reactive; samples),
                  materialized_warm_update=_trial(reactive_materialized_warm_update!, reactive; samples)),
        reactive_highs=(
            cached=_trial(reactive_cached!, reactive_highs; samples),
            warm_update=_trial(reactive_warm_update!, reactive_highs; samples),
            lazy_solution_warm_update=_trial(
                reactive_lazy_solution_warm_update!, reactive_highs; samples),
        ),
        jump=(cached=_trial(jump_cached!, jump; samples),
              warm_update=_trial(jump_warm_update!, jump; samples)),
        moi=(cached=_trial(moi_cached!, moi; samples),
             warm_update=_trial(moi_warm_update!, moi; samples)),
        poi=(cached=_trial(poi_cached!, poi; samples),
             warm_update=_trial(poi_warm_update!, poi; samples)),
        types=type_report(reactive, jump, moi, poi),
    )
    report
end

if abspath(PROGRAM_FILE) == @__FILE__
    display(run_dispatch_comparison())
end
