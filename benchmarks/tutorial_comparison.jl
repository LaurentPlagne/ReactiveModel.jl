"""
Measured warm-update comparisons for the network-flow and reservoir tutorials.

Run from the repository root with:

    julia --startup-file=no --project=benchmarks \
        benchmarks/tutorial_comparison.jl

ReactiveModel uses HiGHSBackend; JuMP uses HiGHS.Optimizer. Each sample
alternates between two parameter points and includes update plus optimization.
"""

using BenchmarkTools
using JuMP
using HiGHS

import ReactiveModel
include(joinpath(@__DIR__, "..", "examples", "examples.jl"))
using .ReactiveModelExamples

const NETWORK_ARCS = [(:north, :centre), (:centre, :south), (:north, :south)]
const NETWORK_NODES = [:north, :centre, :south]
const NETWORK_DATA = transport_const_data(
    NETWORK_ARCS,
    NETWORK_NODES,
    [100.0 100.0 30.0; 120.0 120.0 40.0],
    [1.0 1.0 3.0; 1.2 0.8 3.0],
    [80.0 0.0 -80.0; 100.0 0.0 -100.0],
)

const RESERVOIR_DATA = reservoir_const_data(
    [20.0, 30.0],
    100.0,
    [50.0, 50.0, 50.0],
    [100.0, 100.0, 100.0],
    [1.0, 1.0, 1.0],
    [10.0 5.0 0.0; 15.0 5.0 0.0],
    [20.0 40.0 30.0; 25.0 35.0 45.0],
    [1.0, 1.0, 1.0],
    [5.0, 2.0],
    [0.0, 100.0],
)

mutable struct ReactiveTutorialState{L,P}
    model::L
    first::P
    second::P
    use_second::Bool
end

function reactive_network_state()
    first = TransportParameters(1)
    model = ReactiveModel.LinearProblem(
        NETWORK_DATA,
        first;
        backend=ReactiveModel.HiGHSBackend(output_flag=false),
    )
    ReactiveTutorialState(model, first, TransportParameters(2), false)
end

function reactive_reservoir_state()
    first = ReservoirParameters(1)
    model = ReactiveModel.LinearProblem(
        RESERVOIR_DATA,
        first;
        backend=ReactiveModel.HiGHSBackend(output_flag=false),
    )
    ReactiveTutorialState(model, first, ReservoirParameters(2), false)
end

function reactive_update!(state::ReactiveTutorialState)
    state.use_second = !state.use_second
    ReactiveModel.objective_value(
        state.model,
        state.use_second ? state.second : state.first,
    )
end

mutable struct JuMPNetworkState
    model::Model
    flow::Vector{VariableRef}
    conservation::Vector{ConstraintRef}
    use_second::Bool
end

function jump_network_state()
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, flow[1:length(NETWORK_ARCS)] >= 0)
    for a in eachindex(NETWORK_ARCS)
        set_upper_bound(flow[a], NETWORK_DATA.capacity[1, a])
    end
    conservation = @constraint(
        model,
        [n in eachindex(NETWORK_NODES)],
        sum((NETWORK_ARCS[a][1] == NETWORK_NODES[n] ? 1.0 :
             NETWORK_ARCS[a][2] == NETWORK_NODES[n] ? -1.0 : 0.0) * flow[a]
            for a in eachindex(NETWORK_ARCS)) == NETWORK_DATA.injection[1, n],
    )
    @objective(
        model,
        Min,
        sum(NETWORK_DATA.cost[1, a] * flow[a] for a in eachindex(NETWORK_ARCS)),
    )
    JuMPNetworkState(model, collect(flow), collect(conservation), false)
end

function jump_network_update!(state::JuMPNetworkState)
    state.use_second = !state.use_second
    time_step = state.use_second ? 2 : 1
    for a in eachindex(state.flow)
        set_upper_bound(state.flow[a], NETWORK_DATA.capacity[time_step, a])
        set_objective_coefficient(state.model, state.flow[a], NETWORK_DATA.cost[time_step, a])
    end
    for n in eachindex(state.conservation)
        set_normalized_rhs(state.conservation[n], NETWORK_DATA.injection[time_step, n])
    end
    optimize!(state.model)
    objective_value(state.model)
end

mutable struct JuMPReservoirState
    model::Model
    stock::Vector{VariableRef}
    turbine::Vector{VariableRef}
    water_balance::Vector{ConstraintRef}
    use_second::Bool
end

function jump_reservoir_state()
    data = RESERVOIR_DATA.data
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, 0 <= stock[0:data.horizon] <= data.capacity)
    @variable(model, 0 <= turbine[t=1:data.horizon] <= data.turbine_max[t])
    @variable(model, 0 <= spill[t=1:data.horizon] <= data.spill_max[t])
    @variable(model, future_value)
    set_lower_bound(stock[0], data.initial_stock[1])
    set_upper_bound(stock[0], data.initial_stock[1])
    water_balance = @constraint(
        model,
        [t=1:data.horizon],
        stock[t] - stock[t - 1] + data.duration[t] * (turbine[t] + spill[t]) ==
            data.duration[t] * data.inflow[1, t],
    )
    @constraint(
        model,
        [k=1:data.cuts],
        future_value - data.slope[k] * stock[data.horizon] <= data.intercept[k],
    )
    @objective(
        model,
        Max,
        sum(data.duration[t] * data.price[1, t] * data.yield_factor[t] * turbine[t]
            for t in 1:data.horizon) + future_value,
    )
    JuMPReservoirState(
        model,
        [stock[t] for t in 0:data.horizon],
        collect(turbine),
        collect(water_balance),
        false,
    )
end

function jump_reservoir_update!(state::JuMPReservoirState)
    data = RESERVOIR_DATA.data
    state.use_second = !state.use_second
    scenario = state.use_second ? 2 : 1
    set_lower_bound(state.stock[1], data.initial_stock[scenario])
    set_upper_bound(state.stock[1], data.initial_stock[scenario])
    for t in 1:data.horizon
        set_normalized_rhs(
            state.water_balance[t],
            data.duration[t] * data.inflow[scenario, t],
        )
        set_objective_coefficient(
            state.model,
            state.turbine[t],
            data.duration[t] * data.price[scenario, t] * data.yield_factor[t],
        )
    end
    optimize!(state.model)
    objective_value(state.model)
end

function trial(f, state; samples=100)
    result = @benchmark $f($state) samples=samples evals=1
    estimate = BenchmarkTools.median(result)
    (median_ns=estimate.time, alloc_bytes=estimate.memory, allocations=estimate.allocs)
end

function run_tutorial_comparison(; samples=100)
    reactive_network = reactive_network_state()
    jump_network = jump_network_state()
    optimize!(jump_network.model)
    reactive_network_value = reactive_update!(reactive_network)
    jump_network_value = jump_network_update!(jump_network)
    @assert isapprox(reactive_network_value, jump_network_value; atol=1e-8)

    reactive_reservoir = reactive_reservoir_state()
    jump_reservoir = jump_reservoir_state()
    optimize!(jump_reservoir.model)
    reactive_reservoir_value = reactive_update!(reactive_reservoir)
    jump_reservoir_value = jump_reservoir_update!(jump_reservoir)
    @assert isapprox(reactive_reservoir_value, jump_reservoir_value; atol=1e-8)

    (
        network=(
            reactive=trial(reactive_update!, reactive_network; samples),
            jump=trial(jump_network_update!, jump_network; samples),
        ),
        reservoir=(
            reactive=trial(reactive_update!, reactive_reservoir; samples),
            jump=trial(jump_reservoir_update!, jump_reservoir; samples),
        ),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    display(run_tutorial_comparison())
end
