"""
A small risk-neutral SDDP implementation using ReactiveModel stage LPs.

Run from the repository root:
    julia --startup-file=no --project=. examples/sddp.jl

Three stages, one reservoir, independent two-point inflows after stage one.
Forward passes sample states; backward passes enumerate the finite support.
A four-path policy evaluation and a nonanticipative scenario-tree LP verify
bounds. This teaching example is not an SDDP.jl optimizer integration.
"""
module ReactiveModelSDDPExample
using ReactiveModel
using Random

const COSTS = (1.0, 2.0, 4.0)
const INFLOWS = (0.0, 2.0)
const PROBABILITIES = (0.5, 0.5)
const CAPACITY = 6.0
const DEMAND = 4.0
const INITIAL_STOCK = 4.0

@model HydroStage(data, point) begin
    @variable stock in [0.0, data.capacity]
    @variable hydro in [0.0, data.demand]
    @variable thermal in [0.0, data.demand]
    @variable spill in [0.0, Inf]
    @variable theta in [0.0, data.terminal ? 0.0 : Inf]
    @constraint balance stock + hydro + spill == point.stock + point.inflow
    @constraint supply hydro + thermal == data.demand
    @constraint future[k in 1:data.cut_budget] theta >= data.slopes[k] * stock + data.intercepts[k]
    @objective Min data.cost * thermal + theta
end

# A separate formulation: one decision per observed history, not per full path.
# Sharing each history's stock and control enforces nonanticipativity.
@model ScenarioTree(data, point) begin
    @variable stock[n in data.nodes] in [0.0, data.capacity]
    @variable hydro[n in data.nodes] in [0.0, data.demand]
    @variable thermal[n in data.nodes] in [0.0, data.demand]
    @variable spill[n in data.nodes] in [0.0, Inf]
    @constraint balance[n in data.nodes] begin
        stock[n] + hydro[n] + spill[n] ==
        (n == 1 ? point.stock : stock[data.parent[n]]) + data.inflow[n]
    end
    @constraint supply[n in data.nodes] hydro[n] + thermal[n] == data.demand
    @objective Min sum(data.probability[n] * data.cost[n] * thermal[n] for n in data.nodes)
end

function stage_data(t, cut_budget)
    (cost=COSTS[t], capacity=CAPACITY, demand=DEMAND, terminal=t == 3,
     cut_budget=cut_budget, slopes=zeros(cut_budget), intercepts=zeros(cut_budget))
end

function make_stages(cut_budget)
    data = [stage_data(t, cut_budget) for t in 1:3]
    point = (stock=INITIAL_STOCK, inflow=0.0)
    problems = [LinearProblem(ConstData(HydroStage(d), d), point) for d in data]
    data, problems
end

function exact_value(first_stage, incoming, observed_inflow)
    parent = [0]
    stage = [first_stage]
    inflow = [observed_inflow]
    probability = [1.0]
    for t in first_stage:2
        parents = findall(==(t), stage)
        for n in parents, (xi, weight) in zip(INFLOWS, PROBABILITIES)
            push!(parent, n)
            push!(stage, t + 1)
            push!(inflow, xi)
            push!(probability, probability[n] * weight)
        end
    end
    data = (nodes=eachindex(parent), parent=parent, inflow=inflow,
            probability=probability, cost=[COSTS[t] for t in stage],
            capacity=CAPACITY, demand=DEMAND)
    point = (stock=Float64(incoming),)
    problem = LinearProblem(ConstData(ScenarioTree(data), data), point)
    objective_value(problem, point)
end

function forward_pass(problems, rng)
    states = zeros(3)
    incoming = INITIAL_STOCK
    for t in 1:3
        states[t] = incoming
        inflow = t == 1 ? 0.0 : rand(rng, INFLOWS)
        solution = solve(problems[t], (stock=incoming, inflow=inflow))
        incoming = primal(solution, :stock)
    end
    states
end

function backward_pass!(data, problems, states, iteration)
    for t in 3:-1:2
        incoming = states[t]
        value, slope = 0.0, 0.0
        for (xi, weight) in zip(INFLOWS, PROBABILITIES)
            solution = solve(problems[t], (stock=incoming, inflow=xi))
            value += weight * objective_value(solution)
            # Incoming stock appears on the RHS of the minimization balance.
            slope += weight * dual(solution, :balance)
        end
        # Insert a supporting cut of the expected approximate stage-t value
        # into stage t-1. This mutation requires rebuilding its numerical model.
        previous = data[t - 1]
        previous.slopes[iteration] = slope
        previous.intercepts[iteration] = value - slope * incoming
        rebuild!(problems[t - 1])
    end
    nothing
end

# Exact expectation of the CURRENT feasible policy, not a sampled cost estimate.
function policy_value(problems, t=1, incoming=INITIAL_STOCK)
    outcomes = t == 1 ? ((0.0, 1.0),) : zip(INFLOWS, PROBABILITIES)
    value = 0.0
    for (xi, weight) in outcomes
        solution = solve(problems[t], (stock=incoming, inflow=xi))
        immediate = COSTS[t] * primal(solution, :thermal)
        outgoing = primal(solution, :stock)
        continuation = t == 3 ? 0.0 : policy_value(problems, t + 1, outgoing)
        value += weight * (immediate + continuation)
    end
    value
end

function run_example(; max_iterations=40, seed=20261112, tolerance=1e-8, verbose=false)
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    tolerance >= 0 || throw(ArgumentError("tolerance must be nonnegative"))
    data, problems = make_stages(max_iterations)
    rng = Xoshiro(seed)
    exact = exact_value(1, INITIAL_STOCK, 0.0)
    history = NamedTuple{(:iteration, :lower, :upper),Tuple{Int,Float64,Float64}}[]
    previous_lower = -Inf
    for iteration in 1:max_iterations
        states = forward_pass(problems, rng)
        backward_pass!(data, problems, states, iteration)
        lower = objective_value(problems[1], (stock=INITIAL_STOCK, inflow=0.0))
        upper = policy_value(problems)
        @assert previous_lower <= lower + 1e-7
        @assert lower <= exact + 1e-7
        @assert exact <= upper + 1e-7
        push!(history, (; iteration, lower, upper))
        previous_lower = lower
        upper - lower <= tolerance && break
    end
    final = last(history)
    @assert final.upper - final.lower <= tolerance + 1e-7
    @assert isapprox(final.lower, exact; atol=1e-7)

    # Check learned underestimators against independent scenario-tree LPs.
    for t in 2:3, incoming in 0.0:0.5:CAPACITY
        truth = sum(weight * exact_value(t, incoming, xi)
                    for (xi, weight) in zip(INFLOWS, PROBABILITIES))
        previous = data[t - 1]
        approximation = max(0.0, maximum(previous.slopes .* incoming .+ previous.intercepts))
        @assert approximation <= truth + 1e-7
    end
    # Away from a kink, confirm the cut-dual sign by a finite difference.
    terminal = problems[3]
    p = (stock=1.0, inflow=0.0)
    terminal_solution = solve(terminal, p)
    terminal_slope = dual(terminal_solution, :balance)
    delta = 1e-4
    finite_difference = (objective_value(terminal, (stock=1.0 + delta, inflow=0.0)) -
                         objective_value(terminal, (stock=1.0 - delta, inflow=0.0))) / (2delta)
    @assert isapprox(terminal_slope, -4.0; atol=1e-7)
    @assert isapprox(terminal_slope, finite_difference; atol=1e-7)
    result = (; exact, lower=final.lower, upper=final.upper,
              iterations=final.iteration, seed, history, data)
    if verbose
        println("iteration | lower bound | exact policy cost")
        for row in history
            println(row.iteration, " | ", row.lower, " | ", row.upper)
        end
        println("Scenario-tree optimum: ", exact)
        println("Cut and sensitivity checks passed")
    end
    result
end
end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ReactiveModelSDDPExample.run_example(; verbose=true)
end
