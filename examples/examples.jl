module ReactiveModelExamples

using ..ReactiveModel
import ..ReactiveModel: lower_bound, upper_bound, lhs_terms!, rhs_value, objective_terms!, set_coeff!, @coeff

export DispatchConstData, DispatchParameters, dispatch_const_data, dispatch_definition,
       TransportConstData, TransportParameters, transport_const_data, transport_definition,
       ReservoirConstData, ReservoirParameters, reservoir_const_data, reservoir_definition

"Study data shared by every dispatch solve. Rows are time steps."
struct DispatchConstData
    generator_indices::UnitRange{Int}
    demand::Vector{Float64}
    capacity::Matrix{Float64}
    cost::Matrix{Float64}
end

"Small per-solve selector for the dispatch study."
struct DispatchParameters
    time_step::Int
end

struct DispatchProduction <: AbstractVariable end
lower_bound(::DispatchProduction, c::DispatchConstData, p::DispatchParameters, g) = 0.0
upper_bound(::DispatchProduction, c::DispatchConstData, p::DispatchParameters, g) = c.capacity[p.time_step, g]

struct DispatchBalance <: AbstractConstraint end
function lhs_terms!(::DispatchBalance, equation, c::DispatchConstData, p::DispatchParameters, key)
    for g in c.generator_indices
        @coeff equation production[g] ← 1.0
    end
end
rhs_value(::DispatchBalance, c::DispatchConstData, p::DispatchParameters, key) = c.demand[p.time_step]

struct DispatchObjective <: AbstractObjective end
function objective_terms!(::DispatchObjective, expression, c::DispatchConstData, p::DispatchParameters)
    for g in c.generator_indices
        @coeff expression production[g] ← c.cost[p.time_step, g]
    end
end

"""
    dispatch_definition(number_of_generators)

Build the fixed dispatch formulation. The complete time series lives in
`ConstData`; `Parameters.time_step` only selects the row used by one solve.
"""
function dispatch_definition(number_of_generators)
    variables = (production=VariableBlock(DispatchProduction, 1:number_of_generators),)
    constraints = (balance=ConstraintBlock(:eq, DispatchBalance),)
    objective = Objective(:min, DispatchObjective)
    ProblemDefinition(; variables, constraints, objective)
end

function dispatch_const_data(demand, capacity, cost)
    demand = Float64.(collect(demand))
    capacity = Matrix{Float64}(capacity)
    cost = Matrix{Float64}(cost)
    isempty(demand) && throw(ArgumentError("dispatch study has no time step"))
    size(capacity, 1) == length(demand) ||
        throw(DimensionMismatch("capacity rows must match demand"))
    size(cost) == size(capacity) ||
        throw(DimensionMismatch("cost and capacity must have the same shape"))
    number_of_generators = size(capacity, 2)
    number_of_generators > 0 || throw(ArgumentError("dispatch study has no generator"))
    data = DispatchConstData(1:number_of_generators, demand, capacity, cost)
    ConstData(dispatch_definition(number_of_generators), data)
end

"Fixed topology of a directed capacitated transport problem."
struct TransportConstData{A,N}
    arcs::Vector{A}
    nodes::Vector{N}
    arc_positions::Dict{A,Int}
    node_positions::Dict{N,Int}
    outgoing::Dict{N,Vector{A}}
    incoming::Dict{N,Vector{A}}
    capacity::Matrix{Float64}
    cost::Matrix{Float64}
    injection::Matrix{Float64}
end

"Small per-solve selector for the transport study."
struct TransportParameters
    pdt::Int
end

struct TransportFlow <: AbstractVariable end
lower_bound(::TransportFlow, c::TransportConstData, p::TransportParameters, arc) = 0.0
upper_bound(::TransportFlow, c::TransportConstData, p::TransportParameters, arc) = c.capacity[p.pdt, c.arc_positions[arc]]

struct TransportConservation <: AbstractConstraint end
function lhs_terms!(::TransportConservation, equation, c::TransportConstData, p::TransportParameters, node)
    for arc in c.outgoing[node]
        @coeff equation flow[arc] ← 1.0
    end
    for arc in c.incoming[node]
        @coeff equation flow[arc] ← -1.0
    end
end
rhs_value(::TransportConservation, c::TransportConstData, p::TransportParameters, node) =
    c.injection[p.pdt, c.node_positions[node]]

struct TransportObjective <: AbstractObjective end
function objective_terms!(::TransportObjective, expression, c::TransportConstData, p::TransportParameters)
    for (a, arc) in enumerate(c.arcs)
        @coeff expression flow[arc] ← c.cost[p.pdt, a]
    end
end

"""
    transport_definition(arcs, nodes)

Minimize `sum(cost[pdt, a] * flow[a] for a)`, with
`0 <= flow[a] <= capacity[pdt, a]`. For every node `n`, outgoing flow minus
incoming flow equals `injection[pdt, n]`. All arrays belong to `ConstData`.
"""
function transport_definition(arcs, nodes)
    arcs = collect(arcs)
    nodes = collect(nodes)
    variables = (flow=VariableBlock(TransportFlow, arcs),)
    constraints = (conservation=ConstraintBlock(nodes, :eq, TransportConservation),)
    objective = Objective(:min, TransportObjective)
    ProblemDefinition(; variables, constraints, objective)
end

function transport_const_data(arcs, nodes, capacity, cost, injection)
    arcs = collect(arcs)
    nodes = collect(nodes)
    capacity = Matrix{Float64}(capacity)
    cost = Matrix{Float64}(cost)
    injection = Matrix{Float64}(injection)
    isempty(arcs) && throw(ArgumentError("transport study has no arc"))
    isempty(nodes) && throw(ArgumentError("transport study has no node"))
    size(capacity, 2) == length(arcs) ||
        throw(DimensionMismatch("capacity columns must match arcs"))
    size(cost) == size(capacity) ||
        throw(DimensionMismatch("cost and capacity must have the same shape"))
    size(injection) == (size(capacity, 1), length(nodes)) ||
        throw(DimensionMismatch("injection shape must be (time steps, nodes)"))
    data = TransportConstData(
        arcs,
        nodes,
        Dict(arc => i for (i, arc) in enumerate(arcs)),
        Dict(node => i for (i, node) in enumerate(nodes)),
        Dict(node => [arc for arc in arcs if first(arc) == node] for node in nodes),
        Dict(node => [arc for arc in arcs if last(arc) == node] for node in nodes),
        capacity,
        cost,
        injection,
    )
    ConstData(transport_definition(arcs, nodes), data)
end

"Study data shared by every finite-horizon reservoir solve."
struct ReservoirConstData
    horizon::Int
    cuts::Int
    initial_stock::Vector{Float64}
    capacity::Float64
    turbine_max::Vector{Float64}
    spill_max::Vector{Float64}
    duration::Vector{Float64}
    inflow::Matrix{Float64}
    price::Matrix{Float64}
    yield_factor::Vector{Float64}
    slope::Vector{Float64}
    intercept::Vector{Float64}
end

"Small per-solve selector for the reservoir study."
struct ReservoirParameters
    scenario::Int
end

struct ReservoirStock <: AbstractVariable end
lower_bound(::ReservoirStock, c::ReservoirConstData, p::ReservoirParameters, t) = t == 0 ? c.initial_stock[p.scenario] : 0.0
upper_bound(::ReservoirStock, c::ReservoirConstData, p::ReservoirParameters, t) =
    t == 0 ? c.initial_stock[p.scenario] : c.capacity

struct ReservoirTurbine <: AbstractVariable end
lower_bound(::ReservoirTurbine, c::ReservoirConstData, p::ReservoirParameters, t) = 0.0
upper_bound(::ReservoirTurbine, c::ReservoirConstData, p::ReservoirParameters, t) = c.turbine_max[t]

struct ReservoirSpill <: AbstractVariable end
lower_bound(::ReservoirSpill, c::ReservoirConstData, p::ReservoirParameters, t) = 0.0
upper_bound(::ReservoirSpill, c::ReservoirConstData, p::ReservoirParameters, t) = c.spill_max[t]

struct ReservoirFutureValue <: AbstractVariable end
lower_bound(::ReservoirFutureValue, c::ReservoirConstData, p::ReservoirParameters, key) = -Inf
upper_bound(::ReservoirFutureValue, c::ReservoirConstData, p::ReservoirParameters, key) = Inf

struct ReservoirWaterBalance <: AbstractConstraint end
function lhs_terms!(::ReservoirWaterBalance, equation, c::ReservoirConstData, p::ReservoirParameters, t)
    @coeff equation stock[t] ← 1.0
    @coeff equation stock[t - 1] ← -1.0
    @coeff equation turbine[t] ← c.duration[t]
    @coeff equation spill[t] ← c.duration[t]
end
rhs_value(::ReservoirWaterBalance, c::ReservoirConstData, p::ReservoirParameters, t) =
    c.duration[t] * c.inflow[p.scenario, t]

struct ReservoirTerminalCut <: AbstractConstraint end
function lhs_terms!(::ReservoirTerminalCut, equation, c::ReservoirConstData, p::ReservoirParameters, k)
    @coeff equation future_value ← 1.0
    @coeff equation stock[c.horizon] ← -c.slope[k]
end
rhs_value(::ReservoirTerminalCut, c::ReservoirConstData, p::ReservoirParameters, k) = c.intercept[k]

struct ReservoirObjective <: AbstractObjective end
function objective_terms!(::ReservoirObjective, expression, c::ReservoirConstData, p::ReservoirParameters)
    for t in 1:c.horizon
        @coeff expression turbine[t] ← c.duration[t] * c.price[p.scenario, t] * c.yield_factor[t]
    end
    @coeff expression future_value ← 1.0
end

"""
    reservoir_definition(horizon, cuts)

The water balance is
`stock[t] - stock[t-1] + duration[t]*(turbine[t] + spill[t]) =
 duration[t]*inflow[t]`.
The terminal cuts are
`future_value <= slope[k] * stock[horizon] + intercept[k]`.
The objective maximizes turbine revenue plus `future_value`.
"""
function reservoir_definition(horizon, cuts)
    variables = (
        stock=VariableBlock(ReservoirStock, 0:horizon),
        turbine=VariableBlock(ReservoirTurbine, 1:horizon),
        spill=VariableBlock(ReservoirSpill, 1:horizon),
        future_value=VariableBlock(ReservoirFutureValue),
    )
    constraints = (
        water_balance=ConstraintBlock(1:horizon, :eq, ReservoirWaterBalance),
        terminal_cut=ConstraintBlock(1:cuts, :le, ReservoirTerminalCut),
    )
    objective = Objective(:max, ReservoirObjective)
    ProblemDefinition(; variables, constraints, objective)
end

function reservoir_const_data(initial_stock, capacity, turbine_max, spill_max,
                              duration, inflow, price, yield_factor, slope, intercept)
    initial_stock = Float64.(collect(initial_stock))
    turbine_max = Float64.(collect(turbine_max))
    spill_max = Float64.(collect(spill_max))
    duration = Float64.(collect(duration))
    inflow = Matrix{Float64}(inflow)
    price = Matrix{Float64}(price)
    yield_factor = Float64.(collect(yield_factor))
    slope = Float64.(collect(slope))
    intercept = Float64.(collect(intercept))
    horizon = length(duration)
    cuts = length(slope)
    horizon > 0 || throw(ArgumentError("reservoir study has an empty horizon"))
    cuts > 0 || throw(ArgumentError("reservoir study has no terminal cut"))
    length(turbine_max) == horizon || throw(DimensionMismatch("turbine_max must match horizon"))
    length(spill_max) == horizon || throw(DimensionMismatch("spill_max must match horizon"))
    length(yield_factor) == horizon || throw(DimensionMismatch("yield_factor must match horizon"))
    size(inflow) == (length(initial_stock), horizon) ||
        throw(DimensionMismatch("inflow shape must be (scenarios, horizon)"))
    size(price) == size(inflow) || throw(DimensionMismatch("price and inflow must have the same shape"))
    length(intercept) == cuts || throw(DimensionMismatch("intercept must match slope"))
    data = ReservoirConstData(
        horizon,
        cuts,
        initial_stock,
        Float64(capacity),
        turbine_max,
        spill_max,
        duration,
        inflow,
        price,
        yield_factor,
        slope,
        intercept,
    )
    ConstData(reservoir_definition(horizon, cuts), data)
end

end # module
