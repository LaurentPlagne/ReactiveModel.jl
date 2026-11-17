"""
Frontend comparison on independently specified synthetic cascade LPs.

    julia --startup-file=no --project=benchmarks -t 8 \
        benchmarks/valley_comparison.jl [preset] [nscenarios] [order]

Presets are `small`, `medium`, `large` or `all`. Inputs, topology and
experimental coefficients are specified in `synthetic_valley.jl` and
`PROVENANCE.md`. Frontends share numerical inputs; solver entry points and
adapter behavior can still differ. Validate outputs before interpreting timings.
"""

using Statistics: median
using JuMP
using HiGHS
import MathOptInterface as MOI
import ParametricOptInterface as POI
include(joinpath(@__DIR__, "optional_solvers.jl"))

import ReactiveModel
import ReactiveModel as RM
using ReactiveModel: @model

using ReactiveModel: ClpBackend
const _CLP_LIB = RM._CLP_LIB

include(joinpath(@__DIR__, "synthetic_valley.jl"))
import .SyntheticValley as SV
using .SyntheticValley: ValleyData, Point

# ---------------------------------------------------------------------------
# ReactiveModel formulation (Declarative @model DSL)
# ---------------------------------------------------------------------------

# The bounds of the selected concave region, rescaled onto each cut reservoir's
# usable volume, exactly as `SyntheticValley.row_rhs!` does for the flat model.
region_scale(data::ValleyData, r::Int) =
    data.volume_max[r, data.spec.npostes] / data.volume_max[1, data.spec.npostes]

@model SyntheticValleyDSL(data, sel) begin
    spec = data.spec
    NR, NT, NC = spec.nreservoirs, spec.npostes, spec.ncut_reservoirs
    NU = length(spec.unit_groups)
    R = 1:NR
    T = 1:NT
    U = 1:NU
    E = data.energetic

    # 1. Variables with reactive bounds
    @variable initial_volumes[r in R] in [SV.initial_volume_bound(data, sel, r), SV.initial_volume_bound(data, sel, r)] update_when ((previous, current) -> SV.stocks_changed(previous, current))
    @variable volumes[r in R, t in T] in SV.volume_bounds(data, sel, r, t) update_when false
    @variable debits[u in U, t in T] in SV.flow_bounds(data, sel, u, t) update_when ((previous, current) -> SV.stocks_changed(previous, current))
    @variable group_debits[u in E, g in 1:data.ngroups[u], t in T] in SV.group_flow_bounds(data, sel, u, g, t) update_when ((previous, current) -> SV.stocks_changed(previous, current))
    @variable vB in [-Inf, Inf] update_when false

    # 2. Constraints with selective dirty-block tracking
    @constraint water_balance[r in R, t in T] begin
        (t == 1 ? -1.0 * initial_volumes[r] : -1.0 * volumes[r, t - 1]) +
        1.0 * volumes[r, t] +
        sum(data.duration[t] * debits[u, t] for u in U if spec.unit_upstream[u] == r) -
        sum(data.duration[t] * debits[u, t] for u in U if spec.unit_downstream[u] == r) ==
        SV.inflow_value(data, sel, r, t)
    end update_lhs_when false update_rhs_when ((previous, current) -> SV.scenario_changed(previous, current))

    @constraint final_volume[r in (NC + 1):NR] begin
        1.0 * volumes[r, NT] == SV.final_volume_target(data, sel, r)
    end update_lhs_when false update_rhs_when ((previous, current) -> SV.stocks_changed(previous, current))

    @constraint group_sum[u in E, t in T] begin
        -1.0 * debits[u, t] + sum(1.0 * group_debits[u, g, t] for g in 1:data.ngroups[u]) == 0.0
    end update_when false

    @constraint bellman_cut[s in 1:spec.ncuts] begin
        1.0 * vB - sum(SV.cut_slope(data, s, r) * volumes[r, NT] for r in 1:NC) <= SV.cut_rhs(data, sel, s)
    end update_lhs_when false update_rhs_when ((previous, current) -> SV.region_changed(previous, current))

    @constraint inactive_cut begin
        1.0 * volumes[1, NT] <= Inf
    end update_when false

    @constraint min_final[r in 1:NC] begin
        1.0 * volumes[r, NT] >= (data.volume_max[r, NT] / data.volume_max[1, NT]) * SV.min_final_stock(data, sel)
    end update_lhs_when false update_rhs_when ((previous, current) -> SV.region_changed(previous, current))

    @constraint max_final[r in 1:NC] begin
        1.0 * volumes[r, NT] <= (data.volume_max[r, NT] / data.volume_max[1, NT]) * SV.max_final_stock(data, sel)
    end update_lhs_when false update_rhs_when ((previous, current) -> SV.region_changed(previous, current))

    @constraint inactive_triangle[k in 1:(NC * NC)] begin
        1.0 * volumes[mod1(k, NC), NT] <= Inf
    end update_when false

    # 3. Objective with reactive market profits
    @objective Max begin
        1.0 * vB + sum(SV.profit(data, sel, u, g, t) * group_debits[u, g, t] for t in T, u in E, g in 1:data.ngroups[u])
    end update_when ((previous, current) -> SV.profits_changed(previous, current))
end

"""
    reactive_const_data(data)

The `ConstData` of the synthetic valley, constructed declaratively via the `@model` DSL.
"""
function reactive_const_data(data::ValleyData)
    RM.ConstData(SyntheticValleyDSL(data), data)
end

# ---------------------------------------------------------------------------
# Flat front-ends: JuMP, direct MOI, POI
# ---------------------------------------------------------------------------

const SOLVER_OPTIONS = (("threads", 1), ("parallel", "off"),
                        ("solver", "simplex"), ("simplex_max_concurrency", 1))

"""
    apply_solver_options!(model)

Set the HiGHS options that pin it to a single-threaded simplex. TinyHiGHS has no
such raw attributes, so each is applied only where it is supported.
"""
function apply_solver_options!(model)
    for (name, value) in SOLVER_OPTIONS
        attribute = MOI.RawOptimizerAttribute(name)
        MOI.supports(JuMP.backend(model), attribute) || continue
        JuMP.set_optimizer_attribute(model, name, value)
    end
    return nothing
end

function apply_solver_options!(optimizer::MOI.AbstractOptimizer)
    for (name, value) in SOLVER_OPTIONS
        attribute = MOI.RawOptimizerAttribute(name)
        MOI.supports(optimizer, attribute) || continue
        MOI.set(optimizer, attribute, value)
    end
    return nothing
end

"""
    OPTIMIZERS

The two solvers the MOI-based front-ends can be built on.

**TinyHiGHS has no warm start through MathOptInterface.** Its wrapper rebuilds
the LP and a fresh `SimplexEngine` on every `MOI.optimize!`, so a JuMP, MOI or
POI model backed by it solves cold every time, whatever the modeling layer does
to update it incrementally. Only a direct backend — ReactiveModel's
`TinyHiGHSBackend` — can retain its basis. That asymmetry is the
point of measuring the full matrix, not a defect of the measurement.
"""
const OPTIMIZERS = BENCH_TINY ? (highs=HiGHS.Optimizer, tiny=TinyHiGHS.Optimizer) : (highs=HiGHS.Optimizer,)

"""
    moving_entries(data, layout, points)

Scan the measured point sequence once and report which bounds, right-hand sides
and objective coefficients actually move. The POI models need this before they
are built; it also documents how much of the LP a warm update really touches.
"""
function moving_entries(data::ValleyData, l::SV.Layout, points)
    lower = fill(NaN, l.ncols); upper = fill(NaN, l.ncols)
    rhs = fill(NaN, l.nrows); costs = fill(NaN, length(l.objective_columns))
    lower_moves = falses(l.ncols); upper_moves = falses(l.ncols)
    rhs_moves = falses(l.nrows); cost_moves = falses(length(costs))
    nl = similar(lower); nu = similar(upper); nr = similar(rhs); nc = similar(costs)
    for (index, point) in enumerate(points)
        SV.column_bounds!(nl, nu, data, l, point)
        SV.row_rhs!(nr, data, l, point)
        SV.objective_costs!(nc, data, l, point)
        if index > 1
            for j in 1:l.ncols
                isfinite(nl[j]) == isfinite(lower[j]) &&
                isfinite(nu[j]) == isfinite(upper[j]) ||
                    error("bound finiteness changes at column $j")
                lower_moves[j] |= !isequal(nl[j], lower[j])
                upper_moves[j] |= !isequal(nu[j], upper[j])
            end
            for j in 1:l.nrows
                isfinite(nr[j]) == isfinite(rhs[j]) ||
                    error("right-hand side finiteness changes at row $j")
                rhs_moves[j] |= !isequal(nr[j], rhs[j])
            end
            for i in eachindex(costs)
                cost_moves[i] |= !isequal(nc[i], costs[i])
            end
        end
        copyto!(lower, nl); copyto!(upper, nu)
        copyto!(rhs, nr); copyto!(costs, nc)
    end
    (; lower_moves, upper_moves, rhs_moves, cost_moves)
end

jump_set(op, value) = op === :eq ? MOI.EqualTo(value) :
    op === :le ? MOI.LessThan(value) :
    op === :ge ? MOI.GreaterThan(value) : error("unsupported operator $op")

"""
    jump_reference(data, layout, point; optimizer, direct)

`direct=true` builds a `JuMP.direct_model`, which drops JuMP's caching layer and
talks straight to the optimizer. It is the closest JuMP gets to bypassing
MathOptInterface — it does not bypass it, since MOI is the only language JuMP
speaks. It therefore cannot restore a warm start that the optimizer's own MOI
wrapper throws away; what it does measure is the cost of the cache layer.
"""
function jump_reference(data::ValleyData, l::SV.Layout, point::Point;
        optimizer=HiGHS.Optimizer, direct::Bool=false)
    model = direct ? JuMP.direct_model(optimizer()) : JuMP.Model(optimizer)
    JuMP.set_silent(model)
    apply_solver_options!(model)
    x = JuMP.@variable(model, [1:l.ncols])
    lower = Vector{Float64}(undef, l.ncols); upper = Vector{Float64}(undef, l.ncols)
    SV.column_bounds!(lower, upper, data, l, point)
    for j in 1:l.ncols
        isfinite(lower[j]) && JuMP.set_lower_bound(x[j], lower[j])
        isfinite(upper[j]) && JuMP.set_upper_bound(x[j], upper[j])
    end
    rhs = Vector{Float64}(undef, l.nrows)
    SV.row_rhs!(rhs, data, l, point)
    rowrefs = Vector{JuMP.ConstraintRef}(undef, l.nrows)
    for j in 1:l.nrows
        indices, coefficients, op = l.rows[j]
        expression = JuMP.AffExpr(0.0)
        for k in eachindex(indices)
            JuMP.add_to_expression!(expression, coefficients[k], x[indices[k]])
        end
        rowrefs[j] = JuMP.add_constraint(model,
            JuMP.ScalarConstraint(expression, jump_set(op, rhs[j])))
    end
    costs = Vector{Float64}(undef, length(l.objective_columns))
    SV.objective_costs!(costs, data, l, point)
    objective = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(objective, 1.0, x[l.vB])
    for i in eachindex(costs)
        JuMP.add_to_expression!(objective, costs[i], x[l.objective_columns[i]])
    end
    JuMP.set_objective_sense(model, MOI.MAX_SENSE)
    JuMP.set_objective_function(model, objective)
    (; model, x, rowrefs, lower, upper, rhs, costs,
        scratch_lower=similar(lower), scratch_upper=similar(upper),
        scratch_rhs=similar(rhs), scratch_costs=similar(costs),
        previous=Ref(point), data, l)
end

function jump_update!(state, point::Point)
    old = state.previous[]
    l, data = state.l, state.data
    if SV.stocks_changed(old, point)
        SV.column_bounds!(state.scratch_lower, state.scratch_upper, data, l, point)
        for j in 1:l.ncols
            if !isequal(state.lower[j], state.scratch_lower[j])
                JuMP.set_lower_bound(state.x[j], state.scratch_lower[j])
                state.lower[j] = state.scratch_lower[j]
            end
            if !isequal(state.upper[j], state.scratch_upper[j])
                JuMP.set_upper_bound(state.x[j], state.scratch_upper[j])
                state.upper[j] = state.scratch_upper[j]
            end
        end
    end
    if SV.rhs_changed(old, point)
        SV.row_rhs!(state.scratch_rhs, data, l, point)
        for j in 1:l.nrows
            if !isequal(state.rhs[j], state.scratch_rhs[j])
                JuMP.set_normalized_rhs(state.rowrefs[j], state.scratch_rhs[j])
                state.rhs[j] = state.scratch_rhs[j]
            end
        end
    end
    if SV.profits_changed(old, point)
        SV.objective_costs!(state.scratch_costs, data, l, point)
        for i in eachindex(state.costs)
            if !isequal(state.costs[i], state.scratch_costs[i])
                JuMP.set_objective_coefficient(state.model,
                    state.x[l.objective_columns[i]], state.scratch_costs[i])
                state.costs[i] = state.scratch_costs[i]
            end
        end
    end
    state.previous[] = point
    return nothing
end

function jump_solve!(state, point::Point)
    jump_update!(state, point)
    JuMP.optimize!(state.model)
    JuMP.termination_status(state.model) == MOI.OPTIMAL ||
        error("JuMP failed at $point: $(JuMP.termination_status(state.model))")
    JuMP.objective_value(state.model)
end

function moi_reference(data::ValleyData, l::SV.Layout, point::Point;
        optimizer_type=HiGHS.Optimizer)
    optimizer = optimizer_type()
    MOI.supports(optimizer, MOI.RawOptimizerAttribute("output_flag")) &&
        MOI.set(optimizer, MOI.RawOptimizerAttribute("output_flag"), false)
    MOI.supports(optimizer, MOI.Silent()) && MOI.set(optimizer, MOI.Silent(), true)
    apply_solver_options!(optimizer)
    x = MOI.add_variables(optimizer, l.ncols)
    lower = Vector{Float64}(undef, l.ncols); upper = Vector{Float64}(undef, l.ncols)
    SV.column_bounds!(lower, upper, data, l, point)
    lower_refs = Vector{Union{Nothing,MOI.ConstraintIndex}}(nothing, l.ncols)
    upper_refs = Vector{Union{Nothing,MOI.ConstraintIndex}}(nothing, l.ncols)
    for j in 1:l.ncols
        isfinite(lower[j]) &&
            (lower_refs[j] = MOI.add_constraint(optimizer, x[j], MOI.GreaterThan(lower[j])))
        isfinite(upper[j]) &&
            (upper_refs[j] = MOI.add_constraint(optimizer, x[j], MOI.LessThan(upper[j])))
    end
    rhs = Vector{Float64}(undef, l.nrows)
    SV.row_rhs!(rhs, data, l, point)
    rowrefs = Vector{MOI.ConstraintIndex}(undef, l.nrows)
    for j in 1:l.nrows
        indices, coefficients, op = l.rows[j]
        function_ = MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(coefficients[k], x[indices[k]])
             for k in eachindex(indices)], 0.0)
        rowrefs[j] = MOI.add_constraint(optimizer, function_, jump_set(op, rhs[j]))
    end
    costs = Vector{Float64}(undef, length(l.objective_columns))
    SV.objective_costs!(costs, data, l, point)
    terms = [MOI.ScalarAffineTerm(1.0, x[l.vB])]
    for i in eachindex(costs)
        push!(terms, MOI.ScalarAffineTerm(costs[i], x[l.objective_columns[i]]))
    end
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(terms, 0.0))
    (; optimizer, x, lower_refs, upper_refs, rowrefs, lower, upper, rhs, costs,
        scratch_lower=similar(lower), scratch_upper=similar(upper),
        scratch_rhs=similar(rhs), scratch_costs=similar(costs),
        previous=Ref(point), data, l)
end

function moi_update!(state, point::Point)
    old = state.previous[]
    l, data = state.l, state.data
    if SV.stocks_changed(old, point)
        SV.column_bounds!(state.scratch_lower, state.scratch_upper, data, l, point)
        for j in 1:l.ncols
            if state.lower_refs[j] !== nothing &&
                    !isequal(state.lower[j], state.scratch_lower[j])
                MOI.set(state.optimizer, MOI.ConstraintSet(), state.lower_refs[j],
                    MOI.GreaterThan(state.scratch_lower[j]))
                state.lower[j] = state.scratch_lower[j]
            end
            if state.upper_refs[j] !== nothing &&
                    !isequal(state.upper[j], state.scratch_upper[j])
                MOI.set(state.optimizer, MOI.ConstraintSet(), state.upper_refs[j],
                    MOI.LessThan(state.scratch_upper[j]))
                state.upper[j] = state.scratch_upper[j]
            end
        end
    end
    if SV.rhs_changed(old, point)
        SV.row_rhs!(state.scratch_rhs, data, l, point)
        for j in 1:l.nrows
            if !isequal(state.rhs[j], state.scratch_rhs[j])
                reference = state.rowrefs[j]
                set = MOI.get(state.optimizer, MOI.ConstraintSet(), reference)
                MOI.set(state.optimizer, MOI.ConstraintSet(), reference,
                    typeof(set)(state.scratch_rhs[j]))
                state.rhs[j] = state.scratch_rhs[j]
            end
        end
    end
    if SV.profits_changed(old, point)
        SV.objective_costs!(state.scratch_costs, data, l, point)
        for i in eachindex(state.costs)
            if !isequal(state.costs[i], state.scratch_costs[i])
                MOI.modify(state.optimizer,
                    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                    MOI.ScalarCoefficientChange(state.x[l.objective_columns[i]],
                        state.scratch_costs[i]))
                state.costs[i] = state.scratch_costs[i]
            end
        end
    end
    state.previous[] = point
    return nothing
end

function moi_solve!(state, point::Point)
    moi_update!(state, point)
    MOI.optimize!(state.optimizer)
    MOI.get(state.optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL ||
        error("direct MOI failed at $point")
    MOI.get(state.optimizer, MOI.ObjectiveValue())
end

"""
    poi_reference(data, layout, moving, point; parametric_bounds)

A `MOI.Parameter` may only appear inside an affine function, never as a variable
bound. A transition LP moves variable bounds from point to point, so a faithful
POI model has to restate every moving bound as a row, `x - lb >= 0` and
`x - ub <= 0`, with the parameter on the left. That is a strictly larger LP.

`parametric_bounds=true` builds that pure parametric model.
`parametric_bounds=false` builds the fairer best-effort hybrid: POI parameters
for the moving right-hand sides and objective costs, ordinary JuMP bound updates
for the bounds.
"""
function poi_reference(data::ValleyData, l::SV.Layout, moving, point::Point;
        parametric_bounds::Bool, optimizer=HiGHS.Optimizer)
    model = JuMP.Model(() -> POI.Optimizer(optimizer()))
    JuMP.set_silent(model)
    apply_solver_options!(model)
    x = JuMP.@variable(model, [1:l.ncols])
    lower = Vector{Float64}(undef, l.ncols); upper = Vector{Float64}(undef, l.ncols)
    SV.column_bounds!(lower, upper, data, l, point)
    lower_parameters = Vector{Union{Nothing,JuMP.VariableRef}}(nothing, l.ncols)
    upper_parameters = Vector{Union{Nothing,JuMP.VariableRef}}(nothing, l.ncols)
    for j in 1:l.ncols
        if isfinite(lower[j])
            if parametric_bounds && moving.lower_moves[j]
                parameter = JuMP.@variable(model, set = MOI.Parameter(lower[j]))
                JuMP.@constraint(model, x[j] - parameter >= 0.0)
                lower_parameters[j] = parameter
            else
                JuMP.set_lower_bound(x[j], lower[j])
            end
        end
        if isfinite(upper[j])
            if parametric_bounds && moving.upper_moves[j]
                parameter = JuMP.@variable(model, set = MOI.Parameter(upper[j]))
                JuMP.@constraint(model, x[j] - parameter <= 0.0)
                upper_parameters[j] = parameter
            else
                JuMP.set_upper_bound(x[j], upper[j])
            end
        end
    end
    rhs = Vector{Float64}(undef, l.nrows)
    SV.row_rhs!(rhs, data, l, point)
    rowrefs = Vector{JuMP.ConstraintRef}(undef, l.nrows)
    rhs_parameters = Vector{Union{Nothing,JuMP.VariableRef}}(nothing, l.nrows)
    for j in 1:l.nrows
        indices, coefficients, op = l.rows[j]
        expression = JuMP.AffExpr(0.0)
        for k in eachindex(indices)
            JuMP.add_to_expression!(expression, coefficients[k], x[indices[k]])
        end
        constant = rhs[j]
        if moving.rhs_moves[j]
            isfinite(rhs[j]) || error("a moving right-hand side must be finite; row $j")
            parameter = JuMP.@variable(model, set = MOI.Parameter(rhs[j]))
            JuMP.add_to_expression!(expression, -1.0, parameter)
            rhs_parameters[j] = parameter
            constant = 0.0
        end
        rowrefs[j] = JuMP.add_constraint(model,
            JuMP.ScalarConstraint(expression, jump_set(op, constant)))
    end
    costs = Vector{Float64}(undef, length(l.objective_columns))
    SV.objective_costs!(costs, data, l, point)
    cost_parameters = Vector{Union{Nothing,JuMP.VariableRef}}(nothing, length(costs))
    objective = JuMP.QuadExpr(JuMP.AffExpr(0.0))
    JuMP.add_to_expression!(objective, 1.0, x[l.vB])
    for i in eachindex(costs)
        variable = x[l.objective_columns[i]]
        if moving.cost_moves[i]
            parameter = JuMP.@variable(model, set = MOI.Parameter(costs[i]))
            JuMP.add_to_expression!(objective, 1.0, parameter, variable)
            cost_parameters[i] = parameter
        else
            JuMP.add_to_expression!(objective, costs[i], variable)
        end
    end
    JuMP.set_objective_sense(model, MOI.MAX_SENSE)
    JuMP.set_objective_function(model, objective)
    (; model, x, rowrefs, lower_parameters, upper_parameters, rhs_parameters,
        cost_parameters, lower, upper, rhs, costs,
        scratch_lower=similar(lower), scratch_upper=similar(upper),
        scratch_rhs=similar(rhs), scratch_costs=similar(costs),
        previous=Ref(point), data, l, parametric_bounds)
end

function poi_update!(state, point::Point)
    old = state.previous[]
    l, data = state.l, state.data
    if SV.stocks_changed(old, point)
        SV.column_bounds!(state.scratch_lower, state.scratch_upper, data, l, point)
        for j in 1:l.ncols
            if isfinite(state.lower[j]) &&
                    !isequal(state.lower[j], state.scratch_lower[j])
                parameter = state.lower_parameters[j]
                parameter === nothing ?
                    JuMP.set_lower_bound(state.x[j], state.scratch_lower[j]) :
                    JuMP.set_parameter_value(parameter, state.scratch_lower[j])
                state.lower[j] = state.scratch_lower[j]
            end
            if isfinite(state.upper[j]) &&
                    !isequal(state.upper[j], state.scratch_upper[j])
                parameter = state.upper_parameters[j]
                parameter === nothing ?
                    JuMP.set_upper_bound(state.x[j], state.scratch_upper[j]) :
                    JuMP.set_parameter_value(parameter, state.scratch_upper[j])
                state.upper[j] = state.scratch_upper[j]
            end
        end
    end
    if SV.rhs_changed(old, point)
        SV.row_rhs!(state.scratch_rhs, data, l, point)
        for j in 1:l.nrows
            if !isequal(state.rhs[j], state.scratch_rhs[j])
                parameter = state.rhs_parameters[j]
                parameter === nothing ?
                    JuMP.set_normalized_rhs(state.rowrefs[j], state.scratch_rhs[j]) :
                    JuMP.set_parameter_value(parameter, state.scratch_rhs[j])
                state.rhs[j] = state.scratch_rhs[j]
            end
        end
    end
    if SV.profits_changed(old, point)
        SV.objective_costs!(state.scratch_costs, data, l, point)
        for i in eachindex(state.costs)
            if !isequal(state.costs[i], state.scratch_costs[i])
                parameter = state.cost_parameters[i]
                parameter === nothing ?
                    JuMP.set_objective_coefficient(state.model,
                        state.x[l.objective_columns[i]], state.scratch_costs[i]) :
                    JuMP.set_parameter_value(parameter, state.scratch_costs[i])
                state.costs[i] = state.scratch_costs[i]
            end
        end
    end
    state.previous[] = point
    return nothing
end

function poi_solve!(state, point::Point)
    poi_update!(state, point)
    JuMP.optimize!(state.model)
    JuMP.termination_status(state.model) == MOI.OPTIMAL ||
        error("POI failed at $point: $(JuMP.termination_status(state.model))")
    JuMP.objective_value(state.model)
end

include(joinpath(@__DIR__, "jump_user_model.jl"))

# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

"""
`PATHS` is the default set. `MATRIX_PATHS` crosses the four modeling layers with
the two solvers, which is what separates the cost of a modeling layer from the
cost of a backend.
"""
const PATHS = benchmark_paths((:rm_highs, :rm_tiny, :rm_clp, :jump, :moi, :poi, :poi_hybrid))
const MATRIX_PATHS = benchmark_paths((:rm_highs, :rm_tiny, :jump, :jump_direct, :jump_tiny,
                      :jump_direct_tiny, :moi, :moi_tiny, :poi, :poi_tiny))

"""
    build_models(data, points)

One independent set of front-ends, each with its own solver instance and its own
warm basis. Building a set is expensive and is never inside a measured region.
"""
function build_models(data::ValleyData, points)
    l = SV.layout(data)
    point = first(points)
    moving = moving_entries(data, l, points)
    const_data = reactive_const_data(data)
    (; l, moving, data,
       rm_highs=RM.LinearProblem(const_data, point;
           backend=RM.HiGHSBackend(output_flag=false, threads=1, parallel="off",
               solver="simplex", simplex_max_concurrency=1)),
       rm_tiny=BENCH_TINY ? RM.LinearProblem(const_data, point; backend=RM.TinyHiGHSBackend()) : nothing,
       rm_clp=RM.LinearProblem(const_data, point; backend=ClpBackend()),
       jump=jump_reference(data, l, point),
       jump_user=jump_user_model(data, point),
       jump_user_pred=jump_user_model(data, point),
       jump_user_tiny=BENCH_TINY ? jump_user_model(data, point; optimizer=OPTIMIZERS.tiny) : nothing,
       jump_tiny=BENCH_TINY ? jump_reference(data, l, point; optimizer=OPTIMIZERS.tiny) : nothing,
       jump_direct=jump_reference(data, l, point; direct=true),
       jump_direct_tiny=BENCH_TINY ? jump_reference(data, l, point; optimizer=OPTIMIZERS.tiny,
           direct=true) : nothing,
       moi=moi_reference(data, l, point),
       moi_tiny=BENCH_TINY ? moi_reference(data, l, point; optimizer_type=OPTIMIZERS.tiny) : nothing,
       poi=poi_reference(data, l, moving, point; parametric_bounds=true),
       poi_tiny=BENCH_TINY ? poi_reference(data, l, moving, point; parametric_bounds=true,
           optimizer=OPTIMIZERS.tiny) : nothing,
       poi_hybrid=poi_reference(data, l, moving, point; parametric_bounds=false))
end

"""
    solve_point!(models, path, point; read_dual)

Solve one point through one front-end and return the objective, optionally
adding the dual of the first water-balance row so that the dual read is inside
the measured region. JuMP, MOI and POI report MOI's opposite dual sign, which is
why the sign differs between branches.
"""
function solve_point!(models, path::Symbol, point::Point; read_dual::Bool)
    if path === :rm_highs || path === :rm_tiny || path === :rm_clp
        lp = path === :rm_highs ? models.rm_highs :
             path === :rm_tiny ? models.rm_tiny : models.rm_clp
        read_dual || return RM.objective_value(lp, point)
        solution = RM.solve(lp, point)
        return RM.objective_value(solution) +
            RM.dual(solution, :water_balance, (1, 1))
    elseif path === :jump_user || path === :jump_user_pred ||
            path === :jump_user_tiny
        state = path === :jump_user_tiny ? models.jump_user_tiny :
                path === :jump_user_pred ? models.jump_user_pred : models.jump_user
        value = jump_user_solve!(state, point; predicated = path === :jump_user_pred)
        return read_dual ?
            value - JuMP.dual(state.water_balance[1, 1]) : value
    elseif path === :jump || path === :jump_tiny ||
            path === :jump_direct || path === :jump_direct_tiny
        state = path === :jump ? models.jump :
                path === :jump_tiny ? models.jump_tiny :
                path === :jump_direct ? models.jump_direct : models.jump_direct_tiny
        value = jump_solve!(state, point)
        return read_dual ?
            value - JuMP.dual(state.rowrefs[models.l.water_balance[1, 1]]) : value
    elseif path === :moi || path === :moi_tiny
        state = path === :moi ? models.moi : models.moi_tiny
        value = moi_solve!(state, point)
        return read_dual ? value - MOI.get(state.optimizer, MOI.ConstraintDual(),
            state.rowrefs[models.l.water_balance[1, 1]]) : value
    elseif path === :poi || path === :poi_tiny || path === :poi_hybrid
        state = path === :poi ? models.poi :
                path === :poi_tiny ? models.poi_tiny : models.poi_hybrid
        value = poi_solve!(state, point)
        return read_dual ?
            value - JuMP.dual(state.rowrefs[models.l.water_balance[1, 1]]) : value
    end
    error("unknown path: $path")
end

run_path!(models, path::Symbol, points; read_dual::Bool) =
    sum(solve_point!(models, path, point; read_dual) for point in points)

"""
    check_parity!(models, points; atol, tiny_atol)

Every front-end must return the same objective and the same first water-balance
dual at every point. TinyHiGHS is a different simplex and is held to a looser
tolerance; at a degenerate optimum it may also select a different valid dual, so
a dual mismatch there is reported rather than raised.
"""
function check_parity!(models, points; atol=1e-7, rtol=1e-9, tiny_atol=1e-5,
        clp_atol=1e-6)
    tiny_objective_error = 0.0
    tiny_dual_error = 0.0
    clp_objective_error = 0.0
    clp_dual_error = 0.0
    row = models.l.water_balance[1, 1]
    for point in points
        solution = RM.solve(models.rm_highs, point)
        reference = RM.objective_value(solution)
        reference_dual = RM.dual(solution, :water_balance, (1, 1))
        values = (jump=jump_solve!(models.jump, point),
                  moi=moi_solve!(models.moi, point),
                  poi=poi_solve!(models.poi, point),
                  poi_hybrid=poi_solve!(models.poi_hybrid, point))
        all(isapprox(reference, value; atol, rtol) for value in values) ||
            error("objective mismatch at $point: $reference versus $values")
        duals = (jump=JuMP.dual(models.jump.rowrefs[row]),
                 moi=MOI.get(models.moi.optimizer, MOI.ConstraintDual(),
                     models.moi.rowrefs[row]),
                 poi=JuMP.dual(models.poi.rowrefs[row]),
                 poi_hybrid=JuMP.dual(models.poi_hybrid.rowrefs[row]))
        all(isapprox(-reference_dual, value; atol=1e-6, rtol=1e-6) for value in duals) ||
            error("dual mismatch at $point: $reference_dual versus $duals")
        if BENCH_TINY
            tiny_solution = RM.solve(models.rm_tiny, point)
            tiny_objective_error = max(tiny_objective_error,
                abs(reference - RM.objective_value(tiny_solution)))
            tiny_dual_error = max(tiny_dual_error,
                abs(reference_dual - RM.dual(tiny_solution, :water_balance, (1, 1))))
        end
        clp_solution = RM.solve(models.rm_clp, point)
        clp_objective_error = max(clp_objective_error,
            abs(reference - RM.objective_value(clp_solution)))
        clp_dual_error = max(clp_dual_error,
            abs(reference_dual - RM.dual(clp_solution, :water_balance, (1, 1))))
    end
    tiny_objective_error <= tiny_atol * max(1.0, abs(tiny_objective_error)) ||
        error("TinyHiGHS objective error $tiny_objective_error exceeds $tiny_atol")
    clp_objective_error <= clp_atol * max(1.0, abs(clp_objective_error)) ||
        error("Clp objective error $clp_objective_error exceeds $clp_atol")
    (; points=length(points), tiny_objective_error, tiny_dual_error,
       clp_objective_error, clp_dual_error)
end

"""
    timed_path!(models, path, points; read_dual, passes)

Median of `passes` warmed passes over the whole point sequence. Model
construction, compilation and the first pass are outside the measured region.
"""
function timed_path!(models, path::Symbol, points; read_dual::Bool, passes::Int=5)
    checksum = run_path!(models, path, points; read_dual)
    measurements = map(1:passes) do _
        GC.gc()
        @timed run_path!(models, path, points; read_dual)
    end
    all(m -> isfinite(m.value) && isapprox(checksum, m.value; rtol=1e-6),
        measurements) || error("unstable checksum on path $path")
    seconds = median(m.time for m in measurements)
    bytes = median(m.bytes for m in measurements)
    gc = median(m.gctime for m in measurements)
    (; path, read_dual, per_solve_us=1e6 * seconds / length(points),
       per_solve_bytes=bytes / length(points),
       gc_us_per_solve=1e6 * gc / length(points), gc_share=gc / seconds, checksum)
end

"""
    timed_parallel_path!(worker_models, path, points; read_dual, passes)

Wall time of the whole sequence split across `length(worker_models)` threads,
one persistent model set and one solver instance per worker. The reported
microseconds are wall time divided by the number of points, not single-LP
latency.
"""
function timed_parallel_path!(worker_models, path::Symbol, points;
        read_dual::Bool, passes::Int=5)
    nworkers = length(worker_models)
    nworkers > 1 || error("start Julia with at least two threads")
    nworkers <= Threads.nthreads() ||
        error("$nworkers workers exceed Threads.nthreads()=$(Threads.nthreads())")
    chunks = [points[worker:nworkers:end] for worker in 1:nworkers]
    partial = zeros(nworkers)
    function pass!()
        Threads.@threads :static for worker in 1:nworkers
            partial[worker] = run_path!(worker_models[worker], path, chunks[worker];
                read_dual)
        end
        sum(partial)
    end
    checksum = pass!()
    measurements = map(1:passes) do _
        GC.gc()
        @timed pass!()
    end
    all(m -> isfinite(m.value) && isapprox(checksum, m.value; rtol=1e-6),
        measurements) || error("unstable parallel checksum on path $path")
    seconds = median(m.time for m in measurements)
    bytes = median(m.bytes for m in measurements)
    gc = median(m.gctime for m in measurements)
    (; path, read_dual, nworkers, per_solve_us=1e6 * seconds / length(points),
       per_solve_bytes=bytes / length(points),
       gc_us_per_solve=1e6 * gc / length(points), gc_share=gc / seconds, checksum)
end

"""
    highs_iterations(backend)

Simplex iterations of the **last** solve on this HiGHS instance. With a warm
basis this is often zero or one; from a cold start it is the full re-solve.
"""
function highs_iterations(backend)
    value = Ref{Cint}(0)
    HiGHS.Highs_getIntInfoValue(backend.inner, "simplex_iteration_count", value)
    Int(value[])
end

"""
    warm_start_report(models, points; passes)

What the warm basis is actually worth on this point sequence. The same
persistent model solves the same points twice: once normally, and once with
`Highs_clearSolver` called before each solve, which throws the basis away while
leaving the model untouched. The gap between the two is the warm-start gain,
and the iteration counts say how close consecutive LPs really are.

This measures the *solver*, so it is reported for `rm_highs` only; every other
HiGHS front-end sees the same LPs and the same basis.
"""
function warm_start_report(models, points; passes::Int=3)
    lp = models.rm_highs
    function sweep(cold::Bool)
        iterations = 0
        checksum = 0.0
        for point in points
            cold && HiGHS.Highs_clearSolver(lp.backend.inner)
            checksum += RM.objective_value(lp, point)
            iterations += highs_iterations(lp.backend)
        end
        (iterations, checksum)
    end
    function timed(cold::Bool)
        sweep(cold)
        iterations, checksum = sweep(cold)
        seconds = median(map(1:passes) do _
            GC.gc()
            start = time_ns()
            sweep(cold)
            (time_ns() - start) / 1e9
        end)
        (; cold, iterations_per_solve=iterations / length(points),
           us_per_solve=1e6 * seconds / length(points), checksum)
    end
    warm, cold = timed(false), timed(true)
    isapprox(warm.checksum, cold.checksum; rtol=1e-9) ||
        error("clearing the basis changed the objective: $(warm.checksum) vs $(cold.checksum)")
    (; warm, cold, speedup=cold.us_per_solve / warm.us_per_solve)
end

# ---------------------------------------------------------------------------
# Production traversal
#
# Optional traversal with multiple concave regions:
#
#   scenarios are split into contiguous chunks, one chunk per thread
#     each chunk owns one TransitionProblemSet, which holds
#       one persistent LP and one warm basis per concave region
#     for each scenario of the chunk:
#       for each initial stock, visited successively by the grid refiner:
#         scan the reachable concave regions, each on its own model,
#         and keep the best objective
#
# Two consequences a single-model benchmark cannot show. First, a region's basis
# is never perturbed by another region: between two solves on the same model,
# only the initial stock, and occasionally the scenario, has moved. Second, one
# stock point costs as many LPs as there are reachable regions.
# ---------------------------------------------------------------------------

"""
    ProductionSchedule

`points[i]` is a (scenario, stock) pair — its `region` field is unused — and
`regions[i]` lists the concave regions to solve it on, in scan order.
"""
struct ProductionSchedule{NR}
    points::Vector{Point{NR}}
    regions::Vector{Vector{Int}}
    nregions::Int
    skipped::Int
end

Base.length(schedule::ProductionSchedule) = length(schedule.points)
solves(schedule::ProductionSchedule) = sum(length, schedule.regions)

"""
    production_points(data; scenarios, stocks_per_scenario)

The (scenario, stock) sequence in production order: scenario outermost, stock
levels visited in bisection order inside it.
"""
function production_points(data::ValleyData, scenarios;
        stocks_per_scenario::Int=9)
    NC = data.spec.ncut_reservoirs
    top = ntuple(r -> data.volume_max[r, data.spec.npostes], NC)
    fractions = SV.bisection_fractions(stocks_per_scenario)
    [Point{NC}(scenario, ntuple(r -> fraction * top[r], NC), 1)
     for scenario in scenarios for fraction in fractions]
end

"""
    production_schedule(models, data, points; nregions)

Attach to each point the regions the production solver would scan: the ones the
reachable-stock box selects, minus any that turn out to be infeasible from that
stock. The infeasible ones are removed **here**, outside every measured region,
by solving them once.

Production does pay for its infeasible probes; this harness does not, because
an infeasible solve raises through the front-ends and would force a model
rebuild in the middle of a timed loop. So the measured cost is the cost of the
feasible part of the scan, identically for all front-ends.
"""
function production_schedule(region_models, data::ValleyData, points;
        nregions::Int)
    kept = eltype(points)[]
    scanned = Vector{Int}[]
    skipped = 0
    for point in points
        candidates = SV.reachable_regions(data, point)
        isempty(candidates) && (candidates = 1:nregions)
        feasible = Int[]
        for region in candidates
            probe = typeof(point)(point.scenario, point.stocks, region)
            try
                RM.objective_value(region_models[region].rm_highs, probe)
                push!(feasible, region)
            catch
                # Infeasible from this stock: production skips it too, after
                # paying for the LP. Rebuild the model, which the failed solve
                # invalidated, before it is used again.
                region_models[region] = build_models(data, [probe])
            end
        end
        if isempty(feasible)
            # No reachable region is feasible from this stock. The traversal detects
            # this with `compute_feasible_volume_domain` and falls back on a
            # penalized problem; the harness simply drops the point.
            skipped += 1
            continue
        end
        push!(kept, point)
        push!(scanned, feasible)
    end
    isempty(kept) && error("every stock point was infeasible")
    ProductionSchedule(kept, scanned, nregions, skipped)
end

"""
    build_region_models(data, schedule_points, nregions)

One complete set of front-ends **per concave region**, each warmed only on its
own region, exactly as the production solver keeps one transition problem per
region. Building these is expensive and is never inside a measured region.
"""
function build_region_models(data::ValleyData, points, nregions::Int)
    map(1:nregions) do region
        pinned = [typeof(point)(point.scenario, point.stocks, region)
                  for point in points]
        build_models(data, pinned)
    end
end

"""
    production_run!(region_models, path, schedule)

One full pass of the production traversal through one front-end.
"""
function production_run!(region_models, path::Symbol, schedule::ProductionSchedule)
    total = 0.0
    for index in eachindex(schedule.points)
        point = schedule.points[index]
        best = -Inf
        for region in schedule.regions[index]
            value = solve_point!(region_models[region], path,
                typeof(point)(point.scenario, point.stocks, region);
                read_dual=false)
            value > best && (best = value)
        end
        total += best
    end
    total
end

"""
    timed_production!(region_models, path, schedule; passes)

Per-LP figures, so the cost of scanning several regions per stock point does not
hide inside a per-point average: the divisor is the number of solves, not the
number of stock points.
"""
function timed_production!(region_models, path::Symbol, schedule::ProductionSchedule;
        passes::Int=5)
    checksum = production_run!(region_models, path, schedule)
    measurements = map(1:passes) do _
        GC.gc()
        @timed production_run!(region_models, path, schedule)
    end
    all(m -> isfinite(m.value) && isapprox(checksum, m.value; rtol=1e-6),
        measurements) || error("unstable production checksum on path $path")
    seconds = median(m.time for m in measurements)
    bytes = median(m.bytes for m in measurements)
    nsolves = solves(schedule)
    gc = median(m.gctime for m in measurements)
    (; path, nsolves, per_solve_us=1e6 * seconds / nsolves,
       per_solve_bytes=bytes / nsolves,
       gc_us_per_solve=1e6 * gc / nsolves, gc_share=gc / seconds, checksum)
end

"""
    timed_production_parallel!(worker_models, worker_schedules, path; passes)

Scenarios are split into contiguous chunks, one chunk per worker, and each
worker owns its own per-region model sets. The parallel traversal uses a
worker never sees another worker's scenarios, so its warm bases follow one
continuous trajectory.
"""
function timed_production_parallel!(worker_models, worker_schedules, path::Symbol;
        passes::Int=5)
    nworkers = length(worker_models)
    nworkers > 1 || error("start Julia with at least two threads")
    nworkers <= Threads.nthreads() ||
        error("$nworkers workers exceed Threads.nthreads()=$(Threads.nthreads())")
    partial = zeros(nworkers)
    function pass!()
        Threads.@threads :static for worker in 1:nworkers
            partial[worker] = production_run!(worker_models[worker], path,
                worker_schedules[worker])
        end
        sum(partial)
    end
    checksum = pass!()
    measurements = map(1:passes) do _
        GC.gc()
        @timed pass!()
    end
    all(m -> isfinite(m.value) && isapprox(checksum, m.value; rtol=1e-6),
        measurements) || error("unstable parallel production checksum on $path")
    seconds = median(m.time for m in measurements)
    bytes = median(m.bytes for m in measurements)
    nsolves = sum(solves, worker_schedules)
    gc = median(m.gctime for m in measurements)
    (; path, nworkers, nsolves, per_solve_us=1e6 * seconds / nsolves,
       per_solve_bytes=bytes / nsolves,
       gc_us_per_solve=1e6 * gc / nsolves, gc_share=gc / seconds, checksum)
end

"""
    scenario_chunks(nscenarios, nworkers)

Contiguous scenario chunks, one persistent workspace per chunk.
"""
function scenario_chunks(nscenarios::Int, nworkers::Int)
    1 <= nworkers <= nscenarios || throw(ArgumentError("workers must be between 1 and the scenario count"))
    edges = round.(Int, range(0, nscenarios; length=nworkers + 1))
    [(edges[w] + 1):edges[w + 1] for w in 1:nworkers]
end

"""
    run_production(spec; nscenarios, nregions, stocks_per_scenario, nworkers)

The production traversal, serial and in parallel, for every front-end.
"""
function run_production(spec::SV.ValleySpec; nscenarios::Int=spec.nscenarios,
        nregions::Int=3, stocks_per_scenario::Int=9,
        nworkers::Int=min(Threads.nthreads(), 8))
    data = SV.valley_data(spec; nregions)
    points = production_points(data, 1:nscenarios; stocks_per_scenario)
    models = build_region_models(data, points, nregions)
    schedule = production_schedule(models, data, points; nregions)
    serial = [timed_production!(models, path, schedule) for path in PATHS]
    parallel = if nworkers > 1
        chunks = scenario_chunks(nscenarios, nworkers)
        worker_points = [production_points(data, chunk; stocks_per_scenario)
                         for chunk in chunks]
        worker_models = [build_region_models(data, worker_points[w], nregions)
                         for w in 1:nworkers]
        worker_schedules = [production_schedule(worker_models[w], data,
            worker_points[w]; nregions) for w in 1:nworkers]
        [timed_production_parallel!(worker_models, worker_schedules, path)
         for path in PATHS]
    else
        nothing
    end
    (; shape=SV.describe(data), nregions, npoints=length(schedule),
       skipped=schedule.skipped, nsolves=solves(schedule),
       regions_per_point=solves(schedule) / length(schedule),
       serial, parallel)
end

const PRESETS = Dict("small" => SV.VALLEY_SMALL,
                     "medium" => SV.VALLEY_MEDIUM,
                     "large" => SV.VALLEY_LARGE)

"""
    run_preset(spec; nscenarios, nworkers, parity_points)

Parity check, serial timings and, when more than one thread is available,
parallel timings for one preset.
"""
function run_preset(spec::SV.ValleySpec; nscenarios::Int=spec.nscenarios,
        nworkers::Int=min(Threads.nthreads(), 8), parity_points::Int=24,
        order::Symbol=:scenario_fast, stocks_per_scenario::Int=3)
    data = SV.valley_data(spec)
    points = SV.valley_points(data; nscenarios, order, stocks_per_scenario)
    models = build_models(data, points)
    parity = check_parity!(models, points[1:min(parity_points, length(points))])
    moving = models.moving
    warm_start = warm_start_report(models, points)
    serial = [timed_path!(models, path, points; read_dual)
              for read_dual in (false, true) for path in PATHS]
    parallel = if nworkers > 1
        workers = [build_models(data, points) for _ in 1:nworkers]
        [timed_parallel_path!(workers, path, points; read_dual)
         for read_dual in (false, true) for path in PATHS]
    else
        nothing
    end
    (; shape=SV.describe(data), order, npoints=length(points), parity,
       moving=(; lower=count(moving.lower_moves), upper=count(moving.upper_moves),
               rhs=count(moving.rhs_moves), costs=count(moving.cost_moves)),
       warm_start, serial, parallel)
end

"""
    main(args)

`valley_comparison.jl [preset] [nscenarios] [order] [nregions]`, where `preset`
is `small`, `medium`, `large` or `all`, and `order` is `scenario_fast`,
`stock_fast`, `refinement`, `all`, or `production`. The `production` order runs
the multi-region traversal — contiguous scenario chunks, one persistent LP and
one warm basis per concave region, a linear scan of the reachable regions at
every stock point — and then `nregions` (default 3) sets how many concave
regions the future Bellman function is split into.
"""
function main(args=ARGS)
    preset = isempty(args) ? "all" : lowercase(args[1])
    names = preset == "all" ? ("small", "medium", "large") : (preset,)
    requested = length(args) >= 3 ? lowercase(args[3]) : "scenario_fast"
    orders = requested == "all" ?
        (:scenario_fast, :stock_fast, :refinement) : (Symbol(requested),)
    nregions = length(args) >= 4 ? parse(Int, args[4]) : 3
    for name in names
        haskey(PRESETS, name) || error("unknown preset $name")
        spec = PRESETS[name]
        nscenarios = length(args) >= 2 ? parse(Int, args[2]) : spec.nscenarios
        if requested == "production"
            result = run_production(spec; nscenarios, nregions)
            println("### ", spec.name, "  production traversal, ",
                result.nregions, " concave regions")
            println("  shape      = ", result.shape)
            println("  stock points = ", result.npoints, "  LP solves = ",
                result.nsolves, "  regions/point = ",
                round(result.regions_per_point; digits=2))
            for row in result.serial
                println("  serial   ", row)
            end
            if result.parallel !== nothing
                for row in result.parallel
                    println("  parallel ", row)
                end
            end
            continue
        end
        for order in orders
            result = run_preset(spec; nscenarios, order)
            println("### ", spec.name, "  order=", order)
            println("  shape      = ", result.shape)
            println("  points     = ", result.npoints, "  parity = ", result.parity)
            println("  moving     = ", result.moving)
            println("  warm start = ", result.warm_start)
            for row in result.serial
                println("  serial   ", row)
            end
            result.parallel === nothing && continue
            for row in result.parallel
                println("  parallel ", row)
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
