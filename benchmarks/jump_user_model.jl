"""
The synthetic transition LP expressed with named, indexed JuMP containers.

`valley_comparison.jl` also provides a flat-array JuMP formulation with cached
numerical buffers. Both retain the solver model between points. This file keeps
reservoir/unit/group keys in the formulation and updates coefficients through
the corresponding JuMP references.

Two update policies are measured separately:

* `jump_user_update!` reapplies the parameter-dependent quantities at each point;
* `jump_user_update_predicated!` skips blocks according to the same parameter
  dependency predicates used by ReactiveModel.

The difference measures these implementations and update policies, not an
intrinsic property of all JuMP models. Group flows use dependent index sets,
which produce a SparseAxisArray. Explicit affine rows and set-form constraints
preserve the synthetic layout, including rows with an infinite upper bound.
"""

using JuMP
import MathOptInterface as MOI

"""
    as_row(v)

Return an affine expression containing the variable with coefficient one.
The benchmark uses explicit affine rows for final-stock and inactive constraints
to match its numerical layout and row-result queries.
"""
as_row(v::VariableRef) = (e = AffExpr(0.0); add_to_expression!(e, 1.0, v); e)

"""
    jump_user_model(data, point; optimizer)

Build the transition LP with named, indexed JuMP containers.

    V0[r]                                           initial volume, pinned to the stock
    V[r, t]  in [Vmin(r,t), Vmax(r,t)]              volume at the end of each poste
    D[u, t]  in [Qmin(u,t), Qmax(u,t)]              unit flow
    q[u, g, t] in [qmin, qmax]                      group flow, energetic units only
    vB                                              free, the Bellman value

    water_balance[r, t] :  V[r,t] - V[r,t-1] + dt*(outflows) - dt*(inflows) == apports(r,t)
    final_volume[r]     :  V[r, T] == target(r)                  reservoirs outside the cuts
    group_sum[u, t]     :  sum(q[u,:,t]) - D[u,t] == 0           energetic units
    bellman_cut[s]      :  vB - sum(slope(s,r) * V[r,T]) <= intercept(s)
    min_final[r]        :  V[r, T] >= region lower bound
    max_final[r]        :  V[r, T] <= region upper bound
    inactive_cut, inactive_triangle[k] : kept in the matrix, neutralized by an
                                         infinite right-hand side

    max  vB + sum(profit(u,g,t) * q[u,g,t])
"""

function jump_user_model(data::ValleyData, point::Point; optimizer=HiGHS.Optimizer)
    spec = data.spec
    R, T = 1:spec.nreservoirs, 1:spec.npostes
    U, NC = 1:length(spec.unit_groups), spec.ncut_reservoirs
    E = data.energetic
    last_T = spec.npostes

    model = JuMP.Model(optimizer)
    JuMP.set_silent(model)
    apply_solver_options!(model)

    @variable(model, V0[r in R])
    for r in R
        level = SV.initial_volume_bound(data, point, r)
        set_lower_bound(V0[r], level)
        set_upper_bound(V0[r], level)
    end

    @variable(model, V[r in R, t in T])
    for r in R, t in T
        lower, upper = SV.volume_bounds(data, point, r, t)
        set_lower_bound(V[r, t], lower)
        set_upper_bound(V[r, t], upper)
    end

    @variable(model, D[u in U, t in T])
    for u in U, t in T
        lower, upper = SV.flow_bounds(data, point, u, t)
        set_lower_bound(D[u, t], lower)
        set_upper_bound(D[u, t], upper)
    end

    # Dependent index sets: the group count depends on the unit, and only
    # energetic units have groups at all.
    @variable(model, q[u in E, g in 1:data.ngroups[u], t in T])
    for u in E, g in 1:data.ngroups[u], t in T
        lower, upper = SV.group_flow_bounds(data, point, u, g, t)
        set_lower_bound(q[u, g, t], lower)
        set_upper_bound(q[u, g, t], upper)
    end

    @variable(model, vB)

    upstream(r) = [u for u in U if spec.unit_upstream[u] == r]
    downstream(r) = [u for u in U if spec.unit_downstream[u] == r]

    @constraint(model, water_balance[r in R, t in T],
        V[r, t] - (t == 1 ? V0[r] : V[r, t - 1])
        + sum(data.duration[t] * D[u, t] for u in upstream(r); init = AffExpr(0.0))
        - sum(data.duration[t] * D[u, t] for u in downstream(r); init = AffExpr(0.0))
        == SV.inflow_value(data, point, r, t))

    outside_cuts = (NC + 1):spec.nreservoirs
    @constraint(model, final_volume[r in outside_cuts],
        V[r, last_T] == SV.final_volume_target(data, point, r))

    @constraint(model, group_sum[u in E, t in T],
        sum(q[u, g, t] for g in 1:data.ngroups[u]) - D[u, t] == 0.0)

    # Same trap on the cuts: a cut outside the selected concave region has an
    # infinite right-hand side, so the set form is needed there too.
    @constraint(model, bellman_cut[s in 1:spec.ncuts],
        vB - sum(SV.cut_slope(data, s, r) * V[r, last_T] for r in 1:NC)
        in MOI.LessThan(SV.cut_rhs(data, point, s)))

    # `@constraint(model, V[1, last_T] <= Inf)` is rejected: JuMP moves the bound
    # into the affine expression's constant and refuses a non-finite one. A row
    # retained in the synthetic matrix and neutralized with an infinite
    # right-hand side therefore has no natural spelling in the macro syntax, and
    # a user has to reach for the set form.
    @constraint(model, inactive_cut, as_row(V[1, last_T]) in MOI.LessThan(Inf))

    @constraint(model, min_final[r in 1:NC],
        as_row(V[r, last_T]) >= region_scale(data, r) * SV.min_final_stock(data, point))
    @constraint(model, max_final[r in 1:NC],
        as_row(V[r, last_T]) <= region_scale(data, r) * SV.max_final_stock(data, point))

    @constraint(model, inactive_triangle[k in 1:(NC * NC)],
        as_row(V[mod1(k, NC), last_T]) in MOI.LessThan(Inf))

    @objective(model, Max,
        vB + sum(SV.profit(data, point, u, g, t) * q[u, g, t]
                 for u in E for g in 1:data.ngroups[u] for t in T))

    (; model, V0, V, D, q, vB, water_balance, final_volume, group_sum,
       bellman_cut, min_final, max_final, data, R, T, U, E, NC, last_T,
       previous=Ref(point))
end

"""
    jump_user_update!(state, point)

Reapply the parameter-dependent bounds, RHS values and objective coefficients
through named JuMP references. Infinite RHS cuts stay inactive.
"""
function jump_user_update!(state, point::Point)
    data, R, T, U, E = state.data, state.R, state.T, state.U, state.E

    for r in R
        level = SV.initial_volume_bound(data, point, r)
        set_lower_bound(state.V0[r], level)
        set_upper_bound(state.V0[r], level)
    end
    for u in U, t in T
        lower, upper = SV.flow_bounds(data, point, u, t)
        set_lower_bound(state.D[u, t], lower)
        set_upper_bound(state.D[u, t], upper)
    end
    for u in E, g in 1:data.ngroups[u], t in T
        lower, upper = SV.group_flow_bounds(data, point, u, g, t)
        set_lower_bound(state.q[u, g, t], lower)
        set_upper_bound(state.q[u, g, t], upper)
    end

    for r in R, t in T
        set_normalized_rhs(state.water_balance[r, t],
            SV.inflow_value(data, point, r, t))
    end
    for r in (state.NC + 1):length(R)
        set_normalized_rhs(state.final_volume[r],
            SV.final_volume_target(data, point, r))
    end
    for s in 1:data.spec.ncuts
        value = SV.cut_rhs(data, point, s)
        isfinite(value) && set_normalized_rhs(state.bellman_cut[s], value)
    end
    for r in 1:state.NC
        scale = region_scale(data, r)
        set_normalized_rhs(state.min_final[r], scale * SV.min_final_stock(data, point))
        set_normalized_rhs(state.max_final[r], scale * SV.max_final_stock(data, point))
    end

    for u in E, g in 1:data.ngroups[u], t in T
        set_objective_coefficient(state.model, state.q[u, g, t],
            SV.profit(data, point, u, g, t))
    end
    state.previous[] = point
    return nothing
end

"""
    jump_user_update_predicated!(state, point)

Apply named updates only for blocks selected by the parameter dependency
predicates. These are block-level guards; the flat-array implementation also
compares individual numerical values before setter calls.
"""
function jump_user_update_predicated!(state, point::Point)
    data, R, T, U, E = state.data, state.R, state.T, state.U, state.E
    old = state.previous[]

    if SV.stocks_changed(old, point)
        for r in R
            level = SV.initial_volume_bound(data, point, r)
            set_lower_bound(state.V0[r], level)
            set_upper_bound(state.V0[r], level)
        end
        if data.spec.head_dependent_bounds
            for u in U, t in T
                lower, upper = SV.flow_bounds(data, point, u, t)
                set_lower_bound(state.D[u, t], lower)
                set_upper_bound(state.D[u, t], upper)
            end
            for u in E, g in 1:data.ngroups[u], t in T
                lower, upper = SV.group_flow_bounds(data, point, u, g, t)
                set_lower_bound(state.q[u, g, t], lower)
                set_upper_bound(state.q[u, g, t], upper)
            end
        end
    end

    if SV.scenario_changed(old, point)
        for r in R, t in T
            set_normalized_rhs(state.water_balance[r, t],
                SV.inflow_value(data, point, r, t))
        end
    end
    if SV.region_changed(old, point)
        for s in 1:data.spec.ncuts
            value = SV.cut_rhs(data, point, s)
            isfinite(value) && set_normalized_rhs(state.bellman_cut[s], value)
        end
        for r in 1:state.NC
            scale = region_scale(data, r)
            set_normalized_rhs(state.min_final[r], scale * SV.min_final_stock(data, point))
            set_normalized_rhs(state.max_final[r], scale * SV.max_final_stock(data, point))
        end
    end

    if SV.profits_changed(old, point)
        for u in E, g in 1:data.ngroups[u], t in T
            set_objective_coefficient(state.model, state.q[u, g, t],
                SV.profit(data, point, u, g, t))
        end
    end
    state.previous[] = point
    return nothing
end

function jump_user_solve!(state, point::Point; predicated::Bool=false)
    predicated ? jump_user_update_predicated!(state, point) :
                 jump_user_update!(state, point)
    optimize!(state.model)
    termination_status(state.model) == MOI.OPTIMAL ||
        error("JuMP (user model) failed at $point: $(termination_status(state.model))")
    objective_value(state.model)
end
