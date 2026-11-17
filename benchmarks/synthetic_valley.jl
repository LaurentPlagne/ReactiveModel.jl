"""
    SyntheticValley

A deterministic hydro-cascade transition LP for numerical experiments.
The public presets use explicitly chosen round dimensions, simple chain/bypass
graphs and seeds 11, 22 and 33. They are not calibrated to a particular site.
See `PROVENANCE.md` for the choices and scope of the experiment.

Everything is generated from an integer seed and a `ValleySpec`, allowing
benchmarks to be published, shared, and reproduced by any user.

# What is modeled

A multi-period transition problem with storage, discharge, and future value
cuts (generalizing the single-reservoir model in Tutorial 3 to multi-reservoir cascades).
The model is solved repeatedly across parameter points that vary:
* the scenario (inflows, electricity market prices),
* the initial stock levels of each reservoir,
* the concave piecewise-linear region of the future value function.

Between points, the mathematical structure of the LP never changes: the same
columns, rows, and non-zero sparsity pattern. Only variable bounds, right-hand
sides, and objective coefficients move.

# Column layout

| block | keys | count |
|---|---|---|
| `initial_volumes` | `r` | `nreservoirs` |
| `volumes` | `(r, t)` | `nreservoirs * npostes` |
| `debits` | `(u, t)` | `nunits * npostes` |
| `group_debits` | `(u, g, t)`, energetic `u` only | `sum(ngroups) * npostes` |
| `vB` | - | 1 |

# Row layout

| block | keys | count | op |
|---|---|---|---|
| `water_balance` | `(r, t)` | `nreservoirs * npostes` | `=` |
| `final_volume` | `r`, reservoirs outside the cuts | `nreservoirs - ncut_reservoirs` | `=` |
| `group_sum` | `(u, t)`, energetic `u` | `nenergetic * npostes` | `=` |
| `bellman_cut` | `s` | `ncuts` | `<=` |
| `inactive_cut` | - | 1 | `<=`, rhs `Inf` |
| `min_final_stock` | cut reservoir | `ncut_reservoirs` | `>=` |
| `max_final_stock` | cut reservoir | `ncut_reservoirs` | `<=` |
| `inactive_triangle` | cut reservoir | `ncut_reservoirs^2` | `<=`, rhs `Inf` |
"""
module SyntheticValley

export ValleySpec, ValleyData, Point, valley_data, valley_points, inflow_band,
    bisection_fractions, concave_regions, reachable_regions,
    VALLEY_SMALL, VALLEY_MEDIUM, VALLEY_LARGE, describe

# ---------------------------------------------------------------------------
# Specification and presets
# ---------------------------------------------------------------------------

"""
    ValleySpec

Shape of a synthetic valley. `unit_groups[u] == 0` marks a non-energetic unit:
it has a flow column but no group columns and no group-sum row.

`unit_upstream[u]` is the reservoir the unit draws from, `unit_downstream[u]`
the reservoir it discharges into, or `0` to leave the valley.
"""
Base.@kwdef struct ValleySpec
    name::String
    npostes::Int
    nscenarios::Int
    ncuts::Int
    unit_upstream::Vector{Int}
    unit_downstream::Vector{Int}
    unit_groups::Vector{Int}
    nreservoirs::Int = maximum(vcat(unit_upstream, unit_downstream))
    ncut_reservoirs::Int = 1
    bellman_slopes::Int = 0
    release_fraction::Float64 = 0.25
    inflow_fraction::Tuple{Float64,Float64} = (0.01, 0.15)
    head_dependent_bounds::Bool = false
    seed::Int = 11
end

"""
    inflow_band(spec)

The inflow band actually used: the upper end is clamped just below
`release_fraction` so that incoming flow can be discharged downstream.
"""
function inflow_band(spec::ValleySpec)
    lower, upper = spec.inflow_fraction
    ceiling = 0.95 * spec.release_fraction
    (min(lower, ceiling), min(upper, ceiling))
end

nreservoirs(spec::ValleySpec) = spec.nreservoirs
nunits(spec::ValleySpec) = length(spec.unit_groups)
energetic_units(spec::ValleySpec) = findall(>(0), spec.unit_groups)

"""
`VALLEY_SMALL`: one reservoir, two units, eight time steps and sixteen cut rows.
"""
const VALLEY_SMALL = ValleySpec(
    name="ValleySmall",
    npostes=8, nscenarios=64, ncuts=16,
    unit_upstream=[1, 1], unit_downstream=[0, 0], unit_groups=[1, 2],
    nreservoirs=1,
    release_fraction=0.5, inflow_fraction=(0.02, 0.1), seed=11,
)

"""
`VALLEY_MEDIUM`: a three-reservoir chain with two outlet bypasses,
twelve time steps and twenty-four cut rows.
"""
const VALLEY_MEDIUM = ValleySpec(
    name="ValleyMedium",
    npostes=12, nscenarios=64, ncuts=24,
    unit_upstream=[1, 1, 2, 2, 3], unit_downstream=[2, 0, 3, 0, 0],
    unit_groups=[2, 0, 2, 0, 2],
    nreservoirs=3, bellman_slopes=6,
    release_fraction=0.5, inflow_fraction=(0.02, 0.1), seed=22,
)

"""
`VALLEY_LARGE`: a five-reservoir chain with four parallel bypasses,
sixteen time steps and sixty-four cut rows on two cut reservoirs.
"""
const VALLEY_LARGE = ValleySpec(
    name="ValleyLarge",
    npostes=16, nscenarios=64, ncuts=64,
    unit_upstream=[1, 2, 3, 4, 5, 1, 2, 3, 4],
    unit_downstream=[2, 3, 4, 5, 0, 2, 3, 4, 5],
    unit_groups=[2, 2, 2, 2, 2, 0, 0, 0, 0],
    nreservoirs=5,
    ncut_reservoirs=2,
    release_fraction=0.5, inflow_fraction=(0.02, 0.1), seed=33,
    bellman_slopes=8,
)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

"""
    ValleyData

The constant data of one study: every array here is read-only for the whole life
of a model. This is the payload of a ReactiveModel `ConstData`, and nothing in it
depends on the point being solved.
"""
struct ValleyData
    spec::ValleySpec
    energetic::Vector{Int}
    ngroups::Vector{Int}
    duration::Vector{Float64}                 # [poste], hours
    inflow::Array{Float64,3}                  # [reservoir, poste, scenario]
    volume_min::Matrix{Float64}               # [reservoir, poste]
    volume_max::Matrix{Float64}
    flow_min::Matrix{Float64}                 # [unit, poste]
    flow_max::Matrix{Float64}
    group_min::Array{Float64,3}               # [unit, group, poste]
    group_max::Array{Float64,3}
    price::Matrix{Float64}                    # [poste, scenario]
    efficiency::Array{Float64,3}              # [unit, group, poste], MWh per m3/s-hour
    reference_head::Vector{Float64}           # [reservoir]
    cut_stock::Vector{Float64}                # [ncuts + 1], increasing
    cut_value::Vector{Float64}                # [ncuts + 1], concave
    regions::Vector{Tuple{Int,Int}}           # concave regions, cut-point indices
    initial_volume::Vector{Float64}           # [reservoir], the study's input level
    # Precomputed for the reachable-region filter, so that the filter itself
    # costs O(number of regions) per stock point and allocates nothing.
    max_release::Vector{Float64}              # [reservoir], over the whole period
    total_inflow::Matrix{Float64}             # [reservoir, scenario]
end

"""
    Point{NC}

The light, `isbits` description of one resolution: which scenario, the initial
stock of each **piloted lake**, which concave region of the Bellman function.
Nothing else may live here; the study arrays belong to `ValleyData`.

`NC` is `ncut_reservoirs`, not the number of reservoirs. A reservoir without a
Bellman cut keeps a fixed initial volume; only the cut reservoirs are varied
across stock points. This is an experimental choice about update locality.
"""
struct Point{NC}
    scenario::Int
    stocks::NTuple{NC,Float64}
    region::Int
end

"""
The marginal value of water at an empty reservoir, as a fraction of the best
energy value in the study, and the ratio between that top slope and the slope at
a full reservoir. The public experiment chooses a factor of four and a top
value of 0.8 to keep both energy costs and terminal storage relevant.
"""
const TOP_WATER_VALUE = 0.8
const WATER_VALUE_SPREAD = 4.0

"""
Natural inflows never exceed this fraction of the release capacity below the
reservoir receiving them, and a unit is never less available than
`MINIMUM_AVAILABILITY` of its rated flow. Together these keep every generated
point feasible, including a reservoir that starts full.
"""
const DRAIN_MARGIN = 0.6
const MINIMUM_AVAILABILITY = 0.7

# A small deterministic generator, so that the data depends only on the seed and
# never on the Julia version's global RNG stream.
@inline _hash01(seed::Int, keys::Int...) = begin
    h = UInt64(seed) * 0x9E3779B97F4A7C15
    for k in keys
        h = (h ⊻ (UInt64(k) + 0x9E3779B97F4A7C15 + (h << 6) + (h >> 2))) *
            0x5851F42D4C957F2D
    end
    Float64(h >> 11) / Float64(UInt64(1) << 53)
end

"""
    valley_data(spec; nregions)

Build the constant data of a synthetic valley. The numbers are plausible orders
of magnitude for a hydro valley (Mm3 of storage, m3/s of flow, EUR/MWh of price)
but they are generated, not measured, and describe no real site.

`spec.bellman_slopes` sets how many distinct slopes the function has, which is
not the number of cut rows; see `ValleySpec`.

`nregions` splits the future Bellman function into that many concave regions.
The published timing campaign uses one region. `nregions > 1` provides an
optional experiment with one persistent LP per region and a scan of reachable
regions at each stock point.
"""
function valley_data(spec::ValleySpec; nregions::Int=1)
    NR, NU, NT, NS = nreservoirs(spec), nunits(spec), spec.npostes, spec.nscenarios
    energetic = energetic_units(spec)
    ngroups = spec.unit_groups
    maxgroups = maximum(ngroups)
    seed = spec.seed

    # Postes of unequal length, as in a study whose time step mixes hours and
    # multi-hour blocks.
    duration = [t <= NT ÷ 3 ? 1.0 : (t <= 2NT ÷ 3 ? 4.0 : 8.0) for t in 1:NT]

    horizon = sum(duration)

    volume_min = zeros(NR, NT)
    volume_max = Matrix{Float64}(undef, NR, NT)
    for r in 1:NR, t in 1:NT
        volume_max[r, t] = 20.0 + 40.0 * _hash01(seed, 4, r)
    end

    # --- Sizing, top of the valley downwards -------------------------------
    #
    # A reservoir moves by a fraction of its usable volume over one period, and
    # that fraction is what decides whether the stock bounds bind. Sizing flows
    # against the reservoir controls the relevance of storage bounds.
    #
    # Two feasibility conditions have to hold at every point, or the benchmark
    # would be measuring infeasible solves:
    #
    #  * whatever enters a reservoir must be able to leave it, including the
    #    releases of the units above it. `unit_upstream[u] < unit_downstream[u]`
    #    in every preset, so one pass in increasing reservoir index suffices;
    #  * that must hold poste by poste and not only in total, because a full
    #    reservoir has no headroom to buffer a badly timed inflow. Inflows are
    #    therefore shaped like the release capacity below them and capped at
    #    `DRAIN_MARGIN` of it.
    release_rate = zeros(NR)           # volume per hour the units below r can pass
    drawing = [count(==(r), spec.unit_upstream) for r in 1:NR]
    for r in 1:NR
        own = spec.release_fraction * volume_max[r, NT] / horizon
        upstream = sum(release_rate[q] for q in 1:(r - 1)
                       if any(u -> spec.unit_upstream[u] == q &&
                                   spec.unit_downstream[u] == r, 1:NU); init=0.0)
        # r must pass its own inflow, at most DRAIN_MARGIN of its own release,
        # plus everything the units above it send down.
        release_rate[r] = max(own, 1.05 * upstream / (1.0 - DRAIN_MARGIN))
    end

    flow_min = zeros(NU, NT)
    flow_max = Matrix{Float64}(undef, NU, NT)
    group_min = zeros(NU, maxgroups, NT)
    group_max = zeros(NU, maxgroups, NT)
    efficiency = zeros(NU, maxgroups, NT)
    for u in 1:NU
        upstream = spec.unit_upstream[u]
        share = release_rate[upstream] / max(1, drawing[upstream])
        for t in 1:NT
            # Availability varies over the postes but never reaches zero: a hole
            # in the release capacity of one reservoir, facing a unit above it
            # that is still running, would make a full reservoir infeasible.
            available = MINIMUM_AVAILABILITY +
                (1.0 - MINIMUM_AVAILABILITY) * _hash01(seed, 6, u, t)
            flow_max[u, t] = share * available
        end
        ngroups[u] == 0 && continue
        for g in 1:ngroups[u], t in 1:NT
            group_max[u, g, t] = flow_max[u, t] / ngroups[u]
            # Decreasing marginal efficiency across groups: the objective is
            # concave in the turbined flow, as a real unit's power curve is.
            efficiency[u, g, t] = (0.30 + 0.25 * _hash01(seed, 7, u)) *
                (1.0 - 0.06 * (g - 1))
        end
    end

    # Inflows follow the shape of the release capacity below them, so that water
    # never arrives at a poste where nothing can pass it. The requested total is
    # honoured whenever it stays under `DRAIN_MARGIN` of that capacity, and
    # capped at it otherwise to keep the generated stock points feasible.
    inflow = Array{Float64,3}(undef, NR, NT, NS)
    releasable = Matrix{Float64}(undef, NR, NT)
    for r in 1:NR, t in 1:NT
        releasable[r, t] = duration[t] *
            sum(flow_max[u, t] for u in 1:NU if spec.unit_upstream[u] == r; init=0.0)
    end
    low_inflow, high_inflow = inflow_band(spec)
    for s in 1:NS, r in 1:NR
        capacity = sum(releasable[r, t] for t in 1:NT)
        target = (low_inflow + (high_inflow - low_inflow) * _hash01(seed, 2, r, s)) *
                 volume_max[r, NT]
        scale = capacity > 0 ? min(target / capacity, DRAIN_MARGIN) : 0.0
        for t in 1:NT
            inflow[r, t, s] = scale * releasable[r, t]
        end
    end

    price = Matrix{Float64}(undef, NT, NS)
    for s in 1:NS, t in 1:NT
        price[t, s] = 25.0 + 55.0 * _hash01(seed, 8, t, s) +
            15.0 * sin(2pi * t / NT)
    end

    reference_head = [60.0 + 120.0 * _hash01(seed, 9, r) for r in 1:NR]

    # A concave, piecewise-linear future Bellman function over reservoir 1:
    # increasing values with strictly decreasing slopes.
    ncuts = spec.ncuts
    top = volume_max[1, NT]
    cut_stock = collect(range(0.0, top; length=ncuts + 1))
    cut_value = Vector{Float64}(undef, ncuts + 1)
    cut_value[1] = 0.0
    # The water value has to be on the same scale as the energy value, or the
    # Bellman term decides the whole LP on its own. One unit of volume released
    # at poste t is worth `price * efficiency`, because the water balance carries
    # the duration. The top of the Bellman function is set just below the best
    # such value, and the slope decays to about a third of it across the usable
    # volume. Its exact spread is an explicit synthetic experiment parameter.
    best_energy_value = maximum(price) * maximum(efficiency)
    # `bellman_slopes` distinct slopes spread over `ncuts` cut rows: consecutive
    # rows of a group are collinear, so they are redundant exactly as the facets
    # of one cell of a coarse 2D grid are. The stock axis keeps its `ncuts + 1`
    # breakpoints, so the row count and the matrix shape are unchanged.
    nslopes = spec.bellman_slopes <= 0 ? ncuts : min(spec.bellman_slopes, ncuts)
    top_slope = TOP_WATER_VALUE * best_energy_value
    decay = WATER_VALUE_SPREAD^(-1 / max(nslopes - 1, 1))
    for i in 1:ncuts
        group = div((i - 1) * nslopes, ncuts)          # 0 .. nslopes-1
        slope = top_slope * decay^group
        cut_value[i + 1] = cut_value[i] + slope * (cut_stock[i + 1] - cut_stock[i])
    end
    regions = concave_regions(ncuts, nregions)

    # Reservoirs that carry no Bellman cut are not piloted: their initial volume
    # is a study input and never moves from point to point.
    initial_volume = [0.5 * volume_max[r, NT] for r in 1:NR]

    max_release = [sum(duration[t] * flow_max[u, t]
                       for u in 1:NU if spec.unit_upstream[u] == r for t in 1:NT;
                       init=0.0) for r in 1:NR]
    total_inflow = [sum(inflow[r, t, s] for t in 1:NT) for r in 1:NR, s in 1:NS]

    ValleyData(spec, energetic, ngroups, duration, inflow, volume_min, volume_max,
        flow_min, flow_max, group_min, group_max, price, efficiency,
        reference_head, cut_stock, cut_value, regions, initial_volume,
        max_release, total_inflow)
end

"""
    concave_regions(ncuts, nregions)

Split the cut indices into `nregions` regions of roughly equal width.
Neighbouring regions overlap by one stock-grid interval. A single region spans
all `ncuts + 1` breakpoints.
"""
function concave_regions(ncuts::Int, nregions::Int)
    nregions >= 1 || throw(ArgumentError("need at least one concave region"))
    nregions == 1 && return [(1, ncuts + 1)]
    nregions <= ncuts ||
        throw(ArgumentError("cannot split $ncuts cuts into $nregions regions"))
    edges = round.(Int, range(1, ncuts + 1; length=nregions + 1))
    [(edges[k], k == nregions ? edges[end] : edges[k + 1] + 1) for k in 1:nregions]
end

"""
    bisection_fractions(n)

The stock levels a bisection grid refiner visits, in visit order: the two ends,
then the midpoint of the widest remaining interval, and so on. Consecutive
levels get closer and closer, so consecutive LPs differ by less and less — the
regime in which a warm basis is worth the most.
"""
function bisection_fractions(n::Int)
    n >= 2 || throw(ArgumentError("need at least the two endpoints"))
    grid = [0.0, 1.0]
    visited = [0.0, 1.0]
    while length(visited) < n
        widest = argmax(grid[i + 1] - grid[i] for i in 1:(length(grid) - 1))
        midpoint = (grid[widest] + grid[widest + 1]) / 2
        insert!(grid, widest + 1, midpoint)
        push!(visited, midpoint)
    end
    visited
end

"""
    valley_points(data; nscenarios, stock_fractions, order, stocks_per_scenario)

The measured sequence. `order` decides **what changes between two consecutive
points**, which is what a warm basis actually sees:

* `:scenario_fast` — a fixed stock, scanning every scenario. Every step moves
  all the inflows and all the profits: the largest update per point.
* `:stock_fast` — a fixed scenario, scanning the stock levels. Most steps move
  one bound only, and the scenario changes once per stock sweep.
* `:refinement` — a fixed scenario, visiting `stocks_per_scenario` levels in
  bisection order. This is the order of the production 1D grid refiner, and the
  one where consecutive LPs are closest.

`stock_fractions` is used by the first two orders; `:refinement` generates its
own levels with [`bisection_fractions`](@ref).
"""
function valley_points(data::ValleyData; nscenarios::Int=data.spec.nscenarios,
        stock_fractions=(0.25, 0.5, 0.75), order::Symbol=:scenario_fast,
        stocks_per_scenario::Int=length(stock_fractions))
    NC = data.spec.ncut_reservoirs
    top = ntuple(r -> data.volume_max[r, data.spec.npostes], NC)
    point(scenario, fraction) =
        Point{NC}(scenario, ntuple(r -> fraction * top[r], NC), 1)
    if order === :scenario_fast
        return [point(scenario, fraction)
                for fraction in stock_fractions for scenario in 1:nscenarios]
    elseif order === :stock_fast
        return [point(scenario, fraction)
                for scenario in 1:nscenarios for fraction in stock_fractions]
    elseif order === :refinement
        fractions = bisection_fractions(stocks_per_scenario)
        return [point(scenario, fraction)
                for scenario in 1:nscenarios for fraction in fractions]
    end
    throw(ArgumentError("unknown point order: $order"))
end

# ---------------------------------------------------------------------------
# Point-dependent quantities
#
# These are the single source of truth for every front-end. The ReactiveModel
# blocks, the JuMP model, the MathOptInterface model and the two POI models all
# call exactly these functions, so a timing difference is a difference between
# modeling layers and never a difference in the numbers being pushed.
# ---------------------------------------------------------------------------

"""
The head of a deep reservoir varies by a modest fraction of its reference value
over the whole usable volume. Making profits swing by a factor of two with the
stock, as a naive scaling does, couples the objective to the initial stock far
more tightly than any real site and destroys the warm basis at every stock step.
"""
const HEAD_VARIATION = 0.15

"""
    head(data, point, reservoir)

The head at the start of the period: the reservoir's reference head, modulated by
[`HEAD_VARIATION`](@ref) according to how full it is. Profits and, optionally,
flow limits depend on it, which is what makes the objective move with the
initial stock.
"""
@inline initial_level(data::ValleyData, point::Point, r::Int) =
    r <= data.spec.ncut_reservoirs ? point.stocks[r] : data.initial_volume[r]

@inline head(data::ValleyData, point::Point, r::Int) =
    data.reference_head[r] * (1.0 - HEAD_VARIATION +
        HEAD_VARIATION * initial_level(data, point, r) /
            data.volume_max[r, data.spec.npostes])

@inline head_factor(data::ValleyData, point::Point, u::Int) =
    head(data, point, data.spec.unit_upstream[u]) /
    data.reference_head[data.spec.unit_upstream[u]]

@inline initial_volume_bound(data::ValleyData, point::Point, r::Int) =
    initial_level(data, point, r)

@inline volume_bounds(data::ValleyData, ::Point, r::Int, t::Int) =
    (data.volume_min[r, t], data.volume_max[r, t])

@inline function flow_bounds(data::ValleyData, point::Point, u::Int, t::Int)
    scale = data.spec.head_dependent_bounds ? head_factor(data, point, u) : 1.0
    (data.flow_min[u, t] * scale, data.flow_max[u, t] * scale)
end

@inline function group_flow_bounds(data::ValleyData, point::Point, u::Int, g::Int, t::Int)
    scale = data.spec.head_dependent_bounds ? head_factor(data, point, u) : 1.0
    (data.group_min[u, g, t] * scale, data.group_max[u, g, t] * scale)
end

@inline inflow_value(data::ValleyData, point::Point, r::Int, t::Int) =
    data.inflow[r, t, point.scenario]

"""
    final_volume_target(data, point, r)

A reservoir that carries no Bellman cut has its final volume fixed, at the level
it started the period from. That level is a study input, so this right-hand side
is constant across points.
"""
@inline final_volume_target(data::ValleyData, ::Point, r::Int) =
    data.initial_volume[r]

@inline profit(data::ValleyData, point::Point, u::Int, g::Int, t::Int) =
    data.price[t, point.scenario] * data.efficiency[u, g, t] *
    data.duration[t] * head_factor(data, point, u)

cut_slope(data::ValleyData, s::Int) =
    (data.cut_value[s + 1] - data.cut_value[s]) /
    (data.cut_stock[s + 1] - data.cut_stock[s])

"""
    cut_slope(data, s, r)

Slope of cut `s` with respect to the final volume of cut-reservoir `r`. The
first cut reservoir carries the tabulated slope; the others carry a fixed
fraction of it, which keeps the row concave in the aggregate storage while
giving it the arity of a 2D cut.
"""
cut_slope(data::ValleyData, s::Int, r::Int) =
    r == 1 ? cut_slope(data, s) : cut_slope(data, s) * (0.65 / r)

"""
    cut_rhs(data, point, s)

The intercept of Bellman cut `s`, or `Inf` when the cut lies outside the concave
region selected by the point. Returning `Inf` neutralizes the row instead of
removing it, preserving a fixed matrix layout across points.
"""
function cut_rhs(data::ValleyData, point::Point, s::Int)
    a, b = data.regions[point.region]
    a <= s < b || return Inf
    data.cut_value[s] - cut_slope(data, s) * data.cut_stock[s]
end

min_final_stock(data::ValleyData, point::Point) =
    data.cut_stock[data.regions[point.region][1]]
max_final_stock(data::ValleyData, point::Point) =
    data.cut_stock[data.regions[point.region][2]]

"""
    reachable_stock_box(data, point)

The interval the final stock of the first cut reservoir can reach over the
period, from mass balance and the flow limits: it cannot fall by more than the
whole turbinable volume, nor rise by more than the total inflow.
"""
@inline function reachable_stock_box(data::ValleyData, point::Point)
    top = data.volume_max[1, data.spec.npostes]
    initial = point.stocks[1]
    lower = max(0.0, initial - data.max_release[1])
    upper = min(top, initial + data.total_inflow[1, point.scenario])
    (lower, upper)
end

"""
    reachable_regions(data, point) -> UnitRange{Int}

Return the synthetic concave regions whose final-stock support meets the
reachable box. Both supports are ordered intervals, so the reachable set is
contiguous. The benchmark driver falls back to all regions when this selection
is empty.
"""
function reachable_regions(data::ValleyData, point::Point)
    lower, upper = reachable_stock_box(data, point)
    n = length(data.regions)
    first_region = findfirst(1:n) do k
        data.cut_stock[data.regions[k][2]] >= lower
    end
    last_region = findlast(1:n) do k
        data.cut_stock[data.regions[k][1]] <= upper
    end
    (first_region === nothing || last_region === nothing ||
        first_region > last_region) && return 1:0
    first_region:last_region
end

# Change predicates, shared by every front-end so that all of them skip exactly
# the same work between two points.
@inline stocks_changed(old::Point, new::Point) = old.stocks != new.stocks
@inline scenario_changed(old::Point, new::Point) = old.scenario != new.scenario
@inline region_changed(old::Point, new::Point) = old.region != new.region
@inline flow_bounds_changed(data::ValleyData, old::Point, new::Point) =
    data.spec.head_dependent_bounds && stocks_changed(old, new)
@inline profits_changed(old::Point, new::Point) =
    scenario_changed(old, new) || stocks_changed(old, new)

"""
    rhs_changed(old, new)

Whether any right-hand side can have moved: inflows follow the scenario, and the
Bellman cuts and the final-stock window follow the concave region. The flat
front-ends refresh their whole right-hand-side vector on this predicate, so it
must be the union of what the ReactiveModel blocks react to — otherwise they
would silently keep a stale value that ReactiveModel updates.
"""
@inline rhs_changed(old::Point, new::Point) =
    scenario_changed(old, new) || region_changed(old, new)

# ---------------------------------------------------------------------------
# Flat layout, shared by the JuMP, MOI and POI front-ends
# ---------------------------------------------------------------------------

"""
    Layout

Column and row numbering of the flat LP, in the same block order as the
ReactiveModel formulation. `rows` holds one `(indices, coefficients, op, key)`
tuple per row; `op` is `:eq`, `:le` or `:ge`.
"""
struct Layout
    ncols::Int
    nrows::Int
    initial_volume::Vector{Int}                       # [r]
    volume::Matrix{Int}                               # [r, t]
    debit::Matrix{Int}                                # [u, t]
    group_debit::Array{Int,3}                         # [u, g, t], 0 when absent
    vB::Int
    rows::Vector{Tuple{Vector{Int},Vector{Float64},Symbol}}
    water_balance::Matrix{Int}                        # [r, t] -> row
    final_volume::Vector{Int}                         # [r] -> row, 0 when absent
    group_sum::Matrix{Int}                            # [u, t] -> row, 0 when absent
    bellman_cut::Vector{Int}                          # [s] -> row
    inactive_cut::Int
    min_final::Vector{Int}                            # [cut reservoir] -> row
    max_final::Vector{Int}                            # [cut reservoir] -> row
    inactive_triangle::Vector{Int}
    objective_columns::Vector{Int}                    # columns carrying a cost
    objective_keys::Vector{Tuple{Int,Int,Int}}        # (u, g, t) for each of them
end

function layout(data::ValleyData)
    spec = data.spec
    NR, NU, NT = nreservoirs(spec), nunits(spec), spec.npostes
    maxgroups = maximum(data.ngroups)
    column = 0
    initial_volume = [column += 1 for _ in 1:NR]
    volume = Matrix{Int}(undef, NR, NT)
    for t in 1:NT, r in 1:NR
        volume[r, t] = (column += 1)
    end
    debit = Matrix{Int}(undef, NU, NT)
    for t in 1:NT, u in 1:NU
        debit[u, t] = (column += 1)
    end
    group_debit = zeros(Int, NU, maxgroups, NT)
    for t in 1:NT, u in data.energetic, g in 1:data.ngroups[u]
        group_debit[u, g, t] = (column += 1)
    end
    vB = (column += 1)
    ncols = column

    rows = Tuple{Vector{Int},Vector{Float64},Symbol}[]
    water_balance = Matrix{Int}(undef, NR, NT)
    for t in 1:NT, r in 1:NR
        indices = Int[]; coefficients = Float64[]
        push!(indices, t == 1 ? initial_volume[r] : volume[r, t - 1])
        push!(coefficients, -1.0)
        push!(indices, volume[r, t]); push!(coefficients, 1.0)
        for u in 1:NU
            spec.unit_upstream[u] == r || continue
            push!(indices, debit[u, t]); push!(coefficients, data.duration[t])
        end
        for u in 1:NU
            spec.unit_downstream[u] == r || continue
            push!(indices, debit[u, t]); push!(coefficients, -data.duration[t])
        end
        push!(rows, (indices, coefficients, :eq))
        water_balance[r, t] = length(rows)
    end
    final_volume = zeros(Int, NR)
    for r in (spec.ncut_reservoirs + 1):NR
        push!(rows, ([volume[r, NT]], [1.0], :eq))
        final_volume[r] = length(rows)
    end
    group_sum = zeros(Int, NU, NT)
    for t in 1:NT, u in data.energetic
        indices = [debit[u, t]]; coefficients = [-1.0]
        for g in 1:data.ngroups[u]
            push!(indices, group_debit[u, g, t]); push!(coefficients, 1.0)
        end
        push!(rows, (indices, coefficients, :eq))
        group_sum[u, t] = length(rows)
    end
    NC = spec.ncut_reservoirs
    bellman_cut = Vector{Int}(undef, spec.ncuts)
    for s in 1:spec.ncuts
        indices = [vB]; coefficients = [1.0]
        for r in 1:NC
            push!(indices, volume[r, NT])
            push!(coefficients, -cut_slope(data, s, r))
        end
        push!(rows, (indices, coefficients, :le))
        bellman_cut[s] = length(rows)
    end
    push!(rows, ([volume[1, NT]], [1.0], :le)); inactive_cut = length(rows)
    min_final = Vector{Int}(undef, NC); max_final = Vector{Int}(undef, NC)
    for r in 1:NC
        push!(rows, ([volume[r, NT]], [1.0], :ge)); min_final[r] = length(rows)
    end
    for r in 1:NC
        push!(rows, ([volume[r, NT]], [1.0], :le)); max_final[r] = length(rows)
    end
    # Include NC^2 inactive rows as an explicit fixed-layout stress case.
    inactive_triangle = Int[]
    for k in 1:(NC * NC)
        push!(rows, ([volume[mod1(k, NC), NT]], [1.0], :le))
        push!(inactive_triangle, length(rows))
    end

    objective_columns = Int[]
    objective_keys = Tuple{Int,Int,Int}[]
    for t in 1:NT, u in data.energetic, g in 1:data.ngroups[u]
        push!(objective_columns, group_debit[u, g, t])
        push!(objective_keys, (u, g, t))
    end

    Layout(ncols, length(rows), initial_volume, volume, debit, group_debit, vB,
        rows, water_balance, final_volume, group_sum, bellman_cut, inactive_cut,
        min_final, max_final, inactive_triangle, objective_columns, objective_keys)
end

"""
    column_bounds!(lower, upper, data, layout, point)

Fill the flat bound vectors for `point`. `vB` is free.
"""
function column_bounds!(lower, upper, data::ValleyData, l::Layout, point::Point)
    spec = data.spec
    for r in 1:nreservoirs(spec)
        j = l.initial_volume[r]
        lower[j] = upper[j] = initial_volume_bound(data, point, r)
        for t in 1:spec.npostes
            lower[l.volume[r, t]], upper[l.volume[r, t]] =
                volume_bounds(data, point, r, t)
        end
    end
    for u in 1:nunits(spec), t in 1:spec.npostes
        lower[l.debit[u, t]], upper[l.debit[u, t]] = flow_bounds(data, point, u, t)
    end
    for u in data.energetic, g in 1:data.ngroups[u], t in 1:spec.npostes
        j = l.group_debit[u, g, t]
        lower[j], upper[j] = group_flow_bounds(data, point, u, g, t)
    end
    lower[l.vB] = -Inf; upper[l.vB] = Inf
    return nothing
end

"""
    row_rhs!(rhs, data, layout, point)

Fill the flat right-hand-side vector for `point`.
"""
function row_rhs!(rhs, data::ValleyData, l::Layout, point::Point)
    spec = data.spec
    for r in 1:nreservoirs(spec), t in 1:spec.npostes
        rhs[l.water_balance[r, t]] = inflow_value(data, point, r, t)
    end
    for r in (spec.ncut_reservoirs + 1):nreservoirs(spec)
        rhs[l.final_volume[r]] = final_volume_target(data, point, r)
    end
    for u in data.energetic, t in 1:spec.npostes
        rhs[l.group_sum[u, t]] = 0.0
    end
    for s in 1:spec.ncuts
        rhs[l.bellman_cut[s]] = cut_rhs(data, point, s)
    end
    rhs[l.inactive_cut] = Inf
    for r in eachindex(l.min_final)
        # The other cut reservoirs share the region's bounds, rescaled onto their
        # own usable volume.
        scale = data.volume_max[r, spec.npostes] / data.volume_max[1, spec.npostes]
        rhs[l.min_final[r]] = scale * min_final_stock(data, point)
        rhs[l.max_final[r]] = scale * max_final_stock(data, point)
    end
    for row in l.inactive_triangle
        rhs[row] = Inf
    end
    return nothing
end

"""
    objective_costs!(costs, data, layout, point)

Fill the cost of every column carrying one, in `layout.objective_columns` order.
The `vB` column always has a cost of one and is handled separately.
"""
function objective_costs!(costs, data::ValleyData, l::Layout, point::Point)
    for i in eachindex(l.objective_keys)
        u, g, t = l.objective_keys[i]
        costs[i] = profit(data, point, u, g, t)
    end
    return nothing
end

"""
    describe(data)

Column, row and nonzero counts for the generated LP.
"""
function describe(data::ValleyData)
    l = layout(data)
    (; name=data.spec.name, ncols=l.ncols, nrows=l.nrows,
       nnz=sum(length(first(row)) for row in l.rows),
       nreservoirs=nreservoirs(data.spec), nunits=nunits(data.spec),
       nenergetic=length(data.energetic), npostes=data.spec.npostes,
       ncuts=data.spec.ncuts, nscenarios=data.spec.nscenarios,
       head_dependent_bounds=data.spec.head_dependent_bounds)
end

end # module
