# Synthetic benchmark inputs

The benchmark constructs cascade LPs from explicit graphs, integer dimensions,
seeds and coefficient formulas. It generates the study arrays in code and
requires no operating-site dataset.

| Choice | Small | Medium | Large |
| :--- | ---: | ---: | ---: |
| Reservoirs | 1 | 3 | 5 |
| Time steps | 8 | 12 | 16 |
| Units | 2 | 5 | 9 |
| Cut rows | 16 | 24 | 64 |
| Cut reservoirs | 1 | 1 | 2 |
| Distinct cut slopes | 16 | 6 | 8 |
| Available scenarios | 64 | 64 | 64 |
| Seed | 11 | 22 | 33 |

Small has two outlets from one reservoir. Medium has a three-node chain and
two outlet bypasses. Large has a five-node chain, one outlet and four parallel
bypasses. Groups have one or two pieces; non-generating bypass units have zero
groups.

The base discharge budget is half the generated storage per horizon. Downstream
capacity can be increased to absorb upstream releases. Natural inflow ranges
from 0.02 to 0.1 times storage, subject to a feasibility cap. Terminal water
value starts at 0.8 of the best generated energy value, with slopes decreasing
by a factor of four. These ratios are explicit experiment parameters.

Storage, availability, prices and efficiencies follow the deterministic hash
and formulas in `synthetic_valley.jl`. Durations are blocks of one, four and
eight hours; capacity and coefficient scales are illustrative. The feasibility
construction sizes downstream releases to absorb upstream flow.

The LP includes fixed inactive rows to exercise infinite bounds, changing stock
fixings, scenario-dependent RHS/costs and an optional head-dependent bound.
The default benchmark sequence selects one concave region. Change the point
sequence explicitly to study other update patterns.

This workload covers repeated solves with fixed families and cuts. Rebuilds,
cut insertion and full dynamic-programming or SDDP costs require separate
experiments. Retain the input parameters, source and environment with any
results you produce.
