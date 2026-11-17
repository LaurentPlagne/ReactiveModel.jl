# Diagnosing adapter and scheduler costs

The [synthetic benchmark workflow](parallel-results.md) includes an explicit
JuMP/Clp inner-cache reset alongside ordinary updates. Both paths can be checked
and timed on the same generated points, making the adapter choice explicit.

The reset detaches the inner MOI cache before setters. It still reloads the
complete native LP during optimization. Access to the nested cache is specific
to the adapter layout used by the script; it is not a general JuMP API recommendation.
See the [comparison source](benchmarks-jump.md).

For a new investigation, `benchmarks/diagnose_clp_cache.jl OUTPUT 11` compares
ordinary updates, inner reset and outer reset on the independent inputs. It
retains raw samples, numerical checks and qualitative profiles outside timing.
`benchmarks/diagnose_parallel_clp.jl OUTPUT 11` compares scheduler configurations
within one process. Run them with the benchmark project, eight Julia threads
for the scheduler probe, and separate new/empty output directories.

These probes produce new experiments. They do not validate a general statement
that one scheduler is faster, or that GC explains a difference. Keep thread
counts, spin policy, workload, solver history, environment and external load
under control; distinguish placement observations from a causal explanation.
