# Contributing

ReactiveModel is experimental. Changes should state the mathematical behavior,
usage contract and validation scope they affect. Start from a source checkout
with Julia 1.12 or later; package checks have run on Julia 1.12.7 and 1.13.2.
Recorded compilation and performance measurements use Julia 1.13.2.

## Validate a change

From the repository root:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test(; allow_reresolve=false)'
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```

The documentation build executes examples and checks exported docstrings and
internal references. For numerical changes, check the affected incremental
path against a fresh problem, include failure/recovery where relevant, and run
the main suite. Optional integrations require their own solver environment;
an absent solver is not a successful integration check.

The prepared GitHub Actions workflow runs package tests on Julia 1.12 and the
current stable Julia release, on Linux, macOS and Windows. Its documentation job
executes the manual with four Julia threads and checks local links. These jobs
resolve fresh dependency environments from `Project.toml`; the checked-in
manifests remain the records used by the corresponding local measurements.
The workflow will run when this repository is hosted on GitHub; configuring it
locally does not establish that all platform jobs pass.

Run the same checks locally, using the Julia version you want to validate:

```sh
julia --startup-file=no ci/check.jl test
JULIA_NUM_THREADS=4 julia --startup-file=no ci/check.jl docs
```

The script copies the required source into a temporary directory, without a
manifest or adjacent solver checkout, then resolves and tests that copy. It
does not change dependency environments in the working tree. A populated depot
may supply cached packages and artifacts; this is not an empty-depot test.
Recorded trimming commands remain a separate validation in the trimming guide.

Performance investigations use the separate `benchmarks` environment. See
[the protocol](docs/src/benchmarks.md). Retain raw samples, numerical checks,
source revisions and package versions with any reported timing. Compilation
changes should also follow [the trimming guide](docs/src/guide/trimming.md).

## Write documentation

- Lead with the user task, then a complete example and its expected result.
- Execute numerical examples with Documenter `@example` blocks and assertions.
  Label schematic API snippets and installation commands as such.
- State the required environment, working directory and optional dependencies.
- Distinguish an API contract, an observed measurement and a hypothesis.
- Name the exact workload and environment for a performance claim. Do not infer
  native allocations from Julia allocation counts, or GC causality from timing.
- Describe other tools through their documented capabilities and cite primary
  documentation. Comparisons should identify implementations and settings,
  rather than characterize an entire library by one benchmark path.
- Keep historical measurements dated. Correct context without rewriting raw
  evidence or presenting old results as measurements of current code.

## Report a problem

A source snapshot for publication review can be prepared with
`python3 ci/export_publication.py /absolute/new-or-empty-directory`, choosing
an output outside the checkout. It includes current code, documentation,
synthetic benchmark generators and CI, with a SHA-256 inventory. It excludes
Git metadata, private audit reports and generated build/cache files. Review the
content and reproducibility inputs before creating a public repository from
that snapshot; exporting files does not certify their provenance or rewrite
the original history.

Include a small runnable model and data, the expected and observed numerical
results or error, Julia and package versions, backend/options, source revision
and local changes. For cached or incremental behavior, include the sequence of
points and any study-data mutations and rebuilds.
