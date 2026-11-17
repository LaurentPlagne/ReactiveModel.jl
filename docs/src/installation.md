# Installation and reproduction

The instructions below start from a ReactiveModel source checkout. Run shell
commands from its root directory. A public registry installation is not assumed.
You need Julia 1.12 or later according to `Project.toml`. Package checks have
passed on Julia 1.12.7 and 1.13.2 on macOS aarch64; the recorded compilation and
performance measurements use Julia 1.13.2.

## Run the example

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. examples/getting_started.jl
```

The script asserts costs of `4600.0` and `5300.0` and second-step production
`[80.0, 30.0]`. It exits successfully only if those checks pass. Follow
[the first tutorial](tutorials/first-model.md) for the mathematical explanation.

HiGHS and Clp are installed through Julia binary packages. Neither MOI nor an
adjacent solver checkout is needed for the package. Optional solver extensions
have separate prerequisites; see [backends](backends.md).

To use this checkout from another project, replace the path below with its
absolute location:

```julia
using Pkg
Pkg.develop(path="/absolute/path/to/ReactiveModel.jl")
Pkg.instantiate()
using ReactiveModel
```

This is a Julia setup snippet; it changes the active consumer project's
environment. Keep benchmark dependencies in their separate project.

## Run package checks

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

The main suite checks HiGHS and Clp. Tests requiring an absent optional solver
are reported as skipped/broken; that does not validate those integrations.

## Build and view the manual

```sh
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```

Documenter executes the `@example` blocks, including numerical assertions.
Errors in examples, unresolved references and missing exported docstrings fail
the build. The HTML output is `docs/build/index.html`; open it in a browser.
An optional local server is:

```sh
python3 -m http.server 8000 --directory docs/build
```

Then visit `http://localhost:8000`. No site deployment is performed by the build.
Source and edit links are disabled while no public repository URL is configured;
the output contains no machine-specific `file://` source links.

## Further examples and evidence

The [SDDP tutorial](tutorials/sddp.md) includes a complete small training loop:

```sh
julia --startup-file=no --project=. examples/sddp.jl
```

The tutorial's forward/backward passes, dual-based cuts and bound checks are
also executed in the documentation build. The [synthetic benchmark workflow](parallel-results.md)
provides commands for generating numerical checks and local measurement records.

The tutorial drawings are editable SVGs tracked with the manual. To regenerate
them, run `python3 docs/draw_tutorials.py` from the repository root; no Python
packages are required. This regeneration is optional when building the manual.
The ownership, solve-path and result-lifetime figures are generated separately
with `python3 docs/draw_mechanisms.py`. They use the same dependency-free SVG
approach. On narrow screens, focus or scroll a mechanism figure horizontally
to read its labels; its text alternative and surrounding explanation give the
same contract.

## Record the environment

Keep the source revision, local changes, Julia version, OS, architecture and
manifest with a reproduction report. `Project.toml` describes allowed versions;
a manifest records the versions selected for an environment. A different Julia
version may require resolution and produce a different manifest.

The existing installation checks reused a populated Julia depot. They do not
establish downloads from an empty depot or support on every platform. For timed
experiments, follow the separate [benchmark protocol](benchmarks.md). For an
executable application, follow [compilation and trimming](guide/trimming.md).
