using Documenter
using ReactiveModel

DocMeta.setdocmeta!(
    ReactiveModel,
    :DocTestSetup,
    :(using ReactiveModel);
    recursive=true,
)

makedocs(
    root=@__DIR__,
    sitename="ReactiveModel.jl",
    modules=[ReactiveModel],
    checkdocs=:exports,
    remotes=nothing,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://laurentplagne.github.io/ReactiveModel.jl/",
        assets=["assets/style.css"],
        edit_link=nothing,
        repolink="https://github.com/LaurentPlagne/ReactiveModel.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Installation and reproduction" => "installation.md",
        "Tutorials" => [
            "Your first model" => "tutorials/first-model.md",
            "Build a model with functions" => "tutorials/first-model-functional.md",
            "Business functions and safe updates" => "tutorials/business-functions.md",
            "Indexed families" => "tutorials/indexed-families.md",
            "A reservoir stage problem" => "tutorials/bellman-recursion.md",
            "A small SDDP training loop" => "tutorials/sddp.md",
        ],
        "Context and measurement" => [
            "Modeling approaches" => "why/ergonomics.md",
            "Measuring repeated solves" => "why/performance.md",
        ],
        "Backends" => "backends.md",
        "Guides" => [
            "Data and dependencies" => "guide/concepts.md",
            "DSL grammar" => "guide/dsl.md",
            "Results and the solve cycle" => "guide/solve-cycle.md",
            "Compilation and trimming" => "guide/trimming.md",
            "Design and compilation costs" => "guide/design.md",
        ],
        "Benchmarks" => [
            "Serial and parallel protocol" => "parallel-results.md",
            "Synthetic input specification" => "benchmark-inputs.md",
            "Timing diagnostics" => "benchmark-diagnostics.md",
            "Compilation protocol" => "compilation-results.md",
            "The JuMP comparison code" => "benchmarks-jump.md",
            "Method and reproduction" => "benchmarks.md",
        ],
        "FAQ" => "faq.md",
        "API reference" => "api.md",
    ],
)
