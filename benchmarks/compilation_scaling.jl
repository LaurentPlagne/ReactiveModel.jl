#=
Client compilation diagnostic, one fresh Julia process per configuration.

    julia --startup-file=no --project=. benchmarks/compilation_scaling.jl OUTPUT FAMILIES [KEYS=1]

Generates a fixed family of independent LP rows with the DSL. Measures client
definition, problem construction, first solve/update and a warmed sweep separately.
Use identical Julia/dependency versions and retain the source for comparisons.
=#
using TOML
load_sample=@timed @eval using ReactiveModel

function scaling_model(n)
    body=Expr(:block)
    for i in 1:n
        x=Symbol(:x, i)
        row=Symbol(:floor, i)
        push!(body.args, :(@variable $x[j in 1:data.width] in [0.0, data.capacity]))
        push!(body.args, :(@constraint $row[j in 1:data.width] point.factor*$x[j] >= point.rhs))
    end
    terms=[:(sum(point.cost*$(Symbol(:x, i))[j] for j in 1:data.width)) for i in 1:n]
    push!(body.args, :(@objective Min $(Expr(:call, :+, terms...))))
    return Expr(
        :macrocall, Symbol("@model"), LineNumberNode(1), :(CompilationSweep(data, point)), body
    )
end

function fields(sample)
    result=Dict{String,Any}(
        "seconds"=>sample.time, "julia_bytes"=>sample.bytes, "gc_seconds"=>sample.gctime
    )
    hasproperty(sample, :compile_time) && (result["compile_seconds"]=sample.compile_time)
    return result
end

function warm_sweep!(lp, n, width)
    total=0.0
    for i in 1:20
        point=(rhs=1.0+0.01*i, factor=2.0, cost=4.0)
        value=objective_value(lp, point)
        @assert isapprox(value, 2*n*width*point.rhs; rtol=1e-8)
        total+=value
    end
    return total
end

# A client aggregate makes propagation into application types visible. Record
# actual MethodError text with RM's compact display disabled, after all timings.
struct CompilationClient{T}
    problem::T
end
_compilation_requires_integer(::Int) = nothing
_compilation_client_failure(client) = _compilation_requires_integer(client)

function full_diagnostic(f, value)
    io=IOBuffer()
    f(IOContext(io, :limit=>false, :reactivemodel_full_types=>true), value)
    return String(take!(io))
end

function main(args)
    output=args[1]
    n=parse(Int, args[2])
    width=length(args)>2 ? parse(Int, args[3]) : 1
    n>0 && width>0 || error("family and key counts must be positive")
    expression=scaling_model(n)
    definition=@timed Core.eval(Main, expression)
    data=(width=width, capacity=10.0)
    p=(rhs=1.0, factor=1.0, cost=2.0)
    q=(rhs=3.0, factor=2.0, cost=4.0)
    # invokelatest crosses the world-age boundary from the generated client model.
    construction=@timed Base.invokelatest() do
        return LinearProblem(ConstData(CompilationSweep(data), data), p; backend=ClpBackend())
    end
    lp=construction.value
    first_solve=@timed Base.invokelatest(objective_value, lp, p)
    @assert isapprox(first_solve.value, 2*n*width; rtol=1e-8)
    first_update=@timed Base.invokelatest(objective_value, lp, q)
    @assert isapprox(first_update.value, 6*n*width; rtol=1e-8)
    Base.invokelatest(warm_sweep!, lp, n, width)
    warm=@timed Base.invokelatest(warm_sweep!, lp, n, width)
    # Run through a typed function barrier; exclude invokelatest dispatch boxing.
    warm_bytes=Base.invokelatest(
        (lp, n, width)->(@allocated warm_sweep!(lp, n, width)), lp, n, width
    )
    @assert warm_bytes==0
    # Measure the full compiler type, independently of the compact object display.
    first_result=@timed Base.invokelatest(solve, lp, q)
    solution=first_result.value
    first_collections=@timed Base.invokelatest(s->(s.variables, s.constraints), solution)
    client=CompilationClient(lp)
    diagnostic_error=try
        _compilation_client_failure(client)
    catch error
        error isa MethodError || rethrow()
        error
    end
    diagnostics=Dict(
        "problem_type"=>full_diagnostic(show, typeof(lp)),
        "workspace_type"=>full_diagnostic(show, typeof(lp.workspace)),
        "solution_type"=>full_diagnostic(show, typeof(solution)),
        "client_type"=>full_diagnostic(show, typeof(client)),
        "method_error"=>full_diagnostic(showerror, diagnostic_error),
    )
    report=Dict(
        "metadata"=>Dict(
            "julia"=>string(VERSION),
            "cpu"=>Sys.cpu_info()[1].model,
            "families"=>n,
            "keys_per_family"=>width,
            "columns"=>n*width,
            "rows"=>n*width,
            "threads"=>Threads.nthreads(),
            "backend"=>"Clp_jll",
            "safe_mode"=>false,
            "problem_type_name_bytes"=>sizeof(diagnostics["problem_type"]),
            "workspace_type_name_bytes"=>sizeof(diagnostics["workspace_type"]),
            "solution_type_name_bytes"=>sizeof(diagnostics["solution_type"]),
            "client_type_name_bytes"=>sizeof(diagnostics["client_type"]),
            "method_error_bytes"=>sizeof(diagnostics["method_error"]),
            "timing_scope"=>"fresh client model; package load recorded separately; warmed sweep is 20 RHS updates",
        ),
        "load"=>fields(load_sample),
        "definition"=>fields(definition),
        "construction"=>fields(construction),
        "first_solve"=>fields(first_solve),
        "first_update"=>fields(first_update),
        "first_result"=>fields(first_result),
        "first_collections"=>fields(first_collections),
        "warm"=>fields(warm),
        "warm_kernel_bytes"=>warm_bytes,
        "diagnostics"=>diagnostics,
    )
    open(output, "w") do io
        return TOML.print(io, report; sorted=true)
    end
    return println(
        "families=",
        n,
        " keys=",
        width,
        " construction=",
        construction.time,
        " first solve=",
        first_solve.time,
        " first update=",
        first_update.time,
        " warm Julia bytes=",
        warm_bytes,
    )
end

main(ARGS)
