# Compare first use of two deterministic business functions in fresh processes:
# julia --startup-file=no --project=. benchmarks/reduction_compilation.jl OUTPUT sum|loop LENGTH
using TOML

function sum_rhs(study, point)
    return sum(
        (point.scale * max(weight - point.threshold, 0.0) for weight in study.weights); init=0.0
    )
end

function loop_rhs(study, point)
    total = 0.0
    for weight in study.weights
        total += point.scale * max(weight - point.threshold, 0.0)
    end
    return total
end

function warm_reduction(f::F, study, point) where {F}
    total = 0.0
    for i in 1:20
        total += f(study, (scale=point.scale + 0.01i, threshold=point.threshold))
    end
    return total
end

reduction_allocations(f::F, study, point) where {F} = @allocated warm_reduction(f, study, point)

output, mode = ARGS[1:2]
width = parse(Int, ARGS[3])
width > 0 || error("length must be positive")
reducer = if mode == "sum"
    sum_rhs
elseif mode == "loop"
    loop_rhs
else
    error("use sum or loop")
end
study = (weights=collect(range(0.0, 1.0; length=width)),)
point = (scale=1.5, threshold=0.2)
# The timed invocation crosses a call boundary so compilation of the helper
# occurs after the timer starts, rather than while compiling its caller.
cold = @timed Base.invokelatest(reducer, study, point)
reference = sum(point.scale * max(w - point.threshold, 0.0) for w in study.weights)
@assert isapprox(cold.value, reference; atol=1e-12, rtol=1e-12)
Base.invokelatest(warm_reduction, reducer, study, point)
warm = @timed Base.invokelatest(warm_reduction, reducer, study, point)
bytes = Base.invokelatest(reduction_allocations, reducer, study, point)
@assert bytes == 0
report = Dict(
    "julia" => string(VERSION),
    "cpu" => Sys.cpu_info()[1].model,
    "mode" => mode,
    "length" => width,
    "cold_seconds" => cold.time,
    "cold_julia_bytes" => cold.bytes,
    "warm_seconds_per_call" => warm.time / 20,
    "warm_kernel_bytes" => bytes,
    "value" => cold.value,
    "scope" => "first helper call through invokelatest; 20-call warmed loop; fresh process",
)
open(output, "w") do io
    return TOML.print(io, report; sorted=true)
end
println(mode, " length=", width, " cold=", cold.time, " warm Julia bytes=", bytes)
