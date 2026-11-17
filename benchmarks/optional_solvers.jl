# Optional solvers must be explicitly requested, even if visible on LOAD_PATH.
const BENCH_TINY = get(ENV, "RM_BENCH_TINYHIGHS", "0") == "1"
if BENCH_TINY
    Base.find_package("TinyHiGHS") === nothing && error(
        "RM_BENCH_TINYHIGHS=1 requires TinyHiGHS in a separate benchmark environment; see benchmarks/README.md")
    import TinyHiGHS
end
benchmark_paths(paths) = Tuple(p for p in paths if BENCH_TINY || !occursin("tiny", String(p)))
