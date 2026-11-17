module BackendAssemblyTests

using Test
using SparseArrays
import ReactiveModel as RM

variable(domain, offset, lower, upper) =
    RM._BackendVariableBlock(domain, offset, lower, upper, zeros(length(lower)))
row(indices, values, rhs, sense) = RM._BackendConstraint(indices, values, rhs, sense)
rows(offset, constraints) =
    RM._BackendConstraintBlock(offset, constraints, zeros(length(constraints)))
objective(indices=Int[], values=Float64[]; sense=:min) =
    RM._BackendObjective(RM._BackendObjectiveCoefficients(indices, values), sense)

@testset "Shared backend assembly respects mathematical offsets and sense" begin
    # Deliberately traverse blocks out of column/row order: offsets define the LP.
    variables = [
        variable(:binary, 3, [0.0], [1.0]),
        variable(:continuous, 0, [-Inf, 0.0], [Inf, 9.0]),
        variable(:integer, 2, [-4.0], [6.0]),
    ]
    constraints = [
        rows(2, [row([1, 4], [-1.0, 2.0], -2.0, :ge)]),
        rows(0, [row([1, 3], [2.0, 3.0], 5.0, :eq),
                 row([2, 4], [4.0, -5.0], 8.0, :le)]),
    ]
    obj = objective([4, 1, 3], [7.0, -2.0, 5.0]; sense=:max)
    data = RM._backend_data(variables, constraints, obj)
    @test data.ncols == 4
    @test data.nrows == 3
    @test data.cost == [-2.0, 0.0, 5.0, 7.0]
    @test data.lower == [-Inf, 0.0, -4.0, 0.0]
    @test data.upper == [Inf, 9.0, 6.0, 1.0]
    @test data.row_lower == [5.0, -Inf, -2.0]
    @test data.row_upper == [5.0, 8.0, Inf]
    @test data.is_mip
    @test data.matrix isa SparseMatrixCSC{Float64,Int}
    @test data.matrix.colptr == [1, 3, 4, 5, 7]
    @test data.matrix.rowval == [1, 3, 2, 1, 2, 3]
    @test data.matrix.nzval == [2.0, -1.0, 4.0, 3.0, -5.0, 2.0]
    @test data.matrix * [2.0, 3.0, 4.0, 5.0] == [16.0, -13.0, 8.0]
    @test sum(data.cost .* [2.0, 3.0, 4.0, 5.0]) == 51.0
    @test isempty(data.tiny)

    native = RM._backend_arrays(variables, constraints, obj)
    @test native.integrality == Cint[0, 0, 1, 1]
    @test native.starts isa Vector{Cint}
    @test native.indices isa Vector{Cint}
    @test native.values isa Vector{Float64}
    @test native.starts == Cint[0, 2, 3, 4, 6]
    @test native.indices == Cint[0, 2, 1, 0, 1, 2]
    @test native.values == [2.0, -1.0, 4.0, 3.0, -5.0, 2.0]

    # Shared costs retain the mathematical sign; each solver handles Max locally.
    min_data = RM._backend_data(variables, constraints,
        objective([4, 1, 3], [7.0, -2.0, 5.0]; sense=:min))
    @test min_data.cost == [-2.0, 0.0, 5.0, 7.0]
    @test obj.coefficients.variable_coefficients == [7.0, -2.0, 5.0]
    data.lower[1] = -10.0
    data.cost[1] = 99.0
    @test variables[2].lower_bounds == [-Inf, 0.0]
    @test min_data.lower == [-Inf, 0.0, -4.0, 0.0]
    @test min_data.cost == [-2.0, 0.0, 5.0, 7.0]

    for T in (Cchar, UInt8)
        @test RM._backend_column_types(variables, T) == T.(['C', 'C', 'I', 'B'])
        senses, rhs = RM._backend_row_data(constraints, T)
        @test senses == T.(['E', 'L', 'G'])
        @test rhs == [5.0, 8.0, -2.0]
    end
end

@testset "Backend coefficient filtering retains each solver's threshold contract" begin
    coefficients = [0.0, -0.0, 1e-12, -1e-12, nextfloat(1e-12),
                    -nextfloat(1e-12), 1e-9, -1e-9, nextfloat(1e-9),
                    -nextfloat(1e-9), 0.5, -0.5]
    variables = [variable(:continuous, 0, zeros(12), ones(12))]
    constraints = [rows(0, [row(collect(1:12), copy(coefficients), 0.0, :eq)])]

    unchanged = RM._backend_data(variables, constraints, objective())
    @test !unchanged.is_mip
    @test nnz(unchanged.matrix) == 12 # Includes the two explicit zero entries.
    @test unchanged.matrix.nzval == coefficients
    @test unchanged.matrix.colptr == collect(1:13)

    # ClpJL and commercial imports exclude equality at 1e-12, including zero.
    small = RM._backend_data(variables, constraints, objective(); drop_tolerance=1e-12)
    @test nnz(small.matrix) == 8
    @test small.matrix.nzval == coefficients[5:12]
    @test small.matrix.colptr == [1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    @test findall(!iszero, vec(Matrix(small.matrix))) == collect(5:12)
    native_small = RM._backend_arrays(variables, constraints, objective(); drop_tolerance=1e-12)
    @test native_small.starts == Cint[0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8]
    @test native_small.indices == zeros(Cint, 8)
    @test native_small.values == coefficients[5:12]

    # TinyHiGHS's default import also excludes equality at 1e-9.
    tiny = RM._backend_data(variables, constraints, objective(); drop_tolerance=1e-9)
    @test nnz(tiny.matrix) == 4
    @test tiny.matrix.nzval == coefficients[9:12]
    @test tiny.matrix.colptr == [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5]
    @test findall(!iszero, vec(Matrix(tiny.matrix))) == collect(9:12)

    # Native HiGHS protects strictly below 1e-9; equality is imported directly.
    protected = RM._backend_data(variables, constraints, objective(); protect_tiny=true)
    @test nnz(protected.matrix) == 12
    @test protected.matrix.nzval == [ones(6); coefficients[7:12]]
    @test protected.tiny == [(1, i, coefficients[i]) for i in 1:6]
    @test constraints[1].constraints[1].variable_coefficients == coefficients
    @test unchanged.matrix.nzval == coefficients

    # Filtering takes place per input entry before sparse duplicate aggregation.
    duplicates = [rows(0, [row([1, 1, 2, 2], [0.75e-12, 0.75e-12, 2e-12, -2e-12],
                              0.0, :eq)])]
    combined = RM._backend_data(variables, duplicates, objective(); drop_tolerance=1e-12)
    @test combined.matrix[1, 1] == 0.0
    @test combined.matrix[1, 2] == 0.0
    @test nnz(combined.matrix) == 1 # Surviving opposite entries sum to an explicit zero.
end

@testset "Row senses come from mathematics, including large and infinite RHS" begin
    constraints = [rows(0, [
        row(Int[], Float64[], 1e25, :eq),
        row(Int[], Float64[], -1e25, :le),
        row(Int[], Float64[], 1e25, :ge),
        row(Int[], Float64[], Inf, :eq),
        row(Int[], Float64[], -Inf, :le),
        row(Int[], Float64[], Inf, :ge),
    ])]
    for T in (Cchar, UInt8)
        senses, rhs = RM._backend_row_data(constraints, T)
        @test senses == T.(['E', 'L', 'G', 'E', 'L', 'G'])
        @test rhs == [1e25, -1e25, 1e25, Inf, -Inf, Inf]
    end
end

@testset "Empty backend data retains CSC dimensions and ABI sentinel" begin
    variables = RM._BackendVariableBlock[]
    constraints = RM._BackendConstraintBlock[]
    empty = RM._backend_data(variables, constraints, objective())
    @test (empty.nrows, empty.ncols) == (0, 0)
    @test size(empty.matrix) == (0, 0)
    @test empty.matrix.colptr == [1]
    @test isempty(empty.matrix.rowval)
    @test isempty(empty.matrix.nzval)
    @test isempty(empty.cost)
    @test isempty(empty.lower)
    @test isempty(empty.upper)
    @test isempty(empty.row_lower)
    @test isempty(empty.row_upper)
    @test isempty(empty.tiny)
    @test !empty.is_mip
    native = RM._backend_arrays(variables, constraints, objective())
    @test isempty(native.integrality)
    @test native.starts == Cint[0]
    @test isempty(native.indices)
    @test isempty(native.values)
    for T in (Cchar, UInt8)
        @test isempty(RM._backend_column_types(variables, T))
        senses, rhs = RM._backend_row_data(constraints, T)
        @test isempty(senses)
        @test isempty(rhs)
    end

    columns_only = RM._backend_data(
        [variable(:continuous, 0, zeros(3), ones(3))], constraints, objective())
    @test size(columns_only.matrix) == (0, 3)
    @test columns_only.matrix.colptr == [1, 1, 1, 1]
    @test columns_only.cost == zeros(3)
    @test !columns_only.is_mip
    empty_integer = RM._backend_data(
        [variable(:integer, 0, Float64[], Float64[])], constraints, objective())
    @test !empty_integer.is_mip
    rows_only = RM._backend_data(variables,
        [rows(0, [row(Int[], Float64[], 2.0, :eq)])], objective())
    @test size(rows_only.matrix) == (1, 0)
    @test rows_only.matrix.colptr == [1]
    @test rows_only.row_lower == [2.0]
    @test rows_only.row_upper == [2.0]
end

end # module
