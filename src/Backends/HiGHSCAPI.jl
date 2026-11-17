# Private bindings for the subset of the HiGHS C API used by ReactiveModel.
# ABI reference: HiGHS v1.15.1, highs/interfaces/highs_c_api.h.
# No HiGHS.jl or MathOptInterface dependency is needed to drive the solver.
module _HiGHSCAPI

import HiGHS_jll: libhighs
import LinearAlgebra
import OpenBLAS32_jll

# The native solver uses LP64 BLAS through libblastrampoline. HiGHS.jl used
# to initialize this for us; the direct backend must do so independently.
function __init__()
    for library in LinearAlgebra.BLAS.lbt_get_config().loaded_libs
        library.interface === :lp64 && return nothing
    end
    LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
    nothing
end

# HiGHS_jll uses 32-bit HighsInt. The backend checks the loaded library's ABI
# before passing any index buffers to it.
const kHighsStatusError = Cint(-1)
const kHighsStatusOk = Cint(0)
const kHighsVarTypeContinuous = Cint(0)
const kHighsVarTypeInteger = Cint(1)
const kHighsObjSenseMinimize = Cint(1)
const kHighsObjSenseMaximize = Cint(-1)
const kHighsMatrixFormatColwise = Cint(1)
const kHighsModelStatusOptimal = Cint(7)
const kHighsModelStatusInfeasible = Cint(8)
const kHighsModelStatusUnknown = Cint(15)

Highs_create() = ccall((:Highs_create, libhighs), Ptr{Cvoid}, ())
Highs_destroy(highs) = ccall((:Highs_destroy, libhighs), Cvoid, (Ptr{Cvoid},), highs)
Highs_getSizeofHighsInt(highs) =
    ccall((:Highs_getSizeofHighsInt, libhighs), Cint, (Ptr{Cvoid},), highs)
Highs_clearModel(highs) = ccall((:Highs_clearModel, libhighs), Cint, (Ptr{Cvoid},), highs)
Highs_clearSolver(highs) = ccall((:Highs_clearSolver, libhighs), Cint, (Ptr{Cvoid},), highs)
Highs_run(highs) = ccall((:Highs_run, libhighs), Cint, (Ptr{Cvoid},), highs)
Highs_getModelStatus(highs) =
    ccall((:Highs_getModelStatus, libhighs), Cint, (Ptr{Cvoid},), highs)
Highs_getObjectiveValue(highs) =
    ccall((:Highs_getObjectiveValue, libhighs), Cdouble, (Ptr{Cvoid},), highs)

Highs_setStringOptionValue(highs, option, value) =
    ccall((:Highs_setStringOptionValue, libhighs), Cint,
          (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)

function Highs_passLp(highs, ncol, nrow, nnz, format, sense, offset,
                      cost, lower, upper, row_lower, row_upper, start, index, value)
    ccall((:Highs_passLp, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cdouble,
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
          highs, ncol, nrow, nnz, format, sense, offset,
          cost, lower, upper, row_lower, row_upper, start, index, value)
end

function Highs_passMip(highs, ncol, nrow, nnz, format, sense, offset,
                       cost, lower, upper, row_lower, row_upper, start, index, value,
                       integrality)
    ccall((:Highs_passMip, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Cint, Cdouble,
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}),
          highs, ncol, nrow, nnz, format, sense, offset,
          cost, lower, upper, row_lower, row_upper, start, index, value, integrality)
end

Highs_changeColBounds(highs, col, lower, upper) =
    ccall((:Highs_changeColBounds, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cdouble, Cdouble), highs, col, lower, upper)
Highs_changeColsBoundsBySet(highs, n, indices, lower, upper) =
    ccall((:Highs_changeColsBoundsBySet, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}),
          highs, n, indices, lower, upper)
Highs_changeColsBoundsByRange(highs, first, last, lower, upper) =
    ccall((:Highs_changeColsBoundsByRange, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          highs, first, last, lower, upper)
Highs_changeRowBounds(highs, row, lower, upper) =
    ccall((:Highs_changeRowBounds, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cdouble, Cdouble), highs, row, lower, upper)
Highs_changeRowsBoundsByRange(highs, first, last, lower, upper) =
    ccall((:Highs_changeRowsBoundsByRange, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          highs, first, last, lower, upper)
Highs_changeCoeff(highs, row, col, value) =
    ccall((:Highs_changeCoeff, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Cint, Cdouble), highs, row, col, value)
Highs_changeColsCostBySet(highs, n, indices, cost) =
    ccall((:Highs_changeColsCostBySet, libhighs), Cint,
          (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}), highs, n, indices, cost)
Highs_changeObjectiveSense(highs, sense) =
    ccall((:Highs_changeObjectiveSense, libhighs), Cint,
          (Ptr{Cvoid}, Cint), highs, sense)
Highs_getSolution(highs, col_value, col_dual, row_value, row_dual) =
    ccall((:Highs_getSolution, libhighs), Cint,
          (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          highs, col_value, col_dual, row_value, row_dual)
Highs_writeModel(highs, filename) =
    ccall((:Highs_writeModel, libhighs), Cint,
          (Ptr{Cvoid}, Cstring), highs, filename)

end # module _HiGHSCAPI
