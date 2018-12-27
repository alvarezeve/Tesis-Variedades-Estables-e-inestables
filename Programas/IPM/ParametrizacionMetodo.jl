module ParametrizacionMetodo


using TaylorSeries
using LinearAlgebra
using GenericLinearAlgebra
using PyPlot
using IntervalArithmetic, IntervalRootFinding
using ValidatedNumerics
using StaticArrays

export Inestable
export InestableCoef
export Estable
export EstableCoef
export Convergencia
export Convergencia3
export CalculoError
export Evaluar
export Error4

include("VariedadesIPM.jl")
include("mapeos.jl")
#include("graficar-mapeos.jl")

end
