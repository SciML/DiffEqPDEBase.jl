__precompile__()

module DiffEqPDEBase

using DiffEqBase

using RecipesBase, VectorizedRoutines.Matlab, ChunkedArrays

abstract PDEProblem <: DEProblem
abstract AbstractPoissonProblem <: PDEProblem
abstract AbstractHeatProblem <: PDEProblem

abstract AbstractFEMSolution <: DESolution

"`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
abstract AbstractMesh
abstract AbstractFEMMesh <: AbstractMesh

abstract PDEAlgorithm <: DEAlgorithm
abstract AbstractFEMAlgorithm <: PDEAlgorithm
abstract AbstractFDMAlgorithm <: PDEAlgorithm

# Specific Problem Algorithms
abstract AbstractHeatFEMAlgorithm <: AbstractFEMAlgorithm
abstract AbstractPoissonFEMAlgorithm <: AbstractFEMAlgorithm

import Base: size, length, start, next, done, eltype

include("mesh_tools/meshes.jl")
include("utils.jl")
include("problems.jl")
include("solutions.jl")
include("plotrecipes.jl")
include("fem_assembly.jl")
include("fem_error.jl")
include("mesh_tools/fem_boundary.jl")

# Types
export HeatProblem, PoissonProblem, FEMSolution, FEMMesh, SimpleFEMMesh, FDMMesh

export AbstractPoissonProblem, AbstractHeatProblem, AbstractFEMSolution, PDEProblem

export PDEAlgorithm, AbstractFEMAlgorithm, AbstractFDMAlgorithm

export AbstractHeatFEMAlgorithm, AbstractPoissonFEMAlgorithm

#FEM Functions
export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
        gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
        meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

#Noise Functions

export getNoise

end # module
