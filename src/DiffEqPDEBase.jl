__precompile__()

module DiffEqPDEBase

using DiffEqBase

using RecipesBase, VectorizedRoutines.Matlab, ChunkedArrays

using Compat

@compat abstract type PDEProblem <: DEProblem end
@compat abstract type AbstractPoissonProblem{islinear,isstochastic,MeshType} <: PDEProblem end
@compat abstract type AbstractHeatProblem{islinear,isstochastic,MeshType} <: PDEProblem end

@compat abstract type AbstractFEMSolution{T,N} <: AbstractTimeseriesSolution{T,N} end

@compat abstract type AbstractMesh end
@compat abstract type AbstractFEMMesh <: AbstractMesh end

@compat abstract type PDEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractFEMAlgorithm <: PDEAlgorithm end
@compat abstract type AbstractFDMAlgorithm <: PDEAlgorithm end

# Specific Problem Algorithms
@compat abstract type AbstractHeatFEMAlgorithm <: AbstractFEMAlgorithm end
@compat abstract type AbstractPoissonFEMAlgorithm <: AbstractFEMAlgorithm end

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
