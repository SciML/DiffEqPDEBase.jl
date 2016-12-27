__precompile__()

module DiffEqPDEBase

using DiffEqBase

using RecipesBase, VectorizedRoutines.Matlab, ChunkedArrays


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
export HeatProblem, PoissonProblem, FEMSolution, FEMmesh, SimpleMesh, FDMMesh

export animate

#FEM Functions
export  assemblematrix, findboundary, setboundary, findbdtype, getL2error, quadpts, getH1error,
        gradu, gradbasis, quadfbasis, fem_squaremesh, CFLμ, CFLν,
        meshgrid, notime_squaremesh, parabolic_squaremesh, quadpts1

#Noise Functions

export getNoise

end # module
