using DiffEqPDEBase, FiniteElementDiffEq, DiffEqProblemLibrary
using Base.Test

const FILEIO_ENABLE = false

println("Quadrature Points Tests")
@time @test include("quadpts_test.jl")
println("Boundary Tests")
FILEIO_ENABLE && @time @test include("boundary_tests.jl")
println("Example Mesh Tests")
FILEIO_ENABLE && @time @test include("mesh_examples_tests.jl")
println("Simple Mesh Tests")
@time @test include("mesh_SimpleMesh_tests.jl")
