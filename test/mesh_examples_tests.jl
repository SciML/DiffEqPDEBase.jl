###Example Meshes

using JLD, DiffEqProblemLibrary

pkgdir = Pkg.dir("DiffEqProblemLibrary")
meshes_location = "premade_meshes.jld"

mesh = load("$pkgdir/$meshes_location","bunny")

mesh = load("$pkgdir/$meshes_location","flowpastcylindermesh")

mesh = load("$pkgdir/$meshes_location","lakemesh")

mesh = load("$pkgdir/$meshes_location","Lshapemesh")
mesh = load("$pkgdir/$meshes_location","Lshapeunstructure")

mesh = load("$pkgdir/$meshes_location","oilpump")

mesh = load("$pkgdir/$meshes_location","wavymesh")

mesh = load("$pkgdir/$meshes_location","wavyperturbmesh")

TEST_PLOT && plot(mesh)

true
