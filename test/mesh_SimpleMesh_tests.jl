using DiffEqPDEBase

pdeprob = prob_poisson_wave

res = solve(pdeprob)

mesh = SimpleFEMMesh(pdeprob.mesh.node,pdeprob.mesh.elem)

true
