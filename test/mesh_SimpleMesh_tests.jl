using DiffEqPDEBase

pdeprob = prob_poisson_wave
mesh = SimpleFEMMesh(pdeprob.mesh.node,pdeprob.mesh.elem)

true
