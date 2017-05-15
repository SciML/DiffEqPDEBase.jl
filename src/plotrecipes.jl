@recipe function f(sol::FEMSolution;plot_analytic=false)
  if sol.tslocation==0 #Plot solution at end
    out = Any[]
    for i = 1:size(sol[end],2)
      push!(out,sol.u[end][:,i])
    end
    if plot_analytic
      for i = 1:size(sol.u,2)
        push!(out,sol.u_analytic[:,i])
      end
    end
  else #use timeseries
    out = Any[]
    for i = 1:sol.prob.numvars
      push!(out,sol.u[sol.tslocation][:,i])
    end
  end
  seriestype --> :surface
  layout --> length(out)
  sol.prob.mesh.node[:,1], sol.prob.mesh.node[:,2], out
end


@recipe function f(mesh::AbstractFEMMesh)
  seriestype --> :surface #:plot_trisurf
  #triangles  --> mesh.elem-1
  mesh.node[:,1], mesh.node[:,2], ones(mesh.node[:,1])
end

# mesh = meshExample_lakemesh()
# PyPlot.plot_trisurf(mesh.node[:,1],mesh.node[:,2],ones(mesh.node[:,2]),triangles=mesh.elem-1)
