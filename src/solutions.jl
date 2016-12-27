"""
`FEMSolution`

Holds the data for the solution to a finite element problem.

### Fields

* `fem_mesh::FEMmesh`: The finite element mesh the problem was solved on.
* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `appxtrue::Bool`: Boolean flag for if u_analytic was an approximation.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
  is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.

"""
type FEMSolution <: AbstractFEMSolution
  fem_mesh::FEMmesh
  u#::Array{Number}
  trueknown::Bool
  u_analytic#::AbstractArrayOrVoid
  errors#::Dict{String,Float64}
  appxtrue::Bool
  t#::AbstractArrayOrVoid
  prob::DEProblem
  save_timeseries::Bool
  tslocation::Int # Used in the iterator for animation
  function FEMSolution(fem_mesh::FEMmesh,u,u_analytic,Du,t,prob;save_timeseries=true)
    errors = Dict(:L2=>getL2error(fem_mesh,prob.analytic,u[end]),:H1=>getH1error(fem_mesh,Du,u[end]),
                  :l∞=> maximum(abs.(u[end]-u_analytic)), :l2=> norm(u[end]-u_analytic,2))
    return(new(fem_mesh,u,true,u_analytic,errors,false,t,prob,true,0))
  end
  function FEMSolution(fem_mesh,u,u_analytic,Du,prob)
    FEMSolution(fem_mesh::FEMmesh,u,u_analytic,Du,[],prob,save_timeseries=false)
  end
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,prob)
    return(FEMSolution(fem_mesh,u,[],prob,save_timeseries=false))
  end
  function FEMSolution(fem_mesh::FEMmesh,u::AbstractArray,t,prob;save_timeseries=true)
    return(new(fem_mesh,u,false,nothing,Dict{String,Float64}(),false,t,prob,save_timeseries,0))
  end
end
