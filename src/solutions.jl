type FEMSolution{uType,tType,UTrueType,uElType,ProbType} <: AbstractFEMSolution
  u::uType
  u_analytic::UTrueType
  errors::Dict{Symbol,uElType}
  t::tType
  prob::ProbType
  tslocation::Int # Used in the iterator for animation
end

function FEMSolution(u,u_analytic,Du,t,prob)
  errors = Dict(:L2=>getL2error(prob.mesh,prob.analytic,u[end]),:H1=>getH1error(prob.mesh,Du,u[end]),
                :l∞=> maximum(abs.(u[end]-u_analytic)), :l2=> norm(u[end]-u_analytic,2))
  FEMSolution(u,u_analytic,errors,t,prob,0)
end
function FEMSolution(u,u_analytic,Du,prob)
  FEMSolution(u,u_analytic,Du,[],prob)
end
function FEMSolution(u::AbstractArray,prob)
  FEMSolution(u,[],prob)
end
function FEMSolution(u::AbstractArray,t,prob)
  FEMSolution(u,nothing,Dict{Symbol,Float64}(),t,prob,0)
end

function FEMSolution(sol::AbstractFEMSolution,u_analytic,errors)
  FEMSolution(sol.u,u_analytic,errors,sol.t,sol.prob,sol.tslocation)
end
