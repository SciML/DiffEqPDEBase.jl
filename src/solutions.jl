type FEMSolution{T,N,uType,UTrueType,EType,tType,ProbType} <: AbstractFEMSolution{T,N}
  u::uType
  u_analytic::UTrueType
  errors::EType
  t::tType
  prob::ProbType
  tslocation::Int # Used in the iterator for animation
end

function FEMSolution(u,u_analytic,Du,t,prob)
  errors = Dict(:L2=>getL2error(prob.mesh,prob.analytic,u[end]),:H1=>getH1error(prob.mesh,Du,u[end]),
                :l∞=> maximum(abs.(u[end]-u_analytic)), :l2=> norm(u[end]-u_analytic,2))
  T = eltype(eltype(u))
  N = length((size(u)..., length(u)))
  FEMSolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(prob)}(u,u_analytic,errors,t,prob,0)
end
function FEMSolution(u,u_analytic,Du,prob)
  FEMSolution(u,u_analytic,Du,[],prob)
end
function FEMSolution(u::AbstractArray,prob)
  FEMSolution(u,[],prob)
end
function FEMSolution(u::AbstractArray,t,prob)
  T = eltype(eltype(u))
  N = length((size(u)..., length(u)))
  errors = Dict{Symbol,Float64}()
  FEMSolution{T,N,typeof(u),Void,typeof(errors),typeof(t),typeof(prob)}(u,nothing,errors,t,prob,0)
end

function FEMSolution(sol::AbstractFEMSolution,u_analytic,errors)
  T = eltype(eltype(sol.u))
  N = length((size(sol.u)..., length(sol.u)))
  FEMSolution{T,N,typeof(sol.u),typeof(u_analytic),typeof(errors),typeof(sol.t),typeof(sol.prob)}(sol.u,u_analytic,errors,sol.t,sol.prob,sol.tslocation)
end
function Base.show(io::IO, A::AbstractFEMSolution)
  print(io,"t: ")
  show(io, A.t)
  println(io)
  print(io,"u: ")
  show(io, A.u)
end
function Base.show(io::IO, m::MIME"text/plain", A::AbstractFEMSolution)
  print(io,"t: ")
  show(io, A.t)
  println(io)
  print(io,"u: ")
  show(io, A.u)
end
