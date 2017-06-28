type HeatProblem{islinear,isstochastic,MeshType,F,F2,F3,F4,F5,F6,F7,DiffType} <: AbstractHeatProblem{islinear,isstochastic,MeshType}
  u0::F5
  Du::F2
  f::F
  gD::F3
  gN::F4
  analytic::F7
  numvars::Int
  σ::F6
  noisetype::Symbol
  D::DiffType
  mesh::MeshType
end

function HeatProblem(analytic,Du,f,mesh;gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
  islinear = numargs(f)==2
  u0 = analytic(0,mesh.node)
  numvars = size(u0,2)
  gD = analytic
  if gN == nothing
    gN=(t,x)->zeros(size(x,1),numvars)
  end
  if σ==nothing
    isstochastic=false
    σ=(t,x)->zeros(size(x,1),numvars)
  else
    isstochastic=true
  end
  if D == nothing
    if numvars == 1
      D = 1.0
    else
      D = ones(1,numvars)
    end
  end
  HeatProblem{islinear,isstochastic,typeof(mesh),typeof(f),typeof(Du),typeof(gD),typeof(gN),typeof(u0),typeof(σ),typeof(analytic),typeof(D)}(u0,Du,f,gD,gN,analytic,numvars,σ,noisetype,D,mesh)
end


function HeatProblem(u0,f,mesh;gD=nothing,gN=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
  if σ==nothing
    isstochastic=false
    σ=(t,x)->zeros(size(x,1))
  else
    isstochastic=true
  end
  islinear = numargs(f)==2
  if islinear
    numvars = 1
    if u0==nothing
      u0 = zeros(size(mesh.node,1),numvars)
    end
    if gD == nothing
      gD=(t,x)->zeros(size(x,1))
    end
    if gN == nothing
      gN=(t,x)->zeros(size(x,1))
    end
    if D == nothing
      D = 1.0
    end
  end
  if !islinear #nonlinear
    if u0==nothing && numvars == nothing
      warn("u0 and numvars must be given. numvars assumed 1.")
      numvars = 1
      u0 = zeros(size(mesh.node,1),numvars)
      if gD == nothing
        gD=(t,x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(t,x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        D = 1.0
      end
    elseif u0==nothing #numvars!=nothing
      u0 = zeros(size(mesh.node,1),numvars) #Default to zero
      if gD == nothing
        gD=(t,x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(t,x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        D = ones(1,numvars)
      end
    elseif numvars==nothing #If u0 is given but numvars is not, we're still okay.
      numvars = size(u0,2)
      if gD == nothing
        gD=(t,x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(t,x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        if numvars == 1
          D = 1.0
        else
          D = ones(1,numvars)
        end
      end
    end
  end
  HeatProblem{islinear,isstochastic,typeof(mesh),typeof(f),Void,typeof(gD),typeof(gN),typeof(u0),typeof(σ),Void,typeof(D)}(u0,nothing,f,gD,gN,nothing,numvars,σ,noisetype,D,mesh)
end

Base.summary(prob::HeatProblem) = string(DiffEqBase.parameterless_type(prob))
function Base.show(io::IO, A::HeatProblem)
  println(io,summary(A))
end
function Base.display(io::IO, A::HeatProblem)
  println(io,summary(A))
end

type PoissonProblem{islinear,isstochastic,MeshType,F1,F2,F3,F4,F5,F6,F7,DiffType} <: AbstractPoissonProblem{islinear,isstochastic,MeshType}
  f::F1
  analytic::F2
  Du::F3
  gD::F4
  gN::F5
  u0::F6
  numvars::Int
  σ::F7
  noisetype::Symbol
  D::DiffType
  mesh::MeshType
end

function PoissonProblem(f,analytic,Du,mesh;gN=nothing,σ=nothing,u0=nothing,noisetype=:White,numvars=nothing,D=nothing)
  gD = analytic
  numvars = size(analytic(mesh.node),2)
  islinear = numargs(f)==1
  if gN == nothing
    gN=(x)->zeros(size(x,1),numvars)
  end
  if u0==nothing
    u0 = zeros(size(mesh.node,1),numvars)
  end
  if D == nothing
    if numvars == 1
      D = 1.0
    else
      D = ones(1,numvars)
    end
  end
  if σ==nothing
    isstochastic=false
    σ=(x)->zeros(size(x,1),numvars)
  else
    isstochastic=true
  end
  PoissonProblem{islinear,isstochastic,typeof(mesh),typeof(f),typeof(analytic),typeof(Du),typeof(gD),typeof(gN),typeof(u0),typeof(σ),typeof(D)}(f,analytic,Du,analytic,gN,u0,numvars,σ,noisetype,D,mesh)
end
function PoissonProblem(f,mesh;gD=nothing,gN=nothing,u0=nothing,σ=nothing,noisetype=:White,numvars=nothing,D=nothing)
  if σ==nothing
    isstochastic=false
    σ=(x)->zeros(size(x,1))
  else
    isstochastic = true
  end
  islinear = numargs(f)==1
  if islinear && u0==nothing
    numvars = 1
    u0 = zeros(size(mesh.node,1),numvars)
    if gD == nothing
      gD=(x)->zeros(size(x,1))
    end
    if gN == nothing
      gN=(x)->zeros(size(x,1))
    end
    if D == nothing
      D = 1.0
    end
  end
  if !islinear #nonlinear
    if u0==nothing && numvars == nothing
      warn("u0 and numvars must be given. numvars assumed 1.")
      numvars = 1
      u0 = zeros(size(mesh.node,1),numvars)
      if gD == nothing
        gD=(x)->zeros(size(x,1))
      end
      if gN == nothing
        gN=(x)->zeros(size(x,1))
      end
      if D == nothing
        D = 1.0
      end
    elseif u0==nothing #numvars!=nothing
      u0 = zeros(size(mesh.node,1),numvars) #Default to zero
      if gD == nothing
        gD=(x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        D = ones(1,numvars)
      end
    elseif numvars==nothing #If u0 is given but numvars is not, we're still okay. Generate from size in function.
      numvars = size(u0,2)
      if gD == nothing
        gD=(x)->zeros(size(x,1),numvars)
      end
      if gN == nothing
        gN=(x)->zeros(size(x,1),numvars)
      end
      if D == nothing
        if numvars == 1
          D = 1.0
        else
          D = ones(1,numvars)
        end
      end
    end
  end


  PoissonProblem{islinear,isstochastic,typeof(mesh),typeof(f),Void,Void,typeof(gD),typeof(gN),typeof(u0),typeof(σ),typeof(D)}(f,nothing,nothing,gD,gN,u0,numvars,σ,noisetype,D,mesh)
end

Base.summary(prob::PoissonProblem) = string(DiffEqBase.parameterless_type(prob))
function Base.show(io::IO, A::PoissonProblem)
  println(io,summary(A))
end
function Base.display(io::IO, A::PoissonProblem)
  println(io,summary(A))
end
