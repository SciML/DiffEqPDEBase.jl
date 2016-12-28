type FEMMesh{T1,T2,tType,tspanType} <: AbstractFEMMesh
  node::T1
  elem::Array{Int,2}
  bdnode::Vector{Int}
  freenode::Vector{Int}
  bdedge::Array{Int,2}
  is_bdnode::BitArray{1}
  is_bdelem::BitArray{1}
  bdflag::Array{Int8,2}
  totaledge::Array{Int,2}
  area::T2
  dirichlet::Array{Int,2}
  neumann::Array{Int,2}
  robin::Array{Int,2}
  N::Int
  NT::Int
  dt::tType
  tspan::tspanType
end

function FEMMesh(node,elem,dt,tspan,bdtype)
  N = size(node,1); NT = size(elem,1);
  totaledge = [elem[:,[2,3]]; elem[:,[3,1]]; elem[:,[1,2]]]

  #Compute the area of each element
  ve = Array{eltype(node)}(size(node[elem[:,3],:])...,3)
  ## Compute vedge, edge as a vector, and area of each element
  ve[:,:,1] = node[elem[:,3],:]-node[elem[:,2],:]
  ve[:,:,2] = node[elem[:,1],:]-node[elem[:,3],:]
  ve[:,:,3] = node[elem[:,2],:]-node[elem[:,1],:]
  area = 0.5*abs.(-ve[:,1,3].*ve[:,2,2]+ve[:,2,3].*ve[:,1,2])

  #Boundary Conditions
  bdnode,bdedge,is_bdnode,is_bdelem = findboundary(elem)
  bdflag = setboundary(node::AbstractArray,elem::AbstractArray,bdtype)
  dirichlet = totaledge[vec(bdflag .== 1),:]
  neumann = totaledge[vec(bdflag .== 2),:]
  robin = totaledge[vec(bdflag .== 3),:]
  is_bdnode = falses(N)
  is_bdnode[dirichlet] = true
  bdnode = find(is_bdnode)
  freenode = find(!is_bdnode)
  FEMMesh(node,elem,bdnode,freenode,bdedge,is_bdnode,is_bdelem,bdflag,totaledge,area,dirichlet,neumann,robin,N,NT,dt,tspan)
end
FEMMesh(node,elem,bdtype)=FEMMesh(node,elem,nothing,nothing,bdtype)

"""
`SimpleFEMMesh`

Holds the information describing a finite element mesh. For information on how (node,elem)
can be interpreted as a mesh describing a geometry, see [Programming of Finite
Element Methods by Long Chen](http://www.math.uci.edu/~chenlong/226/Ch3FEMCode.pdf).

### Fields

* `node`: The nodes in the (node,elem) structure.
* `elem`: The elements in the (node,elem) structure.
"""
type SimpleFEMMesh{T} <: AbstractFEMMesh
  node::T
  elem::Array{Int,2}
end


"""
`CFLμ(dt,dx)``

Computes the CFL-condition ``μ= dt/(dx*dx)``
"""
CFLμ(dt,dx)=dt/(dx*dx)

"""
`CFLν(dt,dx)``

Computes the CFL-condition ``ν= dt/dx``
"""
CFLν(dt,dx)=dt/dx

"""
`fem_squaremesh(square,h)`

Returns the grid in the iFEM form of the two arrays (node,elem)
"""
function fem_squaremesh(square,h)
  x0 = square[1]; x1= square[2];
  y0 = square[3]; y1= square[4];
  x,y = meshgrid(x0:h:x1,y0:h:y1)
  node = [x[:] y[:]];

  ni = size(x,1); # number of rows
  N = size(node,1);
  t2nidxMap = 1:N-ni;
  topNode = ni:ni:N-ni;
  t2nidxMap = deleteat!(collect(t2nidxMap),collect(topNode));
  k = t2nidxMap;
  elem = [k+ni k+ni+1 k ; k+1 k k+ni+1];
  return(node,elem)
end

"""
`notime_squaremesh(square,dx,bdtype)`

Computes the (node,elem) square mesh for the square
with the chosen `dx` and boundary settings.

###Example

```julia
square=[0 1 0 1] #Unit Square
dx=.25
notime_squaremesh(square,dx,"dirichlet")
```
"""
function notime_squaremesh(square,dx,bdtype)
  node,elem = fem_squaremesh(square,dx)
  return(FEMMesh(node,elem,bdtype))
end

"""
`parabolic_squaremesh(square,dx,dt,T,bdtype)`

Computes the `(node,elem) x [0,T]` parabolic square mesh
for the square with the chosen `dx` and boundary settings
and with the constant time intervals `dt`.

###Example

```julia
square=[0 1 0 1] #Unit Square
dx=.25; dt=.25;T=2
parabolic_squaremesh(square,dx,dt,T,:dirichlet)
```

"""
function parabolic_squaremesh(square,dx,dt,tspan,bdtype)
  node,elem = fem_squaremesh(square,dx)
  return(FEMMesh(node,elem,dt,tspan,bdtype))
end
