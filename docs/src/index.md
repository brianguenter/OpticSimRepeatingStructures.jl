# Introduction

The Repeat module contains functions for creating regular repeated patterns. This could be pixels in a display grid, or mirrors in an active optics telescope. Repeated patterns are defined by creating an object that inherits from the abstract type [`OpticSimRepeatingStructures.AbstractBasis`](@ref).

Subtypes supporting the AbstractBasis interface should implement these functions:

Returns the neighbors in ring n surrounding centerpoint, excluding centerpoint
```
neighbors(::Type{B},centerpoint::Tuple{T,T},neighborhoodsize::Int) where{T<:Real,B<:AbstractBasis}
```
Returns the lattice basis vectors that define the lattice
```
basismatrix(a::S) where{S<:AbstractBasis}
```
Returns the vertices of the unit polygon for the basis that tiles the plane 
```
tilevertices(a::S) where{S<:AbstractBasis}
```

A lattice is described by a set of lattice vectors eᵢ which are stored in a [`OpticSimRepeatingStructures.AbstractBasis`](@ref) object. You can create bases in any dimension. Points in the lattice are indexed by integer coordinates. These lattice coordinates can be converted to Cartesian coordinates by indexing the LatticeBasis object. 
``` @example example
using OpticSimRepeatingStructures
a = LatticeBasis((1.0,5.0),(0.0,1.0))
a[3,3]
```

The Lattice points are defined by a weighted sum of the basis vectors:
```
latticepoint = ∑αᵢ*eᵢ
```
where the αᵢ are integer weights.

The [`HexBasis1`](@ref) constructor defines a symmetric basis for hexagonal lattices 
```@example 
using OpticSimRepeatingStructures
basismatrix(HexBasis1())
```
The [`rectangularlattice`](@ref) function creates a rectangular lattice basis. 

There are a few visualization functions for special 2D lattices. [`drawcells`](@ref) draws a set of hexagonal cells. Using [`hexcellsinbox`](@ref) we can draw all the hexagonal cells that fit in a rectangular box:

```@example 
using OpticSimRepeatingStructures
drawcells(HexBasis1(),50,hexcellsinbox(2,2))
```

There is also a function to compute the n rings of a cell x, i.e., the cells which can be reached by taking no more than n steps along the lattice from x:

```@example 
using OpticSimRepeatingStructures
drawcells(HexBasis1(),50,neighbors(HexBasis1,(0,0),2))
```
 
You can also draw all the cells contained within an n ring:
 
```@example 
using OpticSimRepeatingStructures
drawcells(HexBasis1(),50,region(HexBasis1,(0,0),2))
```

# Reference
```@autodocs
Modules = [OpticSimRepeatingStructures]
Order   = [:function, :type, :constant,:macro]
```
