# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module OpticSimRepeatingStructures

export basismatrix

using Base: offset_if_vec
using StaticArrays: SVector, MVector, SMatrix, MMatrix
using DataFrames: DataFrame
import LazySets
using LinearAlgebra: norm
import OpticSim #only LensletAssembly uses OpticSim. This doesn't seem like a great idea. Probably should move LensletAssembly somewhere else or at least remove the dependency.
import OpticSim: surfaceintersection
using OpticSim: virtualpoint, SphericalPolygon, processintersection, point, ParaxialLens
import Unitful

include("Lattice.jl")
include("HexagonalLattice.jl")
include("RectangularLattice.jl")
include("Array.jl")
include("Cluster.jl")

include("HexClusters.jl")
include("HexTilings.jl")
include("Analysis.jl")
include("DisplayGeneration.jl")
include("LensletAssignment.jl")
include("Example.jl")

include("LensletAssembly.jl")


end #module
export OpticSimRepeatingStructures