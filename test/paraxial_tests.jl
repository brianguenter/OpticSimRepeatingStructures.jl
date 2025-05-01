using Test
using TestItems
using TestItemRunner

@testitem "Projection" begin
    using StaticArrays
    import OpticSim

    using Unitful, Unitful.DefaultSymbols
    focallength = 10.0
    lens = OpticSim.ParaxialLensRect(focallength, 100.0, 100.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
    display = OpticSimRepeatingStructures.Display(1000, 1000, 1.0μm, 1.0μm, OpticSim.Geometry.translation(0.0, 0.0, -focallength))
    lenslet = OpticSimRepeatingStructures.LensletAssembly(lens, OpticSim.Geometry.identitytransform(), display)
    displaypoint = SVector(0.0, 0.0, -8.0)
    pupilpoints = SMatrix{3,2}(10.0, 10.0, 10.0, -10.0, -10.0, 20.0)
    project(lenslet, displaypoint, pupilpoints)
end

@testitem "BeamEnergy" begin
    using StaticArrays
    import OpticSim
    using Unitful
    using Unitful.DefaultSymbols

    focallength = 10.0
    lens = OpticSim.ParaxialLensRect(focallength, 1.0, 1.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
    display = OpticSimRepeatingStructures.Display(1000, 1000, 1.0μm, 1.0μm, OpticSim.Geometry.translation(0.0, 0.0, -focallength))
    lenslet = OpticSimRepeatingStructures.LensletAssembly(lens, OpticSim.Geometry.identitytransform(), display)
    displaypoint = SVector(0.0, 0.0, -8.0)
    #pupil is placed so that only 1/4 of it (approximately) is illuminated by lens beam
    pupil = OpticSim.Rectangle(1.0, 1.0, SVector(0.0, 0.0, -1.0), SVector(2.0, 2.0, 40.0))
    energy, centroid = OpticSimRepeatingStructures.beamenergy(lenslet, displaypoint, OpticSim.Geometry.vertices3d(pupil))
    @test isapprox(1 / 16, energy, atol=1e-4)
    @test isapprox([0.75, 0.75, 0.0], centroid)
end

@testitem "Repeat" begin
    using Colors
    using StaticArrays
    using DataFrames

    function hex3RGB()
        clusterelements = SVector((0, 0), (-1, 0), (-1, 1))
        colors = [colorant"red", colorant"green", colorant"blue"]
        names = ["R", "G", "B"]
        eltlattice = HexBasis1()
        clusterbasis = LatticeBasis((-1, 2), (2, -1))
        lattice = LatticeCluster(clusterbasis, eltlattice, clusterelements)
        properties = DataFrames.DataFrame(Color=colors, Name=names)
        return ClusterWithProperties(lattice, properties)
    end

    #spherepoint tests
    @test isapprox(spherepoint(1, π / 2, 0.0), [0.0, 1.0, 0.0])
    @test isapprox(spherepoint(1, 0.0, π / 2), [1.0, 0.0, 0.0])
    @test isapprox(spherepoint(1, 0, 0.0), [0.0, 0.0, 1.0])
    @test isapprox(spherepoint(1, 0.0, π / 4), [sqrt(2) / 2, 0.0, sqrt(2) / 2])


    """ Create a LatticeCluster with three elements at (0,0),(-1,0),(-1,1) coordinates in the HexBasis1 lattice"""
    function hex3cluster()
        clusterelts = SVector((0, 0), (-1, 0), (-1, 1))
        eltlattice = HexBasis1()
        clusterbasis = LatticeBasis((-1, 2), (2, -1))
        return LatticeCluster(clusterbasis, eltlattice, clusterelts)
    end

    @test [-1 2; 2 -1] == basismatrix(clusterbasis(hex3RGB()))

    function basistest(a::OpticSimRepeatingStructures.AbstractLatticeCluster)
        return clusterbasis(a)
    end

    @test basistest(hex3cluster()) == basistest(hex3RGB())

    #LatticeCluster testset
    cluster = hex9()

    #generate many tile coordinates. Compute the cluster index and tile index in that cluster for each tile coordinate.
    for iter in 1:100
        (i, j) = rand.((1:1000, 1:1000))
        coords, tileindex = cluster_coordinates_from_tile_coordinates(cluster, i, j)
        reconstructed = tilecoordinates(cluster, coords..., tileindex)
        @test all((i, j) .== reconstructed)
    end

    function testassignment()
        #test assignment of eyebox numbers to RGB clusters
        rgb_cluster = hex12RGB()
        cluster_coords = map(x -> Tuple(x), eachcol(clustercoordinates(rgb_cluster, 0, 0))) #create the cluster coordinates corresponding to each of the tiles in the cluster
        eyeboxnumbers = (1, 1, 2, 1, 2, 2, 3, 3, 3, 4, 4, 4) #correct eyebox number assignment for the tiles in the cluster
        for (index, coord) in enumerate(cluster_coords)
            boxnum = eyeboxnumbers[index]
            num = OpticSimRepeatingStructures.eyebox_number(coord, rgb_cluster, 4)
            if num != boxnum
                return false
            end
        end
        return true
    end

    @test testassignment()

    #verify that the 0,0 cluster is correct
    for (index, element) in pairs(clusterelements(cluster))
        coords, tileindex = cluster_coordinates_from_tile_coordinates(cluster, element...)
        @test all(coords .== 0)
        @test tileindex == index
    end

    #verify that choosecluster assertion doesn't fire incorrectly
    using Unitful, Unitful.DefaultSymbols
    function generate_clusters()
        freq = Vector{Int64}(undef, 0)
        subdivs = Vector{Tuple{Int64,Int64}}(undef, 0)
        areas = Vector(undef, 0)

        try
            for cycles in 15:30
                system_properties(15mm, (10mm, 9mm), (100°, 70°), 3.5mm, 0.1, cycles)
            end
        catch err #if any errors then failure
            return false
        end
        return true
    end

    @test generate_clusters() #shouldn't get assertion failures for any of the frequencies between 15:30
end