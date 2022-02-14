# CURIS: MPAS-Regent Partitioning

# Background

In the Summer of 2021, Raphael, Sylvia and I were tasked with continuing progress on 3 dimensions of the MPAS-Regent project. While Sylvia bolstered testing automation efforts and Raphael worked to integrate CUDA into the code base, I extended the work Victor Lin and others had poured into the project related to partitioning.

To get up to speed on the techniques we were employing, feel free to check out the following with takeaways noted below:

- [Dependent Partitioning](http://theory.stanford.edu/~aiken/publications/papers/oopsla16.pdf)
    - Potential speed up of up to 30x on 64 nodes, if dependent data is correctly partition. Partitioning strategy defines private, shared, and ghost regions, where private includes data specific to a given partition only, shared includes data pertinent to a given partition and neighboring partitions, and ghost includes read-only data pertinent to a given partition but located in a neighboring partition.
    - Partitions are created by supplying a region to partition function and a coloring object: an abstract data type that describes association of colors (generally small integers), assign any number of colors to each element and must allow random access. In the case of MPAS-Regent, this is seeded by the original MPAS code and a call to the METIS partition function.
    - The framework here improves upon the partitioning capabilities of existing programming models by providing maximal expressivity for defining arbitrary independent partitions, but then provides a carefully chosen set of operations for the derivation of dependent partitions.
- [Language Support for Dynamic, Hierarchical Data Partitioning](http://theory.stanford.edu/~aiken/publications/papers/oopsla13a.pdf)
    - Key takeaways from this paper are that the data may be overlapping if we are to use dependent data partitioning; in other words, the data need not be disjoint and bordering subgraphs can be accessed from another subgraph in advancing an algorithm.
    - This expounds upon the type system and operational semantics, if you are coming from a CS 143 background.

The MPAS-Regent project is a complex piece of code. The dependent partitioning strategy originates from [the circuit example of private, shared, and ghost region divvying](https://github.com/StanfordLegion/legion/blob/stable/language/examples/circuit_sparse.rg#L727-L729). There are crucial differences between the two that must be considered moving forward. 

The first key component of its complexity is that the data is unstructured in the Voronoi mesh format. In this, we have a graph of arbitrary-degree nodes with implicit and explicit relationships to follow. This is crucial in distinguishing it from contiguous dependent partitioning, which can often be achieved through GPU/CUDA parallelization techniques.  

# Summer Research Takeaways

What follows is a distillation of the final conclusions made at the end of the summer. Many credits to Elliott Slaughter, who helped explain the particular reasons our region partitioning tree approach may not be compatible with the circuit example. First, here’s a reminder of the approach. 

## Approach & Strategy

To harness implicit parallelization on GPUs and CPUs, Regent requires **strictly disjoint data partitions**. Such partitions are straightforward when performed on arrays. However, **MPAS is an unstructured graph** composed of cells, edges, and vertices. See below.

![Voronoi Cell Mesh](https://github.com/alexaiken/mpas-regent/blob/master/voronoicellmesh.png)

Therefore, we require a preprocessing step to achieve disjoint partitions of cells, edges, and vertices while preserving the ability to calculate values across partitions.

Seeding partitioning with patches over the earth’s surface, we compute dependent
partition layers using the following strategy, given an original partition ***p**:*

- ***Shared**:* cells in ***p*** bordering another partition (in red at right)
- ***Private***: cells in ***p*** not in ***Shared*** (in gold at right)
- ***Ghost**:* cells bordering ***p*** (in gray at right)

![Dependent Partitioning Layers](https://github.com/alexaiken/mpas-regent/blob/master/dependentpartitionlayers.png)

Let ***f(x)*** be the set of all neighbor cells of a cell ***x***. Then the i**mage *f(p)*** of all the cells of a partition ***p*** includes the ghost cells of ***p***. Using **image** and preimage, image’s inverse, we can efficiently compute edge, vertex, and cell partitions to allow **parallel computation across partitions**. This is how we computed partitions (see below in Code).

Ultimately, we seek a full implementation of MPAS’ atmospheric framework in Regent to illustrate task-based parallelism with dependent data partitioning in a smaller codebase. Coupling this with GPUs, we showcase improvements in machine portability, modern hardware usage, and performance.

## Key Results

In the code, you’ll find that the strategy begins with a METIS partition. Through relationships of pointer attributes, we can follow this original partition to get subgraph partitions of edges, vertices, and cells. From this we found halo cells (i.e. the ghost layer), those neighboring the ghost layer in the original partition patch (i.e. shared), and those inside the original partition not in shared (i.e. private). This strategy takes a single unstructured graph region of the cell_region (cr) and splits it into three categories of cells (private, shared, ghost).

Per conversations with Elliott and Alex, we can know analytically (and soundly) that the partitions are disjoint using intersection set logic checks (zero overlap between). However, this alone is not enough to convince the Regent compiler that the partitions are, indeed, disjoint or won’t interfere during future calculations. Even more so, Regent offers the ability to use `*` in the where clause for tasks, but this will not help and incurs the error of `invalid call missing constraint`, which is Regent communicating its inability to mark these partitions as disjoint. In other words, we found the compiler is still confused about the ghost, shared, and private partitions and thinks they aren’t fully disjoint.

At the close of the summer’s investigation for whether we can implement dependent partitioning, Elliott recommended the crucial next step being **sketching a region tree relevant to the constraints we’re aiming to invoke**. If we can find another way to work around following pointers beyond the idea of calls to `image` and `preimage`, we will also be set up to preserve disjointness better as these calls cannot provide this guarantee.

Upon learning about this notion of a region tree, I passed along my most current code to Elliott, who noticed that `private_1`, `shared_1`, and `ghost_1` are all partitions from a single original cell partition. **This is why the region tree is broken from the beginning and cannot guarantee disjointness** because, from that point on, the compiler can no longer tell whether the partitions are disjoint. In contrast, the circuit example, showcased [here](https://github.com/StanfordLegion/legion/blob/stable/language/examples/circuit_sparse.rg#L727-L729), begins with an overall first step of divvying the original nodes into private and shared nodes from the very beginning without the use of `image`. From Elliott: the “issue comes from the structure of the region tree (three partitions on the same region instead of being on separate regions) rather than because one particular partition is not disjoint”:

- Circuit: [you'll see the first thing it does is split the top-level region into two parts, called all_private and all_shared](https://github.com/StanfordLegion/legion/blob/stable/language/examples/circuit_sparse.rg#L727-L729), Regent knows those regions are disjoint because (a) the partition they come from is disjoint, and (b) the indices used to access them are constants (0 and 1).

This may suggest we should do something similar, where we can use partition to achieve cell regions, at the least, which are fully disjoint for MPAS and agreeable to the compiler. As the summer came to a close, we attempted to get the code up and running to get incorrect results with our current partitioning strategy as highlighting potential for speedup is still an important outcome. To do so, I was able to run MPAS-Regent running incrementally with `-foverride-demand-index-launch 1`, though this is a temporary helper flag that should be removed later upon correcting the region tree upstream.

Potential Next Steps:

- Obtain modified source code (uncommitted to Git as not functioning fully and likely should be overhauled in some way(s)), and examine how partitioning has been implemented. See below for excerpts of most crucial changes.
- An idea: load back to disk and then load after creating correct partitions vs. main tradeoffs are complexity and overhead. Loading partitions from disk will get a lot more complicated if your meshes can't fit in memory in a single node. So I'd think very hard about whether that will be the case or not, since that'll be a major feature you'll need to architect for if you go that way.
- Potential for image calls to be annotated in Legion but not Regent as disjoint.

# Code from the Summer

## Sapling > mesh_loading > `mesh_loading.rg`

This code is an excerpt of only the partition task. It has comments for what is in progress and remains to be done. However, it should be noted that iterating on this code might not be the ideal next step as the partition tree may need to change from the initial step of loading the code.

```cpp
-----------------------------------------------
------- TASK: PARTITION REGIONS  --------
-----------------------------------------------

--input: cell_region, edge_region, vertex_region
--return: cell_partition_fs containing private, shared, and ghost partitions for one and two layers.
task partition_regions(num_partitions : int, 
                        cell_region : region(ispace(int2d), cell_fs), 
                        edge_region : region(ispace(int2d), edge_fs), 
                        vertex_region : region(ispace(int2d), 
                        vertex_fs))
where
    reads writes (cell_region, edge_region, vertex_region)
do
    var color_space = ispace(int1d, num_partitions)

    var p = partition(cell_region.partitionNumber, color_space) -- Original partition based on Metis

    format.println("p[1] volume: {}", p[1].volume)

    -- We compute all edges around p in the process of creating halo layers
    -- which assists in our edge partitioning scheme, as well.

    var e0 = image(edge_region, p, cell_region.edgesOnCell0)
    var e1 = image(edge_region, p, cell_region.edgesOnCell1)
    var e2 = image(edge_region, p, cell_region.edgesOnCell2)
    var e3 = image(edge_region, p, cell_region.edgesOnCell3)
    var e4 = image(edge_region, p, cell_region.edgesOnCell4)
    var e5 = image(edge_region, p, cell_region.edgesOnCell5)
    var e6 = image(edge_region, p, cell_region.edgesOnCell6)
    var e7 = image(edge_region, p, cell_region.edgesOnCell7)
    var e8 = image(edge_region, p, cell_region.edgesOnCell8)
    var e9 = image(edge_region, p, cell_region.edgesOnCell9)
    var e = e0 | e1 | e2 | e3 | e4 | e5 | e6 | e7 | e8 | e9

    -- Note e denotes all edges within parition p (i.e. all edges on a cell in p)

    format.println("e0[1].volume={}", e0[1].volume)
    format.println("e1[1].volume={}", e1[1].volume)
    format.println("e2[1].volume={}", e2[1].volume)
    format.println("e3[1].volume={}", e3[1].volume)
    format.println("e4[1].volume={}", e4[1].volume)
    format.println("e5[1].volume={}", e5[1].volume)
    format.println("e6[1].volume={}", e6[1].volume)
    format.println("e7[1].volume={}", e7[1].volume)
    format.println("e8[1].volume={}", e8[1].volume)
    format.println("e9[1].volume={}", e9[1].volume)
    format.println("e[1].volume={}\n", e[1].volume)

    --for k = 0, constants.NUM_PARTITIONS do
    --    format.println("Partition {}", k)
    --    for i in e0[k] do
    --        if i.y == 0 then
    --            format.println("{}, {}", k, i)
    --        end
    --    end
    --end

    for i = 0, constants.nEdges do
        edge_region[{i, 0}].cellOne = rect2d { int2d { edge_region[{i, 0}].cellsOnEdge[0], 0 }, int2d { edge_region[{i, 0}].cellsOnEdge[0], constants.nVertLevels - 1 } }
        edge_region[{i, 0}].cellTwo = rect2d { int2d { edge_region[{i, 0}].cellsOnEdge[1], 0 }, int2d { edge_region[{i, 0}].cellsOnEdge[1], constants.nVertLevels - 1 } }

        -- Now, similarly declare vertices one and two for each edge
        edge_region[{i, 0}].vertexOne = rect2d { int2d { edge_region[{i, 0}].verticesOnEdge[0], 0 }, int2d { edge_region[{i, 0}].verticesOnEdge[0], constants.nVertLevels - 1 } }
        edge_region[{i, 0}].vertexTwo = rect2d { int2d { edge_region[{i, 0}].verticesOnEdge[1], 0 }, int2d { edge_region[{i, 0}].verticesOnEdge[1], constants.nVertLevels - 1 } }
    end

    -- With vertexOne and vertexTwo established, grab the vertex partitions
    var v = image(vertex_region, e, edge_region.vertexOne) | image(vertex_region, e, edge_region.vertexTwo)
    -- Note: v denotes the partition of all vertices within e

    format.println("Upon assigning vertexOne and vertexTwo, vertex-partitioning:")
    format.println("\tv[1].volume={}", v[1].volume)

    --Note: the following code should theoretically be able to replace the 10 images, but seems to give the wrong volume
    var ep_one = preimage(edge_region, p, edge_region.cellOne) -- This should get all edges with cellOne in p, partitioned by cellOne (which cell "originated" it)
    var ep_two = preimage(edge_region, p, edge_region.cellTwo) -- This should get all edges with cellTwo in p, partitioned by cellTwo (which cell it "points to")

    format.println("An experiment to edge-partition with preimage")
    format.println("\tep_one[1].volume={}", ep_one[1].volume)
    format.println("\tep_two[1].volume={}\n", ep_two[1].volume)

    var cp_one = image(cell_region, e, edge_region.cellOne) -- This gets the cellOne on all edges in e
    var cp_two = image(cell_region, e, edge_region.cellTwo) -- This gets the cellTwo on all edges in e
    -- Combined, they should contain all cells who have an edge in e, and therefore the halo.

    var ghost_1_and_p = cp_one | cp_two
    format.println("(cp_one | cp_two)[1] volume: {}", ghost_1_and_p[1].volume)
    var ghost_1 = (cp_one | cp_two) - p

    -- Edge-partitioning from ghost akin to above
    var e0_g = image(edge_region, ghost_1, cell_region.edgesOnCell0)
    var e1_g = image(edge_region, ghost_1, cell_region.edgesOnCell1)
    var e2_g = image(edge_region, ghost_1, cell_region.edgesOnCell2)
    var e3_g = image(edge_region, ghost_1, cell_region.edgesOnCell3)
    var e4_g = image(edge_region, ghost_1, cell_region.edgesOnCell4)
    var e5_g = image(edge_region, ghost_1, cell_region.edgesOnCell5)
    var e6_g = image(edge_region, ghost_1, cell_region.edgesOnCell6)
    var e7_g = image(edge_region, ghost_1, cell_region.edgesOnCell7)
    var e8_g = image(edge_region, ghost_1, cell_region.edgesOnCell8)
    var e9_g = image(edge_region, ghost_1, cell_region.edgesOnCell9)
    var e_g1 = e0_g | e1_g | e2_g | e3_g | e4_g | e5_g | e6_g | e7_g | e8_g | e9_g

    -- Apply same extension to edge partition for vertices
    var v_g1 = image(vertex_region, e_g1, edge_region.vertexOne) | image(vertex_region, e_g1, edge_region.vertexTwo)

    format.println("raw edge_ghost1[1] volume (without -e): {}", e_g1[1].volume)
    format.println("raw vertex_ghost[1] volume (without -v): {}", v_g1[1].volume)

    var e_g1_isolate = e_g1 - e -- Isolate only the 1-layer ghost edge partition
    var v_g1_isolate = v_g1 - v -- Isolate only the 1-layer ghost vertex partition

    --format.println("(isolated) edge_ghost1[1] volume: {}", e_g1_isolate[1].volume)
    --format.println("(isolated) vertex_ghost1[1] volume: {}", v_g1_isolate[1].volume)

    -- Calculate second halo
    var gep_out = preimage(edge_region, ghost_1_and_p, edge_region.cellOne)
    var gep_in = preimage(edge_region, ghost_1_and_p, edge_region.cellTwo)
    var gcp_out = image(cell_region, gep_out, edge_region.cellTwo)
    var gcp_in = image(cell_region, gep_in, edge_region.cellOne)
    var ghost_2 = (gcp_in | gcp_out) - p -- First and second halo layers
    var gep_2 = (gep_out | gep_in) - e -- First and second halo layers 

    -- Compute all cells reachable from ghost_1. shared_1 is intersection of that set with p
    var s1ep_out = preimage(edge_region, ghost_1, edge_region.cellOne)
    var s1ep_in = preimage(edge_region, ghost_1, edge_region.cellTwo)
    var s1cp_out = image(cell_region, s1ep_out, edge_region.cellTwo)
    var s1cp_in = image(cell_region, s1ep_in, edge_region.cellOne)
    var shared_1 = p & (s1cp_out | s1cp_in) -- Cells in p bordering ghost_1
    var private_1 = p - shared_1 -- all cells in p that are not in shared_1

    -- For vertices, shared_1 is the intersection between vertices reachable by ghost and v
    var v_s1 = v & (image(vertex_region, e_g1_isolate, edge_region.vertexOne) | image(vertex_region, e_g1_isolate, edge_region.vertexTwo))

    -- Similarly, edges-shared is the intersection between the ghost edges and original edges
    var e_s1 = e & (s1ep_out | s1ep_in) -- Edges in e bordering ghost_1
    --format.println("edge_shared1[1] volume: {}", e_s1[1].volume)
    --format.println("vertex_shared1[1] volume: {}", v_s1[1].volume)

    -- Finally, private is everything not shared from e
    var e_p1 = (e - e_s1) - e_g1 -- Edges in e not shared or in ghost
    var v_p1 = (v - v_s1) - v_g1 -- Vertices in v not in shared or ghost
    --format.println("edge_private1[1] volume: {}\n", e_p1[1].volume)
    --format.println("vertex_private1[1] volume: {}\n", v_p1[1].volume)

    -- TODO: Extend shared/ghost/private to 2-layers for vertices/eges
    -- shared_2 contains shared_1 and all cells in p bordering shared_1
    var s2ep_out = preimage(edge_region, shared_1, edge_region.cellOne)
    var s2ep_in = preimage(edge_region, shared_1, edge_region.cellTwo)
    var s2cp_out = image(cell_region, s2ep_out, edge_region.cellTwo)
    var s2cp_in = image(cell_region, s2ep_in, edge_region.cellOne)
    var shared_2 = dynamic_cast(partition(disjoint, cell_region, color_space), (shared_1 | (private_1 & (s2cp_out | s2cp_in)))) -- Cells in p bordering ghost_1
    var private_2 = private_1 - shared_2 -- all cells in private_1 that are not in shared_2

    format.println("Cell Partition Sizes (P,S,G):")
    for i = 0, constants.NUM_PARTITIONS do
        format.println("\t({},{},{})", private_1[i].volume, shared_1[i].volume, ghost_1[i].volume)
    end

    format.println("Edge Partition Sizes (P,S,G):")
    for i = 0, constants.NUM_PARTITIONS do
        format.println("\t({},{},{})", e_p1[i].volume, e_s1[i].volume, e_g1_isolate[i].volume)
    end

    format.println("Vertex Partition Sizes (P,S,G):")
    for i = 0, constants.NUM_PARTITIONS do
        format.println("\t({},{},{})", v_p1[i].volume, v_s1[i].volume, v_g1_isolate[i].volume)
    end

    format.println("")
    format.println("\nDisjoint? Check intersection:")
    var interps = private_1 & shared_1
    var interps2 = private_1 & ghost_1
    var interps3 = ghost_1 & shared_1
    format.println("C-Private1[1] & C-Shared1[1]: {}", interps[1].volume)
    format.println("C-Private1[1] & C-Ghost1[1]: {}", interps2[1].volume)
    format.println("C-Ghost[1] & C-Shared1[1]: {}", interps3[1].volume)

    format.println("Cell-Partitioning Results:")
    format.println("\tC-Private1[1] volume:\t{}", private_1[1].volume)
    --format.println("\tC-Private2[1] volume:\t{}", private_2[1].volume)
    format.println("\tC-Shared1[1] volume:\t{}", shared_1[1].volume)
    --format.println("\tC-Shared2[1] volume:\t{}", shared_2[1].volume)
    format.println("\tC-Ghost1[1] volume:\t{}", ghost_1[1].volume)
    --format.println("\tC-Ghost2[1] volume:\t{}", ghost_2[1].volume)

    format.println("Edge-Partitioning Results:")
    format.println("\tE-Private1[1] volume:\t{}", e_p1[1].volume)
    format.println("\tE-Shared1[1] volume:\t{}", e_s1[1].volume)
    format.println("\tE-Ghost1[1] volume:\t{}", e_g1_isolate[1].volume)

    format.println("Vertex-Partitioning Results:")
    format.println("\tE-Private1[1] volume:\t{}", v_p1[1].volume)
    format.println("\tE-Shared1[1] volume:\t{}", v_s1[1].volume)
    format.println("\tE-Ghost1[1] volume:\t{}\n", v_g1_isolate[1].volume)

    --[edge_partition_fs(edge_region)] { e_p1, e_s1, e_g1_isolate} ,
    --[vertex_partition_fs(vertex_region)] { v_p1, v_s1, v_g1_isolate}

    return [cell_partition_fs(cell_region)] {
        private_1, shared_1, ghost_1, private_2, shared_2, ghost_2
    }
end
```

## Sapling > `main.rg`

The following is an excerpt from `main.rg` and lays the ground work for a toggle-able flag-dependent partitioning run. Note that when flagged, the partitioning takes time up front and shouldn’t always be enabled. (This constant should be stored in the `constants.rg` file.)

```cpp
task main()
  -------------------------------------------
  ----- DEFINE INDEX SPACES AND REGIONS -----
  -------------------------------------------

  -- Define index spaces
  var cell_id_space = ispace(int2d, {constants.nCells, constants.nVertLevels + 1})
  var edge_id_space = ispace(int2d, {constants.nEdges, constants.nVertLevels + 1})
  var vertex_id_space = ispace(int2d, {constants.nVertices, constants.nVertLevels + 1})
  var vertical_id_space = ispace(int1d, constants.nVertLevels + 1)
  var ozn_id_space = ispace(int2d, {constants.nCells, constants.nOznLevels + 1})
  var aerosol_id_space = ispace(int2d, {constants.nCells, constants.nAerLevels + 1})

  -- Define regions
  var cell_region = region(cell_id_space, cell_fs)
  var edge_region = region(edge_id_space, edge_fs)
  var vertex_region = region(vertex_id_space, vertex_fs)
  var vertical_region = region(vertical_id_space, vertical_fs)
  var phys_tbls = region(ispace(int1d, 1), phys_tbls_fs)
  fill(phys_tbls.tmin, 0.0)
  fill(phys_tbls.tmax, 0.0)
  for i = 0, constants.plenest do
    phys_tbls[0].estbl[i] = 0.0
  end
  var ozn_region = region(ozn_id_space, ozn_fs)
  var aerosol_region = region(aerosol_id_space, aerosol_fs)

  format.println("Calling load mesh...")
  load_mesh(cell_region, edge_region, vertex_region, constants.FILE_NAME, constants.GRAPH_FILE_NAME)
  format.println("Done calling load mesh...\n")

  -- var cell_partition_fs = partition_regions(constants.NUM_PARTITIONS, cell_region, edge_region, vertex_region)

  -- -- -- PLACEHOLDER PARTITIONING since cannot return multiple values 
  format.println("BEGINNING PARTITIONING")
  var color_space = ispace(int1d, constants.NUM_PARTITIONS)

  var p = partition(cell_region.partitionNumber, color_space) -- Original partition based on Metis

  -- We compute all edges around p in the process of creating halo layers
  -- which assists in our edge partitioning scheme, as well.
  var e0 = image(edge_region, p, cell_region.edgesOnCell0)
  var e1 = image(edge_region, p, cell_region.edgesOnCell1)
  var e2 = image(edge_region, p, cell_region.edgesOnCell2)
  var e3 = image(edge_region, p, cell_region.edgesOnCell3)
  var e4 = image(edge_region, p, cell_region.edgesOnCell4)
  var e5 = image(edge_region, p, cell_region.edgesOnCell5)
  var e6 = image(edge_region, p, cell_region.edgesOnCell6)
  var e7 = image(edge_region, p, cell_region.edgesOnCell7)
  var e8 = image(edge_region, p, cell_region.edgesOnCell8)
  var e9 = image(edge_region, p, cell_region.edgesOnCell9)
  var e = e0 | e1 | e2 | e3 | e4 | e5 | e6 | e7 | e8 | e9

  -- Note e denotes all edges within parition p (i.e. all edges on a cell in p)

  for i = 0, constants.nEdges do
      edge_region[{i, 0}].cellOne = rect2d { int2d { edge_region[{i, 0}].cellsOnEdge[0], 0 }, int2d { edge_region[{i, 0}].cellsOnEdge[0], constants.nVertLevels - 1 } }
      edge_region[{i, 0}].cellTwo = rect2d { int2d { edge_region[{i, 0}].cellsOnEdge[1], 0 }, int2d { edge_region[{i, 0}].cellsOnEdge[1], constants.nVertLevels - 1 } }

      -- Now, similarly declare vertices one and two for each edge
      edge_region[{i, 0}].vertexOne = rect2d { int2d { edge_region[{i, 0}].verticesOnEdge[0], 0 }, int2d { edge_region[{i, 0}].verticesOnEdge[0], constants.nVertLevels - 1 } }
      edge_region[{i, 0}].vertexTwo = rect2d { int2d { edge_region[{i, 0}].verticesOnEdge[1], 0 }, int2d { edge_region[{i, 0}].verticesOnEdge[1], constants.nVertLevels - 1 } }
  end

  var v = image(vertex_region, e, edge_region.vertexOne) | image(vertex_region, e, edge_region.vertexTwo)

  var cp_one = image(cell_region, e, edge_region.cellOne) -- This gets the cellOne on all edges in e
  var cp_two = image(cell_region, e, edge_region.cellTwo) -- This gets the cellTwo on all edges in e
  -- Combined, they should contain all cells who have an edge in e, and therefore the halo.

  var ghost_1_and_p = cp_one | cp_two
  var ghost_1 = (cp_one | cp_two) - p

  -- Edge-partitioning from ghost akin to above
  var e0_g = image(edge_region, ghost_1, cell_region.edgesOnCell0)
  var e1_g = image(edge_region, ghost_1, cell_region.edgesOnCell1)
  var e2_g = image(edge_region, ghost_1, cell_region.edgesOnCell2)
  var e3_g = image(edge_region, ghost_1, cell_region.edgesOnCell3)
  var e4_g = image(edge_region, ghost_1, cell_region.edgesOnCell4)
  var e5_g = image(edge_region, ghost_1, cell_region.edgesOnCell5)
  var e6_g = image(edge_region, ghost_1, cell_region.edgesOnCell6)
  var e7_g = image(edge_region, ghost_1, cell_region.edgesOnCell7)
  var e8_g = image(edge_region, ghost_1, cell_region.edgesOnCell8)
  var e9_g = image(edge_region, ghost_1, cell_region.edgesOnCell9)
  var e_g1 = e0_g | e1_g | e2_g | e3_g | e4_g | e5_g | e6_g | e7_g | e8_g | e9_g

  -- Apply same extension to edge partition for vertices
  var v_g1 = image(vertex_region, e_g1, edge_region.vertexOne) | image(vertex_region, e_g1, edge_region.vertexTwo)
  var e_g1_isolate = e_g1 - e -- Isolate only the 1-layer ghost edge partition
  var v_g1_isolate = v_g1 - v -- Isolate only the 1-layer ghost vertex partition
  
  -- Calculate second halo
  var gep_out = preimage(edge_region, ghost_1_and_p, edge_region.cellOne)
  var gep_in = preimage(edge_region, ghost_1_and_p, edge_region.cellTwo)
  var gcp_out = image(cell_region, gep_out, edge_region.cellTwo)
  var gcp_in = image(cell_region, gep_in, edge_region.cellOne)
  var ghost_2 = (gcp_in | gcp_out) - p -- First and second halo layers
  var gep_2 = (gep_out | gep_in) - e -- First and second halo layers 

  -- Compute all cells reachable from ghost_1. shared_1 is intersection of that set with p
  var s1ep_out = preimage(edge_region, ghost_1, edge_region.cellOne)
  var s1ep_in = preimage(edge_region, ghost_1, edge_region.cellTwo)
  var s1cp_out = image(cell_region, s1ep_out, edge_region.cellTwo)
  var s1cp_in = image(cell_region, s1ep_in, edge_region.cellOne)
  var shared_1 = p & (s1cp_out | s1cp_in) -- Cells in p bordering ghost_1
  var private_1 = p - shared_1 -- all cells in p that are not in shared_1
  -- For vertices, shared_1 is the intersection between vertices reachable by ghost and v
  var v_s1 = v & (image(vertex_region, e_g1_isolate, edge_region.vertexOne) | image(vertex_region, e_g1_isolate, edge_region.vertexTwo))
  -- Similarly, edges-shared is the intersection between the ghost edges and original edges
  var e_s1 = e & (s1ep_out | s1ep_in) -- Edges in e bordering ghost_1
  -- Finally, private is everything not shared from e
  var e_p1 = (e - e_s1) - e_g1 -- Edges in e not shared or in ghost
  var v_p1 = (v - v_s1) - v_g1 -- Vertices in v not in shared or ghost

  -- shared_2 contains shared_1 and all cells in p bordering shared_1
  var s2ep_out = preimage(edge_region, shared_1, edge_region.cellOne)
  var s2ep_in = preimage(edge_region, shared_1, edge_region.cellTwo)
  var s2cp_out = image(cell_region, s2ep_out, edge_region.cellTwo)
  var s2cp_in = image(cell_region, s2ep_in, edge_region.cellOne)
  var shared_2 = dynamic_cast(partition(disjoint, cell_region, color_space), (shared_1 | (private_1 & (s2cp_out | s2cp_in)))) -- Cells in p bordering ghost_1
  var private_2 = private_1 - shared_2 -- all cells in private_1 that are not in shared_2
  
  var cell_partition_fs = [cell_partition_fs(cell_region)] {
        private_1, shared_1, ghost_1, private_2, shared_2, ghost_2
    }
  var edge_partition_fs = [edge_partition_fs(edge_region)] { e_p1, e_s1, e_g1_isolate}
  var vertex_partition_fs = [vertex_partition_fs(vertex_region)] { v_p1, v_s1, v_g1_isolate}
  format.println("FINISHING PARTITIONING")
  -- -- -- PLACEHOLDER PARTITIONING

  --var intersect_check = private_1 & ghost_1
  --for i = 0, constants.NUM_PARTITIONS do
  --  format.println("intersect vol check {}", intersect_check[i].volume)
  --end

  --var cross_prod = cross_product(private_1, ghost_1)
  --for i = 0, constants.nCells do 
  --    format.println("cp check {}", cross_prod[i][i].volume) 
  --end

  fill(cell_region.isShared, false)
  for i = 0, constants.NUM_PARTITIONS do
    mark_shared_cells(cell_partition_fs.shared_1[i])
    mark_shared_cells(cell_partition_fs.shared_2[i])
  end

  -- Notes from before:
  --  the circuit example is helpful for showcasing some of the following logic for private/shared
  --  ultimately we seek __demand(__index_launch) to run over the whole MPAS model
  --  the following is conditioned on a global boolean flag that is toggled to use the partitions to run in parallel

  var start_time = time.clock()

  if constants.parallel_with_partitions then
    -- execute with manually partitioned ghost/private/shared system
    -- note: address "-- loop optimization failed: argument 1 interferes with itself, remove cell_region"
    __demand(__index_launch)
    for i = 0, constants.NUM_PARTITIONS do
      --format.println("Calling init_atm_case_jw...")
      -- loop optimization failed: argument 1 interferes with itself, remove cell_region
      init_atm_case_jw(cell_partition_fs.private_1[i], 
                      cell_partition_fs.shared_1[i], 
                      cell_partition_fs.ghost_1[i], 
                      edge_region, -- edge_partition_fs.edge_private_1[i] later 
                      vertex_region, -- vertex_partition_fs.vertex_private_1[i] later 
                      vertical_region)
      --format.println("Done calling init_atm_case_jw...\n")
    end 

    __demand(__index_launch)
    for i = 0, constants.NUM_PARTITIONS do 
      --format.println("Calling atm_core_init...")
      atm_core_init(cell_partition_fs.private_1[i], 
                    cell_partition_fs.shared_1[i], 
                    cell_partition_fs.ghost_1[i], 
                    edge_region, -- edge_partition_fs.edge_private_1[i] later
                    vertex_region, -- vertex_partition_fs.vertex_private_1[i] later 
                    vertical_region, 
                    phys_tbls)
      --format.println("Done calling atm_core_init...\n")
    end
  
    -- Note reordered todo NUM_TIMESTEPS as outer loop
    for j = 0, constants.NUM_TIMESTEPS do
      __demand(__index_launch)
      for i = 0, constants.NUM_PARTITIONS do
        --format.println("Calling atm_do_timestep...iteration {} \n", j)
        atm_do_timestep(cell_partition_fs.private_1[i], 
                        cell_partition_fs.shared_1[i], 
                        cell_partition_fs.ghost_1[i], 
                        edge_region, vertex_region, vertical_region, phys_tbls, j)
      end
    end

  else
    -- execute the sequential version
    for i = 0, 1 do
      --format.println("Calling init_atm_case_jw...")
      init_atm_case_jw_sequential(cell_region, cell_partition_fs.private_1[i], cell_partition_fs.shared_1[i], cell_partition_fs.ghost_1[i], edge_region, vertex_region, vertical_region)
      --format.println("Done calling init_atm_case_jw...\n")

      --format.println("Calling atm_core_init...")
      atm_core_init_sequential(cell_region, cell_partition_fs.private_1[i], cell_partition_fs.shared_1[i], cell_partition_fs.ghost_1[i], edge_region, vertex_region, vertical_region, phys_tbls)
      --format.println("Done calling atm_core_init...\n")

      for j = 0, constants.NUM_TIMESTEPS do
        --format.println("Calling atm_do_timestep...iteration {} \n", j)
        atm_do_timestep_sequential(cell_region, cell_partition_fs.private_1[i], cell_partition_fs.shared_1[i], cell_partition_fs.ghost_1[i], edge_region, vertex_region, vertical_region, phys_tbls, j)
      end
    end
  end
  
  dummy_wait(cell_region)
  var end_time = time.clock()
  var cpu_time_used = ((double) (end_time - start_time))
  format.println("DONE! Time elapsed: {}", cpu_time_used)

  -- check out circuit example for private/shared for each loop here!  
  -- for i = 0, 1 do
  -- __demand(__index_launch) -- TODO: EB - investigate this more with Alex/Elliott
  --for i = 0, constants.NUM_PARTITIONS do
    --format.println("Calling init_atm_case_jw...")
    -- loop optimization failed: argument 1 interferes with itself, remove cell_region
    --init_atm_case_jw(cell_region, cell_partition_fs.private_1[i], cell_partition_fs.shared_1[i], cell_partition_fs.ghost_1[i], edge_region, vertex_region, vertical_region)
    --format.println("Done calling init_atm_case_jw...\n")

    --format.println("Calling atm_core_init...")
    --atm_core_init(cell_region, cell_partition_fs.private_1[i], cell_partition_fs.shared_1[i], cell_partition_fs.ghost_1[i], edge_region, vertex_region, vertical_region, phys_tbls)
    --format.println("Done calling atm_core_init...\n")

    --for j = 0, constants.NUM_TIMESTEPS do
      --format.println("Calling atm_do_timestep...iteration {} \n", j)
      --atm_do_timestep(cell_region, cell_partition_fs.private_1[i], cell_partition_fs.shared_1[i], cell_partition_fs.ghost_1[i], edge_region, vertex_region, vertical_region, phys_tbls, j)
    --end
  --end

  --atm_compute_output_diagnostics(cell_region)

  --write_output_plotting(cell_region, edge_region, vertex_region)

end
regentlib.start(main)
```

# Contact Info

The above summary is, by no means, comprehensive. If there is any missing component of understanding, happenings from the summer, and ideas for next steps you’d like to reach me about, here is my contact info!

> Eric Bear, ericbear@stanford.edu, (720)-469-2671
>
