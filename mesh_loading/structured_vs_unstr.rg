

--Question: unstructured region vs. int1d index space for the cell IDs


----OPTION 1----
var cell_id_region = region(ispace(ptr, mesh_vars.nCells), int)
var vertex_id_region = region(ispace(ptr, mesh_vars.nVertices), int)
var edge_id_region = region(ispace(ptr, mesh_vars.nEdges), int)

--hexagon centers
var x_i_region = region(ispace(ptr, mesh_vars.nCells), x_i)
--triangle centers
var x_v_region = region(ispace(ptr, mesh_vars.nVertices), x_v)
--edge midpoints
var x_e_region = region(ispace(ptr, mesh_vars.nEdges), x_e)
--distance bw hexagon centers (x_i locations)
var d_e_region = region(ispace(ptr, mesh_vars.nEdges), d_e)
--distance bw triangle centers
var l_e_region = region(ispace(ptr, mesh_vars.nEdges), l_e)


--SET REGIONS
--edges on a cell defining P_i
var EC_region = region(ispace(ptr, mesh_vars.nCells), edges_on_cell(maxEdges, edge_id_region)) 
-- edges on a vertex defining D_v
var EV_region = region(ispace(ptr, mesh_vars.nVertices), edges_on_vertex(edge_id_region))  
-- cells on an edge (primal mesh cells sharing an edge) + edges on the cell pair
var CE_region = region(ispace(ptr, mesh_vars.nEdges), cells_on_edge(maxEdges, cell_id_region, edge_id_region))
-- cells on a vertex (dual mesh cells s)
var CV_region = region(ispace(ptr, mesh_vars.nEdges), cells_on_vertex(cell_id_region))
-- vertices on edge (dual mesh cells sharing a primal edge)
var VE_region = region(ispace(ptr, mesh_vars.nEdges), vertices_on_edge(vertex_id_region))
-- vertices on cell (dual mesh cells centered on vertices of a primal cell)
var VI_region = region(ispace(ptr, mesh_vars.nCells), vertices_on_cell(maxEdges, vertex_id_region))  
-- edge pair associated w a vertex and a cell
var EVC_region = region(ispace(ptr, mesh_vars.nCells*mesh_vars.nEdgesMax), edge_pair(edge_id_region))


--edges that outline P_i
fspace edges_on_cell(maxEdges : int,
                     re : region(int)) {
  edges : ptr(int, re)[6]
}

-edge on primal (hexagonal) mesh made up of two points in the x_v region (centers of dual cells)
fspace edge(rxv : region(x_v),
            rv : region(int),
            rc : region(int),
            re : region(int),
            maxEdges : int,
            rce : region(cells_on_edge(maxEdges, rc, re)),
            rve : region(vertices_on_edge(rv))) {
  angle : double,
  edgeID : int,
  {point1, point2} : ptr(x_v, rxv),
  cells : ptr(cells_on_edge(maxEdges, rc, re), rce),
  vertices : ptr(vertices_on_edge(rv), rve)
}





----OPTION 2----
-- Define index spaces for cell IDs, vertex IDs and edge IDs
var cell_id_space = ispace(int1d, nCells)
var vertex_id_space = ispace(int1d, nVertices)
var edge_id_space = ispace(int1d, nEdges)


-- hexagon centers
var x_i_region = region(cell_id_space, x_i)
--triangle centers
var x_v_region = region(vertex_id_space, x_v)
--edge midpoints
var x_e_region = region(edge_id_space, x_e)
--distance bw hexagon centers (x_i locations)
var d_e_region = region(edge_id_space, d_e)
--distance bw triangle centers
var l_e_region = region(edge_id_space, l_e)


--edges on a cell defining P_i
var EC_region = region(cell_id_space, edges_on_cell)
-- edges on a vertex defining D_v
var EV_region = region(vertex_id_space, edges_on_vertex)
-- cells on an edge (primal mesh cells sharing an edge) + edges on the cell pair
var CE_region = region(edge_id_space, cells_on_edge)
-- cells on a vertex (dual mesh cells s)
var CV_region = region(vertex_id_space, cells_on_vertex)
-- vertices on edge (dual mesh cells sharing a primal edge)
var VE_region = region(edge_id_space, vertices_on_edge)
-- vertices on cell (dual mesh cells centered on vertices of a primal cell)
var VI_region = region(cell_id_space, vertices_on_cell)


--vertex ids for the vertices outlining P_i
fspace vertices_on_cell {
    vertices : int[maxEdges]
}

--edge on primal (hexagonal) mesh made up of two points in the x_v region (centers of dual cells)
fspace edge {
    edgeID : int,
    center : x_e,
    angle : double,
    nEdgesOnEdge : int,
    dvEdge : l_e,
    dv1Edge : double,
    dv2Edge : double,
    dcEdge : d_e,
    cells : cells_on_edge,
    vertices : vertices_on_edge
}