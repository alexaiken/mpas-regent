import "regent"

local c = terralib.includec("stdlib.h")
local cio = terralib.includec("stdio.h")

struct mesh_dimensions {
  nEdgesMax : int,
  nCells : int,
  nVertices : int,
  nEdges : int
}


-----------------------------------------------
------- FIELD SPACES FOR MESH ELEMENTS --------
-----------------------------------------------


--location of center of primal-mesh cells
fspace x_i {
  {lat, lon} : double
}

--location of center of dual-mesh cells
fspace x_v {
  {lat, lon} : double
}

--location of edge points where velocity is defined
fspace x_e {
  {lat, lon} : double
}

--distance bw neighboring primal-mesh centers, i.e. x_i locations
fspace d_e {
  distance : double
}

--distance between neighboring dual-mesh centers, i.e. x_v  locations
fspace l_e {
  distance : double
}


--edges that outline P_i
fspace edges_on_cell(maxEdges : int,
                     re : region(int)) {
  edges : ptr(int, re)[6]
}

--edges that intersect at vertex of D_v
fspace edges_on_vertex(re : region(int)) {
  edges : ptr(int, re)[3]
}

--cell IDs of two cells sharing an edge
fspace cells_on_edge(maxEdges : int,
                     rc : region(int),
                     re : region(int)) {
  {cell1, cell2} : ptr(int, rc),
  edges_of_cell_pair : ptr(int, re)[11]
}


--cell ids of the primal cells sharing a vertex
fspace cells_on_vertex(rc : region(int)) {
  {cell1, cell2, cell3} : ptr(int, rc),
}


--vertex ids for triangles centered on the points of an edge
fspace vertices_on_edge(rv : region(int)) {
  {d1, d2} : ptr(int, rv)
}

--vertex ids for the vertices outlining P_i
fspace vertices_on_cell(maxEdges : int,
                        rv : region(int)) {
  vertices : ptr(int, rv)[6]
}

--edge on primal (hexagonal) mesh made up of two points in the x_v region (centers of dual cells)
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

--edge pair associated w vertex v and mesh cell i
fspace edge_pair(re : region(int)) {
  {e1, e2} : ptr(int, re)
}


--a triangle cell on the dual-mesh
fspace D_v(rxv : region(x_v),
           rv : region(int),
           rc : region(int),
           re : region(int),
           rev : region(edges_on_vertex(re)),
           rcv : region(cells_on_vertex(rc))) {
  vertexID : int, --could make a ptr into rv
  -- map to x_v center point
  center : ptr(x_v, rxv),
  -- map to neighbor vertex ids
  neighbors : ptr(int, rv)[3],

  edge_set : ptr(edges_on_vertex(re), rev),
  cells : ptr(cells_on_vertex(rc), rcv)
}

--a hexagon cell on the primal-mesh
fspace P_i(rxi : region(x_i),
           rv : region(int),
           rc : region(int),
           re : region(int),
           maxEdges : int,
           rec : region(edges_on_cell(maxEdges, re)),
           rvc : region(vertices_on_cell(maxEdges, rv)),
           rep : region(edge_pair(re))) {
  cellID : int,
  -- map to x_i center
  center : ptr(x_i, rxi),
  -- map to neighbor cell ids
  neighbors : ptr(int, rc)[6],
  
  edge_set : ptr(edges_on_cell(maxEdges, re), rec),
  vertices : ptr(vertices_on_cell(maxEdges, rv), rvc),
  edge_pairs : ptr(edge_pair(re), rep)[6]
}


-----------------------------------------------
--------------- REGION SETUP ------------------
-----------------------------------------------

task create_element_regions()
  -- will populate from reading in mesh data later
  var mesh_vars : mesh_dimensions
  mesh_vars.nEdgesMax = 6
  mesh_vars.nCells = 100
  mesh_vars.nVertices = 2*mesh_vars.nCells
  mesh_vars.nEdges = 3*mesh_vars.nCells
  
  --FIELD REGIONS

  var myArray : x_i[5]
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
  

  var maxEdges = mesh_vars.nEdgesMax --compiler error if you call this field directly: symbol has no field "nEdgesMax"

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


  var edge_region = region(ispace(ptr, mesh_vars.nEdges), edge(x_v_region, vertex_id_region, cell_id_region, edge_id_region, maxEdges, CE_region, VE_region))

  --CELL REGIONS

  --primal (hexagonal) cells
  var P_i_region = region(ispace(ptr, mesh_vars.nCells), P_i(x_i_region, vertex_id_region, cell_id_region, edge_id_region, maxEdges, EC_region, VI_region, EVC_region))
  --dual (triangular) cells
  var D_v_region = region(ispace(ptr, mesh_vars.nVertices), D_v(x_v_region, vertex_id_region, cell_id_region, edge_id_region, EV_region, CV_region))

end


task main()
  create_element_regions()
  cio.printf("\n***regions compiled successfully***\n\n") 
end

regentlib.start(main)
