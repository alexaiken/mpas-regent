import "regent"

local c = terralib.includec("stdlib.h")
local cio = terralib.includec("stdio.h")

struct mesh_dimensions {
  nEdgesMax : int
  nCells : int
  nVertices : int
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

--cell IDs of two cells sharing an edge
fspace cells_on_edge(rc : region(int)) {
  {cell1, cell2} : ptr(int, rc)
}

--vertex ids for triangles centered on the points of an edge
fspace vertices_on_edge(rv : region(int)) {
  {d1, d2} : ptr(int, rv)
}

--edge on primal (hexagonal) mesh made up of two points in the x_v region (centers of dual cells)
fspace edge(rxv : region(x_v),
            rv : region(int),
            rc : region(int),
            rce : region(cells_on_edge(rc)),
            rve : region(vertices_on_edge(rv))) {
  angle : double,
  {point1, point2} : ptr(x_v, rxv),
  cells : ptr(cells_on_edge(rc), rce),
  vertices : ptr(vertices_on_edge(rv), rve)
}

--edges that outline P_i
fspace p_edge_set(rxv : region(x_v),
                  rv : region(int),
                  rc : region(int),
                  rce : region(cells_on_edge(rc)),
                  rve : region(vertices_on_edge(rv)),
                  re : region(edge(rxv, rv, rc, rce, rve))) {
  {_0, _1, _2, _3, _4, _5, _6} : ptr(edge(rxv, rv, rc, rce, rve), re)
}

--edges that intersect at vertex of D_v
fspace d_edge_set(rxv : region(x_v),
                  rv : region(int),
                  rc : region(int),
                  rce : region(cells_on_edge(rc)),
                  rve : region(vertices_on_edge(rv)),
                  re : region(edge(rxv, rv, rc, rce, rve))) {  
  {_0, _1, _2} : ptr(edge(rxv, rv, rc, rce, rve), re)
}

--cell ids of the primal cells sharing a vertex
fspace cells_on_vertex(rc : region(int)) {
  {cell1, cell2, cell3} : ptr(int, rc),
}

--vertex ids on a primal cell
fspace vertices_on_cell(rv : region(int)) {
  {d1, d2, d3, d4, d5, d6} : ptr(int, rv)
}

--a triangle cell on the dual-mesh
fspace D_v(rxi: region(x_i),
           rxv : region(x_v),
           rv : region(int),
           rc : region(int),
           rce : region(cells_on_edge(rc)),
           rve : region(vertices_on_edge(rv)),
           re : region(edge(rxv, rv, rc, rce, rve)),
           res : region(d_edge_set(rxv, rv, rc, rce, rve, re)),
           rcv : region(cells_on_vertex(rc))) {
  vertexID : int, --could make a ptr into rv
  -- map to x_v center point
  center : ptr(x_v, rxv),
  -- map to neighbor vertex ids
  {neighbor1, neighbor2, neighbor3} : ptr(int, rv),

  edge_set : ptr(d_edge_set(rxv, rv, rc, rce, rve, re), res),
  cells : ptr(cells_on_vertex(rc), rcv)
}

--a hexagon cell on the primal-mesh
fspace P_i(rxi : region(x_i),
           rxv : region(x_v),
           rv : region(int),
           rc : region(int),
           rce : region(cells_on_edge(rc)),
           rve : region(vertices_on_edge(rv)),
           re : region(edge(rxv, rv, rc, rce, rve)),
           res : region(p_edge_set(rxv, rv, rc, rce, rve, re)),
           rvc : region(vertices_on_cell(rv))) {
  cellID : int,
  -- map to x_i center
  center : ptr(x_i, rxi),
  -- map to neighbor cell ids
  {neighbor1, neighbor2, neighbor3, neighbor4, neighbor5, neighbor6} : ptr(int, rc),
  
  edge_set : ptr(p_edge_set(rxv, rv, rc, rce, rve, re), res),
  vertices : ptr(vertices_on_cell(rv), rvc) 
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


  var cell_id_region = region(ispace(ptr, mesh_vars.nCells), int)
  var vertex_id_region = region(ispace(ptr, mesh_vars.nVertices), int)
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
  
  -- cells on an edge (primal mesh cells sharing an edge)
  var CE_region = region(ispace(ptr, mesh_vars.nEdges), cells_on_edge(cell_id_region))
  -- vertices on edge (dual mesh cells sharing a primal edge)
  var VE_region = region(ispace(ptr, mesh_vars.nEdges), vertices_on_edge(vertex_id_region))
  -- region of all edges
  var edge_region = region(ispace(ptr, mesh_vars.nEdges), edge(x_v_region, vertex_id_region, cell_id_region, CE_region, VE_region))
  --edges on a cell defining P_i
  var EC_region = region(ispace(ptr, mesh_vars.nCells), p_edge_set(x_v_region, vertex_id_region, cell_id_region, CE_region, VE_region, edge_region)) 
  -- edges on a vertex defining D_v
  var EV_region = region(ispace(ptr, mesh_vars.nVertices), d_edge_set(x_v_region, vertex_id_region, cell_id_region, CE_region, VE_region, edge_region))  
  -- cells on a vertex (dual mesh cells s)
  var CV_region = region(ispace(ptr, mesh_vars.nEdges), cells_on_vertex(cell_id_region))
  -- vertices on cell (dual mesh cells centered on vertices of a primal cell)
  var VI_region = region(ispace(ptr, mesh_vars.nCells), vertices_on_cell(vertex_id_region))  


  --CELL REGIONS

  --primal (hexagonal) cells
  var P_i_region = region(ispace(ptr, mesh_vars.nCells), P_i(x_i_region, x_v_region, vertex_id_region, cell_id_region, CE_region, VE_region, edge_region, EC_region, VI_region))
  --dual (triangular) cells
  var D_v_region = region(ispace(ptr, mesh_vars.nVertices), D_v(x_i_region, x_v_region, vertex_id_region, cell_id_region, CE_region, VE_region, edge_region, EV_region, CV_region))

  --TODO: add last two set regions
  -- edges of cell pair meeting at an edge e
  -- edge pair associated with vertex v and edge cell i
end


task main()
  create_element_regions()
  cio.printf("\n***regions compiled successfully***\n\n") 
end

regentlib.start(main)
