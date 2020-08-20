import "regent"

require "data_structures"
require "netcdf_tasks"
require "mesh_loading"

local c = regentlib.c
local cio = terralib.includec("stdio.h")
local clib = terralib.includec("stdlib.h")

terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")
local netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

local nCells = 2562
local nEdges = 7680
local nVertices = 5120


local FILE_NAME = "/home/users/arjunk1/regent_project_2020/mpas-regent/mesh_loading/x1.2562.grid.nc"
local GRAPH_FILE_NAME = "/home/users/arjunk1/regent_project_2020/mpas-regent/mesh_loading/x1.2562.graph.info.part.16"
local NUM_PARTITIONS = 16


task main()
  -------------------------------------------
  ----- DEFINE INDEX SPACES AND REGIONS -----
  -------------------------------------------

  -- Define index spaces for cell IDs, vertex IDs and edge IDs
  var cell_id_space = ispace(int1d, nCells)
  var edge_id_space = ispace(int1d, nEdges)
  var vertex_id_space = ispace(int1d, nVertices)

  -- Define regions
  var cell_region = region(cell_id_space, cell_fs)
  var edge_region = region(edge_id_space, edge_fs)
  var vertex_region = region(vertex_id_space, vertex_fs)


  --TODO: Pass in FILE_NAME and GRAPH_FILE_NAME: how to pass strings?
  load_mesh(cell_region, edge_region, vertex_region)
  calculate_evc(cell_region, edge_region, vertex_region)

  --TODO: This doesn't actually return the halos yet. Need to return the halos.
  partition_regions(NUM_PARTITIONS, cell_region, edge_region, vertex_region)
  
  write_output(cell_region, edge_region, vertex_region)

end
regentlib.start(main)
