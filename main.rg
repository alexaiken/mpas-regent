import "regent"

require "data_structures"
require "netcdf_tasks"
require "mesh_loading"
require "init_atm_cases"
require "dynamics_tasks"
require "rk_timestep"
require "atm_core"
local constants = require("constants")

terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")

task main()
  -------------------------------------------
  ----- DEFINE INDEX SPACES AND REGIONS -----
  -------------------------------------------

  -- Define index spaces for cell IDs, vertex IDs and edge IDs
  var cell_id_space = ispace(int2d, {constants.nCells, constants.nVertLevels + 1})
  var edge_id_space = ispace(int2d, {constants.nEdges, constants.nVertLevels + 1})
  var vertex_id_space = ispace(int2d, {constants.nVertices, constants.nVertLevels + 1})
  var vertical_id_space = ispace(int1d, constants.nVertLevels + 1)

  -- Define regions
  var cell_region = region(cell_id_space, cell_fs)
  var edge_region = region(edge_id_space, edge_fs)
  var vertex_region = region(vertex_id_space, vertex_fs)
  var vertical_region = region(vertical_id_space, vertical_fs)


  load_mesh(cell_region, edge_region, vertex_region, constants.FILE_NAME, constants.GRAPH_FILE_NAME)

  --TODO: This doesn't actually return the halos yet (It creates them in the task but I haven't been able to return them). Need to return the halos.
  partition_regions(constants.NUM_PARTITIONS, cell_region, edge_region, vertex_region)

  init_atm_case_jw(cell_region, edge_region, vertex_region, vertical_region)

  atm_core_init(cell_region, edge_region, vertex_region, vertical_region)

  for i = 0, constants.NUM_TIMESTEPS do
    atm_timestep(cell_region, edge_region, vertex_region, vertical_region, i)
  end

  atm_compute_output_diagnostics(cell_region)

  write_output_plotting(cell_region, edge_region, vertex_region)

end
regentlib.start(main)
