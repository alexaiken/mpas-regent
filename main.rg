import "regent"

require "data_structures"
require "netcdf_tasks"
require "mesh_loading"
require "init_atm_cases"
require "dynamics_tasks"
require "rk_timestep"
require "atm_core"
local constants = require("constants")
local format = require("std/format")

terralib.linklibrary("/home/arjunk1/spack/opt/spack/linux-ubuntu20.04-broadwell/gcc-9.3.0/netcdf-c-4.7.4-zgdvh4hxthdhb3mlsviwhgatvbfnslog/lib/libnetcdf.so")

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

  var cell_partition_fs = partition_regions(constants.NUM_PARTITIONS, cell_region, edge_region, vertex_region)

  format.println("Calling init_atm_case_jw...")
  init_atm_case_jw(cell_region, edge_region, vertex_region, vertical_region)
  format.println("Done calling init_atm_case_jw...\n")

  format.println("Calling atm_core_init...")
  atm_core_init(cell_region, edge_region, vertex_region, vertical_region, phys_tbls)
  format.println("Done calling atm_core_init...\n")

  for i = 0, constants.NUM_TIMESTEPS do
<<<<<<< HEAD
    format.println("Calling atm_do_timestep...iteration {} \n", i)
    atm_do_timestep(cell_region, edge_region, vertex_region, vertical_region, phys_tbls, constants.config_dt)
=======
    atm_do_timestep(cell_region, edge_region, vertex_region, vertical_region, cell_partition_fs, phys_tbls, i)
>>>>>>> Mesh loading compiles with partitions
  end

  atm_compute_output_diagnostics(cell_region)

  write_output_plotting(cell_region, edge_region, vertex_region)

end
regentlib.start(main)
