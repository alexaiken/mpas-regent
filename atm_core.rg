import "regent"
require "data_structures"
require "physics/atmphys_init"
require "physics/atmphys_driver"

local constants = require("constants")
local format = require("std/format")

--__demand(__cuda)
task atm_core_init(cr : region(ispace(int2d), cell_fs),
                   er : region(ispace(int2d), edge_fs),
                   vr : region(ispace(int2d), vertex_fs),
                   vert_r : region(ispace(int1d), vertical_fs),
                   phys_tbls : region(ispace(int1d), phys_tbls_fs))
where 
  reads writes (cr, er, vr, vert_r, phys_tbls) 
do
  format.println("Calling atm_core_init...")

  atm_compute_signs(cr, er, vr)

  atm_adv_coef_compression(cr, er)

  --config_coef_3rd_order = 0.25 in namelist
  atm_couple_coef_3rd_order(0.25, cr, er)

  atm_init_coupled_diagnostics(cr, er, vert_r)

  atm_compute_solve_diagnostics(cr, er, vr, false, -1) --last param is hollingsworth

  mpas_reconstruct_2d(cr, er, false, true) --bools are includeHalos and on_a_sphere

  -- if (moist_physics) then
  physics_init(cr, er, phys_tbls)
  -- endif

  atm_compute_mesh_scaling(cr, er, true)

  --config_zd: default 22000.0, config_xnutr: default 0.2. From config
  atm_compute_damping_coefs(22000, 0.2, cr)

end

task atm_do_timestep(cr : region(ispace(int2d), cell_fs),
                     er : region(ispace(int2d), edge_fs),
                     vr : region(ispace(int2d), vertex_fs),
                     vert_r : region(ispace(int1d), vertical_fs),
                     phys_tbls : region(ispace(int1d), phys_tbls_fs),
                     dt : double)
where 
  reads (phys_tbls),
  reads writes (cr, er, vr, vert_r) 
do

  --if(moist_physics) then
  physics_timetracker()
  physics_driver(cr, phys_tbls)
  --end

  atm_timestep(cr, er, vr, vert_r, dt)
end