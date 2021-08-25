import "regent"
require "data_structures"
require "physics/atmphys_init"
require "physics/atmphys_driver"
require "netcdf_tasks"

local constants = require("constants")
local format = require("std/format")

--__demand(__cuda)
task atm_core_init(cr : region(ispace(int2d), cell_fs),
                   cpr : region(ispace(int2d), cell_fs),
                   csr : region(ispace(int2d), cell_fs),
                   cgr : region(ispace(int2d), cell_fs),
                   er : region(ispace(int2d), edge_fs),
                   vr : region(ispace(int2d), vertex_fs),
                   vert_r : region(ispace(int1d), vertical_fs),
                   phys_tbls : region(ispace(int1d), phys_tbls_fs))
where 
  reads writes (cr, er, vr, vert_r, phys_tbls),
  cpr <= cr, csr <= cr, cgr <= cr
do
  
  atm_compute_signs(cr, er, vr)
  
  atm_adv_coef_compression(cr, er)
  
  --config_coef_3rd_order = 0.25 in namelist
  atm_couple_coef_3rd_order(constants.config_coef_3rd_order, cr, er)
  
  atm_init_coupled_diagnostics(cr, er, vert_r)
  -- This is location 20 in debug mode.
  
  atm_compute_solve_diagnostics(cr, er, vr, false, -1) --last param is hollingsworth
  -- This is location 21 in debug mode.
  
  mpas_reconstruct_2d(cr, er, false, true) --bools are includeHalos and on_a_sphere
  -- This is location 22 in debug mode.
  
  -- if (moist_physics) then
  physics_init(cr, er, phys_tbls)
  -- This is location 23 in debug mode.
  -- endif
  
  atm_compute_mesh_scaling(cr, cpr, csr, cgr, er, true)
  -- This is location 24 in debug mode.
  
  --config_zd: default 22000.0, config_xnutr: default 0.2. From config
  atm_compute_damping_coefs(constants.config_zd, constants.config_xnutr, cr)
  -- This is location 25 in debug mode.

  end

task atm_do_timestep(cr : region(ispace(int2d), cell_fs),
                     cpr : region(ispace(int2d), cell_fs),
                     csr : region(ispace(int2d), cell_fs),
                     cgr : region(ispace(int2d), cell_fs),
                     er : region(ispace(int2d), edge_fs),
                     vr : region(ispace(int2d), vertex_fs),
                     vert_r : region(ispace(int1d), vertical_fs),
                     phys_tbls : region(ispace(int1d), phys_tbls_fs),
                     dt : double)
where 
  reads (phys_tbls),
  reads writes (cr, er, vr, vert_r),
  cpr <= cr, csr <= cr, cgr <= cr
do

  format.println("Inside atm_do_timestep...\n")

  --if(moist_physics) then
  physics_timetracker()
  physics_driver(cr, phys_tbls)
  --end
  atm_timestep(cr, cpr, csr, cgr, er, vr, vert_r, dt)
end