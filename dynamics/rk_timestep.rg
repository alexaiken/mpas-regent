import "regent"

require "data_structures"
require "dynamics_tasks"

local constants = require("constants")


local nCells = 2562
local nEdges = 7680
local nVertices = 5120
local maxEdges = 10
local maxEdges2 = 20
local TWO = 2
local FIFTEEN = 15
local vertexDegree = 3
local nVertLevels = 1


local cio = terralib.includec("stdio.h")

task atm_rk_integration_setup(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs))
where reads writes (cr, er) do
  cio.printf("saving state pre-RK loop\n")
  var edge_range = rect2d { int2d{0, 0}, int2d{nEdges - 1, nVertLevels - 1} }
  var cell_range = rect2d { int2d{0, 0}, int2d{nCells - 1, nVertLevels - 1} }

  for i in edge_range do
    er[i].ru_save = er[i].ru
    er[i].u_2 = er[i].u
  end

  for i in cell_range do
    cr[i].rw_save = cr[i].rw
    cr[i].rtheta_p_save = cr[i].rtheta_p
    cr[i].rho_p_save = cr[i].rho_p

    cr[i].w_2 = cr[i].w
    cr[i].theta_m_2 = cr[i].theta_m
    cr[i].rho_zz_2 = cr[i].rho_zz
    cr[i].rho_zz_old_split = cr[i].rho_zz
    --Not sure how to translate scalars
    --scalars_2(:,:,cellStart:cellEnd) = scalars_1(:,:,cellStart:cellEnd)
  end
end

task summarize_timestep()
  cio.printf("summarizing step\n")
end

task atm_srk3(cr : region(ispace(int2d), cell_fs),
              er : region(ispace(int2d), edge_fs),
              vr : region(ispace(int2d), vertex_fs),
              vert_r : region(ispace(int1d), vertical_fs),
              dt : double)
where reads writes (cr, er, vr, vert_r) do

  -- 2 is default value from Registry.xml
  -- var definition from Registry:
  -- "Number of acoustic steps per full RK step"

  var number_of_sub_steps = 2

  -- assume no dynamics transport, if transport was option, this would be conditionally set
  var dynamics_split = constants.config_dynamics_split_steps
  var dt_dynamics = dt


  -- assume 3rd order RK
  var rk_timestep : double[3]
  rk_timestep[0] = dt_dynamics/3
  rk_timestep[1] = dt_dynamics/2
  rk_timestep[2] = dt_dynamics --holds actual time interval in the step forward

  var rk_sub_timestep : double[3]
  rk_sub_timestep[0] = dt_dynamics/3
  rk_sub_timestep[1] = dt_dynamics/number_of_sub_steps
  rk_sub_timestep[2] = dt_dynamics/number_of_sub_steps

  var number_sub_steps : int[3]
  number_sub_steps[0] = max(1, number_of_sub_steps/2)
  number_sub_steps[1] = max(1, number_of_sub_steps/2)
  number_sub_steps[2] = number_of_sub_steps -- index indicates current rk_step val, num_sub_steps[2] is number for the 3rd rk_step


  -- atm_rk_integration_setup: saves state pre-loop
  atm_rk_integration_setup(cr, er)

  atm_compute_moist_coefficients(cr, er)

  ------------DYNAMICS SUB STEP LOOP ------------------------
  -- we ignore because assume no transport split (ie dynamics_split = 1, so we only loop once)


  --compute original vertical coefficients (for initial step)
  atm_compute_vert_imp_coefs(cr, vert_r, rk_sub_timestep[0])



  ------------------------------------------------------------
  ---------------- BEGIN Runge-Kutta loop --------------------
  -----------------------------------------------------------

  for rk_step = 0, 3 do
    cio.printf("\nRK STEP: %d\n", rk_step)
    if rk_step == 1 then
      --recompute vertical coefficients, same ones for step 3
      atm_compute_vert_imp_coefs(cr, vert_r, rk_sub_timestep[rk_step])
    end

    atm_compute_dyn_tend(cr, er, vr, vert_r, rk_sub_timestep[rk_step], dt, constants.config_horiz_mixing, constants.config_mpas_cam_coef, constants.config_mix_full, constants.config_rayleigh_damp_u)

    atm_set_smlstep_pert_variables(cr, er, vert_r)


    -- SKIPPING if(config_apply_lbcs @ line 683)

    -------------------------------------
    ------begin acoustic steps loop -----
    -------------------------------------
    for small_step = 0, number_sub_steps[rk_step] + 1 do

      cio.printf("performing acoustic substeps within a rk step\n")

      atm_advance_acoustic_step(cr, er, vert_r, rk_sub_timestep[rk_step], small_step)

      atm_divergence_damping_3d(cr, er, rk_sub_timestep[rk_step])
    end

    atm_recover_large_step_variables(cr, er, number_sub_steps[rk_step], vert_r, rk_step, dt, cf1, cf2, cf3)

    -- SKIPPING if(config_apply_lbcs @ line 934)


    -- SKIPPING  if (config_scalar_advection .and. (.not. config_split_dynamics_transport) ) @ line 993

    atm_compute_solve_diagnostics(cr, er, vr, false)

    -- SKIPPING  if (config_scalar_advection .and. (.not. config_split_dynamics_transport) ) @ line 1246

    -- SKIPPING if(config_apply_lbcs @ line 1253)

    --for iEdge = 0, nEdges do
    cio.printf("Horizontal normal velocity at Edge %d is %f\n", 0, er[{0, 0}].u)
    --end

  end

  --SKIPPING if (dynamics_substep < dynamics_split) @ 1282

  atm_rk_dynamics_substep_finish(cr, er, 1, dynamics_split) --Third argument is dynamics_substep, which should eventually go in a loop

  --------------DYNAMICS SUB STEP LOOP WOULD END HERE-----------------

  -- SKIPPING if (config_scalar_advection .and. config_split_dynamics_transport) @ 1355

  --TODO: where is this defined in MPAS
  --mpas_reconstruct

  -- SKIPPING physics @ line 1610
  -- SKIPPING if(config_apply_lbcs @ line 1672 & 1714)

  summarize_timestep()

  --for iCell = 0, nCells do
  --  cio.printf("Surface pressure at Cell %d is %f\n", iCell, cr[{iCell, 0}].pressure_p)
  --end



end

--__demand(__cuda)
task atm_timestep(cr : region(ispace(int2d), cell_fs),
                  er : region(ispace(int2d), edge_fs),
                  vr : region(ispace(int2d), vertex_fs),
                  vert_r : region(ispace(int1d), vertical_fs),
                  dt : double)
where reads writes (cr, er, vr, vert_r) do
--MPAS also uses nowTime and itimestep parameters; itimestep only for physics/IAU, and ignoring timekeeping for now

  atm_srk3(cr, er, vr, vert_r, dt)

end
