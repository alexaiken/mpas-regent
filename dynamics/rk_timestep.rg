import "regent"

require "data_structures"
require "dynamics_tasks"

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


task rk_integration_setup()
  cio.printf("saving state pre-RK loop\n")
end

task summarize_timestep()
  cio.printf("summarizing step\n")
end

task atm_srk3(dt : double,
              vr : region(ispace(int1d), vertex_fs),
              er : region(ispace(int1d), edge_fs),
              cr : region(ispace(int1d), cell_fs))

  -- 2 is default value from Registry.xml
  -- var definition from Registry: 
  -- "Number of acoustic steps per full RK step"

  var number_of_sub_steps = 2

  -- assume no dynamics transport, if transport was option, this would be conditionally set
  var dynamics_split = 1
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


  -- rk_integration_setup: saves state pre-loop
  rk_integration_setup()
  
  atm_compute_moist_coefficients() 

  ------------DYNAMICS SUB STEP LOOP ------------------------
  -- we ignore because assume no transport split (ie dynamics_split = 1, so we only loop once)
  

  --compute original vertical coefficients (for initial step)
  atm_compute_vert_imp_coefs()  

  

  ------------------------------------------------------------
  ---------------- BEGIN Runge-Kutta loop --------------------
  -----------------------------------------------------------  

  for rk_step = 1, 4 do
    cio.printf("\nRK STEP: %d\n", rk_step)
    if rk_step == 2 then 
      --recompute vertical coefficients, same ones for step 3
      atm_compute_vert_imp_coefs() 
    end
    
    atm_compute_dyn_tend()
    
    atm_set_smlstep_pert_variables()


    -- SKIPPING if(config_apply_lbcs @ line 683)

    -------------------------------------
    ------begin acoustic steps loop -----
    -------------------------------------
    for small_step = 1, number_sub_steps[rk_step - 1] + 1 do
      
      cio.printf("performing acoustic substeps within a rk step\n")

      atm_advance_acoustic_step()
      atm_divergence_damping_3d()
    end

    atm_recover_large_step_variables()

    -- SKIPPING if(config_apply_lbcs @ line 934)


    -- SKIPPING  if (config_scalar_advection .and. (.not. config_split_dynamics_transport) ) @ line 993

    atm_compute_solve_diagnostics()

    -- SKIPPING  if (config_scalar_advection .and. (.not. config_split_dynamics_transport) ) @ line 1246

    -- SKIPPING if(config_apply_lbcs @ line 1253)

  end

  --SKIPPING if (dynamics_substep < dynamics_split) @ 1282 

  atm_rk_dynamics_substep_finish()

  --------------DYNAMICS SUB STEP LOOP WOULD END HERE-----------------

  -- SKIPPING if (config_scalar_advection .and. config_split_dynamics_transport) @ 1355

  --TODO: where is this defined in MPAS
  --mpas_reconstruct

  -- SKIPPING physics @ line 1610
  -- SKIPPING if(config_apply_lbcs @ line 1672 & 1714)

  summarize_timestep()

end

task atm_timestep(dt : double,
                  vr : region(ispace(int1d), vertex_fs),
                  er : region(ispace(int1d), edge_fs),
                  cr : region(ispace(int1d), cell_fs))
--MPAS also uses nowTime and itimestep parameters; itimestep only for physics/IAU, and ignoring timekeeping for now

  atm_srk3(dt, vr, er, cr)
  

end


