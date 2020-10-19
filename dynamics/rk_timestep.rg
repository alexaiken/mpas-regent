import "regent"

require "data_structures"
require "dynamics_tasks"

local constants = require("constants")
local format = require("std/format")


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
local cmath = terralib.includec("math.h")

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

--Comments for summarize_timestep
--Unsure what associated(block) is. Also found in atm_srk3 but ignored
--nCellsSolve, nEdgesSolve: not sure what these are
--Scalars
--Translated latCell and lonCell as lat and lon
task summarize_timestep(cr : region(ispace(int2d), cell_fs),
                        er : region(ispace(int2d), edge_fs),
                        config_print_detailed_minmax_vel : bool,
                        config_print_global_minmax_vel : bool,
                        config_print_global_minmax_sca : bool)
where reads writes (cr, er) do
  format.println("summarizing timestep")

  --Variables declared at beginning of function
  var localVals : double[5]
  var globalVals : double[5]
  var pi_const = 2.0 * cmath.asin(1.0) --How does this differ from pi?
  var scalar_min : double
  var scalar_max : double
  var global_scalar_min : double
  var global_scalar_max : double

  --This config needs to be configurable at some point
  if (constants.config_print_detailed_minmax_vel) then
    format.println("")

    -- What is this???
    --block => domain % blocklist
    --do while (associated(block))

---- GLOBAL MIN W ----
      scalar_min = 1.0e20
      var indexMax = -1
      var kMax = -1
      var latMax = 0.0
      var lonMax = 0.0
      -- What is nCellsSolve?
      -- Originally: do iCell = 1, nCellsSolve
      for iCell = 0, nCells do --TODO: Replace with nCellsSolve when resolved
        for k = 0, nVertLevels do
          if (cr[{iCell, k}].w < scalar_min) then
            scalar_min = cr[{iCell, k}].w
            indexMax = iCell
            kMax = k
            latMax = cr[{iCell, 0}].lat
            lonMax = cr[{iCell, 0}].lon
          end
        end
      end
      localVals[0] = scalar_min
      localVals[1] = [double](indexMax)
      localVals[2] = [double](kMax)
      localVals[3] = latMax
      localVals[4] = lonMax
      for i = 0, 5 do
        globalVals[i] = localVals[i]
      end

      global_scalar_min = globalVals[0]
      var indexMax_global = [int](globalVals[1])
      var kMax_global = [int](globalVals[2])
      var latMax_global = globalVals[3]
      var lonMax_global = globalVals[4]
      latMax_global *= 180.0 / pi_const
      lonMax_global *= 180.0 / pi_const
      if (lonMax_global > 180.0) then
        lonMax_global -= 360.0
      end
      -- format statement should be '(a,f9.4,a,i4,a,f7.3,a,f8.3,a)'
      --Is this a print statement?
      --call mpas_log_write(' global min w: $r k=$i, $r lat, $r lon', intArgs=(/kMax_global/), &
      --                    realArgs=(/global_scalar_min, latMax_global, lonMax_global/))
      format.println("global min w: {} k={}, {} lat, {} lon", kMax_global, global_scalar_min, latMax_global, lonMax_global)

---- GLOBAL MAX W ----
      scalar_max = -1.0e20
      indexMax = -1
      kMax = -1
      latMax = 0.0
      lonMax = 0.0
      --Original: do iCell = 1, nCellsSolve
      for iCell = 0, nCells do --TODO: replace with nCellsSolve when resolved
        for k = 0, nVertLevels do
          if (cr[{iCell, k}].w > scalar_max) then
            scalar_max = cr[{iCell, k}].w
            indexMax = iCell
            kMax = k
            latMax = cr[{iCell, k}].lat
            lonMax = cr[{iCell, k}].lon
          end
        end
      end
      localVals[0] = scalar_max
      localVals[1] = [double](indexMax)
      localVals[2] = [double](kMax)
      localVals[3] = latMax
      localVals[4] = lonMax
      for i = 0, 5 do
        globalVals[i] = localVals[i]
      end
      global_scalar_max = globalVals[0]
      indexMax_global = [int](globalVals[1])
      kMax_global = [int](globalVals[2])
      latMax_global = globalVals[3]
      lonMax_global = globalVals[4]
      latMax_global *= 180.0 / pi_const
      lonMax_global *= 180.0 / pi_const
      if (lonMax_global > 180.0) then
        lonMax_global -= 360.0
      end
      -- format statement should be '(a,f9.4,a,i4,a,f7.3,a,f8.3,a)'
      --call mpas_log_write(' global max w: $r k=$i, $r lat, $r lon', intArgs=(/kMax_global/), &
      --                    realArgs=(/global_scalar_max, latMax_global, lonMax_global/))
      format.println("global max w: {} k={}, {} lat, {} lon", kMax_global, global_scalar_max, latMax_global, lonMax_global)

---- GLOBAL MIN U ----
      scalar_min = 1.0e20
      indexMax = -1
      kMax = -1
      latMax = 0.0
      lonMax = 0.0
      -- nEdgesSolve?
      -- Original: do iEdge = 1, nEdgesSolve
      for iEdge = 0, nEdges do -- TODO: Replace with nEdgesSolve when resolved
        for k = 0, nVertLevels do
          if (er[{iEdge, k}].u < scalar_min) then
            scalar_min = er[{iEdge, k}].u
            indexMax = iEdge
            kMax = k
            latMax = er[{iEdge, 0}].lat
            lonMax = er[{iEdge, 0}].lon
          end
        end
      end
      localVals[0] = scalar_min
      localVals[1] = [double](indexMax)
      localVals[2] = [double](kMax)
      localVals[3] = latMax
      localVals[4] = lonMax
      for i = 0, 5 do
        globalVals[i] = localVals[i]
      end

      global_scalar_min = globalVals[0]
      indexMax_global = [int](globalVals[1])
      kMax_global = [int](globalVals[2])
      latMax_global = globalVals[3]
      lonMax_global = globalVals[4]
      latMax_global *= 180.0 / pi_const
      lonMax_global *= 180.0 / pi_const
      if (lonMax_global > 180.0) then
        lonMax_global -= 360.0
      end
      -- format statement should be '(a,f9.4,a,i4,a,f7.3,a,f8.3,a)'
      --call mpas_log_write(' global min u: $r k=$i, $r lat, $r lon', intArgs=(/kMax_global/), &
      --                    realArgs=(/global_scalar_max, latMax_global, lonMax_global/))
      format.println("global min u: {} k={}, {} lat, {} lon", kMax_global, global_scalar_max, latMax_global, lonMax_global)

---- GLOBAL MAX U ----
      scalar_max = -1.0e20
      indexMax = -1
      kMax = -1
      latMax = 0.0
      lonMax = 0.0
      --Original: do iEdge = 1, nEdgesSolve
      for iEdge = 0, nEdges do --TODO: replace with nEdgesSolve when resolved
        for k = 0, nVertLevels do
          if (er[{iEdge, k}].u > scalar_max) then
            scalar_max = er[{iEdge, k}].u
            indexMax = iEdge
            kMax = k
            latMax = er[{iEdge, k}].lat
            lonMax = er[{iEdge, k}].lon
          end
        end
      end
      localVals[0] = scalar_max
      localVals[1] = [double](indexMax)
      localVals[2] = [double](kMax)
      localVals[3] = latMax
      localVals[4] = lonMax
      --Original: call mpas_dmpar_maxattributes_real(domain % dminfo, scalar_max, localVals, globalVals)
      for i = 0, 5 do
        globalVals[i] = localVals[i]
      end

      global_scalar_max = globalVals[0]
      indexMax_global = [int](globalVals[1])
      kMax_global = [int](globalVals[2])
      latMax_global = globalVals[3]
      lonMax_global = globalVals[4]
      latMax_global *= 180.0 / pi_const
      lonMax_global *= 180.0 / pi_const
      if (lonMax_global > 180.0) then
        lonMax_global -= 360.0
      end
      -- format statement should be '(a,f9.4,a,i4,a,f7.3,a,f8.3,a)'
      --call mpas_log_write(' global max u: $r k=$i, $r lat, $r lon', intArgs=(/kMax_global/), &
      --                    realArgs=(/global_scalar_max, latMax_global, lonMax_global/))
      format.println("global max u: {} k={}, {} lat, {} lon", kMax_global, global_scalar_max, latMax_global, lonMax_global)

---- GLOBAL MAX WSP ----
      scalar_max = -1.0e20
      indexMax = -1
      kMax = -1
      latMax = 0.0
      lonMax = 0.0
      --Original: do iEdge = 1, nEdgesSolve
      for iEdge = 0, nEdges do -- TODO: Replace with nEdgesSolve when resolved
        for k = 0, nVertLevels do
          var spd = cmath.sqrt(cmath.pow(er[{iEdge, k}].u, 2) + cmath.pow(er[{iEdge, k}].v, 2))
          if (spd > scalar_max) then
            scalar_max = spd
            indexMax = iEdge
            kMax = k
            latMax = er[{iEdge, 0}].lat
            lonMax = er[{iEdge, 0}].lon
          end
        end
      end
      localVals[0] = scalar_max
      localVals[1] = [double](indexMax)
      localVals[2] = [double](kMax)
      localVals[3] = latMax
      localVals[4] = lonMax
      for i = 0, 5 do
        globalVals[i] = localVals[i]
      end

      global_scalar_max = globalVals[0]
      indexMax_global = [int](globalVals[1])
      kMax_global = [int](globalVals[2])
      latMax_global = globalVals[3]
      lonMax_global = globalVals[4]
      latMax_global *= 180.0 / pi_const
      lonMax_global *= 180.0 / pi_const
      if (lonMax_global > 180.0) then
        lonMax_global -= 360.0
      end
      -- format statement should be '(a,f9.4,a,i4,a,f7.3,a,f8.3,a)'
      --call mpas_log_write(' global max wsp: $r k=$i, $r lat, $r lon', intArgs=(/kMax_global/), &
      --                    realArgs=(/global_scalar_max, latMax_global, lonMax_global/))
      format.println("global max wsp: {} k={}, {} lat, {} lon", kMax_global, global_scalar_max, latMax_global, lonMax_global)

      -- Check for NaNs
      --Original: do iCell = 1, nCellsSolve
      for iCell = 0, nCells do -- TODO: Replace with nCellsSolve when resolved
        for k = 0, nVertLevels do
          if (cmath.isnan(cr[{iCell, k}].w) ~= 0) then
            format.println("NaN detected in w field.")
          end
        end
      end

      --Original: do iEdge = 1, nEdgesSolve
      for iEdge = 0, nEdges do -- TODO: replace with nEdgesSolve when resolved
        for k = 0, nVertLevels do
          if (cmath.isnan(er[{iEdge, k}].u) ~= 0) then
            format.println("NaN detected in u field.")
          end
        end
      end

      -- What is this ????
      -- block => block % next
    --end
  end

  if (config_print_global_minmax_vel) then
    format.println("")

    --block => domain % blocklist
    --do while (associated(block))

      scalar_min = 0.0
      scalar_max = 0.0
      --Original: do iCell = 1, nCellsSolve
      for iCell = 0, nCells do -- TODO: replace with nCellsSolve when resolved
        for k = 0, nVertLevels do
          scalar_min = min(scalar_min, cr[{iCell, k}].w)
          scalar_max = max(scalar_max, cr[{iCell, k}].w)
        end
      end
      --call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
      --call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
      --call mpas_log_write('global min, max w $r $r', realArgs=(/global_scalar_min, global_scalar_max/))
      format.println("global min, max w {} {}", global_scalar_min, global_scalar_max)

      scalar_min = 0.0
      scalar_max = 0.0
      --Original: do iEdge = 1, nEdgesSolve
      for iEdge = 0, nEdges do -- TODO: replace with nEdgesSolve when resolved
        for k = 0, nVertLevels do
          scalar_min = min(scalar_min, er[{iEdge, k}].u)
          scalar_max = max(scalar_max, er[{iEdge, k}].u)
        end
      end
      --call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
      --call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
      --call mpas_log_write('global min, max u $r $r', realArgs=(/global_scalar_min, global_scalar_max/))
      format.println("global min, max u {} {}", global_scalar_min, global_scalar_max)

      --block => block % next
    --end
  end

  if (config_print_global_minmax_sca) then
    if (not (config_print_global_minmax_vel or config_print_detailed_minmax_vel)) then
      format.println("")
    end

    --block => domain % blocklist
    --do while (associated(block))

      --This block of code has not been translated yet
      --do iScalar = 1, num_scalars
        --scalar_min = 0.0
        --scalar_max = 0.0
        --do iCell = 1, nCellsSolve
          --do k = 1, nVertLevels
            --scalar_min = min(scalar_min, scalars(iScalar,k,iCell))
            --scalar_max = max(scalar_max, scalars(iScalar,k,iCell))
          --end do
        --end do
        --call mpas_dmpar_min_real(domain % dminfo, scalar_min, global_scalar_min)
        --call mpas_dmpar_max_real(domain % dminfo, scalar_max, global_scalar_max)
        --call mpas_log_write(' global min, max scalar $i $r $r', intArgs=(/iScalar/), realArgs=(/global_scalar_min, global_scalar_max/))
        --format.println("global min, max scalar {} {} {}", iScalar, global_scalar_min, global_scalar_max)
      --end do

      --block => block % next
    --end
  end
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

    atm_recover_large_step_variables()

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

  mpas_reconstruct_2d(cr, er, false, true) --bools are includeHalos and on_a_sphere

  -- SKIPPING physics @ line 1610
  -- SKIPPING if(config_apply_lbcs @ line 1672 & 1714)

  summarize_timestep(cr, er, constants.config_print_detailed_minmax_vel, constants.config_print_global_minmax_vel, constants.config_print_global_minmax_sca)

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
