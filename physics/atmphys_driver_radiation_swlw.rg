import "regent"
require "data_structures"
require "physics/ra_cam"

local constants = require("constants")
local c = regentlib.c
local format = require("std/format")

local nCellsSolve = constants.nCells --simplification
local nAerLevels = constants.nAerLevels
local nAerosols = constants.nAerosols
local nOznLevels = constants.nOznLevels
local levsiz = constants.nOznLevels
local nVertLevels = constants.nVertLevels

local jts = 0
local jte = 0

--Math function imports
--TODO: These cause an error "macros must be called from inside terra code".
--Note: do we still need these?
--min = regentlib.fmin(double, double)
--sin = regentlib.sin(double)
--asin = regentlib.asin(double)
--cos = regentlib.cos(double)
--acos = regentlib.acos(double)

task radconst(julian : double,
              degrad : double,
              dpd : double)
  var obecl : double
  var sinob : double
  var sxlong : double
  var arg : double
  var decdeg : double
  var djul : double
  var rjul : double
  var eccfac : double
  var declin : double
  var solcon : double

  -- obecl : obliquity = 23.5 degree.

  obecl = 23.5 * degrad
  sinob = sin(obecl)

  -- calculate longitude of the sun from vernal equinox:

  if (julian > 80.0) then
    sxlong = dpd * (julian - 80.0)
  else
    sxlong = dpd * (julian + 285.0)
  end
  sxlong *= degrad
  arg = sinob * sin(sxlong)
  declin = asin(arg)
  decdeg = declin / degrad

  -- solar constant eccentricity factor (paltridge and platt 1976)

  djul = julian * 360.0 / 365.0
  rjul = djul * degrad
  eccfac = 1.000110 + 0.034221 * cos(rjul) + 0.001280 * sin(rjul) + 0.000719 * 
           cos(2 * rjul) + 0.000077 * sin(2 * rjul)
  solcon = constants.solcon_0 * eccfac
end

----------------
-- SHORT WAVE --
----------------

task allocate_radiation_sw()
end

task deallocate_radiation_sw()
end

task radiation_sw_from_MPAS()
end

task radiation_sw_to_MPAS()
end

task driver_radiation_sw(cr : region(ispace(int2d), cell_fs),
                         phys_tbls : region(ispace(int1d), phys_tbls_fs))
where 
  reads (phys_tbls),
  reads writes (cr)
do
  radiation_sw_from_MPAS()
  radconst(0, 0, 0)
  var m_psp = region(ispace(int2d, {5, 5}), double) -- TODO: Placeholder
  var m_psn = region(ispace(int2d, {5, 5}), double) -- TODO: Placeholder
  var curr_julday : double = 0.0 --TODO: Placeholder
  var ozncyc : bool = true --TODO: Placeholder
  var paerlev : int = 0 --TODO: Placeholder
  var naer_c : int = 0 --TODO: Placeholder
  camrad(cr, phys_tbls, constants.nOznLevels, m_psp, m_psn, curr_julday, ozncyc, paerlev, naer_c)
  radiation_sw_to_MPAS()
end

---------------
-- LONG WAVE --
---------------

task allocate_radiation_lw()
end

task deallocate_radiation_lw()
end

task vinterp_ozn()
end

__demand(__cuda)
task radiation_lw_from_MPAS(cr : region(ispace(int2d), cell_fs),
                            ozn_r : region(ispace(int2d), ozn_fs),
                            aer_r : region(ispace(int2d), aerosol_fs),
                            radt_lw_scheme : regentlib.string,
                            microp_scheme : regentlib.string,
                            config_microp_re : bool,
                            config_o3climatology : bool,
                            xtime_s : double,
                            degrad : double)
where
  reads (cr.{cldfrac, lat, lon, m_ps, pres_hyd_p, re_cloud, re_ice, re_snow, sfc_albedo, sfc_emiss,
             skintemp, snow, xice, xland},
         aer_r.{aerosols, m_hybi},
         ozn_r.{o3clim, ozmixm, pin}),
  writes (cr.{absnxt_p, abstot_p, cemiss_p, cldfrac_p, coszr_p, emstot_p, f_ice, f_rain, glw_p, gsw_p,
              lwcf_p, lwdnb_p, lwdnbc_p, lwdnt_p, lwdntc_p, lwupb_p, lwupbc_p, lwupt_p, lwuptc_p, m_psn_p, 
              m_psp_p, o3vmr, olrtoa_p, p2d, recloud_p, reice_p, resnow_p, rrecloud_p, rreice_p, rresnow_p, 
              rthratenlw_p, rthratensw_p, sfc_albedo_p, sfc_emiss_p, snow_p, swcf_p, swdnb_p, swdnbc_p, 
              swdnt_p, swdntc_p, swupb_p, swupbc_p, swupt_p, swuptc_p, taucldc_p, taucldi_p, tsk_p, xice_p, 
              xland_p, xlat_p, xlon_p},
          aer_r.{aerosolcn_p, aerosolcp_p, m_hybi_p},
          ozn_r.{o3clim_p, ozmixm_p, pin_p}),
  reads writes (cr.o32d)
do
  var cell_range_1d = rect2d{ {0, 0}, {nCellsSolve - 1, 0} }
  var cell_range_2d = rect2d{ {0, 0}, {nCellsSolve - 1, nVertLevels - 1} }
  var cell_range_extra_2d = rect2d{ {0, 0}, {nCellsSolve - 1, nVertLevels} } -- for loop goes from 0 to nVertLevels + 1 
  var ozn_range_1d = rect2d{ {0, 0}, {0, nOznLevels - 1} }
  var ozn_range_2d = rect2d{ {0, 0}, {nCellsSolve - 1, nOznLevels - 1} }


  -- Removed j loop because it doesn't seem to do anything in current code.
  --for j = jts, jte do
    for iCell in cell_range_1d do
      cr[iCell].sfc_emiss_p = cr[iCell].sfc_emiss
      cr[iCell].tsk_p       = cr[iCell].skintemp
      cr[iCell].snow_p      = cr[iCell].snow
      cr[iCell].xice_p      = cr[iCell].xice
      cr[iCell].xland_p     = cr[iCell].xland
    end
  --end
    for iCell in cell_range_2d do
      cr[iCell].cldfrac_p = cr[iCell].cldfrac
    end
  --end

  -- initialization:
  --for j = jts, jte do
    for iCell in cell_range_2d do
      cr[iCell].f_ice  = 0.0
      cr[iCell].f_rain = 0.0
    end
  --end

  --for j = jts, jte do
    for iCell in cell_range_1d do
      cr[iCell].glw_p      = 0.0
      cr[iCell].lwcf_p     = 0.0
      cr[iCell].lwdnb_p    = 0.0
      cr[iCell].lwdnbc_p   = 0.0
      cr[iCell].lwdnt_p    = 0.0
      cr[iCell].lwdntc_p   = 0.0
      cr[iCell].lwupb_p    = 0.0
      cr[iCell].lwupbc_p   = 0.0
      cr[iCell].lwupt_p    = 0.0
      cr[iCell].lwuptc_p   = 0.0
      cr[iCell].olrtoa_p   = 0.0
    end
  
    for iCell in cell_range_2d do
        cr[iCell].rthratenlw_p = 0.0
    end
  --end

  if ([rawstring](radt_lw_scheme) == "rrtmg_lw") then -- string should be trimmed
    if ([rawstring](microp_scheme) == "mp_thompson" or [rawstring](microp_scheme) == "mp_wsm6") then
      if (config_microp_re) then

        --for j = jts, jte do
          for iCell in cell_range_2d do
            cr[iCell].recloud_p = cr[iCell].re_cloud
            cr[iCell].reice_p   = cr[iCell].re_ice
            cr[iCell].resnow_p  = cr[iCell].re_snow
          end
        --end
      else
        --for j = jts, jte do
          for iCell in cell_range_2d do
            cr[iCell].recloud_p = 0.0
            cr[iCell].reice_p   = 0.0
            cr[iCell].resnow_p  = 0.0
          end
        --end
      end
      --for j = jts, jte do
        for iCell in cell_range_2d do
          cr[iCell].rrecloud_p = 0.0
          cr[iCell].rreice_p   = 0.0
          cr[iCell].rresnow_p  = 0.0
        end
      --end
    end

    -- Commented out the loop since inner part of the loop was commented out.
    --for j = jts, jte do
    --  for k = 0, nVertLevels + 2 do
    --    for i = 0, nCellsSolve do
          --lwdnflx_p[{i, k, j}]  = 0.0
          --lwdnflxc_p[{i, k, j}] = 0.0
          --lwupflx_p[{i, k, j}]  = 0.0
          --lwupflxc_p[{i, k, j}] = 0.0
    --    end
    --  end
    --end
    
    if (config_o3climatology) then
      -- ozone mixing ratio:
      for index in ozn_range_1d do
        ozn_r[index].pin_p = ozn_r[index].pin
      end
      --for j = jts, jte do
        for index in ozn_range_2d do
          ozn_r[index].o3clim_p = ozn_r[index].o3clim
        end
      --end

      -- Only get used in call vinterp_ozn
      --var nlevs = nVertLevels + 1
      --var ncols = nCellsSolve + 1
      --if(.not.allocated(p2d) ) allocate(p2d(its:ite,kts:kte) )
      --if(.not.allocated(o32d)) allocate(o32d(its:ite,kts:kte))
      --for j = jts, jte do
        for iCell in cell_range_2d do
          cr[iCell].o32d = 0.0
          cr[iCell].p2d  = cr[iCell].pres_hyd_p / 100.0
        end
        --vinterp_ozn(cr, ozn_r, 1, ncols, ncols, nlevs)
        for iCell in cell_range_2d do
            cr[iCell].o3vmr = cr[iCell].o32d
        end
      --end
      --if(allocated(p2d))  deallocate(p2d)
      --if(allocated(o32d)) deallocate(o32d)
    else
      for index in ozn_range_1d do
        ozn_r[index].pin_p = 0.0
      end
      --for j = jts, jte do
        for index in ozn_range_2d do
          ozn_r[index].o3clim_p = 0.0
        end
      --end
    end
  end
  if ([rawstring](radt_lw_scheme) == "cam_lw") then
    --for j = jts, jte do
      for iCell in cell_range_1d do
        cr[iCell].xlat_p = cr[iCell].lat / degrad
        cr[iCell].xlon_p = cr[iCell].lon / degrad
        cr[iCell].sfc_albedo_p = cr[iCell].sfc_albedo

        cr[iCell].coszr_p      = 0.0
        cr[iCell].gsw_p        = 0.0
        cr[iCell].swcf_p       = 0.0
        cr[iCell].swdnb_p      = 0.0
        cr[iCell].swdnbc_p     = 0.0
        cr[iCell].swdnt_p      = 0.0
        cr[iCell].swdntc_p     = 0.0
        cr[iCell].swupb_p      = 0.0
        cr[iCell].swupbc_p     = 0.0
        cr[iCell].swupt_p      = 0.0
        cr[iCell].swuptc_p     = 0.0
      end
      for iCell in cell_range_2d do
        cr[iCell].rthratensw_p = 0.0
        cr[iCell].cemiss_p     = 0.0
        cr[iCell].taucldc_p    = 0.0
        cr[iCell].taucldi_p    = 0.0
      end
    --end

    -- On the first time-step of each model run, the local arrays absnxt_p, absnst_p,
    -- and emstot_p are filled with the MPAS arrays abstot, absnxt, and emstot. If it
    -- is a new run, these three arrays will be initialized to zero;If a restart run,
    -- these three arrays will be filled with the restart values.
    if (xtime_s < 1.0e-12) then
      --for j = jts, jte do
        for iCell in cell_range_2d do
          for n = 0, constants.cam_abs_dim1 do -- Originally from 1 to cam_abs_dim1
            cr[iCell].absnxt_p[n] = 0.0
          end
        end
        for iCell in cell_range_extra_2d do
          for n = 0, constants.cam_abs_dim2 do -- Originally from 1 to cam_abs_dim2
            cr[iCell].abstot_p[n] = 0.0
          end
        end
        for iCell in cell_range_extra_2d do
          cr[iCell].emstot_p = 0.0
        end
      --end
    end

    -- ozone mixing ratio:
    for index in ozn_range_1d do
      ozn_r[index].pin_p = ozn_r[index].pin
    end
    --for n = 0, constants.nMonths do -- Moved inside the loop
      --for j = jts, jte do
        for index in ozn_range_2d do
          for n = 0, constants.nMonths do -- Moved inside the loop
            ozn_r[index].ozmixm_p[n] = ozn_r[index].ozmixm[n]
          end
        end
      --end
    --end
    --aerosol mixing ratio:
    var aer_range_1d = rect2d{ {0, 0}, {0, nAerLevels - 1} }
    var aer_range_2d = rect2d{ {0, 0}, {nCellsSolve - 1, nAerLevels - 1} }
    for index in aer_range_1d do
      aer_r[index].m_hybi_p = aer_r[index].m_hybi
    end
    for iCell in cell_range_1d do
      --for j = jts, jte do
        cr[iCell].m_psp_p = cr[iCell].m_ps
        cr[iCell].m_psn_p = cr[iCell].m_ps
      --end
    end
    --for n = 0, nAerosols do -- moved inside for loop
      --for j = jts, jte do
        for index in aer_range_2d do
          for n = 0, nAerosols do -- moved inside for loop
            aer_r[index].aerosolcp_p[n] = aer_r[index].aerosols[n]
            aer_r[index].aerosolcn_p[n] = aer_r[index].aerosols[n]
          end
        end
      --end
    --end
  end
end

task radiation_lw_to_MPAS()
end

task driver_radiation_lw(cr : region(ispace(int2d), cell_fs),
                         phys_tbls : region(ispace(int1d), phys_tbls_fs),
                         radt_lw_scheme : regentlib.string,
                         config_o3climatology : bool,
                         microp_scheme : regentlib.string,
                         config_microp_re : bool)
where 
  reads (phys_tbls),
  reads writes (cr)
do

  format.println("Calling driver_radiation_lw...")

  var radt : double
  --radiation_lw_from_MPAS()

  -- I don't think we are actually doing rrtmg?
  --if ([rawstring](radt_lw_scheme) == "rrtmg_lw") then
  --  var o3input : int = 0
  --  if (config_o3climatology) then
  --    o3input = 2
  --  end
  --  rrtmg_lwrad() -- o3input will be an argument
  --else if
  if (c.strcmp([rawstring](radt_lw_scheme), "cam_lw") == 0) then
    radconst(0, 0, 0) --TODO: Placeholder arguments! I don't know where actual arguments are from
    radt = constants.config_dt / 60.0
    var m_psp = region(ispace(int2d, {5, 5}), double) -- TODO: Placeholder
    var m_psn = region(ispace(int2d, {5, 5}), double) -- TODO: Placeholder
    var curr_julday : double = 0.0 --TODO: Placeholder
    var ozncyc : bool = true --TODO: Placeholder
    var paerlev : int = 0 --TODO: Placeholder
    var naer_c : int = 0 -- TODO: Placeholder
    camrad(cr, phys_tbls, constants.nOznLevels, m_psp, m_psn, curr_julday, ozncyc, paerlev, naer_c)
  end

  --radiation_lw_to_MPAS(cr, radt_lw_scheme, microp_scheme, config_microp_re)
end
