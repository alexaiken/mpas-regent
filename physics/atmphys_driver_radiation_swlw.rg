import "regent"
require "data_structures"
require "physics/ra_cam"

local constants = require("constants")

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
  --sinob = sin(obecl)
      
  -- calculate longitude of the sun from vernal equinox:

  if (julian > 80.0) then
    sxlong = dpd * (julian - 80.0)
  else
    sxlong = dpd * (julian + 285.0)
  end
  sxlong *= degrad
  --arg = sinob * sin(sxlong)
  --declin = asin(arg)
  --decdeg = declin / degrad

  -- solar constant eccentricity factor (paltridge and platt 1976)

  --djul = julian * 360.0 / 365.0
  --rjul = djul * degrad
  --eccfac = 1.000110 + 0.034221 * cos(rjul) + 0.001280 * sin(rjul) + 0.000719 * 
  --      cos(2 * rjul) + 0.000077 * sin(2 * rjul)
  --solcon = constants.solcon_0 * eccfac
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
  var curr_julday : double = 0.0 --TODO: Placeholder
  var ozncyc : bool = true --TODO: Placeholder
  camrad(cr, phys_tbls, constants.nOznLevels, curr_julday, ozncyc)
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

  for j = jts, jte do
    for i = 0, nCellsSolve do
      cr[{i, 0}].sfc_emiss_p = cr[{i, 0}].sfc_emiss
      cr[{i, 0}].tsk_p       = cr[{i, 0}].skintemp
      cr[{i, 0}].snow_p      = cr[{i, 0}].snow
      cr[{i, 0}].xice_p      = cr[{i, 0}].xice
      cr[{i, 0}].xland_p     = cr[{i, 0}].xland
    end
  end
  for j = jts, jte do
    for k = 0, nVertLevels do
      for i = 0, nCellsSolve do
        cr[{i, k}].cldfrac_p = cr[{i, k}].cldfrac
      end
    end
  end

  -- initialization:
  for j = jts, jte do
    for k = 0, nVertLevels do
      for i = 0, nCellsSolve do
        cr[{i, k}].f_ice  = 0.0
        cr[{i, k}].f_rain = 0.0
      end
    end
  end

  for j = jts, jte do
    for i = 0, nCellsSolve do
      cr[{i, 0}].glw_p      = 0.0
      cr[{i, 0}].lwcf_p     = 0.0
      cr[{i, 0}].lwdnb_p    = 0.0
      cr[{i, 0}].lwdnbc_p   = 0.0
      cr[{i, 0}].lwdnt_p    = 0.0
      cr[{i, 0}].lwdntc_p   = 0.0
      cr[{i, 0}].lwupb_p    = 0.0
      cr[{i, 0}].lwupbc_p   = 0.0
      cr[{i, 0}].lwupt_p    = 0.0
      cr[{i, 0}].lwuptc_p   = 0.0
      cr[{i, 0}].olrtoa_p   = 0.0
    end
  
    for k = 0, nVertLevels do
      for i = 0, nCellsSolve do
        cr[{i, k}].rthratenlw_p = 0.0
      end
    end
  end

  if ([rawstring](radt_lw_scheme) == "rrtmg_lw") then -- string should be trimmed
    if ([rawstring](microp_scheme) == "mp_thompson" or [rawstring](microp_scheme) == "mp_wsm6") then
      if (config_microp_re) then

        for j = jts, jte do
          for k = 0, nVertLevels do
            for i = 0, nCellsSolve do
              cr[{i, k}].recloud_p = cr[{i, k}].re_cloud
              cr[{i, k}].reice_p   = cr[{i, k}].re_ice
              cr[{i, k}].resnow_p  = cr[{i, k}].re_snow
            end
          end
        end
      else
        for j = jts, jte do
          for k = 0, nVertLevels do
            for i = 0, nCellsSolve do
              cr[{i, k}].recloud_p = 0.0
              cr[{i, k}].reice_p   = 0.0
              cr[{i, k}].resnow_p  = 0.0
            end
          end
        end
      end
      for j = jts, jte do
        for k = 0, nVertLevels do
          for i = 0, nCellsSolve do
            cr[{i, k}].rrecloud_p = 0.0
            cr[{i, k}].rreice_p   = 0.0
            cr[{i, k}].rresnow_p  = 0.0
          end
        end
      end
    end

    for j = jts, jte do
      for k = 0, nVertLevels + 2 do
        for i = 0, nCellsSolve do
          --lwdnflx_p[{i, k, j}]  = 0.0
          --lwdnflxc_p[{i, k, j}] = 0.0
          --lwupflx_p[{i, k, j}]  = 0.0
          --lwupflxc_p[{i, k, j}] = 0.0
        end
      end
    end

    if (config_o3climatology) then
      -- ozone mixing ratio:
      for k = 0, nOznLevels do
        ozn_r[{0, k}].pin_p = ozn_r[{0, k}].pin
      end
      for j = jts, jte do
        for k = 0, nOznLevels do
          for i = 0, nCellsSolve do
            ozn_r[{i, k}].o3clim_p = ozn_r[{i, k}].o3clim
          end
        end
      end

      var nlevs = nVertLevels + 1
      var ncols = nCellsSolve + 1
      --if(.not.allocated(p2d) ) allocate(p2d(its:ite,kts:kte) )
      --if(.not.allocated(o32d)) allocate(o32d(its:ite,kts:kte))
      for j = jts, jte do
        for i = 0, nCellsSolve do
          for k = 0, nVertLevels do
            cr[{i, k}].o32d = 0.0
            cr[{i, k}].p2d  = cr[{i, k}].pres_hyd_p / 100.0
          end
        end
        --vinterp_ozn(cr, ozn_r, 1, ncols, ncols, nlevs)
        for i = 0, nCellsSolve do
          for k = 0, nVertLevels do
            cr[{i, k}].o3vmr = cr[{i, k}].o32d
          end
        end
      end
      --if(allocated(p2d))  deallocate(p2d)
      --if(allocated(o32d)) deallocate(o32d)
    else
      for k = 0, nOznLevels do
        ozn_r[{0, k}].pin_p = 0.0
      end
      for j = jts, jte do
        for k = 0, nOznLevels do
          for i = 0, nCellsSolve do
            ozn_r[{i, k}].o3clim_p = 0.0
          end
        end
      end
    end
  end
  if ([rawstring](radt_lw_scheme) == "cam_lw") then
    for j = jts, jte do
      for i = 0, nCellsSolve do
        cr[{i, 0}].xlat_p = cr[{i, 0}].lat / degrad
        cr[{i, 0}].xlon_p = cr[{i, 0}].lon / degrad
        cr[{i, 0}].sfc_albedo_p = cr[{i, 0}].sfc_albedo

        cr[{i, 0}].coszr_p      = 0.0
        cr[{i, 0}].gsw_p        = 0.0
        cr[{i, 0}].swcf_p       = 0.0
        cr[{i, 0}].swdnb_p      = 0.0
        cr[{i, 0}].swdnbc_p     = 0.0
        cr[{i, 0}].swdnt_p      = 0.0
        cr[{i, 0}].swdntc_p     = 0.0
        cr[{i, 0}].swupb_p      = 0.0
        cr[{i, 0}].swupbc_p     = 0.0
        cr[{i, 0}].swupt_p      = 0.0
        cr[{i, 0}].swuptc_p     = 0.0
      end
      for k = 0, nVertLevels do
        for i = 0, nCellsSolve do
          cr[{i, k}].rthratensw_p = 0.0
          cr[{i, k}].cemiss_p     = 0.0
          cr[{i, k}].taucldc_p    = 0.0
          cr[{i, k}].taucldi_p    = 0.0
        end
      end
    end

    -- On the first time-step of each model run, the local arrays absnxt_p, absnst_p,
    -- and emstot_p are filled with the MPAS arrays abstot, absnxt, and emstot. If it
    -- is a new run, these three arrays will be initialized to zero;If a restart run,
    -- these three arrays will be filled with the restart values.
    if (xtime_s < 1.0e-12) then
      for j = jts, jte do
        for n = 0, constants.cam_abs_dim1 do -- Originally from 1 to cam_abs_dim1
          for k = 0, nVertLevels do
            for i = 0, nCellsSolve do
              cr[{i, k}].absnxt_p[n] = 0.0
            end
          end
        end
        for n = 0, constants.cam_abs_dim2 do -- Originally from 1 to cam_abs_dim2
          for k = 0, nVertLevels + 1 do
            for i = 0, nCellsSolve do
              cr[{i, k}].abstot_p[n] = 0.0
            end
          end
        end
        for k = 0, nVertLevels + 1 do
          for i = 0, nCellsSolve do
            cr[{i, k}].emstot_p = 0.0
          end
        end
      end
    end

    -- ozone mixing ratio:
    for k = 0, nOznLevels do
      ozn_r[{0, k}].pin_p = ozn_r[{0, k}].pin
    end
    for n = 0, constants.nMonths do
      for j = jts, jte do
        for k = 0, nOznLevels do
          for i = 0, nCellsSolve do
            ozn_r[{i, k}].ozmixm_p[n] = ozn_r[{i, k}].ozmixm[n]
          end
        end
      end
    end
    --aerosol mixing ratio:
    for k = 0, nAerLevels do
      aer_r[{0, k}].m_hybi_p = aer_r[{0, k}].m_hybi
    end
    for i = 0, nCellsSolve do
      for j = jts, jte do
        cr[{i, 0}].m_psp_p = cr[{i, 0}].m_ps
        cr[{i, 0}].m_psn_p = cr[{i, 0}].m_ps
      end
    end
    for n = 0, nAerosols do
      for j = jts, jte do
        for k = 0, nAerLevels do
          for i = 0, nCellsSolve do
            aer_r[{i, k}].aerosolcp_p[n] = aer_r[{i, k}].aerosols[n]
            aer_r[{i, k}].aerosolcn_p[n] = aer_r[{i, k}].aerosols[n]
          end
        end
      end
    end
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
  if ([rawstring](radt_lw_scheme) == "cam_lw") then
    radconst(0, 0, 0) --TODO: Placeholder arguments! I don't know where actual arguments are from
    radt = constants.config_dt / 60.0
    var curr_julday : double = 0.0 --TODO: Placeholder
    var ozncyc : bool = true --TODO: Placeholder
    camrad(cr, phys_tbls, constants.nOznLevels, curr_julday, ozncyc)
  end

  --radiation_lw_to_MPAS(cr, radt_lw_scheme, microp_scheme, config_microp_re)
end