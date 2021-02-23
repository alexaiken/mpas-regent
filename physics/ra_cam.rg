import "regent"
require "data_structures"
require "physics/ra_cam_cld_support"
require "physics/ra_cam_radctl_support"

local constants = require("constants")
local c = regentlib.c
local format = require("std/format")

task param_cldoptics_calc()
  cldefr()
  cldems()
  cldovrlap()
end

------------------------------------------------------------------------- 
-- MPAS-Model/src/core_atmosphere/physics/physics_wrf/module_ra_cam.F
------------------------------------------------------------------------- 
-- 
-- Purpose: 
-- Driver for radiation computation.
-- 
-- Method: 
-- Radiation uses cgs units, so conversions must be done from
-- model fields to radiation fields.
--
-- Author: CCM1,  CMS Contact: J. Truesdale
-- 
-------------------------------------------------------------------------
task radctl(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
            camrad_1d_r : region(ispace(int1d), camrad_1d_fs),
            camrad_2d_r : region(ispace(int2d), camrad_2d_fs),
            lchnk : int,
            ncol : int,
            pcols : int,
            pver : int, pverp : int, pverr : int, pverrp : int,
            julian : double,
            ozmixmj : region(ispace(int3d), double),
            ozmix : region(ispace(int2d), double),
            levsiz : double,
            pin : region(ispace(int1d), double),
            ozncyc : bool,
            aerosoljp : region(ispace(int3d), double),
            aerosoljn : region(ispace(int3d), double),
            m_hybi : region(ispace(int1d), double),
            paerlev : int)
where reads (
        cr.{pmid, pint, t}, 
        phys_tbls, 
        ozmixmj, 
        ozmix, 
        pin,
        m_psp, m_psn,
        aerosoljp, aerosoljn,
        m_hybi
      ),
      writes (
        cr.pmid,
        ozmix
      )
do

  -----------------------------Local variables-----------------------------

  var i : int
  var k : int

  -- indexes of gases in constituent array
  var in2o : int
  var ich4 : int
  var if11 : int
  var if12 : int

  var eccf : double       -- Earth/sun distance factor

  var radctl_1d_r = region(ispace(int1d, pcols), radctl_1d_fs)
  var radctl_2d_pver_r = region(ispace(int2d, {pcols, pver}), radctl_2d_pver_fs)
  var radctl_2d_pverr_r = region(ispace(int2d, {pcols, pverr}), radctl_2d_pverr_fs)
  var radctl_3d_aero_r = region(ispace(int3d, {pcols, nspint, naer_groups}), radctl_3d_aero_fs)
  
  var pnm = region(ispace(int2d, {pcols, pverrp}), double)      -- Model interface pressures (dynes/cm2)
  
  var aerosol = region(ispace(int3d, {pcols, pver, constants.naer_all}), double) -- aerosol mass mixing ratios
  fill(aerosol, 0.0) -- TODO: TEMP
  
  var scales = region(ispace(int1d, constants.naer_all), double)  -- scaling factors for aerosols
  -------------------------------------------------------------------------

  --
  -- Interpolate ozone volume mixing ratio to model levels
  --
  -- WRF: added pin, levsiz, ozmix here
  oznint(julian, ozmixmj, ozmix, levsiz, pcols, ozncyc)

  radozn(cr, radctl_2d_pverr_r, ncol, pcols, pver, pin, levsiz, ozmix)

  --
  -- Set chunk dependent radiation input
  --
  radinp(cr, radctl_2d_pverr_r, radctl_2d_pverrp_r, ncol, pver, pverp)

  --
  -- Solar radiation computation
  --
  if (dosw) then
    --
    -- calculate heating with aerosols
    --
    aqsat(cr, phys_tbls, radctl_2d_pverr_r, ncol, pver)

    for i = 0, ncol do
      for j = 0, pver do
        radctl_2d_pverr_r[{i, j}].rh = qm1(i, j, 0) / qsat(i, j) *
            ((1.0 - epsilo) * qsat(i, j) + epsilo) /
            ((1.0 - epsilo) * qm1(i, j, 0) + epsilo)
      end
    end

    -- REGENT NOTE: Skipping lines 1662-1715 in module_ra_cam.F because radforce always == false

    get_int_scales() -- TODO

    get_aerosol(cr, phys_tbls, lchnk, julian, m_psp, m_psn, aerosoljp, aerosoljn, m_hybi, paerlev, 
                pcols, pver, pverp, pverr, pverrp, aerosol, scales)

    aerosol_indirect() -- TODO
  
    radcswmx() -- TODO

    -- -- tls ---------------------------------------------------------------2
    --
    -- Convert units of shortwave fields needed by rest of model from CGS to MKS
    --
    for i = 0, ncol do
      radctl_1d_r[i].solin = radctl_1d_r[i].solin * 1.e-3
      camrad_1d_r[i].fsds = camrad_1d_r[i]fsds * 1.e-3
      radctl_1d_r[i].fsnirt = radctl_1d_r[i].fsnirt * 1.e-3
      radctl_1d_r[i].fsnrtc = radctl_1d_r[i].fsnrtc * 1.e-3
      radctl_1d_r[i].fsnirtsq = radctl_1d_r[i].fsnirtsq * 1.e-3
      camrad_1d_r[i].fsnt = camrad_1d_r[i].fsnt * 1.e-3
      camrad_1d_r[i].fsns = camrad_1d_r[i].fsns * 1.e-3
      radctl_1d_r[i].fsntc = radctl_1d_r[i].fsntc * 1.e-3
      radctl_1d_r[i].fsnsc = radctl_1d_r[i].fsnsc * 1.e-3
      radctl_1d_r[i].fsdsc = radctl_1d_r[i].fsdsc * 1.e-3
      radctl_1d_r[i].fsntoa = radctl_1d_r[i].fsntoa * 1.e-3
      radctl_1d_r[i].fsntoac = radctl_1d_r[i].fsntoac * 1.e-3
      swcf[i] = radctl_1d_r[i].fsntoa - radctl_1d_r[i].fsntoac
    end
    for i = 0, ncol do
      for j = 0, pver do
        radctl_2d_pver_r[{i, j}].ftem = camrad_2d_r[{i, j}].qrs / cpair
      end
    end

    -- Added upward/downward total and clear sky fluxes
    for k = 0, pverp
      for i = 0, ncol
        camrad_2d_r[{i, k}].fsup  = camrad_2d_r[{i, k}].fsup * 1.e-3
        camrad_2d_r[{i, k}].fsupc = camrad_2d_r[{i, k}].fsupc * 1.e-3
        camrad_2d_r[{i, k}].fsdn  = camrad_2d_r[{i, k}].fsdn * 1.e-3
        camrad_2d_r[{i, k}].fsdnc = camrad_2d_r[{i, k}].fsdnc * 1.e-3
      end
    end
  end

  --
  -- Longwave radiation computation
  --
  if (dolw) then
    get_int_scales() -- TODO

    get_aerosol(cr, phys_tbls, lchnk, julian, m_psp, m_psn, aerosoljp, aerosoljn, m_hybi, paerlev, 
                pcols, pver, pverp, pverr, pverrp, aerosol, scales)

    --
    -- Convert upward longwave flux units to CGS
    --
    for i = 0, ncol do
      lwupcgs(i) = lwups(i)
    end

    --
    -- Do longwave computation. If not implementing greenhouse gas code then
    -- first specify trace gas mixing ratios. If greenhouse gas code then:
    --  o ixtrcg   => indx of advected n2o tracer
    --  o ixtrcg+1 => indx of advected ch4 tracer
    --  o ixtrcg+2 => indx of advected cfc11 tracer
    --  o ixtrcg+3 => indx of advected cfc12 tracer
    --
    if (trace_gas) then
      radclwmx() -- TODO
    else
      trcmix() -- TODO

      radclwmx() -- TODO
    end

    --
    -- Convert units of longwave fields needed by rest of model from CGS to MKS
    --
    for i = 0, ncol do
      flnt(i)  = flnt(i) * 1.e-3
      flut(i)  = flut(i) * 1.e-3
      flutc(i) = flutc(i) * 1.e-3
      flns(i)  = flns(i) * 1.e-3
      flntc(i) = flntc(i) * 1.e-3
      flnsc(i) = flnsc(i) * 1.e-3
      flwds(i) = flwds(i) * 1.e-3
      lwcf(i)  = flutc(i) - flut(i)
    end

    -- Added upward/downward total and clear sky fluxes
    for k = 0, pverp do
      for i = 0, ncol do
        flup(i,k)  = flup(i,k) * 1.e-3
        flupc(i,k) = flupc(i,k) * 1.e-3
        fldn(i,k)  = fldn(i,k) * 1.e-3
        fldnc(i,k) = fldnc(i,k) * 1.e-3
      end
    end
  end
end

task camrad(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
            levsiz : int,
            julian : double,
            ozncyc : bool,
            paerlev : int,
            naer_c : int)
where 
  reads (cr.{pmid, pint, t}, phys_tbls),
  writes (cr.pmid) -- TEMP needed for fill
do
  -----------------------------Local variables-----------------------------

  var lchnk : int = 10     -- TODO: TEMP
  var ncol : int = 10      -- TODO: TEMP
  var pcols : int = 10     -- TODO: TEMP
  var pver : int = 10      -- TODO: TEMP
  var pverp : int = 10     -- TODO: TEMP
  var pverr : int = 10     -- TODO: TEMP
  var pverrp : int = 10    -- TODO: TEMP

  var pcnst : int
  var pnats : int
  var ppcnst : int
  var i : int
  var j : int
  var k : int
  var ii : int
  var kk : int
  var kk1 : int
  var m : int
  var n : int

  var begchunk : int
  var endchunk : int

  var nyrm : int
  var nyrp : int

  var doymodel : double
  var doydatam : double
  var doydatap : double
  var deltat : double
  var fact1 : double
  var fact2 : double

  var xt24 : double
  var tloctm : double
  var hrang : double
  var xxlat : double
  var oldxt24 : double

  var camrad_1d_r = region(ispace(int1d, constants.nCells), camrad_1d_fs)
  fill(camrad_1d_r.m_psjp, 0.0); -- TODO: TEMP
  fill(camrad_1d_r.m_psjn, 0.0); -- TODO: TEMP

  var camrad_2d_r = region(ispace(int2d, {constants.nCells, constants.nVertLevels}, camrad_2d_fs))

  var ozmixmj = region(ispace(int3d, {constants.nCells, levsiz, constants.nMonths}), double)
  fill(ozmixmj, 0.0); -- TODO: TEMP
  
  var ozmix = region(ispace(int2d, {constants.nCells, levsiz}), double)
  fill(ozmix, 0.0); -- TODO: TEMP

  var pin = region(ispace(int1d, levsiz), double)
  fill(pin, 0.0); -- TODO: TEMP

  var aerosoljp = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  var aerosoljn = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  fill(aerosoljp, 0.0); -- TODO: TEMP
  fill(aerosoljn, 0.0); -- TODO: TEMP

  var m_hybi = region(ispace(int1d, paerlev), double)
  fill(m_hybi, 0.0); -- TODO: TEMP

  var abstot = region(isapce(int3d, {constants.nCells, constants.nVertLevels, constants.nVertLevels}), double) -- Total absorptivity
  var absnxt = region(isapce(int3d, {constants.nCells, constants.nVertLevels, 4}), double) -- Total nearest layer absorptivity
  var emstot = region(isapce(int3d, {constants.nCells, constants.nVertLevels + 1}), double) -- Total emissivity

  fill(cr.pmid, 0.0); -- TODO: TEMP

  -----------------------------

  format.println("Calling camrad...")

  param_cldoptics_calc()

  radctl(
    cr, 
    phys_tbls, 
    camrad_1d_r, camrad_2d_r,
    lchnk,
    ncol, 
    pcols, 
    pver, 
    pverp, 
    pverr, 
    pverrp, 
    julian,
    ozmixmj, 
    ozmix, 
    levsiz, 
    pin, 
    ozncyc,
    aerosoljp, aerosoljn, 
    m_hybi, 
    paerlev
  )

  format.println("Camrad done")

end