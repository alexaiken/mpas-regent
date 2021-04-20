import "regent"
require "data_structures"
require "physics_data_structures"
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
            lchnk : int, ncol : int, pcols : int, pver : int, pverp : int, pverr : int, pverrp : int, ppcnst : int, pcnst : int,
            qm1 : region(ispace(int3d), double),
            julian : double,
            ozmixmj : region(ispace(int3d), double), ozmix : region(ispace(int2d), double),
            levsiz : double,
            pin : region(ispace(int1d), double),
            ozncyc : bool,
            m_psp : region(ispace(int1d), double),
            m_psn : region(ispace(int1d), double),
            aerosoljp : region(ispace(int3d), double),
            aerosoljn : region(ispace(int3d), double),
            m_hybi : region(ispace(int1d), double),
            paerlev : int,
            dolw : bool, dosw : bool,
            swcf : region(ispace(int1d), double),
            lwcf : region(ispace(int1d), double),
            flut : region(ispace(int1d), double))
where reads writes (
        camrad_1d_r.{fsds, fsnt, fsns, flnt, flns, flwds, lwups},
        camrad_2d_r.{pmid, pint, t, qrs, fsup, fsupc, fsdn, fsdnc, flup, flupc, fldn, fldnc},
        qm1,
        ozmixmj, ozmix,
        pin,
        m_psp, m_psn,
        aerosoljp, aerosoljn,
        m_hybi,
        flut
      ),
      reads (
        phys_tbls.{tmin, tmax, estbl, idxVOLC}
      ),
      writes (
        camrad_2d_r.pmid,
        ozmix,
        swcf,
        lwcf
      )
do

  -----------------------------Input variables-----------------------------

   var nspint : int = 19           -- Num of spctrl intervals across solar spectrum
   var naer_groups : int = 7       -- Num of aerosol groups for optical diagnostics
          -- current groupings are sul, sslt, all carbons, all dust, background, and all aerosols

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
  fill(radctl_2d_pverr_r.esat, 0.0)
  var radctl_3d_aero_r = region(ispace(int3d, {pcols, nspint, naer_groups}), radctl_3d_aero_fs)
  
  var pnm = region(ispace(int2d, {pcols, pverrp}), double)      -- Model interface pressures (dynes/cm2)
  
  var aerosol = region(ispace(int3d, {pcols, pver, constants.naer_all}), double) -- aerosol mass mixing ratios
  fill(aerosol, 0.0) -- TODO: TEMP
  
  var scales = region(ispace(int1d, constants.naer_all), double)  -- scaling factors for aerosols
  -------------------------------------------------------------------------

  -- TODO actual computation

  oznint(julian, ozmixmj, ozmix, levsiz, pcols, ozncyc)

  radozn(cr, camrad_2d_r, radctl_2d_pverr_r, ncol, pcols, pver, pin, levsiz, ozmix)

  radinp(cr, camrad_2d_r, radctl_2d_pverr_r, ncol, pver, pverp, pnm)

  aqsat(cr, phys_tbls, camrad_2d_r, radctl_2d_pverr_r, ncol, pver)

  get_int_scales() -- TODO

  get_aerosol() -- TODO

  aerosol_indirect() -- TODO

  radcswmx() -- TODO

  trcmix() -- TODO

  radclwmx() -- TODO
end

task camrad(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
            levsiz : int,
            m_psp : region(ispace(int2d), double),
            m_psn : region(ispace(int2d), double),
            julian : double,
            ozncyc : bool,
            paerlev : int,
            naer_c : int
            --abstot_3d : region(isapce(int3d), double)
          )
where 
  reads (phys_tbls)
do
  -----------------------------Local variables-----------------------------

  var lchnk : int = 10
  var ncol : int = 10
  var pcols : int = 10
  var pver : int = 10
  var pverp : int = 10
  var pverr : int = 10
  var pverrp : int = 10

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

  var n_cldadv : int = 3

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
  var m_psjp = region(ispace(int1d, constants.nCells), double) -- MATCH surface pressure
  fill(m_psjp, 0.0) -- TODO temp
  var m_psjn = region(ispace(int1d, constants.nCells), double) -- MATCH surface pressure
  fill(m_psjn, 0.0) -- TODO temp
  var swcftoa = region(ispace(int1d, constants.nCells), double)
  var lwcftoa = region(ispace(int1d, constants.nCells), double)
  var olrtoa = region(ispace(int1d, constants.nCells), double)

  var m_hybi = region(ispace(int1d, paerlev), double)
  var pin = region(ispace(int1d, levsiz), double)
  fill(pin, 0.0)

  var camrad_2d_r = region(ispace(int2d, {constants.nCells, constants.nVertLevels}), camrad_2d_fs)

  var ozmix = region(ispace(int2d, {constants.nCells, levsiz}), double)
  fill(ozmix, 0.0) -- TODO temp

  var q = region(ispace(int3d, {constants.nCells, constants.nVertLevels, n_cldadv}), double)

  var ozmixmj = region(ispace(int3d, {constants.nCells, levsiz, constants.nMonths}), double)
  fill(ozmixmj, 0.0) -- TODO temp

  var aerosoljp = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  fill(aerosoljp, 0.0) -- TODO temp
  var aerosoljn = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  fill(aerosoljn, 0.0) -- TODO temp

  var abstot = region(ispace(int3d, {constants.nCells, constants.nVertLevels, constants.nVertLevels}), double) -- Total absorptivity
  var absnxt = region(ispace(int3d, {constants.nCells, constants.nVertLevels, 4}), double) -- Total nearest layer absorptivity
  var emstot = region(ispace(int2d, {constants.nCells, constants.nVertLevels + 1}), double) -- Total emissivity

  fill(camrad_2d_r.pint, 0.0)
  fill(camrad_2d_r.pmid, 0.0)
  fill(camrad_2d_r.t, 0.0)

  -----------------------------

  format.println("Calling camrad...")

  -- TODO camrad computation missing

  var dolw : bool = true -- TEMP
  var dosw : bool = true -- TEMP
  radctl(
    cr, phys_tbls, camrad_1d_r, camrad_2d_r,
    lchnk, ncol, pcols, pver, pverp, pverr, pverrp, ppcnst, pcnst,
    q,
    julian,
    ozmixmj, ozmix, 
    levsiz, 
    pin, 
    ozncyc,
    m_psjp, m_psjn,
    aerosoljp, aerosoljn, 
    m_hybi, 
    paerlev,
    dolw, dosw,
    swcftoa, lwcftoa, olrtoa
  )

  format.println("Camrad done")
end