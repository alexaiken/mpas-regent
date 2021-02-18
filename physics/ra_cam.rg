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

task radctl(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
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
            m_psp : region(ispace(int1d), double),
            m_psn : region(ispace(int1d), double),
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

  var in2o : int
  var ich4 : int
  var if11 : int
  var if12 : int

  var radctl_1d_r = region(ispace(int1d, pcols), radctl_1d_fs)
  var radctl_2d_pver_r = region(ispace(int2d, {pcols, pver}), radctl_2d_pver_fs)
  var radctl_2d_pverr_r = region(ispace(int2d, {pcols, pverr}), radctl_2d_pverr_fs)
  var radctl_2d_pverrp_r = region(ispace(int2d, {pcols, pverrp}), radctl_2d_pverrp_fs)
  
  var aerosol = region(ispace(int3d, {pcols, pver, constants.naer_all}), double) -- aerosol mass mixing ratios
  var scales = region(ispace(int1d, constants.naer_all), double)  -- scaling factors for aerosols
  -------------------------------------------------------------------------

  oznint(julian, ozmixmj, ozmix, levsiz, pcols, ozncyc)

  radozn(cr, radctl_2d_pverr_r, ncol, pcols, pver, pin, levsiz, ozmix)

  radinp(cr, radctl_2d_pverr_r, radctl_2d_pverrp_r, ncol, pver, pverp)

  aqsat(cr, phys_tbls, radctl_2d_pverr_r, ncol, pver)

  get_aerosol(cr, phys_tbls, lchnk, julian, m_psp, m_psn, aerosoljp, aerosoljn, m_hybi, paerlev, 
              pcols, pver, pverp, pverr, pverrp, aerosol, scales)

  aerosol_indirect()
  radcswmx()
  get_int_scales()
  radclwmx()
  trcmix()
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

  var ozmixmj = region(ispace(int3d, {constants.nCells, levsiz, constants.nMonths}), double)
  fill(ozmixmj, 0.0); -- TODO: TEMP
  
  var ozmix = region(ispace(int2d, {constants.nCells, levsiz}), double)
  fill(ozmix, 0.0); -- TODO: TEMP

  var pin = region(ispace(int1d, levsiz), double)
  fill(pin, 0.0); -- TODO: TEMP

  fill(cr.pmid, 0.0); -- TODO: TEMP

  var m_psjp = region(ispace(int1d, constants.nCells), double)     -- MATCH surface pressure
  var m_psjn = region(ispace(int1d, constants.nCells), double)     -- MATCH surface pressure
  fill(m_psjp, 0.0); -- TODO: TEMP
  fill(m_psjn, 0.0); -- TODO: TEMP

  var aerosoljp = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  var aerosoljn = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  fill(aerosoljp, 0.0); -- TODO: TEMP
  fill(aerosoljn, 0.0); -- TODO: TEMP

  var m_hybi = region(ispace(int1d, paerlev), double)
  fill(m_hybi, 0.0); -- TODO: TEMP

  -----------------------------

  format.println("Calling camrad...")

  param_cldoptics_calc()

  radctl(
    cr, 
    phys_tbls, 
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
    m_psjp, m_psjn,
    aerosoljp, aerosoljn, 
    m_hybi, 
    paerlev
  )

  format.println("Camrad done")

end