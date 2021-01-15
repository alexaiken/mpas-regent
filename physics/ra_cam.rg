import "regent"
require "data_structures"
require "physics/ra_cam_cld_support"
require "physics/ra_cam_radctl_support"

local constants = require("constants")

task param_cldoptics_calc()
  cldefr()
  cldems()
  cldovrlap()
end

task radctl(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
            ncol : int,
            pcols : int,
            pver : int, pverp : int, pverr : int, pverrp : int,
            julian : double,
            ozmixmj : region(ispace(int3d), double),
            ozmix : region(ispace(int2d), double),
            levsiz : double,
            pin : region(ispace(int1d), double),
            ozncyc : bool)
where reads (cr.{pmid, pint, t}, phys_tbls, ozmixmj, ozmix, pin),
      writes (ozmix)
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
  -------------------------------------------------------------------------

  oznint(julian, ozmixmj, ozmix, levsiz, pcols, ozncyc)

  radozn(cr, radctl_2d_pverr_r, ncol, pcols, pver, pin, levsiz, ozmix)

  radinp(cr, radctl_2d_pverr_r, radctl_2d_pverrp_r, ncol, pver, pverp)

  aqsat(cr, phys_tbls, radctl_2d_pverr_r, ncol, pver)

  get_rf_scales()
  get_aerosol()
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
            ozncyc : bool)
where 
  reads (cr.{pmid, pint, t}, phys_tbls)
do
  -----------------------------Local variables-----------------------------

  var lchnk : int
  var ncol : int
  var pcols : int
  var pver : int
  var pverp : int
  var pverr : int
  var pverrp : int

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
  var ozmix = region(ispace(int2d, {constants.nCells, levsiz}), double)
  var pin = region(ispace(int1d, levsiz), double)

  -----------------------------

  param_cldoptics_calc()

  radctl(cr, phys_tbls, ncol, pcols, pver, pverp, pverr, pverrp, julian,
         ozmixmj, ozmix, levsiz, pin, ozncyc)

end