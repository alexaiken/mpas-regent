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

task radctl(radctl_args : radctl_args_fs,
            pin : region(ispace(int1d, levsiz), double))

  -----------------------------Local variables-----------------------------

  var i : int,
  k : int,

  in2o : int,
  ich4 : int,
  if11 : int,
  if12 : int,

  eccf : double,          -- Earth/sun distance factor





  -- passing temp variables for now
  var ozmixmj = region(ispace(int3d, {pcols, levsiz, num_months}), double)
  for i=0, pcols do
    for j=0, levsiz do
      for k=0, num_months do
        ozmixmj[{i, j, k}] = .5
      end
    end
  end
  var ozmix = region(ispace(int2d, {pcols, levsiz}), double) -- ozone mixing ratio
  for i=0, pcols do
    for j=0, levsiz do
      ozmix[{i, j}] = .5
    end
  end
  oznint(radctl_args.julian, 
         ozmixmj, 
         ozmix, 
         radctl_args.levsiz, 
         radctl_args.num_months, 
         radctl_args.pcols, 
         radctl_args.ozncyc)

  -- passing temp variables for now
  var pmid = region(ispace(int2d, {pcols, pver}), double)    -- level pressures (mks)
  for i=0, pcols do
    for j=0, levsiz do
      pmid[{i, j}] = .5
    end
  end
  var o3vmr = region(ispace(int2d, {pcols, pver}), double)
  for i=0, pcols do
    for j=0, levsiz do
      o3vmr[{i, j}] = .5
    end
  end
  radozn(radctl_args.ncol, 
         radctl_args.pcols, 
         radctl_args.pver, 
         pmid, 
         pin, 
         radctl_args.levsiz, 
         ozmix, 
         o3vmr)

  radinp()
  aqsat()
  get_rf_scales()
  get_aerosol()
  aerosol_indirect()
  radcswmx()
  get_int_scales()
  radclwmx()
  trcmix()
end

task camrad(cr : region(ispace(int2d), cell_fs))

  ------------------
  -- START LOCALS --

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

  -- END LOCALS --
  ----------------

  param_cldoptics_calc()

  var radctl_args = radctl_args_fs {}
  radctl(radctl_args)
end