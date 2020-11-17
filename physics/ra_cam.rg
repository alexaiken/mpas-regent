import "regent"
require "data_structures"
require "physics/ra_cam_cld_support"
require "physics/ra_cam_radctl_support"

task param_cldoptics_calc()
  cldefr()
  cldems()
  cldovrlap()
end

task radctl()

  -- passing temp variables for now
  var julian : double = 1.0
  var levsiz : int = 1
  var num_months : int = 1
  var pcols : int = 1
  var ozncyc : bool = true
  var ozmixmj = region(ispace(int3d, {pcols, levsiz, num_months}), doublefield)
  var ozmix = region(ispace(int2d, {pcols, levsiz}), doublefield) -- ozone mixing ratio
  oznint(julian, ozmixmj, ozmix, levsiz, num_months, pcols, ozncyc)

  -- passing temp variables for now
  var ncol : int = 1     -- number of atmospheric columns
  var pver : int = 1
  var pmid = region(ispace(int2d, {pcols, pver}), doublefield)    -- level pressures (mks)
  var pin = region(ispace(int1d, levsiz), doublefield)            -- ozone data level pressures (mks)
  var o3vmr = region(ispace(int2d, {pcols, pver}), doublefield)
  radozn(ncol, pcols, pver, pmid, pin, levsiz, ozmix, o3vmr)

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

task camrad()
  param_cldoptics_calc()
  radctl()
end