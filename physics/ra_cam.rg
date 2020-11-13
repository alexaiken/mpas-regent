import "regent"
require "data_structures"
require "ra_cam_cld_support"
require "ra_cam_radctl_support"

task param_cldoptics_calc()
  cldefr()
  cldems()
  cldovrlap()
end

task radctl()
  oznint()
  radozn()
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