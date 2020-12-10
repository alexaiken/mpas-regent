import "regent"
require "data_structures"
require "physics/ra_cam_cld_support"
require "physics/ra_cam_radctl_support"

local constants = require("constants")

fspace ozone_mix_fs {
  ozmixmj : double,     -- monthly ozone mixing ratio
  ozmix : double,                             -- ozone mixing ratio (time interpolated)
}

fspace aerosol_fs {
  aerosoljp : double, -- monthly aerosol concentrations
  aerosoljn : double, -- monthly aerosol concentrations
}

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
  oznint(julian, ozmixmj, ozmix, levsiz, num_months, pcols, ozncyc)

  -- passing temp variables for now
  var ncol : int = 1     -- number of atmospheric columns
  var pver : int = 1
  var pmid = region(ispace(int2d, {pcols, pver}), double)    -- level pressures (mks)
  for i=0, pcols do
    for j=0, levsiz do
      pmid[{i, j}] = .5
    end
  end
  var pin = region(ispace(int1d, levsiz), double)            -- ozone data level pressures (mks)
  for i=0, pcols do
    pin[i] = .5
  end
  var o3vmr = region(ispace(int2d, {pcols, pver}), double)
  for i=0, pcols do
    for j=0, levsiz do
      o3vmr[{i, j}] = .5
    end
  end
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

  var n_cldadv : int
  var paerlev : int
  var levsiz : int
  var naer_c : int

  var nCells : int = cr.bounds.lo.x

  var q = region(ispace(int3d, {nCells, constants.nVertLevels, n_cldadv}), double)
  var ozone_mix_locals = region(ispace(int2d, {nCells, levsiz}), ozone_mix_fs)
  var pin = region(ispace(int1d, levsiz), double)
  var aerosol_locals = region(ispace(int3d, {nCells, paerlev, naer_c}), aerosol_fs)
  var m_hybi = region(ispace(int1d, paerlev), double)
  var abstot = region(ispace(int3d, {nCells, constants.nVertLevels+1, constants.nVertLevels+1}), double)   -- Total absorptivity
  var absnxt = region(ispace(int3d, {nCells, constants.nVertLevels, 4}), double)             -- Total nearest layer absorptivity
  var emstot = region(ispace(int2d, {nCells, constants.nVertLevels+1}), double)              -- Total emissivity

  -- END LOCALS --
  ----------------

  param_cldoptics_calc()
  radctl()
end