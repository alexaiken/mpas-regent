import "regent"

local constants = require("constants")

fspace camrad_1d_fs {
  -- Arrays of size (ite-its)
  coszrs : double,
  landfrac : double,
  landm : double,
  snowh : double,
  icefrac : double,
  lwups : double,

  asdir : double,
  asdif : double,
  aldir : double,
  aldif : double,
  ps : double,

  nmxrgn : int,           -- Number of maximally overlapped regions

  fsns : double,          -- Surface absorbed solar flux
  fsnt : double,          -- Net column abs solar flux at model top
  flns : double,          -- Srf longwave cooling (up-down) flux
  flnt : double,          -- Net outgoing lw flux at model top

  sols : double,          -- Downward solar rad onto surface (sw direct)
  soll : double,          -- Downward solar rad onto surface (lw direct)
  solsd : double,         -- Downward solar rad onto surface (sw diffuse)
  solld : double,         -- Downward solar rad onto surface (lw diffuse)

  fsds : double,          -- Flux Shortwave Downwelling Surface
  flwds : double,         -- Surface down longwave flux

  clat : double,          -- latitude in radians for columns
}

fspace camrad_2d_fs {
  -- Arrays of size (ite-its, kte-kts) aka (nCells, nVertLevels)
  -- NOTE: some are kte-ktw+1
  cld : double,
  pmid : double,
  lnpmid : double,
  pdel : double,
  zm : double,
  t : double,
  pint : double,
  lnpint : double,
  
  cicewp : double,        -- in-cloud cloud ice water path
  cliqwp : double,        -- in-cloud cloud liquid water path
  tauxcl : double,        -- cloud water optical depth
  tauxci : double,        -- loud ice optical depth
  emis : double,          -- cloud emissivity
  rel : double,           -- effective drop radius (microns)
  rei : double,           -- ice effective drop size (microns)
  pmxrgn : double,        -- Maximum values of pressure for each

  fsup : double,          -- Upward total sky solar
  fsupc : double,         -- Upward clear sky solar
  fsdn : double,          -- Downward total sky solar
  fsdnc : double,         -- Downward clear sky solar
  flup : double,          -- Upward total sky longwave
  flupc : double,         -- Upward clear sky longwave
  fldn : double,          -- Downward total sky longwave
  fldnc : double,         -- Downward clear sky longwave

  qrs : double,           -- Solar heating rate
  qrl : double,           -- Longwave cooling rate
}

fspace radctl_1d_fs {
  -- Arrays of size (pcols) --
  solin : double,         -- Solar incident flux
  fsntoa : double,        -- Net solar flux at TOA
  fsntoac : double,       -- Clear sky net solar flux at TOA
  fsnirt : double,        -- Near-IR flux absorbed at toa
  fsnrtc : double,        -- Clear sky near-IR flux absorbed at toa
  fsnirtsq : double,      -- Near-IR flux absorbed at toa >= 0.7 microns
  fsntc : double,         -- Clear sky total column abs solar flux
  fsnsc : double,         -- Clear sky surface abs solar flux
  fsdsc : double,         -- Clear sky surface downwelling solar flux
  flutc : double,         -- Upward Clear Sky flux at top of model
  flntc : double,         -- Clear sky lw flux at model top
  flnsc : double,         -- Clear sky lw flux at srf (up-down)

  lwupcgs : double,       -- Upward longwave flux in cgs units

  frc_day : double,       -- = 1 for daylight, =0 for night colums
}

fspace radctl_2d_pver_fs {
  -- Arrays of size (pcols, pver) --
  ftem : double,          -- temporary array for outfld
  
  n2o : double,      -- nitrous oxide mass mixing ratio
  ch4 : double,      -- methane mass mixing ratio
  cfc11 : double,    -- cfc11 mass mixing ratio
  cfc12 : double,    -- cfc12 mass mixing ratio
}

fspace radctl_2d_pverr_fs {
  -- Arrays of size (pcols, pverr) --
  pbr : double,     -- Model mid-level pressures (dynes/cm2)

  o3vmr : double,   -- Ozone volume mixing ratio
  o3mmr : double,   -- Ozone mass mixing ratio
  rh : double,      -- level relative humidity (fraction)

  esat : double,    -- saturation vapor pressure
  qsat : double,    -- saturation specific humidity
}

fspace radctl_3d_aero_fs {
  -- Arrays of size (pcols, nspint, naer_groups)
  aertau : double, -- Aerosol column optical depth
  aerssa : double, -- Aerosol column averaged single scattering albedo
  aerasm : double, -- Aerosol column averaged asymmetry parameter
  aerfwd : double, -- Aerosol column averaged forward scattering
}
