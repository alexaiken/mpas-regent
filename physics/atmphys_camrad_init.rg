import "regent"
require "data_structures"

local constants = require("constants")
local cmath = terralib.includec("math.h")

-- Math library imports
round = regentlib.round(double)


-- Purpose: 
-- Computes saturation vapor pressure over water and/or over ice using
-- Goff & Gratch (1946) relationships. 
-- <Say what the routine does> 
-- 
-- Method: 
-- T (temperature), and itype are input parameters, while es (saturation
-- vapor pressure) is an output parameter.  The input parameter itype
-- serves two purposes: a value of zero indicates that saturation vapor
-- pressures over water are to be returned (regardless of temperature),
-- while a value of one indicates that saturation vapor pressures over
-- ice should be returned when t is less than freezing degrees.  If itype
-- is negative, its absolute value is interpreted to define a temperature
-- transition region below freezing in which the returned
-- saturation vapor pressure is a weighted average of the respective ice
-- and water value.  That is, in the temperature range 0 => -itype
-- degrees c, the saturation vapor pressures are assumed to be a weighted
-- average of the vapor pressure over supercooled water and ice (all
-- water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
task gffgch(t : double,         -- in
            es : double,        -- inout 
            itype : double)     -- inout
  -- TODO
end

task radaeini(pstdx : double,     -- Standard pressure (dynes/cm^2)
              mwdryx : double,    -- Molecular weight of dry air 
              mwco2x : double)    -- Molecular weight of carbon dioxide

  -- ??? Some local variable setup? Unsure

  -- Set up table of H2O saturation vapor pressures for use in calculation effective path RH.  
  -- Need separate table from table in wv_saturation because:
  -- (1. Path temperatures can fall below minimum of that table; and
  -- (2. Abs/Emissivity tables are derived with RH for water only.

  -- TODO where do min_tp_h2o and max_tp_h2o come from?
  var tmin = round(min_tp_h2o)
  var tmax = round(max_tp_h2o) + 1
  var itype = 0
  for t = tmin, tmax, 1 do
    var tdbl = t
    gffgch(tdbl, estblh2o(t - tmin), itype) -- TODO is estblh2o an array?
  end
end

-- Initialize various constants for radiation scheme; note that
-- the radiation scheme uses cgs units.
task radini(gravx : double, 
            cpairx : double, 
            epsilox : double, 
            stebolx : double, 
            pstdx : double)
  radaeini(pstdx, constants.mwdry, constants.mwco2)
end

-- initialization of saturation vapor pressures:
task esinti(cr : region(ispace(int2d), cell_fs),
            er : region(ispace(int2d), edge_fs))
where reads (cr),
      writes (cr)
do
  -- TODO
end

-- initialization of ozone mixing ratios:
task oznini(cr : region(ispace(int2d), cell_fs),
            er : region(ispace(int2d), edge_fs))
where reads (cr),
      writes (cr)
do
  -- TODO
end

-- initialization of aerosol concentrations:
task aerosol_init(cr : region(ispace(int2d), cell_fs),
                  er : region(ispace(int2d), edge_fs))
where reads (cr),
      writes (cr)
do
  -- TODO
end

-- Counterpart of mpas_atmphys_camrad_init.F:camradinit(...)
task camradinit(cr : region(ispace(int2d), cell_fs),
                er : region(ispace(int2d), edge_fs))
where reads (cr),
      writes (cr)
do
  var pstd : double = 101325.0

  radini(constants.gravity, constants.cp, constants.ep_2, constants.stbolt, pstd*10.0)
  -- TODO correct arguments
  esinti(cr, er)
  oznini(cr, er)
  aerosol_init(cr, er)
end
