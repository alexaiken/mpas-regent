import "regent"
require "data_structures"

local constants = require("constants")
local c = regentlib.c

fabs = regentlib.fabs(double)
floor = regentlib.floor(double)
fmod = regentlib.fmod(double, double)


fspace doublefield {
  x         : double;
}

-- check sanity of time interpolation factors to within 32-bit roundoff
task validfactors(fact1 : double,
                  fact2 : double)
  var delta : double = 1.e-6
  if (fabs(fact1 + fact2 - 1.0) > delta or
      fact1 > 1.0 + delta or fact1 < -1 * delta or
      fact2 > 1.0 + delta or fact2 < -1 * delta) then
    return false
  end
  return true
end

-- Purpose: Determine time interpolation factors (normally for a boundary dataset)
--          for linear interpolation.
--
-- Method:  Assume 365 days per year.  Output variable fact1 will be the weight to
--          apply to data at calendar time "cdayminus", and fact2 the weight to apply
--          to data at time "cdayplus".  Combining these values will produce a result
--          valid at time "cday".  Output arguments fact1 and fact2 will be between
--          0 and 1, and fact1 + fact2 = 1 to roundoff.
task getfactors(cycflag : bool,
                np1 : int,
                cdayminus : double,
                cdayplus : double,
                cday : double)
  var fact1 : double
  var fact2 : double

  -- Determine time interpolation factors.  Account for December-January
  -- interpolation if dataset is being cycled yearly.
  if cycflag and np1 == 1 then
    var deltat : double = cdayplus + constants.daysperyear - cdayminus
    if cday > cdayplus then
      fact1 = (cdayplus + constants.daysperyear - cday) / deltat
      fact2 = (cday - cdayminus) / deltat
    else
      fact1 = (cdayplus - cday) / deltat
      fact2 = (cday + constants.daysperyear - cdayminus) / deltat
    end
  else
    var deltat : double = cdayplus - cdayminus
    fact1 = (cdayplus - cday) / deltat
    fact2 = (cday - cdayminus) / deltat
  end

  if validfactors(fact1, fact2) ~= true then
    c.printf("Bad fact1 and/or fact2=%.3f,%.3f", fact1, fact2)
  end

  return fact1
end

-- Sets ozmix
task oznint(julian : double,
            ozmixmj : region(ispace(int3d), doublefield),
            ozmix : region(ispace(int2d), doublefield),
            levsiz : int,
            num_months : int,
            pcols : int,
            ozncyc : bool)
where
  reads (ozmixmj),
  writes (ozmix)
do

  -- julian starts from 0.0 at 0Z on 1 Jan.
  var intjulian : double = julian + 1.0    -- offset by one day
  -- jan 1st 00z is julian = 1.0 here
  var ijul : int = floor(intjulian)
  -- Note that following will drift. 
  -- Need to use actual month/day info to compute julian.
  intjulian = intjulian - ijul
  ijul = fmod(ijul, constants.daysperyear)
  if (ijul == 0) then 
    ijul = constants.daysperyear 
  end
  intjulian = intjulian + ijul

  var date_oz = array(16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350)
  var np1 : int = 0
  var finddate : bool = false
  -- do m = 1, num_months
  for m = 0, num_months do
    if(date_oz[m] > intjulian and finddate == false) then
      np1 = m
      finddate = true
    end
  end

  var cdayozp : double = date_oz[np1]
  var cdayozm : double
  var np : int
  var nm : int
  if np1 > 0 then
    cdayozm = date_oz[np1 - 1]
    np = np1
    nm = np-1
  else
    cdayozm = date_oz[11]
    np = np1
    nm = 11
  end

  var fact1 : double = getfactors(ozncyc, np1, cdayozm, cdayozp, intjulian)
  var fact2 : double = 1.0 - fact1

  -- Time interpolation
  for k = 0, levsiz do 
    for i = 0, pcols do
      ozmix[{i, k}].x = ozmixmj[{i, k, nm}].x * fact1 + ozmixmj[{i, k, np}].x * fact2
    end
  end

end

task radozn()
end

task radinp()
end

task aqsat()
end

task get_rf_scales()
end

task get_aerosol()
end

task aerosol_indirect()
end

task radcswmx()
end

task get_int_scales()
end

task radclwmx()
end

task trcmix()
end