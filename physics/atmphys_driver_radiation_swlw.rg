import "regent"
require "data_structures"
require "physics/ra_cam"

task radconst()
end

----------------
-- SHORT WAVE --
----------------

task allocate_radiation_sw()
end

task deallocate_radiation_sw()
end

task radiation_sw_from_MPAS()
end

task radiation_sw_to_MPAS()
end

task driver_radiation_sw()
  radiation_sw_from_MPAS()
  radconst()
  camrad()
  radiation_sw_to_MPAS()
end

---------------
-- LONG WAVE --
---------------

task allocate_radiation_lw()
end

task deallocate_radiation_lw()
end

task vinterp_ozn()
end

task radiation_lw_from_MPAS()
  vinterp_ozn()
end

task radiation_lw_to_MPAS()
end

task driver_radiation_lw()
  radiation_lw_from_MPAS()
  radconst()
  camrad()
  radiation_lw_to_MPAS()
end