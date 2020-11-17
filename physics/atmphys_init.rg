import "regent"
require "data_structures"
require "physics/atmphys_camrad_init"

-- Counterpart of mpas_atmphys_init.F:physics_init(...)
-- Currently only implementing radiation
-- Radiation work in progress
-- TODO correct arguments
task physics_init(cr : region(ispace(int2d), cell_fs),
                  er : region(ispace(int2d), edge_fs))
where reads (cr),
      writes (cr)
do

  -- SKIPPING LOTS OF mpas_atmphys_init.F:physics_init(...)

  -- Assuming we're always in the "cam_sw" case
  -- -----
  -- If we're in the "cam_sw" case, then "init_radiation_lw" 
  -- and "init_radiation_sw" are identical, and only need to be
  -- called once, so we can directly call camradinit() 
  camradinit(cr, er)
end
