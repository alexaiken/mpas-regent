import "regent"
require "data_structures"
require "physics/atmphys_driver_cloudiness"
require "physics/atmphys_driver_radiation_swlw"

local format = require("std/format")

task physics_timetracker()
end

task allocate_forall_physics()
end

task deallocate_forall_physics()
end

task MPAS_to_physics()
end

task update_radiation_diagnostics()
end

task physics_driver(cr : region(ispace(int2d), cell_fs),
                    phys_tbls : region(ispace(int1d), phys_tbls_fs))
where 
  reads (phys_tbls),
  reads writes (cr)
do
  
  format.println("Calling physics_driver...")

  allocate_forall_physics()
  MPAS_to_physics()

  allocate_cloudiness()
  driver_cloudiness()

  allocate_radiation_sw()
  driver_radiation_sw(cr, phys_tbls)

  allocate_radiation_lw()

  -- TODO these are temp: find out what they actually are
  var radt_lw_scheme : regentlib.string = "cam_lw"
  var config_o3climatology : bool = true
  var microp_scheme : regentlib.string = "DUMMY STRING"
  var config_microp_re : bool = true
  driver_radiation_lw(
    cr, phys_tbls, radt_lw_scheme, config_o3climatology, microp_scheme,
    config_microp_re
  )

  update_radiation_diagnostics()

  deallocate_cloudiness()
  deallocate_radiation_sw()
  deallocate_radiation_lw()
  deallocate_forall_physics()
end