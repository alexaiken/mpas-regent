import "regent"
require "data_structures"
require "physics/atmphys_driver_cloudiness"
require "physics/atmphys_driver_radiation_swlw"

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

task physics_driver(cr : region(ispace(int2d), cell_fs))
  allocate_forall_physics()
  MPAS_to_physics()

  allocate_cloudiness()
  driver_cloudiness()

  allocate_radiation_sw()
  --driver_radiation_sw(cr)

  allocate_radiation_lw()
  --driver_radiation_lw(cr)

  update_radiation_diagnostics()

  deallocate_cloudiness()
  deallocate_radiation_sw()
  deallocate_radiation_lw()
  deallocate_forall_physics()
end