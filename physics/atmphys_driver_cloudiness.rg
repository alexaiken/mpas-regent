import "regent"
require "data_structures"
require "physics/thompson_cldfra3"

task allocate_cloudiness()
end

task deallocate_cloudiness()
end

task cloudiness_from_MPAS()
end

task cloudiness_to_MPAS()
end

task calc_cldincidence()
end

task calc_cldfraction()
end

task driver_cloudiness()
  cloudiness_from_MPAS()
  
  -- temp
  var i = 0
  if i == 0 then
    calc_cldincidence()
  elseif i == 1 then
    calc_cldfraction()
  else
    calc_cldfra3()
  end

  cloudiness_to_MPAS()
end