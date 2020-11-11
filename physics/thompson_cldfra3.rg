import "regent"
require "data_structures"

task adjust_cloudIce()
end

task adjust_cloudH2O()
end

task adjust_cloudFinal()
end

task find_cloudLayers()
  adjust_cloudIce()
  adjust_cloudH2O()
  adjust_cloudFinal()
end

task calc_cldfra3()
  find_cloudLayers()
end