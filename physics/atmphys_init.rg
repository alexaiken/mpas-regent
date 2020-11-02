import "regent"
require "data_structures"

local constants = require("constants")

-- Dummy task for now, pay no mind to the arguments
task physics_init(cr : region(ispace(int2d), cell_fs),
                  er : region(ispace(int2d), edge_fs))
where reads (cr),
      writes (cr)
do
  -- TODO
end
