import "regent"

local constants = {}
constants.my_const = 123

constants.c = regentlib.c
constants.cio = terralib.includec("stdio.h")
constants.clib = terralib.includec("stdlib.h")

constants.netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

constants.FILE_NAME = "mesh_loading/x1.2562.grid.nc"
constants.GRAPH_FILE_NAME = "mesh_loading/x1.2562.graph.info.part.16"
constants.NUM_PARTITIONS = 16
constants.MAXCHAR = 5
constants.NUM_TIMESTEPS = 10

constants.nCells = 2562
constants.nEdges = 7680
constants.nVertices = 5120
constants.maxEdges = 10
constants.maxEdges2 = 20
constants.TWO = 2
constants.FIFTEEN = 15
constants.vertexDegree = 3
constants.nVertLevels = 1
constants.sphere_radius = 6371229.0 --from ncdump of the grid, default is 1.0, but set to "a" from constants.F in init_atm_core.F
constants.nlat = 721

constants.pii = 3.141592653589793
constants.omega = 7.29212E-5
constants.rgas = terralib.constant(double, 287.0)
constants.rv = 461.6
constants.cp = terralib.constant(double, `(7.0*constants.rgas/2.0))
constants.gravity = 9.80616
constants.rvord = terralib.constant(double, `(constants.rv/constants.rgas))
constants.config_epssm = 0.1

constants.nRelaxZone = 5

constants.config_smdiv = 0.1
constants.config_len_disp = terralib.constant(double, 120000.0)








return constants
