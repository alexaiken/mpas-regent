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
constants.sphere_radius = terralib.constant(double, 6371229.0) --from ncdump of the grid, default is 1.0, but set to "a" from constants.F in init_atm_core.F
constants.nlat = 721

constants.pii = 3.141592653589793
constants.omega = 7.29212E-5
constants.rgas = terralib.constant(double, 287.0)
constants.rv = 461.6
constants.cp = terralib.constant(double, `(7.0*constants.rgas/2.0))
constants.cv = terralib.constant(double, `(constants.cp - constants.rgas))
constants.cvpm = terralib.constant(double, `(-constants.cv / constants.cp))
constants.gravity = 9.80616
constants.rvord = terralib.constant(double, `(constants.rv/constants.rgas))
constants.config_epssm = 0.1
constants.prandtl = terralib.constant(double, 1.0)

constants.nRelaxZone = 5

constants.config_smdiv = 0.1
constants.config_len_disp = terralib.constant(double, 120000.0)

constants.config_v_mom_eddy_visc2 = terralib.constant(double, 0.0)
constants.config_h_theta_eddy_visc2 = terralib.constant(double, 0.0)
constants.config_v_theta_eddy_visc2 = terralib.constant(double, 0.0)
constants.config_h_mom_eddy_visc4 = terralib.constant(double, 0.0)
constants.config_h_theta_eddy_visc4 = terralib.constant(double, 0.0)
constants.config_visc4_2dsmag = 0.05
constants.config_smagorinsky_coef = 0.125
constants.config_del4u_div_factor = terralib.constant(double, 10.0)
constants.config_number_rayleigh_damp_u_levels = 6
constants.config_rayleigh_damp_u_timescale_days = terralib.constant(double, 5.0)
constants.config_coef_3rd_order = 0.25
constants.config_dynamics_split_steps = 1 --Default value is 3, but we will temporarily use 1 for simplicity

--Not sure about these as they seem to be parameters:
constants.config_horiz_mixing = "2d_smagorinsky"
constants.config_mpas_cam_coef = terralib.constant(double, 0.0)
constants.config_mix_full = false
constants.config_rayleigh_damp_u = false
constants.config_print_detailed_minmax_vel = false
constants.config_print_global_minmax_vel = false
constants.config_print_global_minmax_sca = false




return constants
