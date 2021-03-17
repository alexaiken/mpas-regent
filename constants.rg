import "regent"

local constants = {}
constants.my_const = 123

constants.c = regentlib.c
constants.cio = terralib.includec("stdio.h")
constants.clib = terralib.includec("stdlib.h")

constants.netcdf = terralib.includec("/home/arjunk1/spack/opt/spack/linux-ubuntu20.04-broadwell/gcc-9.3.0/netcdf-c-4.7.4-h7i6kmblkfnyttdnctplncjm4fpzaxqz/include/netcdf.h")

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
constants.nVertLevels = 5
constants.sphere_radius = terralib.constant(double, 6371229.0) --from ncdump of the grid, default is 1.0, but set to "a" from constants.F in init_atm_core.F
constants.nlat = 721
constants.seconds_per_day = terralib.constant(double, 86400.0)

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
constants.nScalars = 8

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

-- Physics - Radiation
constants.R_d = constants.rgas              -- gas constant for dry air [J kg-1 K-1]
constants.R_v = 461.6                       -- gas constant for water vapor
constants.ep_2 = terralib.constant(double, `(-constants.R_d / constants.R_v)) --Ratio of mol. wght of H2O to dry air
constants.stbolt = 5.67051e-8
constants.mwdry = 28.966                                        -- molecular weight dry air ~ kg/kmole (shr_const_mwdair)
constants.mwco2 = terralib.constant(double, 44.0)               -- molecular weight co2
constants.tmelt = 273.16                                        -- freezing T of fresh water ~ K
constants.daysperyear = 365
constants.solcon_0 = terralib.constant(double, 1370.0)          -- solar constant [W/m2]
constants.nMonths = 12
constants.cam_abs_dim1 = 4
constants.cam_abs_dim2 = constants.nVertLevels + 1
constants.amd = 28.9644                                         -- Molecular weight of dry air (kg/kmol)
constants.amo = terralib.constant(double, 48.0000)              -- Molecular weight of ozone (g/mol)

constants.nAerLevels = 29
constants.nOznLevels = 59
constants.nAerosols = 12
constants.naer = 10
constants.naer_all = 12

-- Physics - Radiation (parameters)
constants.min_tp_h2o = terralib.constant(double, 160.0)        -- min T_p for pre-calculated abs/emis
constants.max_tp_h2o = 349.999999   -- max T_p for pre-calculated abs/emis
constants.ntemp = 192               -- Number of temperatures in H2O sat. table for Tp
constants.plenest = 250             -- length of saturation vapor pressure table

constants.config_dt = terralib.constant(double, 720.0)          -- Model time step, seconds

-- mpas_atmphys_constants.F -- 

constants.c0 = terralib.constant(double, 0.0)
constants.c1 = terralib.constant(double, 1.0)

constants.P0 = terralib.constant(double, 100000.0)                  -- Reference pressure
constants.t00 = 273.15                                              -- Reference temperature
constants.R_v = 461.6                                               -- Gas constant for water vapor
constants.ep_1 = terralib.constant(double, `(constants.R_v / constants.R_d - 1))
constants.ep_2 = terralib.constant(double, `(constants.R_d / constants.R_v))
constants.cpv = terralib.constant(double, `(4.0 * constants.R_v))
constants.rdg = terralib.constant(double, `(constants.R_d / constants.gravity))
constants.rcp = terralib.constant(double, `(constants.R_d / constants.cp))
constants.rcv = terralib.constant(double, `(constants.R_d / (constants.cp - constants.R_d)))

constants.rho_a = 1.28
constants.rho_r = terralib.constant(double, 1000.0)
constants.rho_s = terralib.constant(double, 100.0)
constants.rho_w = terralib.constant(double, 1000.0)

constants.svp1 = 0.6112
constants.svp2 = 17.67
constants.svp3 = 29.65
constants.svpt0 = 273.15

constants.xlv = 2.50e6                              -- latent heat of vaporization
constants.xlf = 3.50e5                              -- latent heat of fusion
constants.xls = terralib.constant(double, `(constants.xlv + constants.xlf))       -- latent heat of sublimation

constants.xlv0 = 3.15e6        
constants.xlv1 = 2370.
constants.xls0 = 2.905e6
constants.xls1 = 259.532

constants.karman = 0.4                              -- Von Karman constant
constants.eomeg = 7.29210e-5
constants.stbolt = 5.67051e-8

constants.cliq    = terralib.constant(double, 4190.0)
constants.cice    = terralib.constant(double, 2106.0)
constants.epsilon = 1.e-15
constants.psat    = 610.78

-- constants specific to long- and short-wave radiation codes:
constants.solcon_0 = terralib.constant(double, 1370.0)              -- solar constant                                     [W/m2]
constants.degrad   = terralib.constant(double, 3.1415926 / 180.)    -- conversion from degree to radiant                     [-]
constants.dpd      = terralib.constant(double, 360.0 / 365.0)

return constants
