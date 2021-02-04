

task atm_compute_signs(cr : region(ispace(int2d), cell_fs),
                       er : region(ispace(int2d), edge_fs),
                       vr : region(ispace(int2d), vertex_fs))
where 
  reads (cr.{edgesOnCell, nEdgesOnCell, verticesOnCell}, er.{cellsOnEdge, verticesOnEdge, zb, zb3},
         vr.{cellsOnVertex, edgesOnVertex}), 
  writes (cr.{edgesOnCellSign, kiteForCell, zb_cell, zb3_cell}, vr.{edgesOnVertexSign})
do

task atm_adv_coef_compression(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs))
where
  reads (cr.{cellsOnCell, nEdgesOnCell}, er.{cellsOnEdge, dcEdge, deriv_two, dvEdge}),
  writes (er.{advCellsForEdge}),
  reads writes (er.{adv_coefs, adv_coefs_3rd, nAdvCellsForEdge})
do

task atm_compute_damping_coefs(config_zd : double,
                               config_xnutr : double,
                               cr : region(ispace(int2d), cell_fs))
where
  reads (cr.{meshDensity, zgrid}),
  reads writes (cr.{dss})
do

task atm_couple_coef_3rd_order(config_coef_3rd_order : double,
                               cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs))
where
  reads writes (cr.zb3_cell, er.adv_coefs_3rd)
do

task atm_compute_solve_diagnostics(cr : region(ispace(int2d), cell_fs),
                                   er : region(ispace(int2d), edge_fs),
                                   vr : region(ispace(int2d), vertex_fs),
                                   hollingsworth : bool,
                                   rk_step : int)
where
  reads (cr.{edgesOnCell, edgesOnCellSign, h, invAreaCell, kiteForCell, nEdgesOnCell, verticesOnCell}, er.{cellsOnEdge, dcEdge, dvEdge, edgesOnEdge_ECP, nEdgesOnEdge, u, verticesOnEdge, weightsOnEdge}, vr.{edgesOnVertex, edgesOnVertexSign, fVertex, invAreaTriangle, kiteAreasOnVertex}),
  writes (er.{h_edge, pv_edge}),
  reads writes (cr.{divergence, ke}, er.{ke_edge, v}, vr.{ke_vertex, pv_vertex, vorticity})
do

task atm_compute_moist_coefficients(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs))
where
  reads (er.{cellsOnEdge}),
  writes (cr.{cqw, er.cqu}),
  reads writes (cr.{qtot})
do

task atm_compute_vert_imp_coefs(cr : region(ispace(int2d), cell_fs),
                                vert_r : region(ispace(int1d), vertical_fs),
                                dts : double)
where
  reads (cr.{cqw, exner, exner_base, qtot, rho_base, rtheta_base, rtheta_p, theta_m, zz},
         vert_r.{rdzu, rdzw, fzm, fzp}),
  reads writes (cr.{a_tri, alpha_tri, b_tri, c_tri, coftz, cofwr, cofwt, cofwz, gamma_tri}, vert_r.{cofrz})
do

task atm_compute_mesh_scaling(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs),
                              config_h_ScaleWithMesh : bool)
where
  reads (cr.{meshDensity}, er.{cellsOnEdge}),
  writes (cr.{meshScalingRegionalCell}, er.{meshScalingDel2, meshScalingDel4, meshScalingRegionalEdge})
do


task atm_init_coupled_diagnostics(cr : region(ispace(int2d), cell_fs),
                                  er : region(ispace(int2d), edge_fs),
                                  vert_r : region(ispace(int1d), vertical_fs))
where
  reads (cr.{edgesOnCell, edgesOnCellSign, nEdgesOnCell, rho_base, theta, theta_base, theta_m, w,
             zb_cell, zb3_cell, zz},
         er.{cellsOnEdge, u},
         vert_r.{fzm, fzp}),
  writes (cr.{pressure_base, pressure_p}),
  reads writes (cr.{exner, exner_base, rho_p, rho_zz, rtheta_base, rtheta_p, rw, theta_m},
                er.{ru})
do

task atm_compute_output_diagnostics(cr : region(ispace(int2d), cell_fs))
where
  reads (cr.{pressure_base, pressure_p, rho_zz, theta_m, zz}),
  writes (cr.{pressure, rho, theta})
do

task atm_rk_integration_setup(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs))
where
  reads (cr.{rho_p, rho_zz, rtheta_p, rw, theta_m, w}, er.{ru, u}),
  writes (cr.{rho_p_save, rho_zz_2, rho_zz_old_split, rtheta_p_save, rw_save, theta_m_2, w_2}, er.{ru_save, u_2})
do


task atm_compute_dyn_tend_work(cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs),
                               vr : region(ispace(int2d), vertex_fs),
                               vert_r : region(ispace(int1d), vertical_fs),
                               rk_step : int,
                               dt : double,
                               config_horiz_mixing : regentlib.string,
                               config_mpas_cam_coef : double,
                               config_mix_full : bool,
                               config_rayleigh_damp_u : bool)
where
  reads (cr.{cqw, defc_a, defc_b, divergence, edgesOnCell, edgesOnCell_sign, invAreaCell, ke, lat,
             nEdgesOnCell, pressure_p, qtot, rho_base, rho_zz, rho_p_save, rt_diabatic_tend, rw, rw_save, t_init, tend_rho_physics,
             tend_rtheta_physics, theta_m, theta_m_save, uReconstructZonal, uReconstructMeridional, w, zgrid, zz},
         er.{advCellsForEdge, adv_coefs, adv_coefs_3rd, angleEdge, cellsOnEdge, cqu, dcEdge, dvEdge,
             edgesOnEdge, invDcEdge, invDvEdge, lat, meshScalingDel2, meshScalingDel4, nAdvCellsForEdge,
             nEdgesOnEdge, pv_edge, rho_edge, ru, ru_save, tend_ru_physics, u, v, verticesOnEdge,
             weightsOnEdge, zxu},
         vr.{edgesOnVertex, edgesOnVertex_sign, invAreaTriangle, vorticity},
         vert_r.{fzm, fzp, rdzu, rdzw, u_init, v_init}),
  writes (cr.{rthdynten, tend_rho, tend_rtheta_adv}),
  reads writes (cr.{delsq_divergence, delsq_theta, delsq_w, dpdz, flux_arr, h_divergence, kdiff, ru_edge_w,
                    tend_theta, tend_theta_euler, w, tend_w_euler, wdtz, wdwz},
                er.{delsq_u, q, tend_u, tend_u_euler, u_mix, wduz},
                vr.{delsq_vorticity})
do



task atm_set_smlstep_pert_variables_work(cr : region(ispace(int2d), cell_fs),
                                         er : region(ispace(int2d), edge_fs),
                                         vert_r : region(ispace(int1d), vertical_fs))
where
  reads (cr.{bdyMaskCell, edgesOnCell, edgesOnCell_sign, nEdgesOnCell, zb_cell, zb3_cell, zz},
         er.{u_tend}, vert_r.{fzm, fzp}),
  reads writes (cr.{w})
do



task atm_divergence_damping_3d(cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs),
                               dts : double)
where
  reads (cr.{rtheta_pp, rtheta_pp_old, theta_m}, er.{cellsOnEdge, specZoneMaskEdge}),
  reads writes (er.{ru_p})
do


task atm_recover_large_step_variables_work(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    vert_r : region(ispace(int1d), vertical_fs),
                                    ns : int,
                                    rk_step : int,
                                    dt : double)
where
  reads (cr.{bdyMaskCell, edgesOnCell, edgesOnCell_sign, exner_base, nEdgesOnCell, rho_base, rho_p_save, rho_pp,
             rt_diabatic_tend, rtheta_base, rtheta_p_save, rtheta_pp, rw_p, rw_save, zb_cell, zb3_cell, zz},
         er.{cellsOnEdge, ru_p, ru_save},
         vert_r.{cf1, cf2, cf3, fzm, fzp}),
  writes (cr.{pressure_p, theta_m}, er.u),
  reads writes (cr.{exner, rho_p, rho_zz, rtheta_p, rw, w, wwAvg},
                er.{ruAvg, ru})
do

task mpas_reconstruct_2d(cr : region(ispace(int2d), cell_fs),
                         er : region(ispace(int2d), edge_fs),
                         includeHalos : bool,
                         on_a_sphere : bool)
where
  reads (cr.{coeffs_reconstruct, edgesOnCell, lat, lon, nEdgesOnCell}, er.{u}),
  writes (cr.{uReconstructMeridional, uReconstructZonal}),
  reads writes (cr.{uReconstructX, uReconstructY, uReconstructZ})
do

task atm_rk_dynamics_substep_finish(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    dynamics_substep : int,
                                    dynamics_split : int)
where
  reads (cr.{rho_p, rho_zz_2, rho_zz_old_split, rtheta_p, rw, theta_m_2, w_2}, er.{ru, u_2}),
  writes (cr.{rho_p_save, rho_zz, rtheta_p_save, rw_save, theta_m, w}, er.{ru_save, u}),
  reads writes (cr.{wwAvg, wwAvg_split}, er.{ruAvg, ruAvg_split})
do




############ FROM INIT_ATM_CASSE.RG ########
task init_atm_case_jw(cr : region(ispace(int2d), cell_fs),
                      er : region(ispace(int2d), edge_fs),
                      vr : region(ispace(int2d), vertex_fs),
                      vertr : region(ispace(int1d), vertical_fs))
where
  reads (cr.lat,
         er.{cellsOnEdge, edgesOnEdge_ECP, lat, lon, nEdgesOnEdge, verticesOnEdge, weightsOnEdge}, 
         vr.{lat, lon}),
  writes (cr.{qsat, relhum, rho, rtheta_base, rtheta_p, surface_pressure, theta, w},
          er.{fEdge, zxu},
          vr.fVertex,
          vertr.{rdzu, rdzw}),
  reads writes (cr.{areaCell, dss, exner, hx, pressure_base, pressure_p, qv, rho_base, rho_p, rho_zz, 
                    rw, surface_pressure, theta_base, theta_m, x, y, z, zgrid, zz},
                er.{dvEdge, dcEdge, ru, u, v, x, y, z, zb, zb3},
                vr.{areaTriangle, kiteAreasOnVertex, x, y, z},
                vertr.{dzu, fzm, fzp, cf1, cf2, cf3})
do

