import "regent"

local constants = require("constants")

-----------------------------------------------
------- FIELD SPACES FOR MESH ELEMENTS --------
-----------------------------------------------

--a hexagonal/primal cell (aka cell)
fspace cell_fs {
    cellID : int,
    lat : double,
    lon : double,
    x : double,
    y : double,
    z : double,
    meshDensity : double,
    nEdgesOnCell : int,
    areaCell : double,
    partitionNumber: int1d,
    edgesOnCell : int[constants.maxEdges],
    cellsOnCell : int[constants.maxEdges],
    verticesOnCell : int[constants.maxEdges],
    evc : int[3*constants.maxEdges],   --edge pair associated with vertex v and mesh cell i. This is stored as (vertexID, edge1, edge2), and each cell has 10 of those triples arranged sequentially in the array
    scalars : double[constants.nScalars],

    -----------begin halo fields------------------
    neighbor0 : int1d,
    neighbor1 : int1d,
    neighbor2 : int1d,
    neighbor3 : int1d,
    neighbor4 : int1d,
    neighbor5 : int1d,
    neighbor6 : int1d,
    neighbor7 : int1d,
    neighbor8 : int1d,
    neighbor9 : int1d,

    neighbor00 : int1d,
    neighbor01 : int1d,
    neighbor02 : int1d,
    neighbor03 : int1d,
    neighbor04 : int1d,
    neighbor05 : int1d,
    neighbor06 : int1d,
    neighbor07 : int1d,
    neighbor08 : int1d,
    neighbor09 : int1d,

    neighbor10 : int1d,
    neighbor11 : int1d,
    neighbor12 : int1d,
    neighbor13 : int1d,
    neighbor14 : int1d,
    neighbor15 : int1d,
    neighbor16 : int1d,
    neighbor17 : int1d,
    neighbor18 : int1d,
    neighbor19 : int1d,

    neighbor20 : int1d,
    neighbor21 : int1d,
    neighbor22 : int1d,
    neighbor23 : int1d,
    neighbor24 : int1d,
    neighbor25 : int1d,
    neighbor26 : int1d,
    neighbor27 : int1d,
    neighbor28 : int1d,
    neighbor29 : int1d,

    neighbor30 : int1d,
    neighbor31 : int1d,
    neighbor32 : int1d,
    neighbor33 : int1d,
    neighbor34 : int1d,
    neighbor35 : int1d,
    neighbor36 : int1d,
    neighbor37 : int1d,
    neighbor38 : int1d,
    neighbor39 : int1d,

    neighbor40 : int1d,
    neighbor41 : int1d,
    neighbor42 : int1d,
    neighbor43 : int1d,
    neighbor44 : int1d,
    neighbor45 : int1d,
    neighbor46 : int1d,
    neighbor47 : int1d,
    neighbor48 : int1d,
    neighbor49 : int1d,

    neighbor50 : int1d,
    neighbor51 : int1d,
    neighbor52 : int1d,
    neighbor53 : int1d,
    neighbor54 : int1d,
    neighbor55 : int1d,
    neighbor56 : int1d,
    neighbor57 : int1d,
    neighbor58 : int1d,
    neighbor59 : int1d,

    neighbor60 : int1d,
    neighbor61 : int1d,
    neighbor62 : int1d,
    neighbor63 : int1d,
    neighbor64 : int1d,
    neighbor65 : int1d,
    neighbor66 : int1d,
    neighbor67 : int1d,
    neighbor68 : int1d,
    neighbor69 : int1d,

    neighbor70 : int1d,
    neighbor71 : int1d,
    neighbor72 : int1d,
    neighbor73 : int1d,
    neighbor74 : int1d,
    neighbor75 : int1d,
    neighbor76 : int1d,
    neighbor77 : int1d,
    neighbor78 : int1d,
    neighbor79 : int1d,

    neighbor80 : int1d,
    neighbor81 : int1d,
    neighbor82 : int1d,
    neighbor83 : int1d,
    neighbor84 : int1d,
    neighbor85 : int1d,
    neighbor86 : int1d,
    neighbor87 : int1d,
    neighbor88 : int1d,
    neighbor89 : int1d,

    neighbor90 : int1d,
    neighbor91 : int1d,
    neighbor92 : int1d,
    neighbor93 : int1d,
    neighbor94 : int1d,
    neighbor95 : int1d,
    neighbor96 : int1d,
    neighbor97 : int1d,
    neighbor98 : int1d,
    neighbor99 : int1d,


    -----------begin vertical structure ----------
    zgrid : double, -- cell + level dependent


    zz : double, -- cell + level dependent
    dss : double, -- cell + level dependent: "w-damping coefficient"

    zb_cell : double[constants.maxEdges], -- cell + level dependent
    zb3_cell : double[constants.maxEdges], -- cell + level dependent

    -----------begin dynamics fields--------------
    rho : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="Dry air density"
    kiteForCell : int[constants.maxEdges], --Index of kite in kiteAreasOnVertex that lies within a cell for each of verticesOnCell
    edgesOnCellSign : double[constants.maxEdges], --Sign for edges surrounding a cell: positive for positive outward normal velocity

    rtheta_pp : double,  --rho*theta_m/zz perturbation from rtheta_p. dimensions="nVertLevels nCells Time". units="kg K m^{-3}"
    rtheta_pp_old : double,  --"old time level values of rho*theta_m/zz perturbation from rtheta_p, used in 3D divergence damping". "nVertLevels nCells Time"
    rho_pp : double, --"rho/zz perturbation from rho_pp, advanced over acoustic steps" "nVertLevels nCells Time"
    rho_zz : double, --"Dry air density divided by d(zeta)/dz". dimensions="nVertLevels nCells Time" units="kg m^{-3}"

    rw : double,  --"rho*omega/zz carried at w points". dimensions="nVertLevelsP1 nCells Time"
    rw_p : double, --"acoustic perturbation rho*omega/zz carried at w points" "nVertLevelsP1 nCells Time"
    rw_save : double, --"predicted value of rho*omega/zz, saved before acoustic steps" dimensions="nVertLevelsP1 nCells Time" units="kg m^{-2} s^{-1}"
    wwAvg : double, --"time-averaged rho*omega/zz used in scalar transport" "nVertLevelsP1 nCells Time"
    invAreaCell : double, --"Inverse of Voronoi cell area"
    theta_m : double, --"Moist potential temperature: theta*(1+q_v*R_v/R_d)" --nVertLevels nCells Time" alias:tend_rt

    tend_rho : double, --name_in_code="rho_zz" "Tendency of dry density from dynamics" dimensions="nVertLevels nCells Time" units="kg m^{-3} s^{-1}" description=

    coftz : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="s K" description="coefficient for implicit contribution of omega vertical derivative to the theta_m update"
    cofwz : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1} K^{-1}" description="coefficient for implicit contribution of density to the vertical velocity update"
    cofwr : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="coefficient for implicit contribution of density to the vertical velocity update"
    cofwt : double, -- type="real" dimensions="nVertLevels nCells Time" units="m s^{-1} K^{-1}" description="coefficient for implicit contribution of density to the vertical velocity update"

    a_tri : double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="implicit tridiagonal solve coefficients"
    b_tri : double, --Note: not found in Registry.xml
    c_tri : double, --Note: not found in Registry.xml
    alpha_tri : double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="implicit tridiagonal solve coefficients"
    gamma_tri : double,  --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="implicit tridiagonal solve coefficients"

    exner: double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="Exner function"
    exner_base: double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="Base-state Exner function"

    w : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="m s^{-1}" description="Vertical velocity at vertical cell faces" aliases:w_tend,tend_rw,tend_w
    specZoneMaskCell: double, --type="real" dimensions="nCells" default_value="0.0" units="-" description="0/1 mask on cells, defined as 1 for cells in the limited-area specified zone"/>
    edgesOnCell_sign: double[constants.maxEdges], --TODO: duplicate of edgesOnCellSign; type="real" dimensions="maxEdges nCells" units="-" description="Sign for edges surrounding a cell: positive for positive outward normal velocity"
    cqw : double, -- type = "reel" dimensions="nVertLevels nCells Time" units="unitless" description="rho_d/rho_m at w points"
    qtot : double, -- defined in timestep for use in compute_dyn_tend and other subroutines
    rho_base : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="Base state dry air density" alias:rb
    rtheta_base : double, -- type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3}" description="reference state rho*theta/zz"
    rtheta_p : double, -- type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3}" description="rho*theta_m/zz perturbation from the reference state value"
    relhum : double, -- type="real" dimensions="nVertLevels nCells Time" units="percent" description="Relative humidity"
    qsat : double, -- defined in init_atm_case_jw for use in moisture calculations dimensions="nVertLevels nCells"
    hx : double, -- type="real" dimensions="nVertLevelsP1 nCells" units="m" description="terrain influence in vertical coordinate, $h_s(x,y,\zeta)$ in Klemp (MWR 2011)"
   surface_pressure : double, -- type="real" dimensions="nCells Time" units="Pa" description="Diagnosed surface pressure"

    qv : double, -- member of scalars superarray
    ----vars first seen in atm_compute_solve_diagnostics_work--
    h : double, --type="real"     dimensions="nVertLevels nCells Time"
    ke : double, --type="real"     dimensions="nVertLevels nCells Time"
    divergence : double, --type="real"     dimensions="nVertLevels nCells Time"


    ---vars first seen in atm_compute_mesh_scaling --
    meshScalingRegionalCell : double, --type="real" dimensions="nCells" units="unitless" description="Cell-centered Scaling coefficient for relaxation zone"

    ---vars first seen in atm_init_coupled_diagnostics --
    pressure_base : double,  --type="real" dimensions="nVertLevels nCells Time" units="Pa"  description="Base state pressure"/>
    pressure_p : double,  --type="real" dimensions="nVertLevels nCells Time" units="Pa"  description="Base state pressure"/> alias:pp
    theta : double, --type="real" dimensions="nVertLevels nCells Time" units="K" description="Potential temperature"/>
    rho_p : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="rho/zz perturbation from the reference state value, advanced over acoustic steps"/>
    theta_base : double, --type="real" dimensions="nVertLevels nCells Time" units="K" description="Base state potential temperature"/>

    --vars first seen in atm_compute_output_diagnostics --
    pressure : double, --type="real" dimensions="nVertLevels nCells Time" units="Pa" description="Pressure"

    --vars first seen in atm_set_smlstep_pert_variables --
    bdyMaskCell : int, --type="integer" dimensions="nCells" default_value="0" units="-" description="Limited-area specified/relaxation zone index for cells"

    --vars first seen in atm_compute_dyn_tend_work --
    defc_a : double[constants.maxEdges], --type="real" dimensions="maxEdges nCells" units="unitless" description="Coefficients for computing the off-diagonal components of the horizontal deformation"
    defc_b : double[constants.maxEdges], --type="real" dimensions="maxEdges nCells" units="unitless" description="Coefficients for computing the diagonal components of the horizontal deformation"
    kdiff : double, --type="real" dimensions="nVertLevels nCells Time" units="m^2 s^{-1}" description="Smagorinsky horizontal eddy viscosity"
    h_divergence : double, --type="real" dimensions="nVertLevels nCells Time" units="???" description="???"
    tend_rho_physics : double, --Note: not found in Registry.xml
    dpdz : double, --Note: not found in Registry.xml
    delsq_divergence : double, --Note: not found in Registry.xml
    delsq_w : double, --Note: not found in Registry.xml
    tend_w_euler : double, --Note: not found in Registry.xml
    tend_theta : double, --type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3} s^{-1}" description="tendency of coupled potential temperature rho*theta_m/zz from dynamics and physics, updated each RK step"
    theta_m_save : double, --Note: not found in Registry.xml
    delsq_theta : double, --Note: not found in Registry.xml
    tend_theta_euler : double, --Note: not found in Registry.xml
    tend_rtheta_adv : double, --type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3} s^{-1}" description="flux divergence for rho*theta_m/zz, used in the Tiedtke convective parameterization"
    rthdynten : double, --type="real" dimensions="nVertLevels nCells Time" units="K s^{-1}" description="tendency of temperature due to horizontal and vertical advections" packages="cu_grell_freitas_in;cu_tiedtke_in"
    rt_diabatic_tend : double, --type="real" dimensions="nVertLevels nCells Time" units="kg K s^{-1}" description="Tendency of coupled potential temperature from physics"
    t_init : double, --type="real" dimensions="nVertLevels nCells" units="K" description="theta reference profile"
    tend_rtheta_physics : double, --Note: not found in Registry.xml

    -- The following are added to the region for parallelization purposes
    ru_edge_w : double, --Note: not found in Registry.xml
    flux_arr : double,
    wdwz : double,

    wdtz : double,

    --vars first seen in atm_rk_integration_setup--
    rtheta_p_save : double, --type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3}" description="predicted value rtheta_p, saved before acoustic steps"
    rho_p_save : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="predicted value rho_p, saved before acoustic steps" alias:rr_save
    w_2 : double, --Note: not found in Registry.xml. Predicted value of w, saved before acoustic steps
    theta_m_2 : double, --Note: not found in Registry.xml. Predicted value of theta_m, saved before acoustic steps
    rho_zz_2 : double, --Note: not found in Registry.xml. Predicted value of rho_zz, saved before acoustic steps
    rho_zz_old_split : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="rho/zz"

    --vars first seen in atm_rk_dynamics_substep_finish--
    wwAvg_split : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="kg m^{-2} s^{-1}" description="time-averaged rho*omega/zz used in scalar transport"

    --vars first seen in mpas_reconstruct--
    uReconstructX : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Cartesian x-component of reconstructed horizontal velocity at cell centers"
    uReconstructY : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Cartesian y-component of reconstructed horizontal velocity at cell centers"
    uReconstructZ : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Cartesian z-component of reconstructed horizontal velocity at cell centers"
    coeffs_reconstruct : double[constants.maxEdges][3], --type="real" dimensions="R3 maxEdges nCells" units="unitless" description="Coefficients to reconstruct velocity vectors at cell centers"
    uReconstructZonal : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Zonal component of reconstructed horizontal velocity at cell centers" alias:ur_cell
    uReconstructMeridional : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Meridional component of reconstructed horizontal velocity at cell centers" alias:vr_cell

    ---------------------------------------------
    -----------begin physics fields--------------
    ---------------------------------------------

    -- Physics - Radiation

    --vars first seen in radiation_lw_from_MPAS
    sfc_emiss : double, --type="real" dimensions="nCells Time" units="unitless" description="surface emissivity"
    skintemp : double, --type="real" dimensions="nCells Time" units="K" description="ground or water surface temperature"
    snow : double, --type="real" dimensions="nCells Time" units="kg m^{-2}" description="snow water equivalent"
    xice : double, --type="real" dimensions="nCells Time" units="unitless" description="fractional area coverage of sea-ice"
    xland : double, --type="real" dimensions="nCells Time" units="unitless" description="land-ocean mask (1=land including sea-ice ; 2=ocean)"
    cldfrac : double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="horizontal cloud fraction"
    re_cloud : double, -- ?
    re_ice : double, -- ?
    re_snow : double, -- ?
    sfc_albedo : double, --type="real" dimensions="nCells Time" units="unitless" description="surface albedo"
    m_ps : double, --type="real" dimensions="nCells Time" units="Pa" description="Surface pressure from match on MPAS grid"

    -- vars first seen in radiation_lw_to_MPAS
    -- The following fields are all of dimensions="nCells Time" units="W m^{-2}"
    glw : double, --all-sky downward surface longwave radiation
    lwcf : double, --top-of-atmosphere cloud longwave radiative forcing
    lwdnb : double, --all-sky downward surface longwave radiation flux
    lwdnbc : double, --clear-sky downward surface longwave radiation flux
    lwdnt : double, --all-sky downward top-of-the-atmosphere longwave radiation flux
    lwdntc : double, --clear-sky downward top-of-the-atmosphere longwave radiation flux
    lwupb : double, --all-sky upward surface longwave radiation flux
    lwupbc : double, --clear-sky upward surface longwave radiation flux
    lwupt : double, --all-sky upward top-of-the-atmosphere longwave radiation flux
    lwuptc : double, --clear-sky upward top-of-the-atmosphere longwave radiation flux
    olrtoa : double, --all-sky top-of-atmosphere outgoing longwave radiation flux

    rthratenlw : double, --type="real" dimensions="nVertLevels nCells Time" units="K s^{-1}" description="tendency of potential temperature due to short wave radiation"
    rre_cloud : double, --type="real" dimensions="nVertLevels nCells Time" units="microns" description="effective radius of cloud water droplets calculated in RRTMG radiation"
    rre_ice : double, --type="real" dimensions="nVertLevels nCells Time" units="microns" description="effective radius of cloud ice crystals calculated in RRTMG radiation"
    rre_snow : double, --type="real" dimensions="nVertLevels nCells Time" units="microns" description="effective radius of snow crystals calculated in RRTMG radiation"

    o32d : double,
    p2d : double,
    o3vmr : double,

    gsw : double, --type="real" dimensions="nCells Time" units="W m^{-2}" description="net surface shortwave radiation flux"
    swcf : double, --type="real" dimensions="nCells Time" units="W m^{-2}" description="top-of-atmosphere cloud shortwave radiative forcing"
    coszr : double, --type="real" dimensions="nCells Time" units="unitless" description="cosine of zenith solar angle"

    -- Temporary versions to be used in radiation. Currently included to match original style of radiation code. TODO: Are these necessary?
    cldfrac_p : double,
    absnxt_p : double[constants.cam_abs_dim1],
    abstot_p : double[constants.cam_abs_dim2],

    ------ From mpas_atmphys_interface.F, allocate_forall_physics() ------

    psfc_p : double,        --allocate(psfc_p(ims:ime,jms:jme))
    ptop_p : double,        --allocate(ptop_p(ims:ime,jms:jme))
    
    u_p : double,         --allocate(u_p(ims:ime,kms:kme,jms:jme))
    v_p : double,         --allocate(v_p(ims:ime,kms:kme,jms:jme))
    fzm_p : double,       --allocate(fzm_p(ims:ime,kms:kme,jms:jme))
    fzp_p : double,       --allocate(fzp_p(ims:ime,kms:kme,jms:jme))
    zz_p : double,        --allocate(zz_p(ims:ime,kms:kme,jms:jme))
    pres_p : double,      --allocate(pres_p(ims:ime,kms:kme,jms:jme))
    pi_p : double,        --allocate(pi_p(ims:ime,kms:kme,jms:jme))
    z_p : double,         --allocate(z_p(ims:ime,kms:kme,jms:jme))
    zmid_p : double,      --allocate(zmid_p(ims:ime,kms:kme,jms:jme))
    dz_p : double,        --allocate(dz_p(ims:ime,kms:kme,jms:jme))
    t_p : double,         --allocate(t_p(ims:ime,kms:kme,jms:jme))
    th_p : double,        --allocate(th_p(ims:ime,kms:kme,jms:jme))
    al_p : double,        --allocate(al_p(ims:ime,kms:kme,jms:jme))
    rh_p : double,        --allocate(rh_p(ims:ime,kms:kme,jms:jme))
    znu_p : double,       --allocate(znu_p(ims:ime,kms:kme,jms:jme))

    w_p : double,         --allocate(w_p(ims:ime,kms:kme,jms:jme))
    pres2_p : double,     --allocate(pres2_p(ims:ime,kms:kme,jms:jme))
    t2_p : double,        --allocate(t2_p(ims:ime,kms:kme,jms:jme))

    qv_p : double,        --allocate(qv_p(ims:ime,kms:kme,jms:jme))
    qc_p : double,        --allocate(qc_p(ims:ime,kms:kme,jms:jme))
    qr_p : double,        --allocate(qr_p(ims:ime,kms:kme,jms:jme))
    qi_p : double,        --allocate(qi_p(ims:ime,kms:kme,jms:jme))
    qs_p : double,        --allocate(qs_p(ims:ime,kms:kme,jms:jme))
    qg_p : double,        --allocate(qg_p(ims:ime,kms:kme,jms:jme))

    ni_p : double,        --allocate(ni_p(ims:ime,kms:kme,jms:jme))

    -- arrays used for calculating the hydrostatic pressure and exner function
    psfc_hyd_p : double,    --allocate(psfc_hyd_p(ims:ime,jms:jme))
    psfc_hydd_p : double,   --allocate(psfc_hydd_p(ims:ime,jms:jme))
    pres_hyd_p : double,    --allocate(pres_hyd_p(ims:ime,kms:kme,jms:jme))
    pres_hydd_p : double,   --allocate(pres_hydd_p(ims:ime,kms:kme,jms:jme))
    pres2_hyd_p : double,   --allocate(pres2_hyd_p(ims:ime,kms:kme,jms:jme))
    pres2_hydd_p : double,  --allocate(pres2_hydd_p(ims:ime,kms:kme,jms:jme))
    znu_hyd_p : double,     --allocate(znu_hyd_p(ims:ime,kms:kme,jms:jme))

    ------ From mpas_atmphys_driver_radiation_lw.F, allocate_radiation_lw() ------

    f_ice : double,             --allocate(f_ice(ims:ime,kms:kme,jms:jme))
    f_rain : double,            --allocate(f_rain(ims:ime,kms:kme,jms:jme))

    sfc_emiss_p : double,   --allocate(sfc_emiss_p(ims:ime,jms:jme))
    snow_p : double,        --allocate(snow_p(ims:ime,jms:jme))
    tsk_p : double,         --allocate(tsk_p(ims:ime,jms:jme))
    xice_p : double,        --allocate(xice_p(ims:ime,jms:jme))
    xland_p : double,       --allocate(xland_p(ims:ime,jms:jme))

    glw_p : double,         --allocate(glw_p(ims:ime,jms:jme))
    lwcf_p : double,        --allocate(lwcf_p(ims:ime,jms:jme))
    lwdnb_p : double,       --allocate(lwdnb_p(ims:ime,jms:jme))
    lwdnbc_p : double,      --allocate(lwdnbc_p(ims:ime,jms:jme))
    lwdnt_p : double,       --allocate(lwdnt_p(ims:ime,jms:jme))
    lwdntc_p : double,      --allocate(lwdntc_p(ims:ime,jms:jme))
    lwupb_p : double,       --allocate(lwupb_p(ims:ime,jms:jme))
    lwupbc_p : double,      --allocate(lwupbc_p(ims:ime,jms:jme))
    lwupt_p : double,       --allocate(lwupt_p(ims:ime,jms:jme))
    lwuptc_p : double,      --allocate(lwuptc_p(ims:ime,jms:jme))
    olrtoa_p : double,      --allocate(olrtoa_p(ims:ime,jms:jme))

    rthratenlw_p : double,      --allocate(rthratenlw_p(ims:ime,kms:kme,jms:jme))

    -- case("rrtmg_lw")
    recloud_p : double,         --allocate(recloud_p(ims:ime,kms:kme,jms:jme))
    reice_p : double,           --allocate(reice_p(ims:ime,kms:kme,jms:jme))
    resnow_p : double,          --allocate(resnow_p(ims:ime,kms:kme,jms:jme))
    rrecloud_p : double,        --allocate(rrecloud_p(ims:ime,kms:kme,jms:jme))
    rreice_p : double,          --allocate(rreice_p(ims:ime,kms:kme,jms:jme))
    rresnow_p : double,         --allocate(rresnow_p(ims:ime,kms:kme,jms:jme))

    -- case("cam_lw")
    xlat_p : double,        --allocate(xlat_p(ims:ime,jms:jme))
    xlon_p : double,        --allocate(xlon_p(ims:ime,jms:jme))
    gsw_p : double,         --allocate(gsw_p(ims:ime,jms:jme))
    swcf_p : double,        --allocate(swcf_p(ims:ime,jms:jme))
    swdnb_p : double,       --allocate(swdnb_p(ims:ime,jms:jme))
    swdnbc_p : double,      --allocate(swdnbc_p(ims:ime,jms:jme))
    swdnt_p : double,       --allocate(swdnt_p(ims:ime,jms:jme))
    swdntc_p : double,      --allocate(swdntc_p(ims:ime,jms:jme))
    swupb_p : double,       --allocate(swupb_p(ims:ime,jms:jme))
    swupbc_p : double,      --allocate(swupbc_p(ims:ime,jms:jme))
    swupt_p : double,       --allocate(swupt_p(ims:ime,jms:jme))
    swuptc_p : double,      --allocate(swuptc_p(ims:ime,jms:jme))
    coszr_p : double,       --allocate(coszr_p(ims:ime,jms:jme))
    sfc_albedo_p : double,  --allocate(sfc_albedo_p(ims:ime,jms:jme))
    rthratensw_p : double,      --allocate(rthratensw_p(ims:ime,kms:kme,jms:jme))
    --
    cemiss_p : double,          --allocate(cemiss_p(ims:ime,kms:kme,jms:jme))
    taucldc_p : double,         --allocate(taucldc_p(ims:ime,kms:kme,jms:jme))
    taucldi_p : double,         --allocate(taucldi_p(ims:ime,kms:kme,jms:jme))
    --
    m_psn_p : double,       --allocate(m_psn_p(ims:ime,jms:jme))
    m_psp_p : double,       --allocate(m_psp_p(ims:ime,jms:jme))

    -- !allocate these arrays on the first time step, only:
    emstot_p : double,          --allocate(emstot_p(ims:ime,kms:kme,jms:jme))

    ------ camrad() 1d locals ------

    coszrs : double,
    landfrac : double,
    landm : double,
    snowh : double,
    icefrac : double,
    lwups : double,
    asdir : double,
    asdif : double,
    aldir : double,
    aldif : double,
    ps : double,
    nmxrgn : int,       -- Number of maximally overlapped regions

    fsns : double,      -- Surface absorbed solar flux
    fsnt : double,      -- Net column abs solar flux at model top
    flns : double,      -- Srf longwave cooling (up-down) flux
    flnt : double,      -- Net outgoing lw flux at model top

    swcftoa : double,   -- Top of the atmosphere solar cloud forcing
    lwcftoa : double,   -- Top of the atmosphere longwave cloud forcing

    sols : double,      -- Downward solar rad onto surface (sw direct)
    soll : double,      -- Downward solar rad onto surface (lw direct)
    solsd : double,     -- Downward solar rad onto surface (sw diffuse)
    solld : double,     -- Downward solar rad onto surface (lw diffuse)
    fsds : double,      -- Flux Shortwave Downwelling Surface
    flwds : double,     -- Surface down longwave flux
    m_psjp : double,    -- MATCH surface pressure
    m_psjn : double,
    clat : double,      -- latitude in radians for columns

    ------ camrad() 2d locals ------

    cld : double, 
    pmid : double,      -- level pressures (mks)
    lnpmid : double, 
    pdel : double, 
    zm : double, 
    t : double,
    cicewp : double,    -- in-cloud cloud ice water path
    cliqwp : double,    -- in-cloud cloud liquid water path
    emis : double,      -- cloud emissivity
    rel : double,       -- effective drop radius (microns)
    rei : double,       -- ice effective drop size (microns)

    qrs : double,       -- Solar heating rate
    qrl : double,       -- Longwave cooling rate

    -- extended (1 to kte-kts + 2) --
    pint : double, 
    lnpint : double,
    pmxrgn : double,    -- Maximum values of pressure for each
    
    -- Added outputs of total and clearsky fluxes etc
    fsup : double,      -- Upward total sky solar
    fsupc : double,     -- Upward clear sky solar
    fsdn : double,      -- Downward total sky solar
    fsdnc : double,     -- Downward clear sky solar
    flup : double,      -- Upward total sky longwave
    flupc : double,     -- Upward clear sky longwave
    fldn : double,      -- Downward total sky longwave  
    fldnc : double,     -- Downward clear sky longwave

    -- left shifted (0 to kte-kts+1) --
    tauxcl : double,    -- cloud water optical depth
    tauxci : double,    -- cloud ice optical depth

    -----------end physics fields----------------
    ---------------------------------------------
}


-- A triangluar/dual cell (aka vertex)
fspace vertex_fs {
    vertexID : int,
    lat : double,
    lon : double,
    x : double,
    y : double,
    z : double,
    areaTriangle : double,
    edgesOnVertex : int[constants.vertexDegree],
    cellsOnVertex : int[constants.vertexDegree],
    kiteAreasOnVertex : double[constants.vertexDegree],

    -----------begin dynamics fields--------------
    edgesOnVertexSign : double[constants.vertexDegree], --Sign for edges incident with a vertex: positive for po    sitive inward tengential velocity
    fVertex : double, --type="real" dimensions="nVertices" units="unitless" description="Coriolis parameter at a vertex"
      ----vars first seen in atm_compute_solve_diagnostics_work--
    vorticity : double,  --type="real"     dimensions="nVertLevels nVertices Time"
    invAreaTriangle : double, --type="real" dimensions="nVertices" units="m^{-2}" description="Inverse area of a Delaunay triangle"
    ke_vertex : double, -- vertex and vertical levels
    pv_vertex : double, --type="real"     dimensions="nVertLevels nVertices Time"

    -- vars first seen in atm_compute_dyn_tend_work --
    delsq_vorticity : double, --Note: not found in Registry.xml
    edgesOnVertex_sign : double[constants.vertexDegree], --type="real" dimensions="vertexDegree nVertices" units="-" description="Sign for edges incident with a vertex: positive for positive inward tengential velocity"
}

fspace edge_fs {
    edgeID : int,
    lat : double,
    lon : double,
    x : double,
    y : double,
    z : double,
    nEdgesOnEdge : int,
    dvEdge : double,
    dv1Edge : double,
    dv2Edge : double,
    dcEdge : double,
    angleEdge : double,
    cellsOnEdge : int[constants.TWO],
    verticesOnEdge : int[constants.TWO],
    edgesOnEdge_ECP : int[constants.maxEdges2],
    weightsOnEdge : double[constants.maxEdges2],
    edgesOnEdge : int[constants.maxEdges2],

    -----------begin vertical structure ---------

    zxu : double, -- EDGE DEPENDENT
    zb : double[constants.TWO], -- EDGE DEPENDENT
    zb3 : double[constants.TWO], -- EDGE DEPENDENT

    -----------begin dynamics fields--------------
    advCellsForEdge : int[constants.FIFTEEN], --Cells used to reconstruct a cell-based field at an edge
    nAdvCellsForEdge : int, --Number of cells used to reconstruct a cell-based field at an edge
    adv_coefs : double[constants.FIFTEEN], --Weighting coefficents used for reconstructing cell-based fields at edges
    adv_coefs_3rd : double[constants.FIFTEEN], --Weighting coefficents used for reconstructing cell-based fields at edges

    deriv_two : double[constants.FIFTEEN * constants.TWO], --weights for cell-centered second derivative, normal to edge, for transport scheme, TODO: where is it initialized?

    invDcEdge : double, --"Inverse distance between cells separated by an edge"
    ru_p : double, --"acoustic perturbation horizontal momentum at cell edge  (rho*u/zz)" "nVertLevels nEdges Time"
    cqu: double, --"rho_d/rho_m at cell edge (u points) dimensions="nVertLevels nEdges Time" units="unitless"
    specZoneMaskEdge: double, --"0/1 mask on edges, defined as 1 for edges in the limited-area specified zone" dimensions="nEdges" default_value="0.0"
    h_edge : double, --type="real"dimensions="nVertLevels nEdges Time"

    tend_ru : double, -- NOT IN REGISTRY
    ruAvg : double, -- NOT IN REGISTRY

    fEdge : double, --type="real" dimensions="nEdges" units="unitless" description="Coriolis parameter at an edge"

    -- vars first seen in atm_compute_solve_diagnostics_work --
    ke_edge : double, -- parameterized by both edges and vertical levels
    u : double, --type="real"     dimensions="nVertLevels nEdges Time"
    v : double,  --type="real"     dimensions="nVertLevels nEdges Time"
    pv_edge : double, --type="real"     dimensions="nVertLevels nEdges Time"


    -- vars first seen in atm_compute_mesh_scaling --
    meshScalingDel2 : double, --type="real"     dimensions="nEdges"
    meshScalingDel4 : double, --type="real"     dimensions="nEdges"
    meshScalingRegionalEdge : double,  --type="real" dimensions="nEdges" units="unitless" description="Edge-centered Scaling coefficient for relaxation zone"/>

    -- vars first seen in atm_init_coupled_diagnostics --
    ru : double, --type="real" dimensions="nVertLevels nEdges Time" units="kg m^{-2} s^{-1}" description="horizontal momentum at cell edge (rho*u/zz)"/>

    -- vars first seen in atm_set_smlstep_pert_variables --
    u_tend : double, --Note: not found in Registry.xml

    -- vars first seen in atm_compute_dyn_tend_work --
    tend_u : double, --type="real" dimensions="nVertLevels nEdges Time" units="m s^{-2}" description="Tendency of u from dynamics"
    rho_edge : double, --type="real" dimensions="nVertLevels nEdges Time" units="kg m^{-3}" description="rho_zz averaged from cell centers to the cell edge"
    tend_u_euler : double, --Note: not found in Registry.xml
    invDvEdge : double, --type="real" dimensions="nEdges" units="m^{-1}" description="Inverse distance between vertex endpoints of an edge"
    delsq_u : double, --Note: not found in Registry.xml
    tend_ru_physics : double, --Note: not found in Registry.xml
    ru_save : double, --type="real" dimensions="nVertLevels nEdges Time" units="kg m^{-2} s^{-1}" description="predicted value of horizontal momentum, saved before acoustic steps"
    -- The following are added to the region for parallelization purposes
    wduz : double, --Note: not found in Registry.xml
    q : double, --Note: not found in Registry.xml
    u_mix : double, --Note: not found in Registry.xml

    -- vars first seen in atm_rk_integration_setup --
    u_2 : double, --Note: not found in Registry.xml. Predicted value of u, saved before acoustic steps

    -- vars first seen in atm_rk_dynamics_substep_finish --
    ruAvg_split : double, --type="real" dimensions="nVertLevels nEdges Time" units="kg m^{-2} s^{-1}" description="time-averaged rho*u/zz used in scalar transport"
}


--Field space for variables that are only level-dependent.
fspace vertical_fs {

  rdzw : double, -- level dependent
  dzu : double, -- level dependent
  rdzu : double, -- level dependent
  fzm : double, -- level dependent
  fzp : double, -- level dependent
  cofrz : double, --type="real" dimensions="nVertLevels Time" units="s m^{-1}" description="coefficient for implicit contribution of Omega to density update"

  -- vars first seen in atm_compute_dyn_tend_work --
  u_init : double, --type="real" dimensions="nVertLevels" units="m s^{-1}" description="u reference profile"
  v_init : double, --type="real" dimensions="nVertLevels" units="m s^{-1}" description="v reference profile"

  --vars first seen in atm_recover_large_step_variables_work--
  cf1 : double, --type="real" dimensions="" units="unitless" description="Surface interpolation weight for level k=1 value"
  cf2 : double, --type="real" dimensions="" units="unitless" description="Surface interpolation weight for level k=1 value"
  cf3 : double, --type="real" dimensions="" units="unitless" description="Surface interpolation weight for level k=1 value"
}

fspace phys_tbls_fs
{
  estblh2o : double[constants.ntemp],   -- table of H2O saturation vapor pressures

  -- Table of saturation vapor pressure values es from tmin degrees
  -- to tmax+1 degrees k in one degree increments.  ttrice defines the
  -- transition region where es is a combination of ice & water values
  estbl : double[constants.plenest],    -- table values of saturation vapor pressure

  itype : double,
  pcf : double[5],          -- polynomial coeffs -> es transition water to ice

  -- es table parameters
  tmin : double,            -- Minimum temperature (K)
  tmax : double,            -- Maximum temperature (K)
  ttrice : double,          -- Trans. range from es over h2o to es over ice
  icephs : bool,            -- Ice phase (true or false)

  -- Physical constants required for es calculation
  epsqs : double,
  hlatv : double,           -- latent heat of evaporation ~ J/kg
  hlatf : double,           -- latent heat of fusion ~ J/kg
  rgasv : double,
  cp : double,
  tmelt : double,

  lentbl : int,

  idxVOLC : int,
}

fspace ozn_fs {
  o3clim : double, --type="real" dimensions="nOznLevels nCells Time" units="mol mol^{-1}" description="climatological ozone on prescribed pressure levels at current time"
  pin : double, --type="real"  dimensions="nOznLevels" units="Pa" description="fixed pressure levels at which climatological ozone is defined"
  ozmixm : double[constants.nMonths], --type="real" dimensions="nMonths nOznLevels nCells" units="mol mol^{-1}" description="monthly-mean climatological ozone defined at fixed pressure levels"

  --Temporary versions to be used in radiation code
  o3clim_p : double,
  pin_p : double,
  ozmixm_p : double[constants.nMonths],
}

fspace aerosol_fs {
  aerosols : double[constants.nAerosols], --type="real" dimensions="nAerLevels nCells Time"
  m_hybi : double, --type="real" dimensions="nAerLevels nCells" units="unitless" description="Matched hybi (needs to be re-checked)"

  --Temporary versions to be used in radiation code
  aerosolcp_p : double[constants.nAerosols],
  aerosolcn_p : double[constants.nAerosols],
  m_hybi_p : double,
}

fspace camrad_1d_fs {
  -- Arrays of size (ite-its)
  coszrs : double,
  landfrac : double,
  landm : double,
  snowh : double,
  icefrac : double,
  lwups : double,

  asdir : double,
  asdif : double,
  aldir : double,
  aldif : double,
  ps : double,

  nmxrgn : int,           -- Number of maximally overlapped regions

  fsns : double,          -- Surface absorbed solar flux
  fsnt : double,          -- Net column abs solar flux at model top
  flns : double,          -- Srf longwave cooling (up-down) flux
  flnt : double,          -- Net outgoing lw flux at model top

  swcftoa : double,
  lwcftoa : double,
  olrtoa : double,

  sols : double,          -- Downward solar rad onto surface (sw direct)
  soll : double,          -- Downward solar rad onto surface (lw direct)
  solsd : double,         -- Downward solar rad onto surface (sw diffuse)
  solld : double,         -- Downward solar rad onto surface (lw diffuse)

  fsds : double,          -- Flux Shortwave Downwelling Surface
  flwds : double,         -- Surface down longwave flux

  m_psjp : double,        -- MATCH surface pressure
  m_psjn : double,        -- MATCH surface pressure

  clat : double,          -- latitude in radians for columns
}

fspace camrad_2d_fs {
  -- Arrays of size (ite-its, kte-kts) aka (nCells, nVertLevels)
  -- NOTE: some are kte-ktw+1
  cld : double,
  pmid : double,
  lnpmid : double,
  pdel : double,
  zm : double,
  t : double,
  pint : double,
  lnpint : double,
  
  cicewp : double,        -- in-cloud cloud ice water path
  cliqwp : double,        -- in-cloud cloud liquid water path
  tauxcl : double,        -- cloud water optical depth
  tauxci : double,        -- loud ice optical depth
  emis : double,          -- cloud emissivity
  rel : double,           -- effective drop radius (microns)
  rei : double,           -- ice effective drop size (microns)
  pmxrgn : double,        -- Maximum values of pressure for each

  fsup : double,          -- Upward total sky solar
  fsupc : double,         -- Upward clear sky solar
  fsdn : double,          -- Downward total sky solar
  fsdnc : double,         -- Downward clear sky solar
  flup : double,          -- Upward total sky longwave
  flupc : double,         -- Upward clear sky longwave
  fldn : double,          -- Downward total sky longwave
  fldnc : double,         -- Downward clear sky longwave

  qrs : double,           -- Solar heating rate
  qrl : double,           -- Longwave cooling rate
}

fspace radctl_1d_fs {
  -- Arrays of size (pcols) --
  solin : double,         -- Solar incident flux
  fsntoa : double,        -- Net solar flux at TOA
  fsntoac : double,       -- Clear sky net solar flux at TOA
  fsnirt : double,        -- Near-IR flux absorbed at toa
  fsnrtc : double,        -- Clear sky near-IR flux absorbed at toa
  fsnirtsq : double,      -- Near-IR flux absorbed at toa >= 0.7 microns
  fsntc : double,         -- Clear sky total column abs solar flux
  fsnsc : double,         -- Clear sky surface abs solar flux
  fsdsc : double,         -- Clear sky surface downwelling solar flux
  flutc : double,         -- Upward Clear Sky flux at top of model
  flntc : double,         -- Clear sky lw flux at model top
  flnsc : double,         -- Clear sky lw flux at srf (up-down)

  lwupcgs : double,       -- Upward longwave flux in cgs units

  frc_day : double,       -- = 1 for daylight, =0 for night colums
}

fspace radctl_2d_pver_fs {
  -- Arrays of size (pcols, pver) --
  ftem : double,          -- temporary array for outfld
  
  n2o : double,      -- nitrous oxide mass mixing ratio
  ch4 : double,      -- methane mass mixing ratio
  cfc11 : double,    -- cfc11 mass mixing ratio
  cfc12 : double,    -- cfc12 mass mixing ratio
}

fspace radctl_2d_pverr_fs {
  -- Arrays of size (pcols, pverr) --
  pbr : double,     -- Model mid-level pressures (dynes/cm2)

  o3vmr : double,   -- Ozone volume mixing ratio
  o3mmr : double,   -- Ozone mass mixing ratio
  rh : double,      -- level relative humidity (fraction)

  esat : double,    -- saturation vapor pressure
  qsat : double,    -- saturation specific humidity
}

fspace radctl_3d_aero_fs {
  -- Arrays of size (pcols, nspint, naer_groups)
  aertau : double, -- Aerosol column optical depth
  aerssa : double, -- Aerosol column averaged single scattering albedo
  aerasm : double, -- Aerosol column averaged asymmetry parameter
  aerfwd : double, -- Aerosol column averaged forward scattering
}
