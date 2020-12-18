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
    theta_m : double, --"Moist potential temperature: theta*(1+q_v*R_v/R_d)" --nVertLevels nCells Time"

    tend_rho : double, --name_in_code="rho_zz" "Tendency of dry density from dynamics" dimensions="nVertLevels nCells Time" units="kg m^{-3} s^{-1}" description=
    tend_rt : double, --NOT IN REGISTRY: TO FIGURE OUT
    tend_rw : double, --NOT IN REGISTRY: TO FIGURE OUT

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

    w : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="m s^{-1}" description="Vertical velocity at vertical cell faces"
    specZoneMaskCell: double, --type="real" dimensions="nCells" default_value="0.0" units="-" description="0/1 mask on cells, defined as 1 for cells in the limited-area specified zone"/>
    edgesOnCell_sign: double[constants.maxEdges], --TODO: duplicate of edgesOnCellSign; type="real" dimensions="maxEdges nCells" units="-" description="Sign for edges surrounding a cell: positive for positive outward normal velocity"
    cqw : double, -- type = "reel" dimensions="nVertLevels nCells Time" units="unitless" description="rho_d/rho_m at w points"
    qtot : double, -- defined in timestep for use in compute_dyn_tend and other subroutines
    rho_base : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="Base state dry air density"
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
    pressure_p : double,  --type="real" dimensions="nVertLevels nCells Time" units="Pa"  description="Base state pressure"/>
    theta : double, --type="real" dimensions="nVertLevels nCells Time" units="K" description="Potential temperature"/>
    rho_p : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="rho/zz perturbation from the reference state value, advanced over acoustic steps"/>
    theta_base : double, --type="real" dimensions="nVertLevels nCells Time" units="K" description="Base state potential temperature"/>

    --vars first seen in atm_compute_output_diagnostics --
    pressure : double, --type="real" dimensions="nVertLevels nCells Time" units="Pa" description="Pressure"

    --vars first seen in atm_set_smlstep_pert_variables --
    bdyMaskCell : int, --type="integer" dimensions="nCells" default_value="0" units="-" description="Limited-area specified/relaxation zone index for cells"
    w_tend : double, --Note: not found in Registry.xml

    --vars first seen in atm_compute_dyn_tend_work --
    defc_a : double[constants.maxEdges], --type="real" dimensions="maxEdges nCells" units="unitless" description="Coefficients for computing the off-diagonal components of the horizontal deformation"
    defc_b : double[constants.maxEdges], --type="real" dimensions="maxEdges nCells" units="unitless" description="Coefficients for computing the diagonal components of the horizontal deformation"
    kdiff : double, --type="real" dimensions="nVertLevels nCells Time" units="m^2 s^{-1}" description="Smagorinsky horizontal eddy viscosity"
    h_divergence : double, --type="real" dimensions="nVertLevels nCells Time" units="???" description="???"
    tend_rho_physics : double, --Note: not found in Registry.xml
    dpdz : double, --Note: not found in Registry.xml
    rb : double, --Note: not found in Registry.xml
    rr_save : double, --Note: not found in Registry.xml
    pp : double, --Note: not found in Registry.xml
    delsq_divergence : double, --Note: not found in Registry.xml
    ur_cell : double, --Note: not found in Registry.xml
    vr_cell : double, --Note: not found in Registry.xml
    tend_w : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="m s^{-2}" description="Tendency of w from dynamics"
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
    rho_p_save : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="predicted value rho_p, saved before acoustic steps" 
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
    uReconstructZonal : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Zonal component of reconstructed horizontal velocity at cell centers"
    uReconstructMeridional : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="Meridional component of reconstructed horizontal velocity at cell centers"

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
    f_ice : double, --Note: not found in Registry.xml
    f_rain : double, --Note: not found in Registry.xml

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
    pres_hyd_p : double,
    o3vmr : double,

    -- Temporary versions to be used in radiation. Currently included to match original style of radiation code. TODO: Are these necessary?
    sfc_emiss_p : double,
    tsk_p : double, --skintemp
    snow_p : double,
    xice_p : double,
    xland_p : double,
    cldfrac_p : double,
    recloud_p : double,
    reice_p : double,
    resnow_p : double,
    rrecloud_p : double,
    rreice_p : double,
    rresnow_p : double,
    sfc_albedo_p : double,
    m_psp_p : double,
    m_psn_p : double,
    glw_p : double,
    lwcf_p : double,
    lwdnb_p : double,
    lwdnbc_p : double,
    lwdnt_p : double,
    lwdntc_p : double,
    lwupb_p : double,
    lwupbc_p : double,
    lwupt_p : double,
    lwuptc_p : double,
    olrtoa_p : double,
    rthratenlw_p : double,
    xlat_p : double,
    xlon_p : double,
    coszr_p : double,
    gsw_p : double,
    swcf_p : double,
    swdnb_p : double,
    swdnbc_p : double,
    swdnt_p : double,
    swdntc_p : double,
    swupb_p : double,
    swupbc_p : double,
    swupt_p : double,
    swuptc_p : double,
    rthratensw_p : double,
    cemiss_p : double,
    taucldc_p : double,
    taucldi_p : double,
    absnxt_p : double[constants.cam_abs_dim1],
    abstot_p : double[constants.cam_abs_dim2],
    emstot_p : double,
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