import "regent"


local nCells = 2562
local nEdges = 7680
local nVertices = 5120
local maxEdges = 10
local maxEdges2 = 20
local TWO = 2
local FIFTEEN = 15
local vertexDegree = 3
local nVertLevels = 1


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
    edgesOnCell : int[maxEdges],
    cellsOnCell : int[maxEdges],
    verticesOnCell : int[maxEdges],
    evc : int[3*maxEdges],   --edge pair associated with vertex v and mesh cell i. This is stored as (vertexID, edge1, edge2), and each cell has 10 of those triples arranged sequentially in the array

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
    rdzw : double, -- level dependent
    dzu : double, -- level dependent
    rdzu : double, -- level dependent
    fzm : double, -- level dependent
    fzp : double, -- level dependent

    zz : double, -- cell + level dependent
    dss : double, -- cell + level dependent: "w-damping coefficient"

    zb_cell : double[maxEdges], -- cell + level dependent
    zb3_cell : double[maxEdges], -- cell + level dependent

    rs : double, --level dependent
    ts : double, --level dependent

    -----------begin dynamics fields--------------
    kiteForCell : int[maxEdges], --Index of kite in kiteAreasOnVertex that lies within a cell for each of verticesOnCell
    edgesOnCellSign : double[maxEdges], --Sign for edges surrounding a cell: positive for positive outward normal velocity

    rtheta_pp : double,  --rho*theta_m/zz perturbation from rtheta_p. dimensions="nVertLevels nCells Time". units="kg K m^{-3}"
    rtheta_pp_old : double,  --"old time level values of rho*theta_m/zz perturbation from rtheta_p, used in 3D divergence damping". "nVertLevels nCells Time"
    rho_pp : double, --"rho/zz perturbation from rho_pp, advanced over acoustic steps" "nVertLevels nCells Time"
    rho_zz : double, --"Dry air density divided by d(zeta)/dz". dimensions="nVertLevels nCells Time" units="kg m^{-3}"

    rw : double,  --"rho*omega/zz carried at w points". dimensions="nVertLevelsP1 nCells Time"
    rw_p : double, --"acoustic perturbation rho*omega/zz carried at w points" "nVertLevelsP1 nCells Time"
    rw_save: double, --"predicted value of rho*omega/zz, saved before acoustic steps" dimensions="nVertLevelsP1 nCells Time" units="kg m^{-2} s^{-1}"
    wwAvg : double, --"time-averaged rho*omega/zz used in scalar transport" "nVertLevelsP1 nCells Time"
    invAreaCell : double, --"Inverse of Voronoi cell area"
    theta_m : double, --"Moist potential temperature: theta*(1+q_v*R_v/R_d)" --nVertLevels nCells Time"

    tend_rho : double, --name_in_code="rho_zz" "Tendency of dry density from dynamics" dimensions="nVertLevels nCells Time" units="kg m^{-3} s^{-1}" description=
    tend_rt : double, --NOT IN REGISTRY: TO FIGURE OUT
    tend_rw : double, --NOT IN REGISTRY: TO FIGURE OUT

    cofrz : double, --type="real" dimensions="nVertLevels Time" units="s m^{-1}" description="coefficient for implicit contribution of Omega to density update"
    coftz : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="s K" description="coefficient for implicit contribution of omega vertical derivative to the theta_m update"
    cofwz : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1} K^{-1}" description="coefficient for implicit contribution of density to the vertical velocity update"
    cofwr : double, --type="real" dimensions="nVertLevels nCells Time" units="m s^{-1}" description="coefficient for implicit contribution of density to the vertical velocity update"
    cofwt : double, -- type="real" dimensions="nVertLevels nCells Time" units="m s^{-1} K^{-1}" description="coefficient for implicit contribution of density to the vertical velocity update"

    a_tri : double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="implicit tridiagonal solve coefficients"
    alpha_tri : double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="implicit tridiagonal solve coefficients"
    gamma_tri : double,  --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="implicit tridiagonal solve coefficients"

    exner: double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="Exner function"
    exner_base: double, --type="real" dimensions="nVertLevels nCells Time" units="unitless" description="Base-state Exner function"
    
    w : double, --type="real" dimensions="nVertLevelsP1 nCells Time" units="m s^{-1}" description="Vertical velocity at vertical cell faces"
    specZoneMaskCell: double, --type="real" dimensions="nCells" default_value="0.0" units="-" description="0/1 mask on cells, defined as 1 for cells in the limited-area specified zone"/>
    edgesOnCell_sign: double[maxEdges], --TODO: duplicate of edgesOnCellSign; type="real" dimensions="maxEdges nCells" units="-" description="Sign for edges surrounding a cell: positive for positive outward normal velocity"
    cqw : double, -- type = "reel" dimensions="nVertLevels nCells Time" units="unitless" description="rho_d/rho_m at w points"
    qtot : double, -- defined in timestep for use in compute_dyn_tend and other subroutines
    rho_base : double, --type="real" dimensions="nVertLevels nCells Time" units="kg m^{-3}" description="Base state dry air density"
    rtheta_base : double, -- type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3}" description="reference state rho*theta/zz"
    rtheta_p : double -- type="real" dimensions="nVertLevels nCells Time" units="kg K m^{-3}" description="rho*theta_m/zz perturbation from the reference state value"
    relhum : double -- type="real" dimensions="nVertLevels nCells Time" units="percent" description="Relative humidity"
    qsat : double -- defined in init_atm_case_jw for use in moisture calculations dimensions="nVertLevels nCells" 


    ----vars first seen in atm_compute_solve_diagnostics_work--
    h : double, --type="real"     dimensions="nVertLevels nCells Time"
    ke : double, --type="real"     dimensions="nVertLevels nCells Time"
    divergence : double, --type="real"     dimensions="nVertLevels nCells Time"
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
    edgesOnVertex : int[vertexDegree],
    cellsOnVertex : int[vertexDegree],
    kiteAreasOnVertex : double[vertexDegree],

    -----------begin dynamics fields--------------
    edgesOnVertexSign : double[vertexDegree], --Sign for edges incident with a vertex: positive for po    sitive inward tengential velocity

      ----vars first seen in atm_compute_solve_diagnostics_work--
    vorticity : double,  --type="real"     dimensions="nVertLevels nVertices Time"
    invAreaTriangle : double, --type="real" dimensions="nVertices" units="m^{-2}" description="Inverse area of a Delaunay triangle"
    ke_vertex : double, -- vertex and vertical levels
    pv_vertex : double, --type="real"     dimensions="nVertLevels nVertices Time"
    fVertex: double,   --type="real"     dimensions="nVertices"/>
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
    cellsOnEdge : int[TWO],
    verticesOnEdge : int[TWO],
    edgesOnEdge_ECP : int[maxEdges2],
    weightsOnEdge : double[maxEdges2],

    -----------begin vertical structure ---------

    zxu : double, -- EDGE DEPENDENT
    zb : double[TWO], -- EDGE DEPENDENT
    zb3 : double[TWO], -- EDGE DEPENDENT

    -----------begin dynamics fields--------------
    advCellsForEdge : int[FIFTEEN], --Cells used to reconstruct a cell-based field at an edge
    nAdvCellsForEdge : int, --Number of cells used to reconstruct a cell-based field at an edge
    adv_coefs : double[FIFTEEN], --Weighting coefficents used for reconstructing cell-based fields at edges
    adv_coefs_3rd : double[FIFTEEN], --Weighting coefficents used for reconstructing cell-based fields at edges

    deriv_two : double[FIFTEEN][TWO], --weights for cell-centered second derivative, normal to edge, for transport scheme, TODO: where is it initialized?

    invDcEdge : double, --"Inverse distance between cells separated by an edge"
    ru_p : double, --"acoustic perturbation horizontal momentum at cell edge  (rho*u/zz)" "nVertLevels nEdges Time"
    cqu: double, --"rho_d/rho_m at cell edge (u points) dimensions="nVertLevels nEdges Time" units="unitless"
    specZoneMaskEdge: double, --"0/1 mask on edges, defined as 1 for edges in the limited-area specified zone" dimensions="nEdges" default_value="0.0"
    h_edge : double, --type="real"dimensions="nVertLevels nEdges Time"

    tend_ru : double, -- NOT IN REGISTRY
    ruAvg : double, -- NOT IN REGISTRY


    ----vars first seen in atm_compute_solve_diagnostics_work--
    ke_edge : double, -- parameterized by both edges and vertical levels
    u : double, --type="real"     dimensions="nVertLevels nEdges Time"
    v : double,  --type="real"     dimensions="nVertLevels nEdges Time"
    pv_edge : double, --type="real"     dimensions="nVertLevels nEdges Time"
}
