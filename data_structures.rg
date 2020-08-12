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
    neighbor11 : int1d,
    neighbor22 : int1d,
    neighbor33 : int1d,
    neighbor44 : int1d,
    neighbor55 : int1d,
    neighbor66 : int1d,
    neighbor77 : int1d,
    neighbor88 : int1d,
    neighbor99 : int1d,


    -----------begin vertical structure ----------
    zgrid : double, -- cell + level dependent
    rdzw : double, -- level dependent
    dzu : double, -- level dependent
    rdzu : double, -- level dependent
    fzm : double, -- level dependent
    fzp : double, -- level dependent

    zz : double, -- cell + level dependent

    zb_cell : double[maxEdges], -- cell + level dependent
    zb3_cell : double[maxEdges], -- cell + level dependent

    -----------begin dynamics fields--------------
    kiteForCell : int[maxEdges], --Index of kite in kiteAreasOnVertex that lies within a cell for eac    h of verticesOnCell
    edgesOnCellSign : double[maxEdges], --Sign for edges surrounding a cell: positive for positive outwa    rd normal velocity

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
    edgesOnVertexSign : double[vertexDegree] --Sign for edges incident with a vertex: positive for po    sitive inward tengential velocity
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

}
