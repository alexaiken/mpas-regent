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


    -----------begin dynamics fields--------------
    advCellsForEdge : int[FIFTEEN], --Cells used to reconstruct a cell-based field at an edge
    nAdvCellsForEdge : int, --Number of cells used to reconstruct a cell-based field at an edge
    adv_coefs : double[FIFTEEN], --Weighting coefficents used for reconstructing cell-based fields at     edges
}
