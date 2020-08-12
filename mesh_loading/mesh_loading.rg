import "regent"

require "data_structures"
require "netcdf_tasks"

local c = regentlib.c
local cio = terralib.includec("stdio.h")
local clib = terralib.includec("stdlib.h")

terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")
local netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

local nCells = 2562
local nEdges = 7680
local nVertices = 5120
local maxEdges = 10
local maxEdges2 = 20
local TWO = 2
local vertexDegree = 3
local nVertLevels = 1

local FILE_NAME = "x1.2562.grid.nc"
local GRAPH_FILE_NAME = "x1.2562.graph.info.part.16"
local MAXCHAR = 5
local NUM_PARTITIONS = 16


--Terra function to read the cell partitions from graph.info file. Returns an array where each element is the partition number of that cell index.
terra read_file(file_name: &int8) : int[nCells]
    var file = c.fopen(file_name, "r")
    regentlib.assert(file ~= nil, "failed to open graph.info file")
    var str : int8[MAXCHAR]
    var partition_array : int[nCells]
    var i = 0
    while c.fgets(str, MAXCHAR, file) ~= nil do
        partition_array[i] = c.atoi(str)
        i = i+1
    end
    return partition_array
end


task main()

    -------------------------------------------
    ----- READ VARIABLES FROM NETCDF FILE -----
    -------------------------------------------
    cio.printf("Starting to read file... \n")
    var ncid : int

    -- Open the file and store the NCID
    open_file(&ncid, FILE_NAME)

    -- Define the variable IDs
    var latCell_varid : int
    var lonCell_varid : int
    var meshDensity_varid : int
    var xCell_varid : int
    var yCell_varid : int
    var zCell_varid : int
    var indexToCellID_varid : int
    var latEdge_varid : int
    var lonEdge_varid : int
    var xEdge_varid : int
    var yEdge_varid : int
    var zEdge_varid : int
    var indexToEdgeID_varid : int
    var latVertex_varid : int
    var lonVertex_varid : int
    var xVertex_varid : int
    var yVertex_varid : int
    var zVertex_varid : int
    var indexToVertexID_varid : int
    var cellsOnEdge_varid : int
    var nEdgesOnCell_varid : int
    var nEdgesOnEdge_varid : int
    var edgesOnCell_varid : int
    var edgesOnEdge_varid : int
    var weightsOnEdge_varid : int
    var dvEdge_varid : int
    var dv1Edge_varid : int
    var dv2Edge_varid : int
    var dcEdge_varid : int
    var angleEdge_varid : int
    var areaCell_varid : int
    var areaTriangle_varid : int
    var cellsOnCell_varid : int
    var verticesOnCell_varid : int
    var verticesOnEdge_varid : int
    var edgesOnVertex_varid : int
    var cellsOnVertex_varid : int
    var kiteAreasOnVertex_varid : int

    -- Define and malloc the data structures to store the variable values
    var latCell_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var lonCell_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var meshDensity_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var xCell_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var yCell_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var zCell_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var indexToCellID_in : &int = [&int](c.malloc([sizeof(int)] * nCells))
    var latEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var lonEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var xEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var yEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var zEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var indexToEdgeID_in : &int = [&int](c.malloc([sizeof(int)] * nEdges))
    var latVertex_in : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var lonVertex_in : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var xVertex_in : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var yVertex_in : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var zVertex_in : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var indexToVertexID_in : &int = [&int](c.malloc([sizeof(int)] * nVertices))
    var cellsOnEdge_in : &int = [&int](c.malloc([sizeof(int)] * nEdges*TWO))
    var nEdgesOnCell_in : &int = [&int](c.malloc([sizeof(int)] * nCells))
    var nEdgesOnEdge_in : &int = [&int](c.malloc([sizeof(int)] * nEdges))
    var edgesOnCell_in : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
    var edgesOnEdge_in : &int = [&int](c.malloc([sizeof(int)] * nEdges*maxEdges2))
    var weightsOnEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges*maxEdges2))
    var dvEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var dv1Edge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var dv2Edge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var dcEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var angleEdge_in : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var areaCell_in : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var areaTriangle_in : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var cellsOnCell_in : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
    var verticesOnCell_in : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
    var verticesOnEdge_in : &int = [&int](c.malloc([sizeof(int)] * nEdges*TWO))
    var edgesOnVertex_in : &int = [&int](c.malloc([sizeof(int)] * nVertices*vertexDegree))
    var cellsOnVertex_in : &int = [&int](c.malloc([sizeof(int)] * nVertices*vertexDegree))
    var kiteAreasOnVertex_in : &double = [&double](c.malloc([sizeof(double)] * nVertices*vertexDegree))


    -- Get the variable IDs of all the variables
    get_varid(ncid, "latCell", &latCell_varid)
    get_varid(ncid, "lonCell", &lonCell_varid)
    get_varid(ncid, "meshDensity", &meshDensity_varid)
    get_varid(ncid, "xCell", &xCell_varid)
    get_varid(ncid, "yCell", &yCell_varid)
    get_varid(ncid, "zCell", &zCell_varid)
    get_varid(ncid, "indexToCellID", &indexToCellID_varid)
    get_varid(ncid, "latEdge", &latEdge_varid)
    get_varid(ncid, "lonEdge", &lonEdge_varid)
    get_varid(ncid, "xEdge", &xEdge_varid)
    get_varid(ncid, "yEdge", &yEdge_varid)
    get_varid(ncid, "zEdge", &zEdge_varid)
    get_varid(ncid, "indexToEdgeID", &indexToEdgeID_varid)
    get_varid(ncid, "latVertex", &latVertex_varid)
    get_varid(ncid, "lonVertex", &lonVertex_varid)
    get_varid(ncid, "xVertex", &xVertex_varid)
    get_varid(ncid, "yVertex", &yVertex_varid)
    get_varid(ncid, "zVertex", &zVertex_varid)
    get_varid(ncid, "indexToVertexID", &indexToVertexID_varid)
    get_varid(ncid, "cellsOnEdge", &cellsOnEdge_varid)
    get_varid(ncid, "cellsOnEdge", &cellsOnEdge_varid)
    get_varid(ncid, "nEdgesOnCell", &nEdgesOnCell_varid)
    get_varid(ncid, "nEdgesOnEdge", &nEdgesOnEdge_varid)
    get_varid(ncid, "edgesOnCell", &edgesOnCell_varid)
    get_varid(ncid, "edgesOnEdge", &edgesOnEdge_varid)
    get_varid(ncid, "weightsOnEdge", &weightsOnEdge_varid)
    get_varid(ncid, "dvEdge", &dvEdge_varid)
    get_varid(ncid, "dv1Edge", &dv1Edge_varid)
    get_varid(ncid, "dv2Edge", &dv2Edge_varid)
    get_varid(ncid, "dcEdge", &dcEdge_varid)
    get_varid(ncid, "angleEdge", &angleEdge_varid)
    get_varid(ncid, "areaCell", &areaCell_varid)
    get_varid(ncid, "areaTriangle", &areaTriangle_varid)
    get_varid(ncid, "cellsOnCell", &cellsOnCell_varid)
    get_varid(ncid, "verticesOnCell", &verticesOnCell_varid)
    get_varid(ncid, "verticesOnEdge", &verticesOnEdge_varid)
    get_varid(ncid, "edgesOnVertex", &edgesOnVertex_varid)
    get_varid(ncid, "cellsOnVertex", &cellsOnVertex_varid)
    get_varid(ncid, "kiteAreasOnVertex", &kiteAreasOnVertex_varid)

    -- Get the variable values, given the variable IDs
    get_var_double(ncid, latCell_varid, latCell_in)
    get_var_double(ncid, lonCell_varid, lonCell_in)
    get_var_double(ncid, meshDensity_varid, meshDensity_in)
    get_var_double(ncid, xCell_varid, xCell_in)
    get_var_double(ncid, yCell_varid, yCell_in)
    get_var_double(ncid, zCell_varid, zCell_in)
    get_var_int(ncid, indexToCellID_varid, indexToCellID_in)
    get_var_double(ncid, latEdge_varid, latEdge_in)
    get_var_double(ncid, lonEdge_varid, lonEdge_in)
    get_var_double(ncid, xEdge_varid, xEdge_in)
    get_var_double(ncid, yEdge_varid, yEdge_in)
    get_var_double(ncid, zEdge_varid, zEdge_in)
    get_var_int(ncid, indexToEdgeID_varid, indexToEdgeID_in)
    get_var_double(ncid, latVertex_varid, latVertex_in)
    get_var_double(ncid, lonVertex_varid, lonVertex_in)
    get_var_double(ncid, xVertex_varid, xVertex_in)
    get_var_double(ncid, yVertex_varid, yVertex_in)
    get_var_double(ncid, zVertex_varid, zVertex_in)
    get_var_int(ncid, indexToVertexID_varid, indexToVertexID_in)
    get_var_int(ncid, cellsOnEdge_varid, cellsOnEdge_in)
    get_var_int(ncid, nEdgesOnCell_varid, nEdgesOnCell_in)
    get_var_int(ncid, nEdgesOnEdge_varid, nEdgesOnEdge_in)
    get_var_int(ncid, edgesOnCell_varid, edgesOnCell_in)
    get_var_int(ncid, edgesOnEdge_varid, edgesOnEdge_in)
    get_var_double(ncid, weightsOnEdge_varid, weightsOnEdge_in)
    get_var_double(ncid, dvEdge_varid, dvEdge_in)
    get_var_double(ncid, dv1Edge_varid, dv1Edge_in)
    get_var_double(ncid, dv2Edge_varid, dv2Edge_in)
    get_var_double(ncid, dcEdge_varid, dcEdge_in)
    get_var_double(ncid, angleEdge_varid, angleEdge_in)
    get_var_double(ncid, areaCell_varid, areaCell_in)
    get_var_double(ncid, areaTriangle_varid, areaTriangle_in)
    get_var_int(ncid, cellsOnCell_varid, cellsOnCell_in)
    get_var_int(ncid, verticesOnCell_varid, verticesOnCell_in)
    get_var_int(ncid, verticesOnEdge_varid, verticesOnEdge_in)
    get_var_int(ncid, edgesOnVertex_varid, edgesOnVertex_in)
    get_var_int(ncid, cellsOnVertex_varid, cellsOnVertex_in)
    get_var_double(ncid, kiteAreasOnVertex_varid, kiteAreasOnVertex_in)

    -------------------------------------------
    ----- DEFINE INDEX SPACES AND REGIONS -----
    -------------------------------------------

    -- Define index spaces for cell IDs, vertex IDs and edge IDs
    var cell_id_space = ispace(int2d, {nCells, nVertLevels + 1})
    var edge_id_space = ispace(int1d, nEdges)
    var vertex_id_space = ispace(int2d, {nVertices, nVertLevels + 1})

    -- Define regions
    var cell_region = region(cell_id_space, cell_fs)
    var edge_region = region(edge_id_space, edge_fs)
    var vertex_region = region(vertex_id_space, vertex_fs)

    var partition_array = read_file(GRAPH_FILE_NAME)

    ----------------------------------
    ----- COPY DATA INTO REGIONS -----
    ----------------------------------

    -- Copy data into cell region
    for i = 0, nCells do
        cell_region[{i, 0}].cellID = indexToCellID_in[i]
        cell_region[{i, 0}].lat = latCell_in[i]
        cell_region[{i, 0}].lon = lonCell_in[i]
        cell_region[{i, 0}].x = xCell_in[i]
        cell_region[{i, 0}].y = yCell_in[i]
        cell_region[{i, 0}].z = zCell_in[i]
        cell_region[{i, 0}].meshDensity = meshDensity_in[i]
        cell_region[{i, 0}].nEdgesOnCell = nEdgesOnCell_in[i]
        cell_region[{i, 0}].areaCell = areaCell_in[i]
        cell_region[{i, 0}].partitionNumber = partition_array[i]

        --cio.printf("Cell : Cell ID %d, partitionNumber %d\n", cell_region[{i, 0}].cellID, cell_region[{i, 0}].partitionNumber)

        for j = 0, maxEdges do
            cell_region[{i, 0}].edgesOnCell[j] = edgesOnCell_in[i*maxEdges + j] --cell_region[{i, 0}].edgesOnCell is a int[maxEdges]
            cell_region[{i, 0}].verticesOnCell[j] = verticesOnCell_in[i*maxEdges + j] --cell_region[{i, 0}].verticesOnCell is a int[maxEdges]
            cell_region[{i, 0}].cellsOnCell[j] = cellsOnCell_in[i*maxEdges + j] --cell_region[{i, 0}].cellsOnCell is a int[maxEdges]
            --cio.printf("edgesOnCell : Cell %d, Edge %d: edge index is %d\n", i, j, cell_region[{i, 0}].edgesOnCell[j])
            --cio.printf("verticesOnCell : Cell %d, Vertex %d: Vertex index is %d\n", i, j, cell_region[{i, 0}].verticesOnCell[j])
            --cio.printf("cellsOnCell : InnerCell %d, OuterCell %d: Cell index is %d\n", i, j, cell_region[{i, 0}].cellsOnCell[j])
        end
        --cio.printf("Cell : Cell ID %d, nEdgesOnCell is %d\n", cell_region[{i, 0}].cellID, cell_region[{i, 0}].nEdgesOnCell)
    end

    -- Copy data into edge region
    for i = 0, nEdges do
        edge_region[{i, 0}].edgeID = indexToEdgeID_in[i]
        edge_region[{i, 0}].lat = latEdge_in[i]
        edge_region[{i, 0}].lon = lonEdge_in[i]
        edge_region[{i, 0}].x = xEdge_in[i]
        edge_region[{i, 0}].y = yEdge_in[i]
        edge_region[{i, 0}].z = zEdge_in[i]
        edge_region[{i, 0}].nEdgesOnEdge = nEdgesOnEdge_in[i]
        edge_region[{i, 0}].angleEdge = angleEdge_in[i]
        edge_region[{i, 0}].dvEdge = dvEdge_in[i]
        edge_region[{i, 0}].dv1Edge = dv1Edge_in[i]
        edge_region[{i, 0}].dv2Edge = dv2Edge_in[i]
        edge_region[{i, 0}].dcEdge = dcEdge_in[i]


        for j = 0, TWO do
            edge_region[{i, 0}].cellsOnEdge[j] = cellsOnEdge_in[i*TWO + j]
            edge_region[{i, 0}].verticesOnEdge[j] = verticesOnEdge_in[i*TWO + j]
            --cio.printf("cellsOnEdge : Edge %d, Cell %d is %d\n", i, j, edge_region[{i, 0}].cellsOnEdge[j])
            --cio.printf("VerticesOnEdge : Edge %d: Vertex %d is $d\n", i, j, edge_region[{i, 0}].verticesOnEdge[j])
        end

        for j = 0, maxEdges2 do
            edge_region[{i, 0}].edgesOnEdge_ECP[j] = edgesOnEdge_in[i*maxEdges2 + j]
            edge_region[{i, 0}].weightsOnEdge[j] = weightsOnEdge_in[i*maxEdges2 + j]
            --cio.printf("edgesOnEdge_ECP : InnerEdge %d, OuterEdge %d is %d\n", i, j, edge_region[{i, 0}].edgesOnEdge_ECP[j])
            --cio.printf("weightsOnEdge : Edge %d: Weight %d is $f\n", i, j, edge_region[{i, 0}].weightsOnEdge[j])
        end
        --cio.printf("Edge: ID is %d, xEdge is %f, yEdge is %f, zEdge is %f \n", i, edge_region[{i, 0}].x, edge_region[{i, 0}].y, edge_region[{i, 0}].z)
    end

    -- Copy data into vertex region
    for i = 0, nVertices do
        vertex_region[i].vertexID = indexToVertexID_in[i]
        vertex_region[i].lat = latVertex_in[i]
        vertex_region[i].lon = lonVertex_in[i]
        vertex_region[i].x = xVertex_in[i]
        vertex_region[i].y = yVertex_in[i]
        vertex_region[i].z = zVertex_in[i]
        vertex_region[i].areaTriangle = areaTriangle_in[i]

        for j = 0, vertexDegree do
            vertex_region[i].edgesOnVertex[j] = edgesOnVertex_in[i*vertexDegree + j]
            vertex_region[i].cellsOnVertex[j] = cellsOnVertex_in[i*vertexDegree + j]
            vertex_region[i].kiteAreasOnVertex[j] = kiteAreasOnVertex_in[i*vertexDegree + j]

            --cio.printf("edgesOnVertex : Vertex %d, Edge %d: Edge index is %d\n", i, j, vertex_region[i].edgesOnVertex[j])
            --cio.printf("cellsOnVertex : Vertex %d, Cell %d: Cell index is %d\n", i, j, vertex_region[i].cellsOnVertex[j])
            --cio.printf("kiteAreasOnVertex : Vertex %d, Kite %d: Kite Area is %f\n", i, j, vertex_region[i].kiteAreasOnVertex[j])
        end
        --cio.printf("Vertex ID is %d, xVertex is %f, yVertex is %f, zVertex is %f \n", i, vertex_region[i].x, vertex_region[i].y, vertex_region[i].z)
    end

    -------------------------
    ----- CALCULATE EVC -----
    -------------------------
    --I know I should do something more intelligent to get the common elements: but for now we do a brute force search to get EVC

    --First, we iterate through the cells and get the edgesOnCell array for each cell
    for i = 0, nCells do
        var curr_edgesOnCell = cell_region[{i, 0}].edgesOnCell

    --Then we iterate through the vertices of that cell
        for j = 0, maxEdges do
            var currVertexID = cell_region[{i, 0}].verticesOnCell[j]
            cell_region[{i,0}].evc[j*3] = currVertexID
            --cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3, cell_region[{i, 0}].evc[j*3])

            if currVertexID == 0 then
                cell_region[{i,0}].evc[j*3 + 1] = 0
                cell_region[{i,0}].evc[j*3 + 2] = 0
                --cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3 + 1, cell_region[{i, 0}].evc[j*3 + 1])
                --cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3 + 2, cell_region[{i, 0}].evc[j*3 + 2])

            --If there is a vertex, we get the edges on that vertex
            elseif currVertexID ~= 0 then
                var curr_edgesOnVertex = vertex_region[currVertexID-1].edgesOnVertex
                var count = 1

                --Then, we get overlapping edges between curr_edgesOnVertex and curr_edgesOnCell to get EVC
                for k = 0, vertexDegree do
                    var currEdgeID = curr_edgesOnVertex[k]
                    for l = 0, maxEdges do
                        if currEdgeID == curr_edgesOnCell[l] and count < 3 then
                            cell_region[{i,0}].evc[j*3 + count] = currEdgeID
                            --cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3 + count, cell_region[{i, 0}].evc[j*3 + count])
                            count = count+1
                        end
                    end
                end
            end
        end
    end

     -- Close the file
	  file_close(ncid)

    -- Free allocated arrays
    c.free(latCell_in)
    c.free(lonCell_in)
    c.free(meshDensity_in)
    c.free(xCell_in)
    c.free(yCell_in)
    c.free(zCell_in)
    c.free(indexToCellID_in)
    c.free(latEdge_in)
    c.free(lonEdge_in)
    c.free(xEdge_in)
    c.free(yEdge_in)
    c.free(zEdge_in)
    c.free(indexToEdgeID_in)
    c.free(latVertex_in)
    c.free(lonVertex_in)
    c.free(xVertex_in)
    c.free(yVertex_in)
    c.free(zVertex_in)
    c.free(indexToVertexID_in)
    c.free(cellsOnEdge_in)
    c.free(nEdgesOnCell_in)
    c.free(nEdgesOnEdge_in)
    c.free(edgesOnCell_in)
    c.free(edgesOnEdge_in)
    c.free(weightsOnEdge_in)
    c.free(dvEdge_in)
    c.free(dv1Edge_in)
    c.free(dv2Edge_in)
    c.free(dcEdge_in)
    c.free(angleEdge_in)
    c.free(areaCell_in)
    c.free(areaTriangle_in)
    c.free(cellsOnCell_in)
    c.free(verticesOnCell_in)
    c.free(verticesOnEdge_in)
    c.free(edgesOnVertex_in)
    c.free(cellsOnVertex_in)
    c.free(kiteAreasOnVertex_in)

    cio.printf("Successfully read file! \n")

    -----------------------------------
    ----- Copy Neighbours -----
    -----------------------------------
    for i = 0, nCells do
      -- cell.cellsOnCell[0] contains an integer with the index of the cell neighbour. I would like cell.neighbour to point to that cell in the region
      -- we subtract 1 because the index spaces are 0-indexed but the cellIDs are 1-indexed
      cell_region[{i, 0}].neighbor0 = cell_region[{i, 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor1 = cell_region[{i, 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor2 = cell_region[{i, 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor3 = cell_region[{i, 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor4 = cell_region[{i, 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor5 = cell_region[{i, 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor6 = cell_region[{i, 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor7 = cell_region[{i, 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor8 = cell_region[{i, 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor9 = cell_region[{i, 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor00 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor01 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor02 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor03 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor04 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor05 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor06 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor07 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor08 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor09 = cell_region[{cell_region[{i, 0}].cellsOnCell[0], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor10 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor11 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor12 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor13 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor14 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor15 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor16 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor17 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor18 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor19 = cell_region[{cell_region[{i, 0}].cellsOnCell[1], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor20 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor21 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor22 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor23 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor24 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor25 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor26 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor27 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor28 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor29 = cell_region[{cell_region[{i, 0}].cellsOnCell[2], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor30 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor31 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor32 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor33 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor34 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor35 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor36 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor37 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor38 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor39 = cell_region[{cell_region[{i, 0}].cellsOnCell[3], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor40 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor41 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor42 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor43 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor44 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor45 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor46 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor47 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor48 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor49 = cell_region[{cell_region[{i, 0}].cellsOnCell[4], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor50 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor51 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor52 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor53 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor54 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor55 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor56 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor57 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor58 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor59 = cell_region[{cell_region[{i, 0}].cellsOnCell[5], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor60 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor61 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor62 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor63 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor64 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor65 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor66 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor67 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor68 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor69 = cell_region[{cell_region[{i, 0}].cellsOnCell[6], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor70 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor71 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor72 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor73 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor74 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor75 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor76 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor77 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor78 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor79 = cell_region[{cell_region[{i, 0}].cellsOnCell[7], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor80 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor81 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor82 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor83 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor84 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor85 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor86 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor87 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor88 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor89 = cell_region[{cell_region[{i, 0}].cellsOnCell[8], 0}].cellsOnCell[9] - 1

      cell_region[{i, 0}].neighbor90 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[0] - 1
      cell_region[{i, 0}].neighbor91 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[1] - 1
      cell_region[{i, 0}].neighbor92 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[2] - 1
      cell_region[{i, 0}].neighbor93 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[3] - 1
      cell_region[{i, 0}].neighbor94 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[4] - 1
      cell_region[{i, 0}].neighbor95 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[5] - 1
      cell_region[{i, 0}].neighbor96 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[6] - 1
      cell_region[{i, 0}].neighbor97 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[7] - 1
      cell_region[{i, 0}].neighbor98 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[8] - 1
      cell_region[{i, 0}].neighbor99 = cell_region[{cell_region[{i, 0}].cellsOnCell[9], 0}].cellsOnCell[9] - 1
    end

    -----------------------------------------------
    ------- PARTITION REGIONS  --------
    -----------------------------------------------

    --Create initial partition based on METIS graph partition
    var color_space = ispace(int1d, NUM_PARTITIONS)
    var cell_partition_initial = partition(complete, cell_region.partitionNumber, color_space)

    --Create dependent partitions based on neighbour fields
    -- Syntax for image: var p = image(parent_region, source_partition, data_region.field)
    -- cell_partition_i = {cell \in cell_region | cell_region[j].neighbour0 matches that cell index and cell_region[j] is in subregion i of pc_initial}
    var cell_partition_neighbor0 = image(cell_region, cell_partition_initial, cell_region.neighbor0)
    var cell_partition_neighbor1 = image(cell_region, cell_partition_initial, cell_region.neighbor1)
    var cell_partition_neighbor2 = image(cell_region, cell_partition_initial, cell_region.neighbor2)
    var cell_partition_neighbor3 = image(cell_region, cell_partition_initial, cell_region.neighbor3)
    var cell_partition_neighbor4 = image(cell_region, cell_partition_initial, cell_region.neighbor4)
    var cell_partition_neighbor5 = image(cell_region, cell_partition_initial, cell_region.neighbor5)
    var cell_partition_neighbor6 = image(cell_region, cell_partition_initial, cell_region.neighbor6)
    var cell_partition_neighbor7 = image(cell_region, cell_partition_initial, cell_region.neighbor7)
    var cell_partition_neighbor8 = image(cell_region, cell_partition_initial, cell_region.neighbor8)
    var cell_partition_neighbor9 = image(cell_region, cell_partition_initial, cell_region.neighbor9)

    -- For second halo
    var cell_partition_neighbor00 = image(cell_region, cell_partition_initial, cell_region.neighbor00)
    var cell_partition_neighbor01 = image(cell_region, cell_partition_initial, cell_region.neighbor01)
    var cell_partition_neighbor02 = image(cell_region, cell_partition_initial, cell_region.neighbor02)
    var cell_partition_neighbor03 = image(cell_region, cell_partition_initial, cell_region.neighbor03)
    var cell_partition_neighbor04 = image(cell_region, cell_partition_initial, cell_region.neighbor04)
    var cell_partition_neighbor05 = image(cell_region, cell_partition_initial, cell_region.neighbor05)
    var cell_partition_neighbor06 = image(cell_region, cell_partition_initial, cell_region.neighbor06)
    var cell_partition_neighbor07 = image(cell_region, cell_partition_initial, cell_region.neighbor07)
    var cell_partition_neighbor08 = image(cell_region, cell_partition_initial, cell_region.neighbor08)
    var cell_partition_neighbor09 = image(cell_region, cell_partition_initial, cell_region.neighbor09)

    var cell_partition_neighbor10 = image(cell_region, cell_partition_initial, cell_region.neighbor10)
    var cell_partition_neighbor11 = image(cell_region, cell_partition_initial, cell_region.neighbor11)
    var cell_partition_neighbor12 = image(cell_region, cell_partition_initial, cell_region.neighbor12)
    var cell_partition_neighbor13 = image(cell_region, cell_partition_initial, cell_region.neighbor13)
    var cell_partition_neighbor14 = image(cell_region, cell_partition_initial, cell_region.neighbor14)
    var cell_partition_neighbor15 = image(cell_region, cell_partition_initial, cell_region.neighbor15)
    var cell_partition_neighbor16 = image(cell_region, cell_partition_initial, cell_region.neighbor16)
    var cell_partition_neighbor17 = image(cell_region, cell_partition_initial, cell_region.neighbor17)
    var cell_partition_neighbor18 = image(cell_region, cell_partition_initial, cell_region.neighbor18)
    var cell_partition_neighbor19 = image(cell_region, cell_partition_initial, cell_region.neighbor19)

    var cell_partition_neighbor20 = image(cell_region, cell_partition_initial, cell_region.neighbor20)
    var cell_partition_neighbor21 = image(cell_region, cell_partition_initial, cell_region.neighbor21)
    var cell_partition_neighbor22 = image(cell_region, cell_partition_initial, cell_region.neighbor22)
    var cell_partition_neighbor23 = image(cell_region, cell_partition_initial, cell_region.neighbor23)
    var cell_partition_neighbor24 = image(cell_region, cell_partition_initial, cell_region.neighbor24)
    var cell_partition_neighbor25 = image(cell_region, cell_partition_initial, cell_region.neighbor25)
    var cell_partition_neighbor26 = image(cell_region, cell_partition_initial, cell_region.neighbor26)
    var cell_partition_neighbor27 = image(cell_region, cell_partition_initial, cell_region.neighbor27)
    var cell_partition_neighbor28 = image(cell_region, cell_partition_initial, cell_region.neighbor28)
    var cell_partition_neighbor29 = image(cell_region, cell_partition_initial, cell_region.neighbor29)

    var cell_partition_neighbor30 = image(cell_region, cell_partition_initial, cell_region.neighbor30)
    var cell_partition_neighbor31 = image(cell_region, cell_partition_initial, cell_region.neighbor31)
    var cell_partition_neighbor32 = image(cell_region, cell_partition_initial, cell_region.neighbor32)
    var cell_partition_neighbor33 = image(cell_region, cell_partition_initial, cell_region.neighbor33)
    var cell_partition_neighbor34 = image(cell_region, cell_partition_initial, cell_region.neighbor34)
    var cell_partition_neighbor35 = image(cell_region, cell_partition_initial, cell_region.neighbor35)
    var cell_partition_neighbor36 = image(cell_region, cell_partition_initial, cell_region.neighbor36)
    var cell_partition_neighbor37 = image(cell_region, cell_partition_initial, cell_region.neighbor37)
    var cell_partition_neighbor38 = image(cell_region, cell_partition_initial, cell_region.neighbor38)
    var cell_partition_neighbor39 = image(cell_region, cell_partition_initial, cell_region.neighbor39)

    var cell_partition_neighbor40 = image(cell_region, cell_partition_initial, cell_region.neighbor40)
    var cell_partition_neighbor41 = image(cell_region, cell_partition_initial, cell_region.neighbor41)
    var cell_partition_neighbor42 = image(cell_region, cell_partition_initial, cell_region.neighbor42)
    var cell_partition_neighbor43 = image(cell_region, cell_partition_initial, cell_region.neighbor43)
    var cell_partition_neighbor44 = image(cell_region, cell_partition_initial, cell_region.neighbor44)
    var cell_partition_neighbor45 = image(cell_region, cell_partition_initial, cell_region.neighbor45)
    var cell_partition_neighbor46 = image(cell_region, cell_partition_initial, cell_region.neighbor46)
    var cell_partition_neighbor47 = image(cell_region, cell_partition_initial, cell_region.neighbor47)
    var cell_partition_neighbor48 = image(cell_region, cell_partition_initial, cell_region.neighbor48)
    var cell_partition_neighbor49 = image(cell_region, cell_partition_initial, cell_region.neighbor49)

    var cell_partition_neighbor50 = image(cell_region, cell_partition_initial, cell_region.neighbor50)
    var cell_partition_neighbor51 = image(cell_region, cell_partition_initial, cell_region.neighbor51)
    var cell_partition_neighbor52 = image(cell_region, cell_partition_initial, cell_region.neighbor52)
    var cell_partition_neighbor53 = image(cell_region, cell_partition_initial, cell_region.neighbor53)
    var cell_partition_neighbor54 = image(cell_region, cell_partition_initial, cell_region.neighbor54)
    var cell_partition_neighbor55 = image(cell_region, cell_partition_initial, cell_region.neighbor55)
    var cell_partition_neighbor56 = image(cell_region, cell_partition_initial, cell_region.neighbor56)
    var cell_partition_neighbor57 = image(cell_region, cell_partition_initial, cell_region.neighbor57)
    var cell_partition_neighbor58 = image(cell_region, cell_partition_initial, cell_region.neighbor58)
    var cell_partition_neighbor59 = image(cell_region, cell_partition_initial, cell_region.neighbor59)

    var cell_partition_neighbor60 = image(cell_region, cell_partition_initial, cell_region.neighbor60)
    var cell_partition_neighbor61 = image(cell_region, cell_partition_initial, cell_region.neighbor61)
    var cell_partition_neighbor62 = image(cell_region, cell_partition_initial, cell_region.neighbor62)
    var cell_partition_neighbor63 = image(cell_region, cell_partition_initial, cell_region.neighbor63)
    var cell_partition_neighbor64 = image(cell_region, cell_partition_initial, cell_region.neighbor64)
    var cell_partition_neighbor65 = image(cell_region, cell_partition_initial, cell_region.neighbor65)
    var cell_partition_neighbor66 = image(cell_region, cell_partition_initial, cell_region.neighbor66)
    var cell_partition_neighbor67 = image(cell_region, cell_partition_initial, cell_region.neighbor67)
    var cell_partition_neighbor68 = image(cell_region, cell_partition_initial, cell_region.neighbor68)
    var cell_partition_neighbor69 = image(cell_region, cell_partition_initial, cell_region.neighbor69)

    var cell_partition_neighbor70 = image(cell_region, cell_partition_initial, cell_region.neighbor70)
    var cell_partition_neighbor71 = image(cell_region, cell_partition_initial, cell_region.neighbor71)
    var cell_partition_neighbor72 = image(cell_region, cell_partition_initial, cell_region.neighbor72)
    var cell_partition_neighbor73 = image(cell_region, cell_partition_initial, cell_region.neighbor73)
    var cell_partition_neighbor74 = image(cell_region, cell_partition_initial, cell_region.neighbor74)
    var cell_partition_neighbor75 = image(cell_region, cell_partition_initial, cell_region.neighbor75)
    var cell_partition_neighbor76 = image(cell_region, cell_partition_initial, cell_region.neighbor76)
    var cell_partition_neighbor77 = image(cell_region, cell_partition_initial, cell_region.neighbor77)
    var cell_partition_neighbor78 = image(cell_region, cell_partition_initial, cell_region.neighbor78)
    var cell_partition_neighbor79 = image(cell_region, cell_partition_initial, cell_region.neighbor79)

    var cell_partition_neighbor80 = image(cell_region, cell_partition_initial, cell_region.neighbor80)
    var cell_partition_neighbor81 = image(cell_region, cell_partition_initial, cell_region.neighbor81)
    var cell_partition_neighbor82 = image(cell_region, cell_partition_initial, cell_region.neighbor82)
    var cell_partition_neighbor83 = image(cell_region, cell_partition_initial, cell_region.neighbor83)
    var cell_partition_neighbor84 = image(cell_region, cell_partition_initial, cell_region.neighbor84)
    var cell_partition_neighbor85 = image(cell_region, cell_partition_initial, cell_region.neighbor85)
    var cell_partition_neighbor86 = image(cell_region, cell_partition_initial, cell_region.neighbor86)
    var cell_partition_neighbor87 = image(cell_region, cell_partition_initial, cell_region.neighbor87)
    var cell_partition_neighbor88 = image(cell_region, cell_partition_initial, cell_region.neighbor88)
    var cell_partition_neighbor89 = image(cell_region, cell_partition_initial, cell_region.neighbor89)

    var cell_partition_neighbor90 = image(cell_region, cell_partition_initial, cell_region.neighbor90)
    var cell_partition_neighbor91 = image(cell_region, cell_partition_initial, cell_region.neighbor91)
    var cell_partition_neighbor92 = image(cell_region, cell_partition_initial, cell_region.neighbor92)
    var cell_partition_neighbor93 = image(cell_region, cell_partition_initial, cell_region.neighbor93)
    var cell_partition_neighbor94 = image(cell_region, cell_partition_initial, cell_region.neighbor94)
    var cell_partition_neighbor95 = image(cell_region, cell_partition_initial, cell_region.neighbor95)
    var cell_partition_neighbor96 = image(cell_region, cell_partition_initial, cell_region.neighbor96)
    var cell_partition_neighbor97 = image(cell_region, cell_partition_initial, cell_region.neighbor97)
    var cell_partition_neighbor98 = image(cell_region, cell_partition_initial, cell_region.neighbor98)
    var cell_partition_neighbor99 = image(cell_region, cell_partition_initial, cell_region.neighbor99)

    var partition_s2_0 = cell_partition_neighbor00 | cell_partition_neighbor01 | cell_partition_neighbor02 | cell_partition_neighbor03 | cell_partition_neighbor04 | cell_partition_neighbor05 | cell_partition_neighbor06 | cell_partition_neighbor07 | cell_partition_neighbor08 | cell_partition_neighbor09
    var partition_s2_1 = cell_partition_neighbor10 | cell_partition_neighbor11 | cell_partition_neighbor12 | cell_partition_neighbor13 | cell_partition_neighbor14 | cell_partition_neighbor15 | cell_partition_neighbor16 | cell_partition_neighbor17 | cell_partition_neighbor18 | cell_partition_neighbor19
    var partition_s2_2 = cell_partition_neighbor20 | cell_partition_neighbor21 | cell_partition_neighbor22 | cell_partition_neighbor23 | cell_partition_neighbor24 | cell_partition_neighbor25 | cell_partition_neighbor26 | cell_partition_neighbor27 | cell_partition_neighbor28 | cell_partition_neighbor29
    var partition_s2_3 = cell_partition_neighbor30 | cell_partition_neighbor31 | cell_partition_neighbor32 | cell_partition_neighbor33 | cell_partition_neighbor34 | cell_partition_neighbor35 | cell_partition_neighbor36 | cell_partition_neighbor37 | cell_partition_neighbor38 | cell_partition_neighbor39
    var partition_s2_4 = cell_partition_neighbor40 | cell_partition_neighbor41 | cell_partition_neighbor42 | cell_partition_neighbor43 | cell_partition_neighbor44 | cell_partition_neighbor45 | cell_partition_neighbor46 | cell_partition_neighbor47 | cell_partition_neighbor48 | cell_partition_neighbor49
    var partition_s2_5 = cell_partition_neighbor50 | cell_partition_neighbor51 | cell_partition_neighbor52 | cell_partition_neighbor53 | cell_partition_neighbor54 | cell_partition_neighbor55 | cell_partition_neighbor56 | cell_partition_neighbor57 | cell_partition_neighbor58 | cell_partition_neighbor59
    var partition_s2_6 = cell_partition_neighbor60 | cell_partition_neighbor61 | cell_partition_neighbor62 | cell_partition_neighbor63 | cell_partition_neighbor64 | cell_partition_neighbor65 | cell_partition_neighbor66 | cell_partition_neighbor67 | cell_partition_neighbor68 | cell_partition_neighbor69
    var partition_s2_7 = cell_partition_neighbor70 | cell_partition_neighbor71 | cell_partition_neighbor72 | cell_partition_neighbor73 | cell_partition_neighbor74 | cell_partition_neighbor75 | cell_partition_neighbor76 | cell_partition_neighbor77 | cell_partition_neighbor78 | cell_partition_neighbor79
    var partition_s2_8 = cell_partition_neighbor80 | cell_partition_neighbor81 | cell_partition_neighbor82 | cell_partition_neighbor83 | cell_partition_neighbor84 | cell_partition_neighbor85 | cell_partition_neighbor86 | cell_partition_neighbor87 | cell_partition_neighbor88 | cell_partition_neighbor89
    var partition_s2_9 = cell_partition_neighbor90 | cell_partition_neighbor91 | cell_partition_neighbor92 | cell_partition_neighbor93 | cell_partition_neighbor94 | cell_partition_neighbor95 | cell_partition_neighbor96 | cell_partition_neighbor97 | cell_partition_neighbor98 | cell_partition_neighbor99
    var partition_s2 = partition_s2_0 | partition_s2_1 | partition_s2_2 | partition_s2_3 | partition_s2_4 | partition_s2_5 | partition_s2_6 | partition_s2_7 | partition_s2_8 | partition_s2_9

    -- Construct halos
    var partition_s_1 = cell_partition_neighbor0 | cell_partition_neighbor1 | cell_partition_neighbor2 | cell_partition_neighbor3 | cell_partition_neighbor4 | cell_partition_neighbor5 | cell_partition_neighbor6 | cell_partition_neighbor7 | cell_partition_neighbor8 | cell_partition_neighbor9
    var partition_halo_1 = partition_s_1 - cell_partition_initial
    var partition_halo_2 = partition_s2  - partition_halo_1 - cell_partition_initial

    --Test code by printing out first neighbours in the original partition
    --var i = 0
    --for p in color_space do
    --    var sub_region = cell_partition_initial[p]
    --    cio.printf("Sub region %d\n", i)
    --    for cell in sub_region do
    --        cio.printf("%d\n", cell.cellsOnCell[0])
    --    end
    --    i=i+1
    --end

    --Print out first neighbour partition to check against the original partition
    --i = 0
    --for p in cell_partition_neighbor0.colors do
    --    var sub_region = cell_partition_neighbor0[p]
    --    cio.printf("Sub region %d\n", i)
    --    for cell in sub_region do
    --        cio.printf("%d\n", cell.cellID)
    --    end
    --    i=i+1
    --end


    ----------------------------------------------------
    ------- TESTING CODE: WRITING NETCDF OUTPUT --------
    ----------------------------------------------------

    -- We create a netcdf file using the data in the regions, to test whether the data was written correctly.
    cio.printf("Starting to write netcdf file..\n")
    var ncid_copy = 65537

    --Create a netcdf file
    file_create("newfile.nc", &ncid_copy)

    --Initialize the file's dimension variables
    var nCells_dimid_copy : int
    var nEdges_dimid_copy : int
    var nVertices_dimid_copy : int
    var maxEdges_dimid_copy : int
    var maxEdges2_dimid_copy : int
    var TWO_dimid_copy : int
    var vertexDegree_dimid_copy : int
    var nVertLevels_dimid_copy : int
    var time_dimid_copy : int

    --Define the dimension variables
    --define_dim(ncid: int, dim_name: &int, dim_size: int, dim_id_ptr: &int)
    define_dim(ncid_copy, "nCells", nCells, &nCells_dimid_copy)
    define_dim(ncid_copy, "nEdges", nEdges, &nEdges_dimid_copy)
    define_dim(ncid_copy, "nVertices", nVertices, &nVertices_dimid_copy)
    define_dim(ncid_copy, "maxEdges", maxEdges, &maxEdges_dimid_copy)
    define_dim(ncid_copy, "maxEdges2", maxEdges2, &maxEdges2_dimid_copy)
    define_dim(ncid_copy, "TWO", TWO, &TWO_dimid_copy)
    define_dim(ncid_copy, "vertexDegree", vertexDegree, &vertexDegree_dimid_copy)
    define_dim(ncid_copy, "nVertLevels", nVertLevels, &nVertLevels_dimid_copy)
    define_dim(ncid_copy, "Time", netcdf.NC_UNLIMITED, &time_dimid_copy)

    --For the 2D variables, the dimIDs need to be put in arrays
    var nEdges_TWO_dimids = array(nEdges_dimid_copy, TWO_dimid_copy)
    var nCells_maxEdges_dimids = array(nCells_dimid_copy, maxEdges_dimid_copy)
    var nEdges_maxEdges2_dimids = array(nEdges_dimid_copy, maxEdges2_dimid_copy)
    var nVertices_vertexDegree_dimids = array(nVertices_dimid_copy, vertexDegree_dimid_copy)

    --Initialize the variable IDs
    var latCell_varid_copy : int
    var lonCell_varid_copy : int
    var meshDensity_varid_copy : int
    var xCell_varid_copy : int
    var yCell_varid_copy : int
    var zCell_varid_copy : int
    var indexToCellID_varid_copy : int
    var latEdge_varid_copy : int
    var lonEdge_varid_copy : int
    var xEdge_varid_copy : int
    var yEdge_varid_copy : int
    var zEdge_varid_copy : int
    var indexToEdgeID_varid_copy : int
    var latVertex_varid_copy : int
    var lonVertex_varid_copy : int
    var xVertex_varid_copy : int
    var yVertex_varid_copy : int
    var zVertex_varid_copy : int
    var indexToVertexID_varid_copy : int
    var cellsOnEdge_varid_copy : int
    var nEdgesOnCell_varid_copy : int
    var nEdgesOnEdge_varid_copy : int
    var edgesOnCell_varid_copy : int
    var edgesOnEdge_varid_copy : int
    var weightsOnEdge_varid_copy : int
    var dvEdge_varid_copy : int
    var dv1Edge_varid_copy : int
    var dv2Edge_varid_copy : int
    var dcEdge_varid_copy : int
    var angleEdge_varid_copy : int
    var areaCell_varid_copy : int
    var areaTriangle_varid_copy : int
    var cellsOnCell_varid_copy : int
    var verticesOnCell_varid_copy : int
    var verticesOnEdge_varid_copy : int
    var edgesOnVertex_varid_copy : int
    var cellsOnVertex_varid_copy : int
    var kiteAreasOnVertex_varid_copy : int

    --Define the variable IDs
    define_var(ncid_copy, "latCell", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &latCell_varid_copy)
    define_var(ncid_copy, "lonCell", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &lonCell_varid_copy)
    define_var(ncid_copy, "meshDensity", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &meshDensity_varid_copy)
    define_var(ncid_copy, "xCell", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &xCell_varid_copy)
    define_var(ncid_copy, "yCell", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &yCell_varid_copy)
    define_var(ncid_copy, "zCell", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &zCell_varid_copy)
    define_var(ncid_copy, "indexToCellID", netcdf.NC_INT, 1, &nCells_dimid_copy, &indexToCellID_varid_copy)
    define_var(ncid_copy, "latEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &latEdge_varid_copy)
    define_var(ncid_copy, "lonEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &lonEdge_varid_copy)
    define_var(ncid_copy, "xEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &xEdge_varid_copy)
    define_var(ncid_copy, "yEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &yEdge_varid_copy)
    define_var(ncid_copy, "zEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &zEdge_varid_copy)
    define_var(ncid_copy, "indexToEdgeID", netcdf.NC_INT, 1, &nEdges_dimid_copy, &indexToEdgeID_varid_copy)
    define_var(ncid_copy, "latVertex", netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &latVertex_varid_copy)
    define_var(ncid_copy, "lonVertex", netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &lonVertex_varid_copy)
    define_var(ncid_copy, "xVertex", netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &xVertex_varid_copy)
    define_var(ncid_copy, "yVertex", netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &yVertex_varid_copy)
    define_var(ncid_copy, "zVertex", netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &zVertex_varid_copy)
    define_var(ncid_copy, "indexToVertexID", netcdf.NC_INT, 1, &nVertices_dimid_copy, &indexToVertexID_varid_copy)
    define_var(ncid_copy, "cellsOnEdge", netcdf.NC_INT, 2, nEdges_TWO_dimids, &cellsOnEdge_varid_copy)
    define_var(ncid_copy, "nEdgesOnCell", netcdf.NC_INT, 1, &nCells_dimid_copy, &nEdgesOnCell_varid_copy)
    define_var(ncid_copy, "nEdgesOnEdge", netcdf.NC_INT, 1, &nEdges_dimid_copy, &nEdgesOnEdge_varid_copy)
    define_var(ncid_copy, "edgesOnCell", netcdf.NC_INT, 2, nCells_maxEdges_dimids, &edgesOnCell_varid_copy)
    define_var(ncid_copy, "edgesOnEdge", netcdf.NC_INT, 2, nEdges_maxEdges2_dimids, &edgesOnEdge_varid_copy)
    define_var(ncid_copy, "weightsOnEdge", netcdf.NC_DOUBLE, 2, nEdges_maxEdges2_dimids, &weightsOnEdge_varid_copy)
    define_var(ncid_copy, "dvEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dvEdge_varid_copy)
    define_var(ncid_copy, "dv1Edge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dv1Edge_varid_copy)
    define_var(ncid_copy, "dv2Edge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dv2Edge_varid_copy)
    define_var(ncid_copy, "dcEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dcEdge_varid_copy)
    define_var(ncid_copy, "angleEdge", netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &angleEdge_varid_copy)
    define_var(ncid_copy, "areaCell", netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &areaCell_varid_copy)
    define_var(ncid_copy, "areaTriangle", netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &areaTriangle_varid_copy)
    define_var(ncid_copy, "cellsOnCell", netcdf.NC_INT, 2, nCells_maxEdges_dimids, &cellsOnCell_varid_copy)
    define_var(ncid_copy, "verticesOnCell", netcdf.NC_INT, 2, nCells_maxEdges_dimids, &verticesOnCell_varid_copy)
    define_var(ncid_copy, "verticesOnEdge", netcdf.NC_INT, 2, nEdges_TWO_dimids, &verticesOnEdge_varid_copy)
    define_var(ncid_copy, "edgesOnVertex", netcdf.NC_INT, 2, nVertices_vertexDegree_dimids, &edgesOnVertex_varid_copy)
    define_var(ncid_copy, "cellsOnVertex", netcdf.NC_INT, 2, nVertices_vertexDegree_dimids, &cellsOnVertex_varid_copy)
    define_var(ncid_copy, "kiteAreasOnVertex", netcdf.NC_DOUBLE, 2, nVertices_vertexDegree_dimids, &kiteAreasOnVertex_varid_copy)

    --This function signals that we're done writing the metadata.
    end_def(ncid_copy)

    --Now define the new arrays to hold the data that will be put in the netcdf files
    var latCell_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var lonCell_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var meshDensity_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var xCell_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var yCell_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var zCell_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var indexToCellID_in_copy : &int = [&int](c.malloc([sizeof(int)] * nCells))
    var latEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var lonEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var xEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var yEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var zEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var indexToEdgeID_in_copy : &int = [&int](c.malloc([sizeof(int)] * nEdges))
    var latVertex_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var lonVertex_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var xVertex_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var yVertex_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var zVertex_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var indexToVertexID_in_copy : &int = [&int](c.malloc([sizeof(int)] * nVertices))
    var cellsOnEdge_in_copy : &int = [&int](c.malloc([sizeof(int)] * nEdges*TWO))
    var nEdgesOnCell_in_copy : &int = [&int](c.malloc([sizeof(int)] * nCells))
    var nEdgesOnEdge_in_copy : &int = [&int](c.malloc([sizeof(int)] * nEdges))
    var edgesOnCell_in_copy : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
    var edgesOnEdge_in_copy : &int = [&int](c.malloc([sizeof(int)] * nEdges*maxEdges2))
    var weightsOnEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges*maxEdges2))
    var dvEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var dv1Edge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var dv2Edge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var dcEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var angleEdge_in_copy : &double = [&double](c.malloc([sizeof(double)] * nEdges))
    var areaCell_in_copy : &double = [&double](c.malloc([sizeof(double)] * nCells))
    var areaTriangle_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices))
    var cellsOnCell_in_copy : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
    var verticesOnCell_in_copy : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
    var verticesOnEdge_in_copy : &int = [&int](c.malloc([sizeof(int)] * nEdges*TWO))
    var edgesOnVertex_in_copy : &int = [&int](c.malloc([sizeof(int)] * nVertices*vertexDegree))
    var cellsOnVertex_in_copy : &int = [&int](c.malloc([sizeof(int)] * nVertices*vertexDegree))
    var kiteAreasOnVertex_in_copy : &double = [&double](c.malloc([sizeof(double)] * nVertices*vertexDegree))

    --Now we copy the data into the arrays so they can be read into the netcdf files
    for i = 0, nCells do
        latCell_in_copy[i] = cell_region[{i, 0}].lat
        lonCell_in_copy[i] = cell_region[{i, 0}].lon
        xCell_in_copy[i] = cell_region[{i, 0}].x
        yCell_in_copy[i] = cell_region[{i, 0}].y
        zCell_in_copy[i] = cell_region[{i, 0}].z
        indexToCellID_in_copy[i] = cell_region[{i, 0}].cellID
        meshDensity_in_copy[i] = cell_region[{i, 0}].meshDensity
        nEdgesOnCell_in_copy[i] = cell_region[{i, 0}].nEdgesOnCell
        areaCell_in_copy[i] = cell_region[{i, 0}].areaCell

        for j = 0, maxEdges do
            edgesOnCell_in_copy[i*maxEdges + j] = cell_region[{i, 0}].edgesOnCell[j]
            verticesOnCell_in_copy[i*maxEdges + j] = cell_region[{i, 0}].verticesOnCell[j]
            cellsOnCell_in_copy[i*maxEdges + j] = cell_region[{i, 0}].cellsOnCell[j]
        end
        --cio.printf("Cell COPY : Cell ID %d, nEdgesOnCell is %d\n", indexToCellID_in_copy[i], nEdgesOnCell_in_copy[i])
    end

    for i = 0, nEdges do
        latEdge_in_copy[i] = edge_region[{i, 0}].lat
        lonEdge_in_copy[i] = edge_region[{i, 0}].lon
        xEdge_in_copy[i] = edge_region[{i, 0}].x
        yEdge_in_copy[i] = edge_region[{i, 0}].y
        zEdge_in_copy[i] = edge_region[{i, 0}].z
        indexToEdgeID_in_copy[i] = edge_region[{i, 0}].edgeID
        nEdgesOnEdge_in_copy[i] = edge_region[{i, 0}].nEdgesOnEdge
        dvEdge_in_copy[i] = edge_region[{i, 0}].dvEdge
        dv1Edge_in_copy[i] = edge_region[{i, 0}].dv1Edge
        dv2Edge_in_copy[i] = edge_region[{i, 0}].dv2Edge
        dcEdge_in_copy[i] = edge_region[{i, 0}].dcEdge
        angleEdge_in_copy[i] = edge_region[{i, 0}].angleEdge

        for j = 0, TWO do
            cellsOnEdge_in_copy[i*TWO + j] = edge_region[{i, 0}].cellsOnEdge[j]
            verticesOnEdge_in_copy[i*TWO + j] = edge_region[{i, 0}].verticesOnEdge[j]
        end

        for j = 0, maxEdges2 do
            edgesOnEdge_in_copy[i*maxEdges2 + j] = edge_region[{i, 0}].edgesOnEdge_ECP[j]
            weightsOnEdge_in_copy[i*maxEdges2 + j] = edge_region[{i, 0}].weightsOnEdge[j]
        end
    end

    for i = 0, nVertices do
        latVertex_in_copy[i] = vertex_region[i].lat
        lonVertex_in_copy[i] = vertex_region[i].lon
        xVertex_in_copy[i] = vertex_region[i].x
        yVertex_in_copy[i] = vertex_region[i].y
        zVertex_in_copy[i] = vertex_region[i].z
        indexToVertexID_in_copy[i] = vertex_region[i].vertexID
        areaTriangle_in_copy[i] = vertex_region[i].areaTriangle

        for j = 0, vertexDegree do
            edgesOnVertex_in_copy[i*vertexDegree + j] = vertex_region[i].edgesOnVertex[j]
            cellsOnVertex_in_copy[i*vertexDegree + j] = vertex_region[i].cellsOnVertex[j]
            kiteAreasOnVertex_in_copy[i*vertexDegree + j] = vertex_region[i].kiteAreasOnVertex[j]
        end
    end


    --Now we put the data into the netcdf file.
    put_var_double(ncid_copy, latCell_varid_copy, latCell_in_copy)
    put_var_double(ncid_copy, lonCell_varid_copy, lonCell_in_copy)
    put_var_double(ncid_copy, meshDensity_varid_copy, meshDensity_in_copy)
    put_var_double(ncid_copy, xCell_varid_copy, xCell_in_copy)
    put_var_double(ncid_copy, yCell_varid_copy, yCell_in_copy)
    put_var_double(ncid_copy, zCell_varid_copy, zCell_in_copy)
    put_var_int(ncid_copy, indexToCellID_varid_copy, indexToCellID_in_copy)
    put_var_int(ncid_copy, nEdgesOnCell_varid_copy, nEdgesOnCell_in_copy)
    put_var_double(ncid_copy, areaCell_varid_copy, areaCell_in_copy)
    put_var_int(ncid_copy, edgesOnCell_varid_copy, edgesOnCell_in_copy)
    put_var_int(ncid_copy, verticesOnCell_varid_copy, verticesOnCell_in_copy)
    put_var_int(ncid_copy, cellsOnCell_varid_copy, cellsOnCell_in_copy)

    put_var_double(ncid_copy, latEdge_varid_copy, latEdge_in_copy)
    put_var_double(ncid_copy, lonEdge_varid_copy, lonEdge_in_copy)
    put_var_double(ncid_copy, xEdge_varid_copy, xEdge_in_copy)
    put_var_double(ncid_copy, yEdge_varid_copy, yEdge_in_copy)
    put_var_double(ncid_copy, zEdge_varid_copy, zEdge_in_copy)
    put_var_int(ncid_copy, indexToEdgeID_varid_copy, indexToEdgeID_in_copy)
    put_var_int(ncid_copy, nEdgesOnEdge_varid_copy, nEdgesOnEdge_in_copy)
    put_var_double(ncid_copy, dvEdge_varid_copy, dvEdge_in_copy)
    put_var_double(ncid_copy, dv1Edge_varid_copy, dv1Edge_in_copy)
    put_var_double(ncid_copy, dv2Edge_varid_copy, dv2Edge_in_copy)
    put_var_double(ncid_copy, dcEdge_varid_copy, dcEdge_in_copy)
    put_var_double(ncid_copy, angleEdge_varid_copy, angleEdge_in_copy)
    put_var_int(ncid_copy, cellsOnEdge_varid_copy, cellsOnEdge_in_copy)
    put_var_int(ncid_copy, verticesOnEdge_varid_copy, verticesOnEdge_in_copy)
    put_var_int(ncid_copy, edgesOnEdge_varid_copy, edgesOnEdge_in_copy)
    put_var_double(ncid_copy, weightsOnEdge_varid_copy, weightsOnEdge_in_copy)

    put_var_double(ncid_copy, latVertex_varid_copy, latVertex_in_copy)
    put_var_double(ncid_copy, lonVertex_varid_copy, lonVertex_in_copy)
    put_var_double(ncid_copy, xVertex_varid_copy, xVertex_in_copy)
    put_var_double(ncid_copy, yVertex_varid_copy, yVertex_in_copy)
    put_var_double(ncid_copy, zVertex_varid_copy, zVertex_in_copy)
    put_var_int(ncid_copy, indexToVertexID_varid_copy, indexToVertexID_in_copy)
    put_var_double(ncid_copy, areaTriangle_varid_copy, areaTriangle_in_copy)
    put_var_int(ncid_copy, edgesOnVertex_varid_copy, edgesOnVertex_in_copy)
    put_var_int(ncid_copy, cellsOnVertex_varid_copy, cellsOnVertex_in_copy)
    put_var_double(ncid_copy, kiteAreasOnVertex_varid_copy, kiteAreasOnVertex_in_copy)

    -- Lastly, we free the allocated memory for the 'copy' arrays
    c.free(latCell_in_copy)
    c.free(lonCell_in_copy)
    c.free(meshDensity_in_copy)
    c.free(xCell_in_copy)
    c.free(yCell_in_copy)
    c.free(zCell_in_copy)
    c.free(indexToCellID_in_copy)
    c.free(latEdge_in_copy)
    c.free(lonEdge_in_copy)
    c.free(xEdge_in_copy)
    c.free(yEdge_in_copy)
    c.free(zEdge_in_copy)
    c.free(indexToEdgeID_in_copy)
    c.free(latVertex_in_copy)
    c.free(lonVertex_in_copy)
    c.free(xVertex_in_copy)
    c.free(yVertex_in_copy)
    c.free(zVertex_in_copy)
    c.free(indexToVertexID_in_copy)
    c.free(cellsOnEdge_in_copy)
    c.free(nEdgesOnCell_in_copy)
    c.free(nEdgesOnEdge_in_copy)
    c.free(edgesOnCell_in_copy)
    c.free(edgesOnEdge_in_copy)
    c.free(weightsOnEdge_in_copy)
    c.free(dvEdge_in_copy)
    c.free(dv1Edge_in_copy)
    c.free(dv2Edge_in_copy)
    c.free(dcEdge_in_copy)
    c.free(angleEdge_in_copy)
    c.free(areaCell_in_copy)
    c.free(areaTriangle_in_copy)
    c.free(cellsOnCell_in_copy)
    c.free(verticesOnCell_in_copy)
    c.free(verticesOnEdge_in_copy)
    c.free(edgesOnVertex_in_copy)
    c.free(cellsOnVertex_in_copy)
    c.free(kiteAreasOnVertex_in_copy)

    -- Close the file
    file_close(ncid_copy)
    cio.printf("Successfully written netcdf file!\n")

end
regentlib.start(main)
