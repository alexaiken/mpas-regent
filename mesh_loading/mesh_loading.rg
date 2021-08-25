import "regent"

require "data_structures"
require "netcdf_tasks"
local constants = require("constants")
local format = require("std/format")

terralib.linklibrary("/home/arjunk1/spack/opt/spack/linux-ubuntu20.04-broadwell/gcc-9.3.0/netcdf-c-4.7.4-zgdvh4hxthdhb3mlsviwhgatvbfnslog/lib/libnetcdf.so")

--Terra function to read the cell partitions from graph.info file. Returns an array where each element is the partition number of that cell index.
terra read_file(file_name: &int8) : int[constants.nCells]
    var file = constants.c.fopen(file_name, "r")
    regentlib.assert(file ~= nil, "failed to open graph.info file")
    var str : int8[constants.MAXCHAR]
    var partition_array : int[constants.nCells]
    var i = 0
    while constants.c.fgets(str, constants.MAXCHAR, file) ~= nil do
        partition_array[i] = constants.c.atoi(str)
        i = i+1
    end
    return partition_array
end


--FILE_NAME, GRAPH_FILE_NAME
--return: cell_region, edge_region, vertex_region
task load_mesh(cell_region : region(ispace(int2d), cell_fs), edge_region : region(ispace(int2d), edge_fs), vertex_region : region(ispace(int2d), vertex_fs), file_name : regentlib.string, graph_file_name : regentlib.string)
where
    reads writes (cell_region, edge_region, vertex_region)
do

    -------------------------------------------
    ----- READ VARIABLES FROM NETCDF FILE -----
    -------------------------------------------
    constants.cio.printf("Starting to read file... (in task load_mesh) \n")
    var ncid : int

    -- Open the file and store the NCID
    open_file(&ncid, [rawstring](file_name))

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
    var latCell_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var lonCell_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var meshDensity_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var xCell_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var yCell_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var zCell_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var indexToCellID_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells))
    var latEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var lonEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var xEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var yEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var zEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var indexToEdgeID_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges))
    var latVertex_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var lonVertex_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var xVertex_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var yVertex_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var zVertex_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var indexToVertexID_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices))
    var cellsOnEdge_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.TWO))
    var nEdgesOnCell_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells))
    var nEdgesOnEdge_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges))
    var edgesOnCell_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var edgesOnEdge_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.maxEdges2))
    var weightsOnEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges*constants.maxEdges2))
    var dvEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dv1Edge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dv2Edge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dcEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var angleEdge_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var areaCell_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var areaTriangle_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var cellsOnCell_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var verticesOnCell_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var verticesOnEdge_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.TWO))
    var edgesOnVertex_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices*constants.vertexDegree))
    var cellsOnVertex_in : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices*constants.vertexDegree))
    var kiteAreasOnVertex_in : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices*constants.vertexDegree))


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



    var partition_array = read_file([rawstring](graph_file_name))

    ----------------------------------
    ----- COPY DATA INTO REGIONS -----
    ----------------------------------

    -- Copy data into cell region
    for i = 0, constants.nCells do
        cell_region[{i, 0}].cellID = indexToCellID_in[i]
        cell_region[{i, 0}].lat = latCell_in[i]
        cell_region[{i, 0}].lon = lonCell_in[i]
        cell_region[{i, 0}].x = xCell_in[i]
        cell_region[{i, 0}].y = yCell_in[i]
        cell_region[{i, 0}].z = zCell_in[i]
        cell_region[{i, 0}].meshDensity = meshDensity_in[i]
        cell_region[{i, 0}].nEdgesOnCell = nEdgesOnCell_in[i]
        cell_region[{i, 0}].areaCell = areaCell_in[i]
        for k = 0, constants.nVertLevels do
            cell_region[{i, k}].partitionNumber = partition_array[i]
        end

        --constants.cio.printf("Cell : Cell ID %d, partitionNumber %d\n", cell_region[{i, 0}].cellID, cell_region[{i, 0}].partitionNumber)

        for j = 0, constants.maxEdges do
            cell_region[{i, 0}].edgesOnCell[j] = edgesOnCell_in[i*constants.maxEdges + j] --cell_region[{i, 0}].edgesOnCell is a int[constants.maxEdges]
            cell_region[{i, 0}].verticesOnCell[j] = verticesOnCell_in[i*constants.maxEdges + j] --cell_region[{i, 0}].verticesOnCell is a int[constants.maxEdges]
            cell_region[{i, 0}].cellsOnCell[j] = cellsOnCell_in[i*constants.maxEdges + j] --cell_region[{i, 0}].cellsOnCell is a int[constants.maxEdges]
            --constants.cio.printf("edgesOnCell : Cell %d, Edge %d: edge index is %d\n", i, j, cell_region[{i, 0}].edgesOnCell[j])
            --constants.cio.printf("verticesOnCell : Cell %d, Vertex %d: Vertex index is %d\n", i, j, cell_region[{i, 0}].verticesOnCell[j])
            --constants.cio.printf("cellsOnCell : InnerCell %d, OuterCell %d: Cell index is %d\n", i, j, cell_region[{i, 0}].cellsOnCell[j])
        end

        --cell_region[{i, 0}].edgesOnCell0 = ptr(cell_region[{i, 0}].edgesOnCell[0])
        cell_region[{i, 0}].edgesOnCell0 = { int2d {cell_region[{i, 0}].edgesOnCell[0], 0}, int2d {cell_region[{i, 0}].edgesOnCell[0], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell1 = { int2d {cell_region[{i, 0}].edgesOnCell[1], 0}, int2d {cell_region[{i, 0}].edgesOnCell[1], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell2 = { int2d {cell_region[{i, 0}].edgesOnCell[2], 0}, int2d {cell_region[{i, 0}].edgesOnCell[2], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell3 = { int2d {cell_region[{i, 0}].edgesOnCell[3], 0}, int2d {cell_region[{i, 0}].edgesOnCell[3], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell4 = { int2d {cell_region[{i, 0}].edgesOnCell[4], 0}, int2d {cell_region[{i, 0}].edgesOnCell[4], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell5 = { int2d {cell_region[{i, 0}].edgesOnCell[5], 0}, int2d {cell_region[{i, 0}].edgesOnCell[5], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell6 = { int2d {cell_region[{i, 0}].edgesOnCell[6], 0}, int2d {cell_region[{i, 0}].edgesOnCell[6], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell7 = { int2d {cell_region[{i, 0}].edgesOnCell[7], 0}, int2d {cell_region[{i, 0}].edgesOnCell[7], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell8 = { int2d {cell_region[{i, 0}].edgesOnCell[8], 0}, int2d {cell_region[{i, 0}].edgesOnCell[8], constants.nVertLevels - 1} }
        cell_region[{i, 0}].edgesOnCell9 = { int2d {cell_region[{i, 0}].edgesOnCell[9], 0}, int2d {cell_region[{i, 0}].edgesOnCell[9], constants.nVertLevels - 1} }

        --constants.cio.printf("Cell : Cell ID %d, nEdgesOnCell is %d\n", cell_region[{i, 0}].cellID, cell_region[{i, 0}].nEdgesOnCell)
    end

    -- Copy data into edge region
    for i = 0, constants.nEdges do
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


        for j = 0, constants.TWO do
            edge_region[{i, 0}].cellsOnEdge[j] = cellsOnEdge_in[i*constants.TWO + j]
            edge_region[{i, 0}].verticesOnEdge[j] = verticesOnEdge_in[i*constants.TWO + j]
            --constants.cio.printf("cellsOnEdge : Edge %d, Cell %d is %d\n", i, j, edge_region[{i, 0}].cellsOnEdge[j])
            --constants.cio.printf("VerticesOnEdge : Edge %d: Vertex %d is $d\n", i, j, edge_region[{i, 0}].verticesOnEdge[j])
        end

        for j = 0, constants.maxEdges2 do
            edge_region[{i, 0}].edgesOnEdge_ECP[j] = edgesOnEdge_in[i*constants.maxEdges2 + j]
            edge_region[{i, 0}].weightsOnEdge[j] = weightsOnEdge_in[i*constants.maxEdges2 + j]
            --constants.cio.printf("edgesOnEdge_ECP : InnerEdge %d, OuterEdge %d is %d\n", i, j, edge_region[{i, 0}].edgesOnEdge_ECP[j])
            --constants.cio.printf("weightsOnEdge : Edge %d: Weight %d is $f\n", i, j, edge_region[{i, 0}].weightsOnEdge[j])
        end
        --constants.cio.printf("Edge: ID is %d, xEdge is %f, yEdge is %f, zEdge is %f \n", i, edge_region[{i, 0}].x, edge_region[{i, 0}].y, edge_region[{i, 0}].z)
    end

    -- Copy data into vertex region
    for i = 0, constants.nVertices do
        vertex_region[{i, 0}].vertexID = indexToVertexID_in[i]
        vertex_region[{i, 0}].lat = latVertex_in[i]
        vertex_region[{i, 0}].lon = lonVertex_in[i]
        vertex_region[{i, 0}].x = xVertex_in[i]
        vertex_region[{i, 0}].y = yVertex_in[i]
        vertex_region[{i, 0}].z = zVertex_in[i]
        vertex_region[{i, 0}].areaTriangle = areaTriangle_in[i]

        for j = 0, constants.vertexDegree do
            vertex_region[{i, 0}].edgesOnVertex[j] = edgesOnVertex_in[i*constants.vertexDegree + j]
            vertex_region[{i, 0}].cellsOnVertex[j] = cellsOnVertex_in[i*constants.vertexDegree + j]
            vertex_region[{i, 0}].kiteAreasOnVertex[j] = kiteAreasOnVertex_in[i*constants.vertexDegree + j]

            --constants.cio.printf("edgesOnVertex : Vertex %d, Edge %d: Edge index is %d\n", i, j, vertex_region[{i, 0}].edgesOnVertex[j])
            --constants.cio.printf("cellsOnVertex : Vertex %d, Cell %d: Cell index is %d\n", i, j, vertex_region[{i, 0}].cellsOnVertex[j])
            --constants.cio.printf("kiteAreasOnVertex : Vertex %d, Kite %d: Kite Area is %f\n", i, j, vertex_region[{i, 0}].kiteAreasOnVertex[j])
        end
        --constants.cio.printf("Vertex ID is %d, xVertex is %f, yVertex is %f, zVertex is %f \n", i, vertex_region[{i, 0}].x, vertex_region[{i, 0}].y, vertex_region[{i, 0}].z)
    end

    -------------------------
    ----- CALCULATE EVC -----
    -------------------------
    --I know I should do something more intelligent to get the common elements: but for now we do a brute force search to get EVC

    --First, we iterate through the cells and get the edgesOnCell array for each cell
    for i = 0, constants.nCells do
        var curr_edgesOnCell = cell_region[{i, 0}].edgesOnCell

    --Then we iterate through the vertices of that cell
        for j = 0, constants.maxEdges do
            var currVertexID = cell_region[{i, 0}].verticesOnCell[j]
            cell_region[{i,0}].evc[j*3] = currVertexID
            --constants.cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3, cell_region[{i, 0}].evc[j*3])

            if currVertexID == 0 then
                cell_region[{i,0}].evc[j*3 + 1] = 0
                cell_region[{i,0}].evc[j*3 + 2] = 0
                --constants.cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3 + 1, cell_region[{i, 0}].evc[j*3 + 1])
                --constants.cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3 + 2, cell_region[{i, 0}].evc[j*3 + 2])

            --If there is a vertex, we get the edges on that vertex
            elseif currVertexID ~= 0 then
                var curr_edgesOnVertex = vertex_region[{currVertexID-1, 0}].edgesOnVertex
                var count = 1

                --Then, we get overlapping edges between curr_edgesOnVertex and curr_edgesOnCell to get EVC
                for k = 0, constants.vertexDegree do
                    var currEdgeID = curr_edgesOnVertex[k]
                    for l = 0, constants.maxEdges do
                        if currEdgeID == curr_edgesOnCell[l] and count < 3 then
                            cell_region[{i,0}].evc[j*3 + count] = currEdgeID
                            --constants.cio.printf("cell_region[%d].evc[%d] = %d\n", i, j*3 + count, cell_region[{i, 0}].evc[j*3 + count])
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
    constants.c.free(latCell_in)
    constants.c.free(lonCell_in)
    constants.c.free(meshDensity_in)
    constants.c.free(xCell_in)
    constants.c.free(yCell_in)
    constants.c.free(zCell_in)
    constants.c.free(indexToCellID_in)
    constants.c.free(latEdge_in)
    constants.c.free(lonEdge_in)
    constants.c.free(xEdge_in)
    constants.c.free(yEdge_in)
    constants.c.free(zEdge_in)
    constants.c.free(indexToEdgeID_in)
    constants.c.free(latVertex_in)
    constants.c.free(lonVertex_in)
    constants.c.free(xVertex_in)
    constants.c.free(yVertex_in)
    constants.c.free(zVertex_in)
    constants.c.free(indexToVertexID_in)
    constants.c.free(cellsOnEdge_in)
    constants.c.free(nEdgesOnCell_in)
    constants.c.free(nEdgesOnEdge_in)
    constants.c.free(edgesOnCell_in)
    constants.c.free(edgesOnEdge_in)
    constants.c.free(weightsOnEdge_in)
    constants.c.free(dvEdge_in)
    constants.c.free(dv1Edge_in)
    constants.c.free(dv2Edge_in)
    constants.c.free(dcEdge_in)
    constants.c.free(angleEdge_in)
    constants.c.free(areaCell_in)
    constants.c.free(areaTriangle_in)
    constants.c.free(cellsOnCell_in)
    constants.c.free(verticesOnCell_in)
    constants.c.free(verticesOnEdge_in)
    constants.c.free(edgesOnVertex_in)
    constants.c.free(cellsOnVertex_in)
    constants.c.free(kiteAreasOnVertex_in)

    constants.cio.printf("Successfully read file! (in task load_mesh) \n")
end


-----------------------------------------------
------- TASK: PARTITION REGIONS  --------
-----------------------------------------------

--input: cell_region, edge_region, vertex_region
--return: cell_partition_fs containing private, shared, and ghost partitions for one and two layers.
task partition_regions(num_partitions : int, cell_region : region(ispace(int2d), cell_fs), edge_region : region(ispace(int2d), edge_fs), vertex_region : region(ispace(int2d), vertex_fs))
where
    reads writes (cell_region, edge_region, vertex_region)
do

    var color_space = ispace(int1d, num_partitions)
    var p = partition(cell_region.partitionNumber, color_space) -- Original partition based on Metis

    format.println("p[1] volume: {}", p[1].volume)

    var e0 = image(edge_region, p, cell_region.edgesOnCell0)
    var e1 = image(edge_region, p, cell_region.edgesOnCell1)
    var e2 = image(edge_region, p, cell_region.edgesOnCell2)
    var e3 = image(edge_region, p, cell_region.edgesOnCell3)
    var e4 = image(edge_region, p, cell_region.edgesOnCell4)
    var e5 = image(edge_region, p, cell_region.edgesOnCell5)
    var e6 = image(edge_region, p, cell_region.edgesOnCell6)
    var e7 = image(edge_region, p, cell_region.edgesOnCell7)
    var e8 = image(edge_region, p, cell_region.edgesOnCell8)
    var e9 = image(edge_region, p, cell_region.edgesOnCell9)
    var e = e0 | e1 | e2 | e3 | e4 | e5 | e6 | e7 | e8 | e9

    format.println("e[1].volume={}", e[1].volume)

    --for k = 0, constants.NUM_PARTITIONS do
    --    format.println("Partition {}", k)
    --    for i in e0[k] do
    --        if i.y == 0 then
    --            format.println("{}, {}", k, i)
    --        end
    --    end
    --end

    -- TODO: Not sure if we need multiple levels at once. If not, we can use the int2d type and declare with only the first half `int2d { edge_region[{i, 0}].cellsOnEdge[1], 0 }`
    for i = 0, constants.nEdges do
        edge_region[{i, 0}].cellOne = rect2d { int2d { edge_region[{i, 0}].cellsOnEdge[0], 0 }, int2d { edge_region[{i, 0}].cellsOnEdge[0], constants.nVertLevels - 1 } }
        edge_region[{i, 0}].cellTwo = rect2d { int2d { edge_region[{i, 0}].cellsOnEdge[1], 0 }, int2d { edge_region[{i, 0}].cellsOnEdge[1], constants.nVertLevels - 1 } }
    end

    --Note: the following code should theoretically be able to replace the 10 images, but seems to give the wrong volume
    --var ep_one = preimage(edge_region, p, edge_region.cellOne) -- This should get all edges with cellOne in p, partitioned by cellOne (which cell "originated" it)
    --var ep_two = preimage(edge_region, p, edge_region.cellTwo) -- This should get all edges with cellTwo in p, partitioned by cellTwo (which cell it "points to")

    var cp_one = image(cell_region, e, edge_region.cellOne) -- This gets the cellOne on all edges in e
    var cp_two = image(cell_region, e, edge_region.cellTwo) -- This gets the cellTwo on all edges in e
    -- Combined, they should contain all cells who have an edge in e, and therefore the halo.

    var ghost_1_and_p = cp_one | cp_two
    format.println("(cp_one | cp_two)[1] volume: {}", ghost_1_and_p[1].volume)
    var ghost_1 = (cp_one | cp_two) - p

    -- Calculate second halo
    var gep_out = preimage(edge_region, ghost_1_and_p, edge_region.cellOne)
    var gep_in = preimage(edge_region, ghost_1_and_p, edge_region.cellTwo)
    var gcp_out = image(cell_region, gep_out, edge_region.cellTwo)
    var gcp_in = image(cell_region, gep_in, edge_region.cellOne)
    var ghost_2 = (gcp_in | gcp_out) - p -- First and second halo layers

    -- Compute all cells reachable from ghost_1. shared_1 is intersection of that set with p
    var s1ep_out = preimage(edge_region, ghost_1, edge_region.cellOne)
    var s1ep_in = preimage(edge_region, ghost_1, edge_region.cellTwo)
    var s1cp_out = image(cell_region, s1ep_out, edge_region.cellTwo)
    var s1cp_in = image(cell_region, s1ep_in, edge_region.cellOne)
    var shared_1 = p & (s1cp_out | s1cp_in) -- Cells in p bordering ghost_1
    var private_1 = p - shared_1 -- all cells in p that are not in shared_1

    -- shared_2 contains shared_1 and all cells in p bordering shared_1
    var s2ep_out = preimage(edge_region, shared_1, edge_region.cellOne)
    var s2ep_in = preimage(edge_region, shared_1, edge_region.cellTwo)
    var s2cp_out = image(cell_region, s2ep_out, edge_region.cellTwo)
    var s2cp_in = image(cell_region, s2ep_in, edge_region.cellOne)
    var shared_2 = dynamic_cast(partition(disjoint, cell_region, color_space), (shared_1 | (private_1 & (s2cp_out | s2cp_in)))) -- Cells in p bordering ghost_1
    var private_2 = private_1 - shared_2 -- all cells in private_1 that are not in shared_2

    format.println("Private1[1] volume: {}", private_1[1].volume)
    format.println("Private2[1] volume: {}", private_2[1].volume)
    format.println("Shared1[1] volume: {}", shared_1[1].volume)
    format.println("Shared2[1] volume: {}", shared_2[1].volume)
    format.println("Ghost1[1] volume: {}", ghost_1[1].volume)
    format.println("Ghost2[1] volume: {}", ghost_2[1].volume)

    return [cell_partition_fs(cell_region)] {
        private_1, shared_1, ghost_1, private_2, shared_2, ghost_2
    }
end

----------------------------------------------------
------- TASK: TESTING CODE: WRITING NETCDF OUTPUT --------
----------------------------------------------------

--input: cell_region, edge_region, vertex_region
task write_output(cell_region : region(ispace(int2d), cell_fs), edge_region : region(ispace(int2d), edge_fs), vertex_region : region(ispace(int2d), vertex_fs))
where
    reads writes (cell_region, edge_region, vertex_region)
do


    ----------------------------------------------------
    ------- TESTING CODE: WRITING NETCDF OUTPUT --------
    ----------------------------------------------------

    -- We create a netcdf file using the data in the regions, to test whether the data was written correctly.
    constants.cio.printf("Starting to write netcdf file..\n")
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
    var nVertLevelsplus_dimid_copy : int
    var time_dimid_copy : int

    --Define the dimension variables
    --define_dim(ncid: int, dim_name: &int, dim_size: int, dim_id_ptr: &int)
    define_dim(ncid_copy, "nCells", constants.nCells, &nCells_dimid_copy)
    define_dim(ncid_copy, "nEdges", constants.nEdges, &nEdges_dimid_copy)
    define_dim(ncid_copy, "nVertices", constants.nVertices, &nVertices_dimid_copy)
    define_dim(ncid_copy, "maxEdges", constants.maxEdges, &maxEdges_dimid_copy)
    define_dim(ncid_copy, "maxEdges2", constants.maxEdges2, &maxEdges2_dimid_copy)
    define_dim(ncid_copy, "TWO", constants.TWO, &TWO_dimid_copy)
    define_dim(ncid_copy, "vertexDegree", constants.vertexDegree, &vertexDegree_dimid_copy)
    define_dim(ncid_copy, "nVertLevels", constants.nVertLevels, &nVertLevels_dimid_copy)
    define_dim(ncid_copy, "nVertLevelsplus", constants.nVertLevels + 1, &nVertLevelsplus_dimid_copy)
    define_dim(ncid_copy, "Time", constants.netcdf.NC_UNLIMITED, &time_dimid_copy)

    --For the 2D variables, the dimIDs need to be put in arrays
    var nEdges_TWO_dimids = array(nEdges_dimid_copy, TWO_dimid_copy)
    var nCells_maxEdges_dimids = array(nCells_dimid_copy, maxEdges_dimid_copy)
    var nEdges_maxEdges2_dimids = array(nEdges_dimid_copy, maxEdges2_dimid_copy)
    var nVertices_vertexDegree_dimids = array(nVertices_dimid_copy, vertexDegree_dimid_copy)
    var nEdges_nVertLevelsplus_dimids = array(nEdges_dimid_copy, nVertLevelsplus_dimid_copy)
    var nCells_nVertLevels_dimids = array(nCells_dimid_copy, nVertLevels_dimid_copy)

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
    define_var(ncid_copy, "latCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &latCell_varid_copy)
    define_var(ncid_copy, "lonCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &lonCell_varid_copy)
    define_var(ncid_copy, "meshDensity", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &meshDensity_varid_copy)
    define_var(ncid_copy, "xCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &xCell_varid_copy)
    define_var(ncid_copy, "yCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &yCell_varid_copy)
    define_var(ncid_copy, "zCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &zCell_varid_copy)
    define_var(ncid_copy, "indexToCellID", constants.netcdf.NC_INT, 1, &nCells_dimid_copy, &indexToCellID_varid_copy)
    define_var(ncid_copy, "latEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &latEdge_varid_copy)
    define_var(ncid_copy, "lonEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &lonEdge_varid_copy)
    define_var(ncid_copy, "xEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &xEdge_varid_copy)
    define_var(ncid_copy, "yEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &yEdge_varid_copy)
    define_var(ncid_copy, "zEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &zEdge_varid_copy)
    define_var(ncid_copy, "indexToEdgeID", constants.netcdf.NC_INT, 1, &nEdges_dimid_copy, &indexToEdgeID_varid_copy)
    define_var(ncid_copy, "latVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &latVertex_varid_copy)
    define_var(ncid_copy, "lonVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &lonVertex_varid_copy)
    define_var(ncid_copy, "xVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &xVertex_varid_copy)
    define_var(ncid_copy, "yVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &yVertex_varid_copy)
    define_var(ncid_copy, "zVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &zVertex_varid_copy)
    define_var(ncid_copy, "indexToVertexID", constants.netcdf.NC_INT, 1, &nVertices_dimid_copy, &indexToVertexID_varid_copy)
    define_var(ncid_copy, "cellsOnEdge", constants.netcdf.NC_INT, 2, nEdges_TWO_dimids, &cellsOnEdge_varid_copy)
    define_var(ncid_copy, "nEdgesOnCell", constants.netcdf.NC_INT, 1, &nCells_dimid_copy, &nEdgesOnCell_varid_copy)
    define_var(ncid_copy, "nEdgesOnEdge", constants.netcdf.NC_INT, 1, &nEdges_dimid_copy, &nEdgesOnEdge_varid_copy)
    define_var(ncid_copy, "edgesOnCell", constants.netcdf.NC_INT, 2, nCells_maxEdges_dimids, &edgesOnCell_varid_copy)
    define_var(ncid_copy, "edgesOnEdge", constants.netcdf.NC_INT, 2, nEdges_maxEdges2_dimids, &edgesOnEdge_varid_copy)
    define_var(ncid_copy, "weightsOnEdge", constants.netcdf.NC_DOUBLE, 2, nEdges_maxEdges2_dimids, &weightsOnEdge_varid_copy)
    define_var(ncid_copy, "dvEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dvEdge_varid_copy)
    define_var(ncid_copy, "dv1Edge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dv1Edge_varid_copy)
    define_var(ncid_copy, "dv2Edge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dv2Edge_varid_copy)
    define_var(ncid_copy, "dcEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dcEdge_varid_copy)
    define_var(ncid_copy, "angleEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &angleEdge_varid_copy)
    define_var(ncid_copy, "areaCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &areaCell_varid_copy)
    define_var(ncid_copy, "areaTriangle", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &areaTriangle_varid_copy)
    define_var(ncid_copy, "cellsOnCell", constants.netcdf.NC_INT, 2, nCells_maxEdges_dimids, &cellsOnCell_varid_copy)
    define_var(ncid_copy, "verticesOnCell", constants.netcdf.NC_INT, 2, nCells_maxEdges_dimids, &verticesOnCell_varid_copy)
    define_var(ncid_copy, "verticesOnEdge", constants.netcdf.NC_INT, 2, nEdges_TWO_dimids, &verticesOnEdge_varid_copy)
    define_var(ncid_copy, "edgesOnVertex", constants.netcdf.NC_INT, 2, nVertices_vertexDegree_dimids, &edgesOnVertex_varid_copy)
    define_var(ncid_copy, "cellsOnVertex", constants.netcdf.NC_INT, 2, nVertices_vertexDegree_dimids, &cellsOnVertex_varid_copy)
    define_var(ncid_copy, "kiteAreasOnVertex", constants.netcdf.NC_DOUBLE, 2, nVertices_vertexDegree_dimids, &kiteAreasOnVertex_varid_copy)

    --This function signals that we're done writing the metadata.
    end_def(ncid_copy)

    --Now define the new arrays to hold the data that will be put in the netcdf files
    var latCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var lonCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var meshDensity_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var xCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var yCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var zCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var indexToCellID_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells))
    var latEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var lonEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var xEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var yEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var zEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var indexToEdgeID_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges))
    var latVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var lonVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var xVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var yVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var zVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var indexToVertexID_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices))
    var cellsOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.TWO))
    var nEdgesOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells))
    var nEdgesOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges))
    var edgesOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var edgesOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.maxEdges2))
    var weightsOnEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges*constants.maxEdges2))
    var dvEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dv1Edge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dv2Edge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dcEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var angleEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var areaCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var areaTriangle_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var cellsOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var verticesOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var verticesOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.TWO))
    var edgesOnVertex_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices*constants.vertexDegree))
    var cellsOnVertex_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices*constants.vertexDegree))
    var kiteAreasOnVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices*constants.vertexDegree))

    --Now we copy the data into the arrays so they can be read into the netcdf files
    for i = 0, constants.nCells do
        latCell_in_copy[i] = cell_region[{i, 0}].lat
        lonCell_in_copy[i] = cell_region[{i, 0}].lon
        xCell_in_copy[i] = cell_region[{i, 0}].x
        yCell_in_copy[i] = cell_region[{i, 0}].y
        zCell_in_copy[i] = cell_region[{i, 0}].z
        indexToCellID_in_copy[i] = cell_region[{i, 0}].cellID
        meshDensity_in_copy[i] = cell_region[{i, 0}].meshDensity
        nEdgesOnCell_in_copy[i] = cell_region[{i, 0}].nEdgesOnCell
        areaCell_in_copy[i] = cell_region[{i, 0}].areaCell

        for j = 0, constants.maxEdges do
            edgesOnCell_in_copy[i*constants.maxEdges + j] = cell_region[{i, 0}].edgesOnCell[j]
            verticesOnCell_in_copy[i*constants.maxEdges + j] = cell_region[{i, 0}].verticesOnCell[j]
            cellsOnCell_in_copy[i*constants.maxEdges + j] = cell_region[{i, 0}].cellsOnCell[j]
        end
        --constants.cio.printf("Cell COPY : Cell ID %d, nEdgesOnCell is %d\n", indexToCellID_in_copy[i], nEdgesOnCell_in_copy[i])
    end

    for i = 0, constants.nEdges do
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

        for j = 0, constants.TWO do
            cellsOnEdge_in_copy[i*constants.TWO + j] = edge_region[{i, 0}].cellsOnEdge[j]
            verticesOnEdge_in_copy[i*constants.TWO + j] = edge_region[{i, 0}].verticesOnEdge[j]
        end

        for j = 0, constants.maxEdges2 do
            edgesOnEdge_in_copy[i*constants.maxEdges2 + j] = edge_region[{i, 0}].edgesOnEdge_ECP[j]
            weightsOnEdge_in_copy[i*constants.maxEdges2 + j] = edge_region[{i, 0}].weightsOnEdge[j]
        end
    end

    for i = 0, constants.nVertices do
        latVertex_in_copy[i] = vertex_region[{i, 0}].lat
        lonVertex_in_copy[i] = vertex_region[{i, 0}].lon
        xVertex_in_copy[i] = vertex_region[{i, 0}].x
        yVertex_in_copy[i] = vertex_region[{i, 0}].y
        zVertex_in_copy[i] = vertex_region[{i, 0}].z
        indexToVertexID_in_copy[i] = vertex_region[{i, 0}].vertexID
        areaTriangle_in_copy[i] = vertex_region[{i, 0}].areaTriangle

        for j = 0, constants.vertexDegree do
            edgesOnVertex_in_copy[i*constants.vertexDegree + j] = vertex_region[{i, 0}].edgesOnVertex[j]
            cellsOnVertex_in_copy[i*constants.vertexDegree + j] = vertex_region[{i, 0}].cellsOnVertex[j]
            kiteAreasOnVertex_in_copy[i*constants.vertexDegree + j] = vertex_region[{i, 0}].kiteAreasOnVertex[j]
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
    constants.c.free(latCell_in_copy)
    constants.c.free(lonCell_in_copy)
    constants.c.free(meshDensity_in_copy)
    constants.c.free(xCell_in_copy)
    constants.c.free(yCell_in_copy)
    constants.c.free(zCell_in_copy)
    constants.c.free(indexToCellID_in_copy)
    constants.c.free(latEdge_in_copy)
    constants.c.free(lonEdge_in_copy)
    constants.c.free(xEdge_in_copy)
    constants.c.free(yEdge_in_copy)
    constants.c.free(zEdge_in_copy)
    constants.c.free(indexToEdgeID_in_copy)
    constants.c.free(latVertex_in_copy)
    constants.c.free(lonVertex_in_copy)
    constants.c.free(xVertex_in_copy)
    constants.c.free(yVertex_in_copy)
    constants.c.free(zVertex_in_copy)
    constants.c.free(indexToVertexID_in_copy)
    constants.c.free(cellsOnEdge_in_copy)
    constants.c.free(nEdgesOnCell_in_copy)
    constants.c.free(nEdgesOnEdge_in_copy)
    constants.c.free(edgesOnCell_in_copy)
    constants.c.free(edgesOnEdge_in_copy)
    constants.c.free(weightsOnEdge_in_copy)
    constants.c.free(dvEdge_in_copy)
    constants.c.free(dv1Edge_in_copy)
    constants.c.free(dv2Edge_in_copy)
    constants.c.free(dcEdge_in_copy)
    constants.c.free(angleEdge_in_copy)
    constants.c.free(areaCell_in_copy)
    constants.c.free(areaTriangle_in_copy)
    constants.c.free(cellsOnCell_in_copy)
    constants.c.free(verticesOnCell_in_copy)
    constants.c.free(verticesOnEdge_in_copy)
    constants.c.free(edgesOnVertex_in_copy)
    constants.c.free(cellsOnVertex_in_copy)
    constants.c.free(kiteAreasOnVertex_in_copy)

    -- Close the file
    file_close(ncid_copy)
    constants.cio.printf("Successfully written netcdf file!\n")

end

--input: cell_region, edge_region, vertex_region
task write_output_plotting(cell_region : region(ispace(int2d), cell_fs), edge_region : region(ispace(int2d), edge_fs), vertex_region : region(ispace(int2d), vertex_fs))
where
    reads writes (cell_region, edge_region, vertex_region)
do


    ----------------------------------------------------
    ------- TESTING CODE: WRITING NETCDF OUTPUT --------
    ----------------------------------------------------

    -- We create a netcdf file using the data in the regions, to test whether the data was written correctly.
    constants.cio.printf("Starting to write netcdf file..\n")
    var ncid_copy = 65538

    --Create a netcdf file
    file_create("timestep_output.nc", &ncid_copy)

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
    var nVertLevelsplus_dimid_copy : int

    --Define the dimension variables
    --define_dim(ncid: int, dim_name: &int, dim_size: int, dim_id_ptr: &int)
    define_dim(ncid_copy, "nCells", constants.nCells, &nCells_dimid_copy)
    define_dim(ncid_copy, "nEdges", constants.nEdges, &nEdges_dimid_copy)
    define_dim(ncid_copy, "nVertices", constants.nVertices, &nVertices_dimid_copy)
    define_dim(ncid_copy, "maxEdges", constants.maxEdges, &maxEdges_dimid_copy)
    define_dim(ncid_copy, "maxEdges2", constants.maxEdges2, &maxEdges2_dimid_copy)
    define_dim(ncid_copy, "TWO", constants.TWO, &TWO_dimid_copy)
    define_dim(ncid_copy, "vertexDegree", constants.vertexDegree, &vertexDegree_dimid_copy)
    define_dim(ncid_copy, "nVertLevels", constants.nVertLevels, &nVertLevels_dimid_copy)
    define_dim(ncid_copy, "nVertLevelsplus", constants.nVertLevels + 1, &nVertLevelsplus_dimid_copy)
    define_dim(ncid_copy, "Time", constants.netcdf.NC_UNLIMITED, &time_dimid_copy)

    --For the 2D variables, the dimIDs need to be put in arrays
    var nEdges_TWO_dimids = array(nEdges_dimid_copy, TWO_dimid_copy)
    var nCells_maxEdges_dimids = array(nCells_dimid_copy, maxEdges_dimid_copy)
    var nEdges_maxEdges2_dimids = array(nEdges_dimid_copy, maxEdges2_dimid_copy)
    var nVertices_vertexDegree_dimids = array(nVertices_dimid_copy, vertexDegree_dimid_copy)
    var nCells_nVertLevels_dimids = array(nCells_dimid_copy, nVertLevels_dimid_copy)
    var nEdges_nVertLevels_dimids = array(nEdges_dimid_copy, nVertLevels_dimid_copy)
    var nCells_nVertLevelsplus_dimids = array(nCells_dimid_copy, nVertLevelsplus_dimid_copy)
    

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

    var u_varid_copy : int
    var v_varid_copy : int

    var w_varid_copy : int
    var pressure_varid_copy : int
    var pressure_p_varid_copy : int
    var rho_varid_copy : int
    var theta_varid_copy : int
    var surface_pressure_varid_copy : int

    --Define the variable IDs
    define_var(ncid_copy, "latCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &latCell_varid_copy)
    define_var(ncid_copy, "lonCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &lonCell_varid_copy)
    define_var(ncid_copy, "meshDensity", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &meshDensity_varid_copy)
    define_var(ncid_copy, "xCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &xCell_varid_copy)
    define_var(ncid_copy, "yCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &yCell_varid_copy)
    define_var(ncid_copy, "zCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &zCell_varid_copy)
    define_var(ncid_copy, "indexToCellID", constants.netcdf.NC_INT, 1, &nCells_dimid_copy, &indexToCellID_varid_copy)
    define_var(ncid_copy, "latEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &latEdge_varid_copy)
    define_var(ncid_copy, "lonEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &lonEdge_varid_copy)
    define_var(ncid_copy, "xEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &xEdge_varid_copy)
    define_var(ncid_copy, "yEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &yEdge_varid_copy)
    define_var(ncid_copy, "zEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &zEdge_varid_copy)
    define_var(ncid_copy, "indexToEdgeID", constants.netcdf.NC_INT, 1, &nEdges_dimid_copy, &indexToEdgeID_varid_copy)
    define_var(ncid_copy, "latVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &latVertex_varid_copy)
    define_var(ncid_copy, "lonVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &lonVertex_varid_copy)
    define_var(ncid_copy, "xVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &xVertex_varid_copy)
    define_var(ncid_copy, "yVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &yVertex_varid_copy)
    define_var(ncid_copy, "zVertex", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &zVertex_varid_copy)
    define_var(ncid_copy, "indexToVertexID", constants.netcdf.NC_INT, 1, &nVertices_dimid_copy, &indexToVertexID_varid_copy)
    define_var(ncid_copy, "cellsOnEdge", constants.netcdf.NC_INT, 2, nEdges_TWO_dimids, &cellsOnEdge_varid_copy)
    define_var(ncid_copy, "nEdgesOnCell", constants.netcdf.NC_INT, 1, &nCells_dimid_copy, &nEdgesOnCell_varid_copy)
    define_var(ncid_copy, "nEdgesOnEdge", constants.netcdf.NC_INT, 1, &nEdges_dimid_copy, &nEdgesOnEdge_varid_copy)
    define_var(ncid_copy, "edgesOnCell", constants.netcdf.NC_INT, 2, nCells_maxEdges_dimids, &edgesOnCell_varid_copy)
    define_var(ncid_copy, "edgesOnEdge", constants.netcdf.NC_INT, 2, nEdges_maxEdges2_dimids, &edgesOnEdge_varid_copy)
    define_var(ncid_copy, "weightsOnEdge", constants.netcdf.NC_DOUBLE, 2, nEdges_maxEdges2_dimids, &weightsOnEdge_varid_copy)
    define_var(ncid_copy, "dvEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dvEdge_varid_copy)
    define_var(ncid_copy, "dv1Edge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dv1Edge_varid_copy)
    define_var(ncid_copy, "dv2Edge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dv2Edge_varid_copy)
    define_var(ncid_copy, "dcEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &dcEdge_varid_copy)
    define_var(ncid_copy, "angleEdge", constants.netcdf.NC_DOUBLE, 1, &nEdges_dimid_copy, &angleEdge_varid_copy)
    define_var(ncid_copy, "areaCell", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &areaCell_varid_copy)
    define_var(ncid_copy, "areaTriangle", constants.netcdf.NC_DOUBLE, 1, &nVertices_dimid_copy, &areaTriangle_varid_copy)
    define_var(ncid_copy, "cellsOnCell", constants.netcdf.NC_INT, 2, nCells_maxEdges_dimids, &cellsOnCell_varid_copy)
    define_var(ncid_copy, "verticesOnCell", constants.netcdf.NC_INT, 2, nCells_maxEdges_dimids, &verticesOnCell_varid_copy)
    define_var(ncid_copy, "verticesOnEdge", constants.netcdf.NC_INT, 2, nEdges_TWO_dimids, &verticesOnEdge_varid_copy)
    define_var(ncid_copy, "edgesOnVertex", constants.netcdf.NC_INT, 2, nVertices_vertexDegree_dimids, &edgesOnVertex_varid_copy)
    define_var(ncid_copy, "cellsOnVertex", constants.netcdf.NC_INT, 2, nVertices_vertexDegree_dimids, &cellsOnVertex_varid_copy)
    define_var(ncid_copy, "kiteAreasOnVertex", constants.netcdf.NC_DOUBLE, 2, nVertices_vertexDegree_dimids, &kiteAreasOnVertex_varid_copy)

    define_var(ncid_copy, "u", constants.netcdf.NC_DOUBLE, 2, nEdges_nVertLevels_dimids, &u_varid_copy)
    define_var(ncid_copy, "v", constants.netcdf.NC_DOUBLE, 2, nEdges_nVertLevels_dimids, &v_varid_copy)

    define_var(ncid_copy, "w", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevelsplus_dimids, &w_varid_copy)
    define_var(ncid_copy, "pressure", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &pressure_varid_copy)

    define_var(ncid_copy, "pressure_p", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &pressure_p_varid_copy)
    define_var(ncid_copy, "rho", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &rho_varid_copy)
    define_var(ncid_copy, "theta", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &theta_varid_copy)
    define_var(ncid_copy, "surface_pressure", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &surface_pressure_varid_copy)


    --This function signals that we're done writing the metadata.
    end_def(ncid_copy)

    --Now define the new arrays to hold the data that will be put in the netcdf files
    var latCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var lonCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var meshDensity_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var xCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var yCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var zCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var indexToCellID_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells))
    var latEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var lonEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var xEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var yEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var zEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var indexToEdgeID_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges))
    var latVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var lonVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var xVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var yVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var zVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var indexToVertexID_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices))
    var cellsOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.TWO))
    var nEdgesOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells))
    var nEdgesOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges))
    var edgesOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var edgesOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.maxEdges2))
    var weightsOnEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges*constants.maxEdges2))
    var dvEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dv1Edge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dv2Edge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var dcEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var angleEdge_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges))
    var areaCell_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var areaTriangle_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices))
    var cellsOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var verticesOnCell_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nCells*constants.maxEdges))
    var verticesOnEdge_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nEdges*constants.TWO))
    var edgesOnVertex_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices*constants.vertexDegree))
    var cellsOnVertex_in_copy : &int = [&int](constants.c.malloc([sizeof(int)] * constants.nVertices*constants.vertexDegree))
    var kiteAreasOnVertex_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nVertices*constants.vertexDegree))

    var u_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges * constants.nVertLevels))
    var v_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges * constants.nVertLevels))

    var w_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * (constants.nVertLevels + 1)))
    var pressure_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))

    var pressure_p_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var rho_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var theta_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var surface_pressure_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))


    --Now we copy the data into the arrays so they can be read into the netcdf files
    for i = 0, constants.nCells do
        latCell_in_copy[i] = cell_region[{i, 0}].lat
        lonCell_in_copy[i] = cell_region[{i, 0}].lon
        xCell_in_copy[i] = cell_region[{i, 0}].x
        yCell_in_copy[i] = cell_region[{i, 0}].y
        zCell_in_copy[i] = cell_region[{i, 0}].z
        indexToCellID_in_copy[i] = cell_region[{i, 0}].cellID
        meshDensity_in_copy[i] = cell_region[{i, 0}].meshDensity
        nEdgesOnCell_in_copy[i] = cell_region[{i, 0}].nEdgesOnCell
        areaCell_in_copy[i] = cell_region[{i, 0}].areaCell
        surface_pressure_in_copy[i]  = cell_region[{i, 0}].surface_pressure

        for j = 0, constants.maxEdges do
            edgesOnCell_in_copy[i*constants.maxEdges + j] = cell_region[{i, 0}].edgesOnCell[j]
            verticesOnCell_in_copy[i*constants.maxEdges + j] = cell_region[{i, 0}].verticesOnCell[j]
            cellsOnCell_in_copy[i*constants.maxEdges + j] = cell_region[{i, 0}].cellsOnCell[j]
        end

        for k = 0, constants.nVertLevels + 1 do
            if k < constants.nVertLevels + 1 then
                pressure_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].pressure
                pressure_p_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].pressure_p
                rho_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].rho
                theta_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].theta
                w_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].w
            else
                w_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].w
            end
        end
        
        --constants.cio.printf("Cell COPY : Cell ID %d, nEdgesOnCell is %d\n", indexToCellID_in_copy[i], nEdgesOnCell_in_copy[i])
    end

    for i = 0, constants.nEdges do
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
        u_in_copy[i]  = edge_region[{i, 0}].u

        v_in_copy[i]  = edge_region[{i, 0}].v

        for j = 0, constants.TWO do
            cellsOnEdge_in_copy[i*constants.TWO + j] = edge_region[{i, 0}].cellsOnEdge[j]
            verticesOnEdge_in_copy[i*constants.TWO + j] = edge_region[{i, 0}].verticesOnEdge[j]
        end

        for j = 0, constants.maxEdges2 do
            edgesOnEdge_in_copy[i*constants.maxEdges2 + j] = edge_region[{i, 0}].edgesOnEdge_ECP[j]
            weightsOnEdge_in_copy[i*constants.maxEdges2 + j] = edge_region[{i, 0}].weightsOnEdge[j]
        end
        for k = 0, constants.nVertLevels do
            u_in_copy[i*constants.nVertLevels + k] = edge_region[{i, k}].u
            v_in_copy[i*constants.nVertLevels + k] = edge_region[{i, k}].v
        end
    end

    for i = 0, constants.nVertices do
        latVertex_in_copy[i] = vertex_region[{i, 0}].lat
        lonVertex_in_copy[i] = vertex_region[{i, 0}].lon
        xVertex_in_copy[i] = vertex_region[{i, 0}].x
        yVertex_in_copy[i] = vertex_region[{i, 0}].y
        zVertex_in_copy[i] = vertex_region[{i, 0}].z
        indexToVertexID_in_copy[i] = vertex_region[{i, 0}].vertexID
        areaTriangle_in_copy[i] = vertex_region[{i, 0}].areaTriangle

        for j = 0, constants.vertexDegree do
            edgesOnVertex_in_copy[i*constants.vertexDegree + j] = vertex_region[{i, 0}].edgesOnVertex[j]
            cellsOnVertex_in_copy[i*constants.vertexDegree + j] = vertex_region[{i, 0}].cellsOnVertex[j]
            kiteAreasOnVertex_in_copy[i*constants.vertexDegree + j] = vertex_region[{i, 0}].kiteAreasOnVertex[j]
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
    put_var_double(ncid_copy, w_varid_copy, w_in_copy)
    put_var_double(ncid_copy, pressure_varid_copy, pressure_in_copy)
    put_var_double(ncid_copy, pressure_p_varid_copy, pressure_p_in_copy)

    put_var_double(ncid_copy, rho_varid_copy, rho_in_copy)
    put_var_double(ncid_copy, theta_varid_copy, theta_in_copy)
    put_var_double(ncid_copy, surface_pressure_varid_copy, surface_pressure_in_copy)



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
    put_var_double(ncid_copy, u_varid_copy, u_in_copy)
    put_var_double(ncid_copy, v_varid_copy, v_in_copy)

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
    constants.c.free(latCell_in_copy)
    constants.c.free(lonCell_in_copy)
    constants.c.free(meshDensity_in_copy)
    constants.c.free(xCell_in_copy)
    constants.c.free(yCell_in_copy)
    constants.c.free(zCell_in_copy)
    constants.c.free(indexToCellID_in_copy)
    constants.c.free(latEdge_in_copy)
    constants.c.free(lonEdge_in_copy)
    constants.c.free(xEdge_in_copy)
    constants.c.free(yEdge_in_copy)
    constants.c.free(zEdge_in_copy)
    constants.c.free(indexToEdgeID_in_copy)
    constants.c.free(latVertex_in_copy)
    constants.c.free(lonVertex_in_copy)
    constants.c.free(xVertex_in_copy)
    constants.c.free(yVertex_in_copy)
    constants.c.free(zVertex_in_copy)
    constants.c.free(indexToVertexID_in_copy)
    constants.c.free(cellsOnEdge_in_copy)
    constants.c.free(nEdgesOnCell_in_copy)
    constants.c.free(nEdgesOnEdge_in_copy)
    constants.c.free(edgesOnCell_in_copy)
    constants.c.free(edgesOnEdge_in_copy)
    constants.c.free(weightsOnEdge_in_copy)
    constants.c.free(dvEdge_in_copy)
    constants.c.free(dv1Edge_in_copy)
    constants.c.free(dv2Edge_in_copy)
    constants.c.free(dcEdge_in_copy)
    constants.c.free(angleEdge_in_copy)
    constants.c.free(areaCell_in_copy)
    constants.c.free(areaTriangle_in_copy)
    constants.c.free(cellsOnCell_in_copy)
    constants.c.free(verticesOnCell_in_copy)
    constants.c.free(verticesOnEdge_in_copy)
    constants.c.free(edgesOnVertex_in_copy)
    constants.c.free(cellsOnVertex_in_copy)
    constants.c.free(kiteAreasOnVertex_in_copy)

    constants.c.free(u_in_copy)
    constants.c.free(v_in_copy)
    constants.c.free(w_in_copy)
    constants.c.free(pressure_in_copy)
    constants.c.free(pressure_p_in_copy)
    constants.c.free(rho_in_copy)
    constants.c.free(theta_in_copy)
    constants.c.free(surface_pressure_in_copy)

    -- Close the file
    file_close(ncid_copy)
    constants.cio.printf("Successfully written netcdf file!\n")

end
