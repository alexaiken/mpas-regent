import "regent"

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


--Terra function to open the netcdf file and store NCID in the variable passed in.
terra open_file(ncid: int, file_name: &int8)
	var retval = netcdf.nc_open(file_name, netcdf.NC_NOWRITE, &ncid)
	return ncid	
end

--Terra function to extract file information (no. of dims, etc).
terra file_inquiry(ncid: int, ndims_in: int , nvars_in: int, ngatts_in: int, unlimdimid_in: int)
	var retval = netcdf.nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in, &unlimdimid_in)
	return ndims_in
end

--Terra function to get the variable ID, given the variable name.
terra get_varid(ncid: int, name: &int8,  varid: int)
	var retval = netcdf.nc_inq_varid(ncid, name, &varid)
	return varid
end

--Terra function to get the variable values, given the variable ID: For variables with type double.
terra get_var_double(ncid: int, varid: int,  var_array_ptr: &double)
	var retval = netcdf.nc_get_var_double(ncid, varid, var_array_ptr)
	return retval
end

--Terra function to get the variable values, given the variable ID: For variables with type int.
terra get_var_int(ncid: int, varid: int,  var_array_ptr: &int)
	var retval = netcdf.nc_get_var_int(ncid, varid, var_array_ptr)
	return retval
end

terra file_close(ncid: int)
	var retval = netcdf.nc_close(ncid)
	return retval
end


task main()
	var ncid : int

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
    


    -- Define the arrays to store the variable values
	var latCell_in : double[nCells]
	var lonCell_in : double[nCells]
	var meshDensity_in : double[nCells]
	var xCell_in : double[nCells]
	var yCell_in : double[nCells]
	var zCell_in : double[nCells]
	var indexToCellID_in : int[nCells]
    var latEdge_in : double[nEdges]
	var lonEdge_in : double[nEdges]
	var xEdge_in : double[nEdges]
	var yEdge_in : double[nEdges]
	var zEdge_in : double[nEdges]
	var indexToEdgeID_in : int[nEdges]
    var latVertex_in : double[nVertices]
	var lonVertex_in : double[nVertices]
	var xVertex_in : double[nVertices]
	var yVertex_in : double[nVertices]
	var zVertex_in : double[nVertices]
	var indexToVertexID_in : int[nVertices]
    var cellsOnEdge_in : int[nEdges][TWO]
	var nEdgesOnCell_in : int[nCells]
	var nEdgesOnEdge_in : int[nEdges]
	var edgesOnCell_in : int[nCells][maxEdges]
	var edgesOnEdge_in : int[nEdges][maxEdges2]
	var weightsOnEdge_in : double[nEdges][maxEdges2]
    var dvEdge_in : double[nEdges]
	var dv1Edge_in : double[nEdges]
	var dv2Edge_in : double[nEdges]
	var dcEdge_in : double[nEdges]
	var angleEdge_in : double[nEdges]
	var areaCell_in : double[nCells]
	var areaTriangle_in : double[nVertices]
	var cellsOnCell_in : int[nCells][maxEdges]
	var verticesOnCell_in : int[nCells][maxEdges]
	var verticesOnEdge_in : int[nEdges][TWO]
	var edgesOnVertex_in : int[nVertices][vertexDegree]
	var cellsOnVertex_in : int[nVertices][vertexDegree]
	var kiteAreasOnVertex_in : double[nVertices][vertexDegree]

	
    -- Open the file and store the NCID
	ncid = open_file(ncid, FILE_NAME)

    -- Get the variable IDs of all the variables
	latCell_varid = get_varid(ncid, "latCell", latCell_varid)	
	lonCell_varid = get_varid(ncid, "lonCell", lonCell_varid)
    meshDensity_varid = get_varid(ncid, "meshDensity", meshDensity_varid)
	xCell_varid = get_varid(ncid, "xCell", xCell_varid) 
	yCell_varid = get_varid(ncid, "yCell", yCell_varid) 
	zCell_varid = get_varid(ncid, "zCell", zCell_varid)
	indexToCellID_varid = get_varid(ncid, "indexToCellID", indexToCellID_varid)
    latEdge_varid = get_varid(ncid, "latEdge", latEdge_varid)
	lonEdge_varid = get_varid(ncid, "lonEdge", lonEdge_varid)
	xEdge_varid = get_varid(ncid, "xEdge", xEdge_varid)
	yEdge_varid = get_varid(ncid, "yEdge", yEdge_varid)
	zEdge_varid = get_varid(ncid, "zEdge", zEdge_varid)
	indexToEdgeID_varid = get_varid(ncid, "indexToEdgeID", indexToEdgeID_varid)
    latVertex_varid = get_varid(ncid, "latVertex", latVertex_varid)
	lonVertex_varid = get_varid(ncid, "lonVertex", lonVertex_varid)
	xVertex_varid = get_varid(ncid, "xVertex", xVertex_varid)
	yVertex_varid = get_varid(ncid, "yVertex", yVertex_varid)
	zVertex_varid = get_varid(ncid, "zVertex", zVertex_varid)
	indexToVertexID_varid = get_varid(ncid, "indexToVertexID", indexToVertexID_varid)
    cellsOnEdge_varid = get_varid(ncid, "cellsOnEdge", cellsOnEdge_varid)
	nEdgesOnCell_varid = get_varid(ncid, "nEdgesOnCell", nEdgesOnCell_varid)
	nEdgesOnEdge_varid = get_varid(ncid, "nEdgesOnEdge", nEdgesOnEdge_varid)
	edgesOnCell_varid = get_varid(ncid, "edgesOnCell", edgesOnCell_varid)
	edgesOnEdge_varid = get_varid(ncid, "edgesOnEdge", edgesOnEdge_varid)
	weightsOnEdge_varid = get_varid(ncid, "weightsOnEdge", weightsOnEdge_varid)
    dvEdge_varid = get_varid(ncid, "dvEdge", dvEdge_varid)
	dv1Edge_varid = get_varid(ncid, "dv1Edge", dv1Edge_varid)
	dv2Edge_varid = get_varid(ncid, "dv2Edge", dv2Edge_varid)
	dcEdge_varid = get_varid(ncid, "dcEdge", dcEdge_varid)
	angleEdge_varid = get_varid(ncid, "angleEdge", angleEdge_varid)
	areaCell_varid = get_varid(ncid, "areaCell", areaCell_varid)
	areaTriangle_varid = get_varid(ncid, "areaTriangle", areaTriangle_varid)
	cellsOnCell_varid = get_varid(ncid, "cellsOnCell", cellsOnCell_varid)
	verticesOnCell_varid = get_varid(ncid, "verticesOnCell", verticesOnCell_varid)
	verticesOnEdge_varid = get_varid(ncid, "verticesOnEdge", verticesOnEdge_varid)
	edgesOnVertex_varid = get_varid(ncid, "edgesOnVertex", edgesOnVertex_varid)
	cellsOnVertex_varid = get_varid(ncid, "cellsOnVertex", cellsOnVertex_varid)
	kiteAreasOnVertex_varid = get_varid(ncid, "kiteAreasOnVertex", kiteAreasOnVertex_varid)

    -- Printing IDs out to test
	cio.printf("weightsOnEdge_varid is %d\n", weightsOnEdge_varid)
    cio.printf("verticesOnCell_varid is %d\n", verticesOnCell_varid)
    cio.printf("verticesOnEdge_varid is %d\n", verticesOnEdge_varid)
    cio.printf("cellsOnVertex_varid is %d\n", cellsOnVertex_varid)
    cio.printf("kiteAreasOnVertex_varid is %d\n", kiteAreasOnVertex_varid)
   
	-- Get the variable values, given the variable IDs
	var retval = get_var_double(ncid, latCell_varid, &latCell_in[0])
	retval = get_var_double(ncid, lonCell_varid, &lonCell_in[0])
    retval = get_var_double(ncid, meshDensity_varid, &meshDensity_in[0])
    retval = get_var_double(ncid, xCell_varid, &xCell_in[0])
    retval = get_var_double(ncid, yCell_varid, &yCell_in[0])
    retval = get_var_double(ncid, zCell_varid, &zCell_in[0])
    retval = get_var_int(ncid, indexToCellID_varid, &indexToCellID_in[0])
    retval = get_var_double(ncid, latEdge_varid, &latEdge_in[0])
    retval = get_var_double(ncid, lonEdge_varid, &lonEdge_in[0])
    retval = get_var_double(ncid, xEdge_varid, &xEdge_in[0])
    retval = get_var_double(ncid, yEdge_varid, &yEdge_in[0])
    retval = get_var_double(ncid, zEdge_varid, &zEdge_in[0])
    retval = get_var_int(ncid, indexToEdgeID_varid, &indexToEdgeID_in[0])
    retval = get_var_double(ncid, latVertex_varid, &latVertex_in[0])
    retval = get_var_double(ncid, lonVertex_varid, &lonVertex_in[0])
    retval = get_var_double(ncid, xVertex_varid, &xVertex_in[0])
    retval = get_var_double(ncid, yVertex_varid, &yVertex_in[0])
    retval = get_var_double(ncid, zVertex_varid, &zVertex_in[0])
    retval = get_var_int(ncid, indexToVertexID_varid, &indexToVertexID_in[0])
    retval = get_var_int(ncid, cellsOnEdge_varid, &cellsOnEdge_in[0][0])
    retval = get_var_int(ncid, nEdgesOnCell_varid, &nEdgesOnCell_in[0])
    retval = get_var_int(ncid, nEdgesOnEdge_varid, &nEdgesOnEdge_in[0])
    retval = get_var_int(ncid, edgesOnCell_varid, &edgesOnCell_in[0][0])
    retval = get_var_int(ncid, edgesOnEdge_varid, &edgesOnEdge_in[0][0])
    --retval = get_var_double(ncid, weightsOnEdge_varid, &weightsOnEdge_in[0][0])
    retval = get_var_double(ncid, dvEdge_varid, &dvEdge_in[0])
    retval = get_var_double(ncid, dv1Edge_varid, &dv1Edge_in[0])
    retval = get_var_double(ncid, dv2Edge_varid, &dv2Edge_in[0])
    retval = get_var_double(ncid, dcEdge_varid, &dcEdge_in[0])
    retval = get_var_double(ncid, angleEdge_varid, &angleEdge_in[0])
    retval = get_var_double(ncid, areaCell_varid, &areaCell_in[0])
    retval = get_var_double(ncid, areaTriangle_varid, &areaTriangle_in[0])
    retval = get_var_int(ncid, cellsOnCell_varid, &cellsOnCell_in[0][0])
    --retval = get_var_int(ncid, verticesOnCell_varid, &verticesOnCell_in[0][0])
    --retval = get_var_int(ncid, verticesOnEdge_varid, &verticesOnEdge_in[0][0])
    retval = get_var_int(ncid, edgesOnVertex_varid, &edgesOnVertex_in[0][0])
    --retval = get_var_int(ncid, cellsOnVertex_varid, &cellsOnVertex_in[0][0])
    --retval = get_var_double(ncid, kiteAreasOnVertex_varid, &kiteAreasOnVertex_in[0][0])

	retval = file_close(ncid)
	
	cio.printf("The 2nd element of the areaTriangle_in array is %f \n", areaTriangle_in[1])
    
end
regentlib.start(main)
