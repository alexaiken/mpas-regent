import "regent"

local cio = terralib.includec("stdio.h")
local clib = terralib.includec("stdlib.h")

terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")
local netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

struct mesh_dimensions {
  nCells : uint64,
  nVertices : uint64,
  nEdges : uint64
}

local nCells = 2562
local nEdges = 7680
local nVertices = 5120
local maxEdges = 10
local maxEdges2 = 20
local TWO = 2
local vertexDegree = 3
local nVertLevels = 1

local FILE_NAME = "x1.2562.grid.nc"

-----------------------------------------------
------- FIELD SPACES FOR MESH ELEMENTS --------
-----------------------------------------------

--Field space with location of center of primal-mesh cells
fspace x_i {
  {lat, lon} : double
}

--location of center of dual-mesh cells
fspace x_v {
  {lat, lon} : double
}

--location of edge points where velocity is defined
fspace x_e {
  {lat, lon} : double
}
--distance bw neighboring primal-mesh centers, i.e. x_i locations
fspace d_e {
  distance : double
}

--distance between neighboring dual-mesh centers, i.e. x_v  locations
fspace l_e {
  distance : double
}


--Terra function to open the netcdf file and store NCID in the variable passed in.
terra open_file(ncid: &int, file_name: &int8)
	var retval = netcdf.nc_open(file_name, netcdf.NC_NOWRITE, ncid)
    if retval == 1 then 
        cio.printf("Error opening file %s \n", file_name)
    end
end

--Terra function to extract file information (no. of dims, etc).
terra file_inquiry(ncid: int, ndims_in: &int , nvars_in: &int, ngatts_in: &int, unlimdimid_in: &int)
	var retval = netcdf.nc_inq(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
    if retval == 1 then 
        cio.printf("Error extracting file information of NCID %d \n", ncid)
    end
end

--Terra function to extract dimension ID 
terra get_dimid(ncid: int, name: &int8, dimid: &int)
	var retval = netcdf.nc_inq_dimid(ncid, name, dimid)
    if retval == 1 then 
        cio.printf("Error extracting dimension ID of dimension %s \n", name)
    end
end

--Terra function to extract dimension value 
terra get_dimlength(ncid: int, dimid: int,  dimlength: &uint64)
	var retval = netcdf.nc_inq_dimlen(ncid, dimid, dimlength)
    if retval == 1 then 
        cio.printf("Error extracting dimension length of dimID %d \n", dimid)
    end
end


--Terra function to get the variable ID, given the variable name.
terra get_varid(ncid: int, name: &int8,  varid: &int)
	var retval = netcdf.nc_inq_varid(ncid, name, varid)
    if retval == 1 then 
        cio.printf("Error extracting variable ID of %s\n", name)
    end
end

--Terra function to get the variable values, given the variable ID: For variables with type double.
terra get_var_double(ncid: int, varid: int,  var_array_ptr: &double)
	var retval = netcdf.nc_get_var_double(ncid, varid, var_array_ptr)
    if retval == 1 then 
        cio.printf("Error extracting variable values of variable ID %d\n", varid)
    end
end

--Terra function to get the variable values, given the variable ID: For variables with type int.
terra get_var_int(ncid: int, varid: int,  var_array_ptr: &int)
	var retval = netcdf.nc_get_var_int(ncid, varid, var_array_ptr)
    if retval == 1 then 
        cio.printf("Error extracting variable values of variable ID %d\n", varid)
    end
end

--Tera function to close file given NCID
terra file_close(ncid: int)
	var retval = netcdf.nc_close(ncid)
    if retval == 1 then 
        cio.printf("Error closing file of NCID %d \n", ncid)
    end
end

task main()
	var ncid : int

	-- Open the file and store the NCID
	open_file(&ncid, FILE_NAME)

    	-- Define dim IDs
    	var nCells_dimid : int
    	var nEdges_dimid : int
    	var nVertices_dimid : int 

	-- Get dimIDs of all of the mesh_dimensions
	get_dimid(ncid, "nCells", &nCells_dimid)	
	get_dimid(ncid, "nVertices", &nVertices_dimid)
	get_dimid(ncid, "nEdges", &nEdges_dimid)

	var mesh_vars : mesh_dimensions
  	get_dimlength(ncid, nCells_dimid, &mesh_vars.nCells)
  	get_dimlength(ncid, nVertices_dimid, &mesh_vars.nVertices)
  	get_dimlength(ncid, nEdges_dimid, &mesh_vars.nEdges)

	---TEST
	cio.printf("nCells_dimid is %d\n nVertices_dimid is %d\n nEdges_dimid is %d\n", nCells_dimid, nVertices_dimid, nEdges_dimid)
	cio.printf("nCells is %u\n nVertices is %u\n nEdges is %u\n", mesh_vars.nCells, mesh_vars.nVertices, mesh_vars.nEdges)

    	-- Define the variable IDs
	var latCell_varid : int
	var lonCell_varid : int
    	var latVertex_varid : int
	var lonVertex_varid : int
	var latEdge_varid : int
	var lonEdge_varid : int
    	var dvEdge_varid : int
    	var dcEdge_varid : int

    	-- Define the data structures to store the variable values
	var latCell_in : double[nCells]
	var lonCell_in : double[nCells]
	var latVertex_in : double[nVertices]
	var lonVertex_in : double[nVertices]
    	var latEdge_in : double[nEdges]
	var lonEdge_in : double[nEdges]
    	var dvEdge_in : double[nEdges]
    	var dcEdge_in : double[nEdges]


    	-- Get the variable IDs of all the variables
	get_varid(ncid, "latCell", &latCell_varid)	
	get_varid(ncid, "lonCell", &lonCell_varid)
    	get_varid(ncid, "latVertex", &latVertex_varid)
	get_varid(ncid, "lonVertex", &lonVertex_varid)
    	get_varid(ncid, "latEdge", &latEdge_varid)
	get_varid(ncid, "lonEdge", &lonEdge_varid)
    	get_varid(ncid, "dvEdge", &dvEdge_varid)
    	get_varid(ncid, "dcEdge", &dcEdge_varid)
	
    	-- Get the variable values, given the variable IDs (To add: error checking)
	get_var_double(ncid, latCell_varid, &latCell_in[0])
	get_var_double(ncid, lonCell_varid, &lonCell_in[0])
    	get_var_double(ncid, latVertex_varid, &latVertex_in[0])
    	get_var_double(ncid, lonVertex_varid, &lonVertex_in[0])
    	get_var_double(ncid, latEdge_varid, &latEdge_in[0])
    	get_var_double(ncid, lonEdge_varid, &lonEdge_in[0])
    	get_var_double(ncid, dvEdge_varid, &dvEdge_in[0])
    	get_var_double(ncid, dcEdge_varid, &dcEdge_in[0])


    	--for i = 0, nEdges do 
	--	for j = 0, TWO do
	--		cio.printf("cellsOnEdge index %d %d is %d\n", i, j, cellsOnEdge_in[i][j])
    	--	end
	--end


    	-- Define region for cell latitude and longitude
    	var cell_id_region = region(ispace(ptr, mesh_vars.nCells), int)
   	var vertex_id_region = region(ispace(ptr, mesh_vars.nVertices), int)
    	var edge_id_region = region(ispace(ptr, mesh_vars.nEdges), int)

	--hexagon centers
    	var x_i_region = region(ispace(ptr, mesh_vars.nCells), x_i)
    	--triangle centers
   	var x_v_region = region(ispace(ptr, mesh_vars.nVertices), x_v)
    	--edge midpoints
    	var x_e_region = region(ispace(ptr, mesh_vars.nEdges), x_e)
    	--distance bw hexagon centers (x_i locations)
    	var d_e_region = region(ispace(ptr, mesh_vars.nEdges), d_e)
    	--distance bw triangle centers
    	var l_e_region = region(ispace(ptr, mesh_vars.nEdges), l_e)

    	-- Copy data from array into region
	var i = 0
	for elem in x_i_region do
		-- elem is a ptr(x_i, x_i_region) (i.e. a pointer to a x_i field space in the x_i region)
		elem.lat = latCell_in[i]
		elem.lon = lonCell_in[i]
		i = i+1
		cio.printf("Cell Lat/Long loop index %d: latCell is %f, lonCell is %f\n", i, elem.lat, elem.lon)
	end

	i = 0
	for elem in x_v_region do
		-- elem is a ptr(x_v, x_v_region) (i.e. a pointer to a x_v field space in the x_v region)
		elem.lat = latVertex_in[i]
		elem.lon = lonVertex_in[i]
		i = i+1
		cio.printf("Vertex Lat/Long loop index %d: latVertex is %f, onVertex is %f\n", i, elem.lat, elem.lon)
	end

	i = 0
	for elem in x_e_region do
		-- elem is a ptr(x_e, x_e_region) (i.e. a pointer to a x_e field space in the x_e_region)
		elem.lat = latEdge_in[i]
		elem.lon = lonEdge_in[i]
		i = i+1
		cio.printf("Edge Lat/long loop index %d: latEdge is %f, lonEdge is %f\n", i, elem.lat, elem.lon)
	end

	i = 0
		for elem in l_e_region do
		-- elem is a ptr(l_e, l_e_region) (i.e. a pointer to a l_e field space in the l_e_region)
		elem.distance = dvEdge_in[i]
		i = i+1
		cio.printf("Vertex Distance Loop %d: dvEdge / l_e is %f\n", i, elem.distance)
	end

	i = 0
		for elem in d_e_region do
		-- elem is a ptr(d_e, d_e_region) (i.e. a pointer to a d_e field space in the d_e_region)
		elem.distance = dcEdge_in[i]
		i = i+1
		cio.printf("Cell Distance Loop %d: dcEdge / d_e is %f\n", i, elem.distance)
	end

    	-- Close the file
	file_close(ncid)
	
end
regentlib.start(main)
