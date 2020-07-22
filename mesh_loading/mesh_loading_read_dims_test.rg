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


local FILE_NAME = "x1.2562.grid.nc"

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
	cio.printf("CAST: nCells is %u\n nVertices is %u\n nEdges is %u\n", terralib.cast(int,mesh_vars.nCells), terralib.cast(int,mesh_vars.nVertices), terralib.cast(int,mesh_vars.nEdges))

	
    -- Define the variable IDs
	var latCell_varid : int
	var lonEdge_varid : int
	var indexToCellID_varid : int

    -- Define the data structures to store the variable values
	var latCell_in : double[mesh_vars.nCells]
	var lonEdge_in : double[nEdges]
	var indexToCellID_in : int[nCells]


    -- Get the variable IDs of all the variables
	get_varid(ncid, "latCell", &latCell_varid)	
	get_varid(ncid, "lonEdge", &lonEdge_varid)
	get_varid(ncid, "indexToCellID", &indexToCellID_varid)
	
    -- Get the variable values, given the variable IDs (To add: error checking)
	get_var_double(ncid, latCell_varid, &latCell_in[0])
	get_var_double(ncid, lonEdge_varid, &lonEdge_in[0])
    get_var_int(ncid, indexToCellID_varid, &indexToCellID_in[0])



    --for i = 0, nEdges do 
	--	cio.printf("lonEdge index %d is %f\n", i, lonEdge_in[i])
	--end


    -- Close the file
	file_close(ncid)
	
end
regentlib.start(main)
