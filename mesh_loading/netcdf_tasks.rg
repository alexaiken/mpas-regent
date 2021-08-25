import "regent"

local constants = require("constants")

terralib.linklibrary("/home/arjunk1/spack/opt/spack/linux-ubuntu20.04-broadwell/gcc-9.3.0/netcdf-c-4.7.4-zgdvh4hxthdhb3mlsviwhgatvbfnslog/lib/libnetcdf.so")


-----------------------------------------------
----- TERRA WRAPPERS FOR NETCDF FUNCTIONS -----
-----------------------------------------------

--Terra function to open the constants.netcdf file and store NCID in the variable passed in.
terra open_file(ncid: &int, file_name: &int8)
    var retval = constants.netcdf.nc_open(file_name, constants.netcdf.NC_NOWRITE, ncid)
    if retval == 1 then
        constants.cio.printf("Error opening file %s \n", file_name)
    end
end

--Terra function to extract file information (no. of dims, etc).
terra file_inquiry(ncid: int, ndims_in: &int , nvars_in: &int, ngatts_in: &int, unlimdimid_in: &int)
    var retval = constants.netcdf.nc_inq(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in)
    if retval == 1 then
        constants.cio.printf("Error extracting file information of NCID %d \n", ncid)
    end
end

--Terra function to extract dimension ID
terra get_dimid(ncid: int, name: &int8, dimid: &int)
    var retval = constants.netcdf.nc_inq_dimid(ncid, name, dimid)
    if retval == 1 then
        constants.cio.printf("Error extracting dimension ID of dimension %s \n", name)
    end
end

--Terra function to extract dimension value
terra get_dimlength(ncid: int, dimid: int,  dimlength: &uint64)
    var retval = constants.netcdf.nc_inq_dimlen(ncid, dimid, dimlength)
    if retval == 1 then
        constants.cio.printf("Error extracting dimension length of dimID %d \n", dimid)
    end
end

--Terra function to get the variable ID, given the variable name.
terra get_varid(ncid: int, name: &int8,  varid: &int)
	var retval = constants.netcdf.nc_inq_varid(ncid, name, varid)
    if retval == 1 then
        constants.cio.printf("Error extracting variable ID of %s\n", name)
    end
end

--Terra function to get the variable values, given the variable ID: For variables with type double.
terra get_var_double(ncid: int, varid: int,  var_array_ptr: &double)
    var retval = constants.netcdf.nc_get_var_double(ncid, varid, var_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error extracting variable values of variable ID %d\n", varid)
    end
end

--Terra function to get the variable values, given the variable ID: For variables with type int.
terra get_var_int(ncid: int, varid: int,  var_array_ptr: &int)
	var retval = constants.netcdf.nc_get_var_int(ncid, varid, var_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error extracting variable values of variable ID %d\n", varid)
    end
end

--Terra function to get the global attribute length
terra get_global_att_len(ncid: int, name: &int8, varlen_ptr: &uint64)
	var retval = constants.netcdf.nc_inq_attlen(ncid, constants.netcdf.NC_GLOBAL, name, varlen_ptr)
    if retval == 1 then
        constants.cio.printf("Error extracting global attribute length of varname %s\n", name)
    end
end

--Terra function to get the global attributes, given the variable ID: For variables with type int.
terra get_global_att_double(ncid: int, name: &int8, att_array_ptr: &double)
	var retval = constants.netcdf.nc_get_att_double(ncid, constants.netcdf.NC_GLOBAL, name, att_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error extracting global attribute values of varname %s\n", name)
    end
end

--Terra function to get the global attributes, given the variable ID: For variables with type int.
terra get_global_att_int(ncid: int, name: &int8, att_array_ptr: &int)
	var retval = constants.netcdf.nc_get_att_int(ncid, constants.netcdf.NC_GLOBAL, name, att_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error extracting global attribute values of varname %s\n", name)
    end
end

--Terra function to get the global attributes, given the variable ID: For variables with type int.
terra get_global_att_text(ncid: int, name: &int8, att_array_ptr: &int8)
	var retval = constants.netcdf.nc_get_att_text(ncid, constants.netcdf.NC_GLOBAL, name, att_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error extracting global attribute values of varname %s\n", name)
    end
end

--Tera function to close file given NCID
terra file_close(ncid: int)
    var retval = constants.netcdf.nc_close(ncid)
    if retval == 1 then
        constants.cio.printf("Error closing file of NCID %d \n", ncid)
    end
end

--Terra function to close file given NCID
terra file_create(file_name: &int8, ncid_ptr: &int)
    var retval = constants.netcdf.nc_create(file_name, constants.netcdf.NC_CLOBBER, ncid_ptr)
    if retval == 1 then
        constants.cio.printf("Error creating file of name %s \n", file_name)
    end
end

--Terra function to define a dimension
terra define_dim(ncid: int, dim_name: &int8, dim_size: int, dim_id_ptr: &int)
    var retval = constants.netcdf.nc_def_dim(ncid, dim_name, dim_size, dim_id_ptr)
    if retval == 1 then
        constants.cio.printf("Error defining dimension of name %s \n", dim_name)
    end
end

--Terra function to define a variable
terra define_var(ncid: int, var_name: &int8, var_type: int, num_dims: int, dim_ids: &int, var_id_ptr: &int)
    var retval = constants.netcdf.nc_def_var(ncid, var_name, var_type, num_dims, dim_ids, var_id_ptr)
    if retval == 1 then
        constants.cio.printf("Error defining variable of name %s \n", var_name)
    end
end

--Terra function to end metadata reading
terra end_def(ncid: int)
	var retval = constants.netcdf.nc_enddef(ncid)
    if retval == 1 then
        constants.cio.printf("Error ending def of ncid %d \n", ncid)
    end
end

--Terra function to read variable
terra put_var_double(ncid: int, varid: int, var_array_ptr: &double)
    var retval = constants.netcdf.nc_put_var_double(ncid, varid, var_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error writing variable def of varid %d \n", varid)
    end
end

--Terra function to read variable
terra put_var_int(ncid: int, varid: int, var_array_ptr: &int)
    var retval = constants.netcdf.nc_put_var_int(ncid, varid, var_array_ptr)
    if retval == 1 then
        constants.cio.printf("Error writing variable def of varid %d \n", varid)
    end
end

task print_var_to_file(cell_region : region(ispace(int2d), cell_fs), edge_region : region(ispace(int2d), edge_fs), file_name: regentlib.string, j : int)
where
    reads writes (cell_region, edge_region)
do

    -- We create a netcdf file using the data in the regions, to test whether the data was written correctly.
    var ncid_copy = 65539

    -- temp location for modifying file names

    file_create(file_name, &ncid_copy)

    --Initialize the file's dimension variables
    var nCells_dimid_copy : int
    var nEdges_dimid_copy : int
    var nVertLevels_dimid_copy : int

    --Define the dimension variable
    define_dim(ncid_copy, "nCells", constants.nCells, &nCells_dimid_copy)
    define_dim(ncid_copy, "nEdges", constants.nEdges, &nEdges_dimid_copy)
    define_dim(ncid_copy, "nVertLevels", constants.nVertLevels, &nVertLevels_dimid_copy)
    var nCells_nVertLevels_dimids = array(nCells_dimid_copy, nVertLevels_dimid_copy)
    var nEdges_nVertLevels_dimids = array(nEdges_dimid_copy, nVertLevels_dimid_copy)

    --Initialize the variable IDs
    var pressure_varid_copy : int
    var rho_varid_copy : int
    var theta_varid_copy : int
    var surface_pressure_varid_copy : int
    var pressure_p_varid_copy : int
    var v_varid_copy : int
    
    --Define the variable IDs
    define_var(ncid_copy, "pressure_p", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &pressure_p_varid_copy)
    define_var(ncid_copy, "pressure", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &pressure_varid_copy)
    define_var(ncid_copy, "rho", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &rho_varid_copy)
    define_var(ncid_copy, "theta", constants.netcdf.NC_DOUBLE, 2, nCells_nVertLevels_dimids, &theta_varid_copy)
    define_var(ncid_copy, "surface_pressure", constants.netcdf.NC_DOUBLE, 1, &nCells_dimid_copy, &surface_pressure_varid_copy)
    define_var(ncid_copy, "v", constants.netcdf.NC_DOUBLE, 2, nEdges_nVertLevels_dimids, &v_varid_copy)

    --This function signals that we're done writing the metadata.
    end_def(ncid_copy)

    --Now define the new arrays to hold the data that will be put in the netcdf files
    var pressure_p_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var pressure_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var rho_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var theta_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells * constants.nVertLevels))
    var surface_pressure_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nCells))
    var v_in_copy : &double = [&double](constants.c.malloc([sizeof(double)] * constants.nEdges * constants.nVertLevels))

    --Now we copy the data into the arrays so they can be read into the netcdf files
    for i = 0, constants.nCells do
        for k = 0, constants.nVertLevels + 1 do
            pressure_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].pressure
            pressure_p_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].pressure_p
            rho_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].rho
            theta_in_copy[i*constants.nVertLevels + k] = cell_region[{i, k}].theta
        end
        surface_pressure_in_copy[i] = cell_region[{i, 0}].surface_pressure
    end

    for i = 0, constants.nEdges do
        for k = 0, constants.nVertLevels + 1 do
            v_in_copy[i*constants.nVertLevels + k] = edge_region[{i, k}].v
        end
    end

    --Now we put the data into the netcdf file.
    put_var_double(ncid_copy, pressure_varid_copy, pressure_in_copy)
    put_var_double(ncid_copy, pressure_p_varid_copy, pressure_p_in_copy)
    put_var_double(ncid_copy, rho_varid_copy, rho_in_copy)
    put_var_double(ncid_copy, theta_varid_copy, theta_in_copy)
    put_var_double(ncid_copy, surface_pressure_varid_copy, surface_pressure_in_copy)
    put_var_double(ncid_copy, v_varid_copy, v_in_copy)

    -- Lastly, we free the allocated memory for the 'copy' arrays
    constants.c.free(pressure_in_copy)
    constants.c.free(pressure_p_in_copy)
    constants.c.free(rho_in_copy)
    constants.c.free(theta_in_copy)
    constants.c.free(surface_pressure_in_copy)
    constants.c.free(v_in_copy)

    -- Close the file
    file_close(ncid_copy)

end
