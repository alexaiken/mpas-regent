import "regent"

local c = regentlib.c
local cio = terralib.includec("stdio.h")
local clib = terralib.includec("stdlib.h")

terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")
local netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

-----------------------------------------------
----- TERRA WRAPPERS FOR NETCDF FUNCTIONS -----
-----------------------------------------------

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

--Terra function to get the global attribute length
terra get_global_att_len(ncid: int, name: &int8, varlen_ptr: &uint64)
	var retval = netcdf.nc_inq_attlen(ncid, netcdf.NC_GLOBAL, name, varlen_ptr)
    if retval == 1 then
        cio.printf("Error extracting global attribute length of varname %s\n", name)
    end
end

--Terra function to get the global attributes, given the variable ID: For variables with type int.
terra get_global_att_double(ncid: int, name: &int8, att_array_ptr: &double)
	var retval = netcdf.nc_get_att_double(ncid, netcdf.NC_GLOBAL, name, att_array_ptr)
    if retval == 1 then
        cio.printf("Error extracting global attribute values of varname %s\n", name)
    end
end

--Terra function to get the global attributes, given the variable ID: For variables with type int.
terra get_global_att_int(ncid: int, name: &int8, att_array_ptr: &int)
	var retval = netcdf.nc_get_att_int(ncid, netcdf.NC_GLOBAL, name, att_array_ptr)
    if retval == 1 then
        cio.printf("Error extracting global attribute values of varname %s\n", name)
    end
end

--Terra function to get the global attributes, given the variable ID: For variables with type int.
terra get_global_att_text(ncid: int, name: &int8, att_array_ptr: &int8)
	var retval = netcdf.nc_get_att_text(ncid, netcdf.NC_GLOBAL, name, att_array_ptr)
    if retval == 1 then
        cio.printf("Error extracting global attribute values of varname %s\n", name)
    end
end

--Tera function to close file given NCID
terra file_close(ncid: int)
    var retval = netcdf.nc_close(ncid)
    if retval == 1 then
        cio.printf("Error closing file of NCID %d \n", ncid)
    end
end

--Terra function to close file given NCID
terra file_create(file_name: &int8, ncid_ptr: &int)
    var retval = netcdf.nc_create(file_name, netcdf.NC_CLOBBER, ncid_ptr)
    if retval == 1 then
        cio.printf("Error creating file of name %s \n", file_name)
    end
end

--Terra function to define a dimension
terra define_dim(ncid: int, dim_name: &int8, dim_size: int, dim_id_ptr: &int)
    var retval = netcdf.nc_def_dim(ncid, dim_name, dim_size, dim_id_ptr)
    if retval == 1 then
        cio.printf("Error defining dimension of name %s \n", dim_name)
    end
end

--Terra function to define a variable
terra define_var(ncid: int, var_name: &int8, var_type: int, num_dims: int, dim_ids: &int, var_id_ptr: &int)
    var retval = netcdf.nc_def_var(ncid, var_name, var_type, num_dims, dim_ids, var_id_ptr)
    if retval == 1 then
        cio.printf("Error defining variable of name %s \n", var_name)
    end
end

--Terra function to end metadata reading
terra end_def(ncid: int)
	var retval = netcdf.nc_enddef(ncid)
    if retval == 1 then
        cio.printf("Error ending def of ncid %d \n", ncid)
    end
end

--Terra function to read variable
terra put_var_double(ncid: int, varid: int, var_array_ptr: &double)
    var retval = netcdf.nc_put_var_double(ncid, varid, var_array_ptr)
    if retval == 1 then
        cio.printf("Error writing variable def of varid %d \n", varid)
    end
end

--Terra function to read variable
terra put_var_int(ncid: int, varid: int, var_array_ptr: &int)
    var retval = netcdf.nc_put_var_int(ncid, varid, var_array_ptr)
    if retval == 1 then
        cio.printf("Error writing variable def of varid %d \n", varid)
    end
end
