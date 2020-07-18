import "regent"

local cio = terralib.includec("stdio.h")
local clib = terralib.includec("stdlib.h


-- Link the .so and header files: If working in sherlock, these should be the same.
terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")
local netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

-- Function open netcdf file 
terra read_file(ncid: int)
        var retval = netcdf.nc_open("x1.2562.grid.nc", netcdf.NC_NOWRITE, &ncid)
        return ncid
end

-- Function to extract file information from netcdf file
terra file_inquiry(ncid: int, ndims_in: int , nvars_in: int, ngatts_in: int, unlimdimid_in: int)
        var retval = netcdf.nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in, &unlimdimid_in)
        return ndims_in
end

task main()
        var latCell_varid : int
        var lonCell_varid : int
        var ncid : int
        var ndims_in : int
        var nvars_in : int
        var ngatts_in : int
        var unlimdimid_in : int

        ncid = read_file(ncid)
        ndims_in = file_inquiry(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in)

        cio.printf("Success reading file: NCID is %d!\n File Dimensions are: %d.\n", ncid, ndims_in)
        
end
regentlib.start(main)
