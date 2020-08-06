import "regent"

require "data_structures"
require "netcdf_tasks"

local c = regentlib.c
local cio = terralib.includec("stdio.h")
local clib = terralib.includec("stdlib.h")

terralib.linklibrary("/share/software/user/open/netcdf/4.4.1.1/lib/libnetcdf.so")
local netcdf = terralib.includec("/share/software/user/open/netcdf/4.4.1.1/include/netcdf.h")

local nCells = 2562
local maxEdges = 10

local FILE_NAME = "x1.2562.grid.nc"
local GRAPH_FILE_NAME = "x1.2562.graph.info.part.16"
local MAXCHAR = 5
local NUM_PARTITIONS = 16


-----------------------------------------------
------- FIELD SPACES FOR MESH ELEMENTS --------
-----------------------------------------------

--a hexagonal/primal cell (aka cell)
fspace cell_fs(rc : region(cell_fs(wild))){
    cellID : int,
    partitionNumber: int1d,
    cellsOnCell : int[maxEdges],
    neighbour : ptr(cell_fs(wild), rc),
}


-----------------------------------------------
------- TERRA WRAPPER FUNCTIONS  --------
-----------------------------------------------

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

  -----------------------------------------------
  ------- READING NETCDF DATA  --------
  -----------------------------------------------

  cio.printf("Starting to read file... \n")
  var ncid : int

  -- Open the file and store the NCID
  open_file(&ncid, FILE_NAME)

  -- Read in variables: can ignore this code
  var indexToCellID_varid : int
  var indexToCellID_in : &int = [&int](c.malloc([sizeof(int)] * nCells))
  get_varid(ncid, "indexToCellID", &indexToCellID_varid)
  get_var_int(ncid, indexToCellID_varid, indexToCellID_in)

  var cellsOnCell_varid : int
  var cellsOnCell_in : &int = [&int](c.malloc([sizeof(int)] * nCells*maxEdges))
  get_varid(ncid, "cellsOnCell", &cellsOnCell_varid)
  get_var_int(ncid, cellsOnCell_varid, cellsOnCell_in)
  var partition_array = read_file(GRAPH_FILE_NAME)

  -----------------------------------------------
  ------- CREATING INDEX SPACES AND REGION  --------
  -----------------------------------------------

  --Creating index space and region
  var cell_id_space = ispace(ptr, nCells)
  var cell_region = region(cell_id_space, cell_fs(wild))

  -----------------------------------------------
  ------- POPULATING REGION  --------
  -----------------------------------------------

  --Population region
  var i = 0
  for cell in cell_region do
    cell.cellID = indexToCellID_in[i]
    cell.partitionNumber = partition_array[i]

    for j = 0, maxEdges do
        cell.cellsOnCell[j] = cellsOnCell_in[i*maxEdges + j] --cell_region[i].cellsOnCell is a int[maxEdges]
    end
    -- cell.cellsOnCell[0] contains an integer with the index of the cell neighbour. I would like cell.neighbour to point to that cell in the region
    cell.neighbour = cell
    i = i+1
  end


  --var partition_is = ispace(int1d, NUM_PARTITIONS)
  --var pc_initial = partition(complete, cell_region.partitionNumber, partition_is)

  --var p = image(parent_region, source_partition, data_region.field)
  --var pc_neighbor0 = image(cell_region, cell_partition, cell_region.neighbour0)

  --Test code by printing out partitions
  --var i = 0
  --for p in partition_is do
  --    var sub_region = cell_partition[p]
  --    cio.printf("Sub region %d\n", i)
  --    for cell in sub_region do
  --        cio.printf("partitionNumber is %d\n", cell.partitionNumber)
  --    end
  --    i=i+1
  --end

  cio.printf("Successfully read file! \n")

end
regentlib.start(main)
