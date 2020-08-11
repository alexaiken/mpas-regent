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
fspace cell_fs{
    cellID : int,
    partitionNumber: int1d,
    cellsOnCell : int[maxEdges],

    curr_neighbor_halo1 : int1d,
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
  var cell_id_space = ispace(int1d, nCells)
  var cell_region = region(cell_id_space, cell_fs)

  -----------------------------------------------
  ------- POPULATING REGION  --------
  -----------------------------------------------

  --Population region
  for i = 0, nCells do
    cell_region[i].cellID = indexToCellID_in[i]
    cell_region[i].partitionNumber = partition_array[i]

    for j = 0, maxEdges do
        cell_region[i].cellsOnCell[j] = cellsOnCell_in[i*maxEdges + j] --cell_region[i].cellsOnCell is a int[maxEdges]
    end

  end

  -----------------------------------------------
  ------- PARTITION REGIONS  --------
  -----------------------------------------------

  --Create initial partition based on METIS graph partition
  var color_space = ispace(int1d, NUM_PARTITIONS)
  var cell_partition_initial = partition(complete, cell_region.partitionNumber, color_space)

  --Define the s_1 partitions, that will eventually be the union of all of partitions created by the neighbour field in a subregion (as defined on call with alex)
  var partition_s_1 : partition(aliased, cell_region, ispace(int1d))

  for j = 0, maxEdges do
    for i = 0, nCells do
      -- populate the current neighbour
      -- we subtract 1 because the index spaces are 0-indexed but the cellIDs are 1-indexed

      cell_region[i].curr_neighbor_halo1 = cell_region[i].cellsOnCell[j] - 1
    end

    --create a partition based on the current neighbour, and add it to our union of neighbour partitions
    var cell_partition_curr_neighbor_halo1 = image(cell_region, cell_partition_initial, cell_region.curr_neighbor_halo1)
    partition_s_1 = partition_s_1 | cell_partition_curr_neighbor_halo1

  end

  var partition_halo_1 = partition_s_1 - cell_partition_initial


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

  --Print our last neighbour partition to check against the original partition
  --i = 0
  --for p in cell_partition_curr_neighbor_halo1.colors do
  --    var sub_region = cell_partition_curr_neighbor_halo1[p]
  --    cio.printf("Sub region %d\n", i)
  --    for cell in sub_region do
  --        cio.printf("%d\n", cell.cellID)
  --    end
  --    i=i+1
  --end

  cio.printf("Successfully read file! \n")

end
regentlib.start(main)
