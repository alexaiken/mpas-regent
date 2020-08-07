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

    neighbor0 : int1d,
    neighbor1 : int1d,
    neighbor2 : int1d,
    neighbor3 : int1d,
    neighbor4 : int1d,
    neighbor5 : int1d,
    neighbor6 : int1d,
    neighbor7 : int1d,
    neighbor8 : int1d,
    neighbor9 : int1d,

    neighbor00 : int1d,
    neighbor11 : int1d,
    neighbor22 : int1d,
    neighbor33 : int1d,
    neighbor44 : int1d,
    neighbor55 : int1d,
    neighbor66 : int1d,
    neighbor77 : int1d,
    neighbor88 : int1d,
    neighbor99 : int1d,
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

  for i = 0, nCells do
    -- cell.cellsOnCell[0] contains an integer with the index of the cell neighbour. I would like cell.neighbour to point to that cell in the region
    -- we subtract 1 because the index spaces are 0-indexed but the cellIDs are 1-indexed
    cell_region[i].neighbor0 = cell_region[i].cellsOnCell[0] - 1
    cell_region[i].neighbor1 = cell_region[i].cellsOnCell[1] - 1
    cell_region[i].neighbor2 = cell_region[i].cellsOnCell[2] - 1
    cell_region[i].neighbor3 = cell_region[i].cellsOnCell[3] - 1
    cell_region[i].neighbor4 = cell_region[i].cellsOnCell[4] - 1
    cell_region[i].neighbor5 = cell_region[i].cellsOnCell[5] - 1
    cell_region[i].neighbor6 = cell_region[i].cellsOnCell[6] - 1
    cell_region[i].neighbor7 = cell_region[i].cellsOnCell[7] - 1
    cell_region[i].neighbor8 = cell_region[i].cellsOnCell[8] - 1
    cell_region[i].neighbor9 = cell_region[i].cellsOnCell[9] - 1

    cell_region[i].neighbor00 = cell_region[cell_region[i].cellsOnCell[0]].cellsOnCell[0] - 1
    cell_region[i].neighbor11 = cell_region[cell_region[i].cellsOnCell[1]].cellsOnCell[1] - 1
    cell_region[i].neighbor22 = cell_region[cell_region[i].cellsOnCell[2]].cellsOnCell[2] - 1
    cell_region[i].neighbor33 = cell_region[cell_region[i].cellsOnCell[3]].cellsOnCell[3] - 1
    cell_region[i].neighbor44 = cell_region[cell_region[i].cellsOnCell[4]].cellsOnCell[4] - 1
    cell_region[i].neighbor55 = cell_region[cell_region[i].cellsOnCell[5]].cellsOnCell[5] - 1
    cell_region[i].neighbor66 = cell_region[cell_region[i].cellsOnCell[6]].cellsOnCell[6] - 1
    cell_region[i].neighbor77 = cell_region[cell_region[i].cellsOnCell[7]].cellsOnCell[7] - 1
    cell_region[i].neighbor88 = cell_region[cell_region[i].cellsOnCell[8]].cellsOnCell[8] - 1
    cell_region[i].neighbor99 = cell_region[cell_region[i].cellsOnCell[9]].cellsOnCell[9] - 1
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

  -- Construct halo1
  var partition_s_1 = cell_partition_neighbor0 | cell_partition_neighbor1 | cell_partition_neighbor2 | cell_partition_neighbor3 | cell_partition_neighbor4 | cell_partition_neighbor5 | cell_partition_neighbor6 | cell_partition_neighbor7 | cell_partition_neighbor8 | cell_partition_neighbor9
  var partition_halo_1 = partition_s_1 - cell_partition_initial

  -- Repeat to create halo2
  var cell_partition_neighbor00 = image(cell_region, cell_partition_initial, cell_region.neighbor00)
  var cell_partition_neighbor11 = image(cell_region, cell_partition_initial, cell_region.neighbor11)
  var cell_partition_neighbor22 = image(cell_region, cell_partition_initial, cell_region.neighbor22)
  var cell_partition_neighbor33 = image(cell_region, cell_partition_initial, cell_region.neighbor33)
  var cell_partition_neighbor44 = image(cell_region, cell_partition_initial, cell_region.neighbor44)
  var cell_partition_neighbor55 = image(cell_region, cell_partition_initial, cell_region.neighbor55)
  var cell_partition_neighbor66 = image(cell_region, cell_partition_initial, cell_region.neighbor66)
  var cell_partition_neighbor77 = image(cell_region, cell_partition_initial, cell_region.neighbor77)
  var cell_partition_neighbor88 = image(cell_region, cell_partition_initial, cell_region.neighbor88)
  var cell_partition_neighbor99 = image(cell_region, cell_partition_initial, cell_region.neighbor99)

  var partition_s_2 = cell_partition_neighbor00 | cell_partition_neighbor11 | cell_partition_neighbor22 | cell_partition_neighbor33 | cell_partition_neighbor44 | cell_partition_neighbor55 | cell_partition_neighbor66 | cell_partition_neighbor77 | cell_partition_neighbor88 | cell_partition_neighbor99
  var partition_halo_2 = partition_s_2  - cell_partition_initial

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

  cio.printf("Successfully read file! \n")

end
regentlib.start(main)
