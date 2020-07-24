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

--edges that outline P_i
fspace edges_on_cell {
	edges : int[maxEdges]
}

--edges that intersect at vertex of D_v
fspace edges_on_vertex {
	edges : int[vertexDegree]
}

--cell IDs of two cells sharing an edge
fspace cells_on_edge {
	{cell1, cell2} : int
}

--cell ids of the primal cells sharing a vertex
fspace cells_on_vertex {
	{cell1, cell2, cell3} : int,
}

--vertex ids for triangles centered on the points of an edge
fspace vertices_on_edge{
	{d1, d2} : int
}

--vertex ids for the vertices outlining P_i
fspace vertices_on_cell {
	vertices : int[maxEdges]
}

--ECP: edges of cell pair meeting at edge e
fspace ecp_edges_of_cell_pair {
	{ec1, ec2} : edges_on_cell
}

--EVC: edge pair associated with vertex v and mesh cell i
fspace evc_edge_pair {
	{e1, e2} : int
}

--edge on primal (hexagonal) mesh made up of two points in the x_v region (centers of dual cells)
fspace edge {
	edgeID : int,
	center : x_e,
	angle : double,
	nEdgesOnEdge : int,
	dvEdge : l_e,
	dv1Edge : double,
	dv2Edge : double,
	dcEdge : d_e,
	cells : cells_on_edge,
	vertices : vertices_on_edge
	--cells : ptr(cells_on_edge, rce)
	--vertices : ptr(vertices_on_edge, rve)
}

--a triangle cell on the dual-mesh (aka vertex)
fspace D_v {
	vertexID : int,
	center : x_v,
	areaTriangle : double,
	edge_set : edges_on_vertex,
	cells : cells_on_vertex
}

--a hexagon cell on the primal-mesh (aka cell)
fspace P_i {
	cellID : int,
	center : x_i,
	area : double,
	meshDensity : double,
	nEdgesOnCell : int,
	--cellsOnCells : cells_on_cell,
	edge_set : edges_on_cell,
	vertices : vertices_on_cell
}



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

--Tera function to close file given NCID
terra file_close(ncid: int)
	var retval = netcdf.nc_close(ncid)
	if retval == 1 then 
		cio.printf("Error closing file of NCID %d \n", ncid)
	end
end

--terra copy_2d_to_1d(array_ptr: &int, index1_lim: int, index2_lim: int)
--   var temp : int[index1_lim*index2_lim]
--    var k = 0
--        for i = 0, index1_lim do 
--           for j = 0, index2_lim do
--                temp[k] = @array_ptr[i][j]
--                k = k+1
--            end
--        end
--    return temp
--end


task main()
	var ncid : int

	-- Open the file and store the NCID
	open_file(&ncid, FILE_NAME)

	---------------------------------------------------------------------------
	----- Code to read dimensions from file directly: temporarily removed -----
	--------------------------------------------------------------------------- 

	-- Define dim IDs
	--var nCells_dimid : int
	--var nEdges_dimid : int
	--var nVertices_dimid : int 

	-- Get dimIDs of all of the mesh_dimensions
	--get_dimid(ncid, "nCells", &nCells_dimid)	
	--get_dimid(ncid, "nVertices", &nVertices_dimid)
	--get_dimid(ncid, "nEdges", &nEdges_dimid)

	var mesh_vars : mesh_dimensions
	--get_dimlength(ncid, nCells_dimid, &mesh_vars.nCells)
	--get_dimlength(ncid, nVertices_dimid, &mesh_vars.nVertices)
	--get_dimlength(ncid, nEdges_dimid, &mesh_vars.nEdges)

	---TEST
	--cio.printf("nCells_dimid is %d\n nVertices_dimid is %d\n nEdges_dimid is %d\n", nCells_dimid, nVertices_dimid, nEdges_dimid)
	--cio.printf("nCells is %u\n nVertices is %u\n nEdges is %u\n", mesh_vars.nCells, mesh_vars.nVertices, mesh_vars.nEdges)

	-- Define the variable IDs
	var latCell_varid : int
	var lonCell_varid : int
	var latVertex_varid : int
	var lonVertex_varid : int
	var latEdge_varid : int
	var lonEdge_varid : int
	var dvEdge_varid : int
	var dcEdge_varid : int
	var edgesOnCell_varid : int
	var edgesOnVertex_varid : int
	var cellsOnEdge_varid : int
	var cellsOnVertex_varid : int
	var verticesOnCell_varid : int
	var verticesOnEdge_varid : int
	var angleEdge_varid : int
	var areaCell_varid : int
	var nEdgesOnCell_varid : int
	var meshDensity_varid : int
	var nEdgesOnEdge_varid : int
	var dv1Edge_varid : int
	var dv2Edge_varid : int
	var areaTriangle_varid : int


	-- Define the data structures to store the variable values
	var latCell_in : double[nCells]
	var lonCell_in : double[nCells]
	var latVertex_in : double[nVertices]
	var lonVertex_in : double[nVertices]
	var latEdge_in : double[nEdges]
	var lonEdge_in : double[nEdges]
	var dvEdge_in : double[nEdges]
	var dcEdge_in : double[nEdges]
	var edgesOnCell_in : int[nCells][maxEdges]
	var edgesOnVertex_in : int[nVertices][vertexDegree]
	var cellsOnEdge_in : int[nEdges][TWO]
	var cellsOnVertex_in : int[nVertices][vertexDegree]
	var verticesOnCell_in : int[nCells][maxEdges]
	var verticesOnEdge_in : int[nEdges][TWO]
	var angleEdge_in : double[nEdges]
	var areaCell_in : double[nCells]
	var nEdgesOnCell_in : int[nCells]
	var meshDensity_in : double[nCells]
	var nEdgesOnEdge_in : int[nEdges]
	var dv1Edge_in : double[nEdges]
	var dv2Edge_in : double[nEdges]
	var areaTriangle_in : double[nVertices]

	-- Get the variable IDs of all the variables
	get_varid(ncid, "latCell", &latCell_varid)	
	get_varid(ncid, "lonCell", &lonCell_varid)
	get_varid(ncid, "latVertex", &latVertex_varid)
	get_varid(ncid, "lonVertex", &lonVertex_varid)
	get_varid(ncid, "latEdge", &latEdge_varid)
	get_varid(ncid, "lonEdge", &lonEdge_varid)
	get_varid(ncid, "dvEdge", &dvEdge_varid)
	get_varid(ncid, "dcEdge", &dcEdge_varid)
	get_varid(ncid, "edgesOnCell", &edgesOnCell_varid)
	get_varid(ncid, "edgesOnVertex", &edgesOnVertex_varid)
	get_varid(ncid, "cellsOnEdge", &cellsOnEdge_varid)
	get_varid(ncid, "cellsOnVertex", &cellsOnVertex_varid)
	get_varid(ncid, "verticesOnCell", &verticesOnCell_varid)
	get_varid(ncid, "verticesOnEdge", &verticesOnEdge_varid)
	get_varid(ncid, "angleEdge", &angleEdge_varid)
	get_varid(ncid, "areaCell", &areaCell_varid)
	get_varid(ncid, "nEdgesOnCell", &nEdgesOnCell_varid)
	get_varid(ncid, "meshDensity", &meshDensity_varid)
	get_varid(ncid, "nEdgesOnEdge", &nEdgesOnEdge_varid)
	get_varid(ncid, "dv1Edge", &dv1Edge_varid)
	get_varid(ncid, "dv2Edge", &dv2Edge_varid)
	get_varid(ncid, "areaTriangle", &areaTriangle_varid)


	-- Get the variable values, given the variable IDs (To add: error checking)
	get_var_double(ncid, latCell_varid, &latCell_in[0])
	get_var_double(ncid, lonCell_varid, &lonCell_in[0])
	get_var_double(ncid, latVertex_varid, &latVertex_in[0])
	get_var_double(ncid, lonVertex_varid, &lonVertex_in[0])
	get_var_double(ncid, latEdge_varid, &latEdge_in[0])
	get_var_double(ncid, lonEdge_varid, &lonEdge_in[0])
	get_var_double(ncid, dvEdge_varid, &dvEdge_in[0])
	get_var_double(ncid, dcEdge_varid, &dcEdge_in[0])
	get_var_int(ncid, edgesOnCell_varid, &edgesOnCell_in[0][0])
	get_var_int(ncid, edgesOnVertex_varid, &edgesOnVertex_in[0][0])
	get_var_int(ncid, cellsOnEdge_varid, &cellsOnEdge_in[0][0])
	get_var_int(ncid, cellsOnVertex_varid, &cellsOnVertex_in[0][0])
	get_var_int(ncid, verticesOnCell_varid, &verticesOnCell_in[0][0])
	get_var_int(ncid, verticesOnEdge_varid, &verticesOnEdge_in[0][0])
	get_var_double(ncid, angleEdge_varid, &angleEdge_in[0])
	get_var_double(ncid, areaCell_varid, &areaCell_in[0])
	get_var_int(ncid, nEdgesOnCell_varid, &nEdgesOnCell_in[0])
	get_var_double(ncid, meshDensity_varid, &meshDensity_in[0])
	get_var_int(ncid, nEdgesOnEdge_varid, &nEdgesOnEdge_in[0])
	get_var_double(ncid, dv1Edge_varid, &dv1Edge_in[0])
	get_var_double(ncid, dv2Edge_varid, &dv2Edge_in[0])
	get_var_double(ncid, areaTriangle_varid, &areaTriangle_in[0])


    -----------------------------------------------------------------
    ----- Converting incorrectly indexed 2d array into 1d array -----
    ----------------------------------------------------------------- 
    -- very ugly/inefficient code, I know: tried to decompose this into a function / find another way around this
	var EC_temp : int[nCells*maxEdges]
	var k = 0
	for i = 0, maxEdges do 
		for j = 0, nCells do
			EC_temp[k] = edgesOnCell_in[i][j]
			--cio.printf("edgesOnCell: Edge index %d, Cell Index %d is value:%d\n", i, j, edgesOnCell_in[i][j])
			--cio.printf("edgesOnCell Array index %d is value: %d \n", k, edgesOnCell_in[i][j])
			k = k+1
		end
	end

	var EV_temp : int[nVertices*vertexDegree]
	k = 0
	for i = 0, vertexDegree do 
		for j = 0, nVertices do
			EV_temp[k] = edgesOnVertex_in[i][j]
			--cios.printf("edgesOnVertex: Edge index %d, Vertex Index %d is value:%d\n", i, j, edgesOnVertex_in[i][j])
			--cio.printf("edgesOnVertex Array index %d is value: %d \n", k, edgesOnVertex_in[i][j])
			k = k+1
		end
	end

	var CE_temp : int[nEdges*TWO]
	k = 0
	for i = 0, TWO do 
		for j = 0, nEdges do
			CE_temp[k] = cellsOnEdge_in[i][j]
			--cios.printf("cellsOnEdge: Cell index %d, Edge Index %d is value:%d\n", i, j, cellsOnEdge_in[i][j])
			--cio.printf("cellsOnEdge Array index %d is value: %d \n", k, cellsOnEdge_in[i][j])
			k = k+1
		end
	end

	--Failed attempt to decompose this routine into a task
	--var CV_temp = copy_2d_to_1d(&cellsOnVertex_in[0][0], vertexDegree, nVertices)

	var CV_temp : int[nVertices*vertexDegree]
	k = 0
	for i = 0, vertexDegree do 
		for j = 0, nVertices do
			CV_temp[k] = cellsOnVertex_in[i][j]
			--cios.printf("cellsOnVertex: Cell index %d, Vertex Index %d is value:%d\n", i, j, cellsOnVertex_in[i][j])
			--cio.printf("cellsOnVertex Array index %d is value: %d \n", k, cellsOnVertex_in[i][j])
			k = k+1
		end
	end

	var VE_temp : int[nEdges*TWO]
	k = 0
	for i = 0, TWO do 
		for j = 0, nEdges do
			VE_temp[k] = verticesOnEdge_in[i][j]
			--cios.printf("verticesOnEdge: Vertex index %d, Edge Index %d is value:%d\n", i, j, verticesOnEdge_in[i][j])
			--cio.printf("verticesOnEdge Array index %d is value: %d \n", k, verticesOnEdge_in[i][j])
			k = k+1
		end
	end

	var VI_temp : int[nCells*maxEdges]
	k = 0
	for i = 0, maxEdges do 
		for j = 0, nCells do
			VI_temp[k] = verticesOnCell_in[i][j]
			--cios.printf("verticesOnCell: Vertex index %d, Cell Index %d is value:%d\n", i, j, verticesOnCell_in[i][j])
			--cio.printf("verticesOnCell Array index %d is value: %d \n", k, verticesOnCell_in[i][j])
			k = k+1
		end
	end



	-- Define index spaces for cell IDs, vertex IDs and edge IDs
	var cell_id_space = ispace(int1d, nCells)
	var vertex_id_space = ispace(int1d, nVertices)
	var edge_id_space = ispace(int1d, nEdges)
	var vertex_cell_space = ispace(int2d, {nVertices, nCells})


	--hexagon centers
	var x_i_region = region(cell_id_space, x_i)
	--triangle centers
	var x_v_region = region(vertex_id_space, x_v)
	--edge midpoints
	var x_e_region = region(edge_id_space, x_e)
	--distance bw hexagon centers (x_i locations)
	var d_e_region = region(edge_id_space, d_e)
	--distance bw triangle centers
	var l_e_region = region(edge_id_space, l_e)


	--edges on a cell defining P_i
	var EC_region = region(cell_id_space, edges_on_cell)
	-- edges on a vertex defining D_v
	var EV_region = region(vertex_id_space, edges_on_vertex)
	-- cells on an edge (primal mesh cells sharing an edge) + edges on the cell pair
	var CE_region = region(edge_id_space, cells_on_edge)
	-- cells on a vertex (dual mesh cells s)
	var CV_region = region(vertex_id_space, cells_on_vertex)
	-- vertices on edge (dual mesh cells sharing a primal edge)
	var VE_region = region(edge_id_space, vertices_on_edge)
	-- vertices on cell (dual mesh cells centered on vertices of a primal cell)
	var VI_region = region(cell_id_space, vertices_on_cell)


	-- edge pair associated w a vertex and a cell
	-- edge pair associated w a vertex and a cell
	var ECP_region = region(edge_id_space, ecp_edges_of_cell_pair)
	var EVC_region = region(vertex_cell_space, evc_edge_pair)

	var edge_region = region(edge_id_space, edge)
	var cell_region = region(cell_id_space, P_i)
	var vertex_region = region(vertex_id_space, D_v)


	-----------------------------------------------
	---------- DATA TRANSFER TO REGIONS -----------
	-----------------------------------------------

	-- Copy data from array into region
	var i = 0
	for elem in x_i_region do
		elem.lat = latCell_in[i]
		elem.lon = lonCell_in[i]
		i = i+1
		--cio.printf("Cell Lat/Long loop index %d: latCell is %f, lonCell is %f\n", i, elem.lat, elem.lon)
	end

	i = 0
	for elem in x_v_region do
		elem.lat = latVertex_in[i]
		elem.lon = lonVertex_in[i]
		i = i+1
		--cio.printf("Vertex Lat/Long loop index %d: latVertex is %f, onVertex is %f\n", i, elem.lat, elem.lon)
	end

	i = 0
	for elem in x_e_region do
		elem.lat = latEdge_in[i]
		elem.lon = lonEdge_in[i]
		i = i+1
		--cio.printf("Edge Lat/long loop index %d: latEdge is %f, lonEdge is %f\n", i, elem.lat, elem.lon)
	end

	i = 0
	for elem in l_e_region do
		elem.distance = dvEdge_in[i]
		i = i+1
		--cio.printf("Vertex Distance Loop %d: dvEdge / l_e is %f\n", i, elem.distance)
	end

	i = 0
	for elem in d_e_region do
		elem.distance = dcEdge_in[i]
		i = i+1
		--cio.printf("Cell Distance Loop %d: dcEdge / d_e is %f\n", i, elem.distance)
	end

	i = 0
	for elem in EC_region do
		for j = 0, maxEdges do
			elem.edges[j] = EC_temp[i*maxEdges + j] --now, elem.edges is a int[maxEdges]   
			--cio.printf("Edge on Cell : Cell %d, Edge %d: edge index is %d\n", i, j, elem.edges[j])
		end           
		i = i+1 
	end

	i = 0
	for elem in EV_region do
		for j = 0, vertexDegree do
			elem.edges[j] = EV_temp[i*vertexDegree + j] --now, elem.edges is a int[vertexDegree]   
			--cio.printf("Edge on Vertex : Vertex %d, Edge %d: edge index is %d\n", i, j, elem.edges[j])
		end           
		i = i+1 
	end

	i = 0
	for elem in CE_region do
		elem.cell1 = CE_temp[i*TWO] --now, elem.edges is a int[vertexDegree]
		elem.cell2 = CE_temp[i*TWO + 1]
		--cio.printf("Cells on Edge : Edge %d: Cell 1 is %d\n, Cell 2 is %d\n", i, elem.cell1, elem.cell2)
		i = i+1 
	end

	i = 0
	for elem in CV_region do
		elem.cell1 = CV_temp[i*vertexDegree]
		elem.cell2 = CV_temp[i*vertexDegree + 1]
		elem.cell3 = CV_temp[i*vertexDegree + 2]
		--cio.printf("Cells on Vertex : Vertex %d: Cell 1 is %d, Cell 2 is %d, Cell 3 is %d\n", i, elem.cell1, elem.cell2, elem.cell3)
		i = i+1
	end

	i = 0
	for elem in VE_region do
		elem.d1 = VE_temp[i*TWO] 
		elem.d2 = VE_temp[i*TWO + 1]
		--cio.printf("VerticesOnEdge : Edge %d: Vertex 1 is %d, Vertex 2 is %d\n", i, elem.d1, elem.d2)
		i = i+1
	end

	i = 0
	for elem in VI_region do
		for j = 0, maxEdges do
			elem.vertices[j] = VI_temp[i*maxEdges + j] --now, elem.vertices is a int[maxEdges]   
			--cio.printf("Vertices on Cell : Cell %d, Vertex %d: Vertex index is %d\n", i, j, elem.vertices[j])
		end  
		i = i+1
	end

	i = 0
	for elem in VI_region do
		for j = 0, maxEdges do
			elem.vertices[j] = VI_temp[i*maxEdges + j] --now, elem.vertices is a int[maxEdges]   
			--cio.printf("Vertices on Cell : Cell %d, Vertex %d: Vertex index is %d\n", i, j, elem.vertices[j])
		end  
		i = i+1
	end

	i = 0 -- Test this code
	for elem in ECP_region do
		var cell1 = CE_region[i].cell1
		var cell2 = CE_region[i].cell2
		elem.ec1 = EC_region[cell1]
		elem.ec2 = EC_region[cell2]
		i = i+1
	end

	i = 0
	for elem in edge_region do
		edge_region[i].edgeID = i
		edge_region[i].center = x_e_region[i]
		edge_region[i].nEdgesOnEdge = nEdgesOnEdge_in[i]
		edge_region[i].dvEdge = l_e_region[i]
		edge_region[i].dv1Edge = dv1Edge_in[i]
		edge_region[i].dv2Edge = dv2Edge_in[i]
		edge_region[i].dcEdge = d_e_region[i]
		edge_region[i].angle = angleEdge_in[i]
		edge_region[i].cells = CE_region[i]
		edge_region[i].vertices = VE_region[i]
		--cio.printf("Edge: ID is %d, dvEdge is %f, dv1Edge is %f, AngleEdge is %f \n", i, edge_region[i].edgeID, edge_region[i].dvEdge.distance, edge_region[i].dv1Edge, edge_region[i].angle)
		i = i+1
	end

	i = 0
	for elem in vertex_region do
		vertex_region[i].vertexID = i
		vertex_region[i].center = x_v_region[i]
		vertex_region[i].areaTriangle = areaTriangle_in[i]
		vertex_region[i].edge_set = EV_region[i]
		vertex_region[i].cells = CV_region[i]
		cio.printf("Vertex ID is %d, Latitude:%f, Longitude: %f, areaTriangle: %f \n", vertex_region[i].vertexID, vertex_region[i].center.lat, vertex_region[i].center.lon, vertex_region[i].areaTriangle)
		i = i+1
	end

	i = 0
	for elem in cell_region do
		cell_region[i].cellID = i
		cell_region[i].center = x_i_region[i]
		cell_region[i].area = areaCell_in[i]
		cell_region[i].nEdgesOnCell = nEdgesOnCell_in[i]
		cell_region[i].meshDensity = meshDensity_in[i]

		cell_region[i].edge_set = EC_region[i]
		cell_region[i].vertices = VI_region[i]
		--cio.printf("Cell : Cell %d, ID is %d, Area is %f, nEdges is %d, meshDensity is %f\n", i, cell_region[i].cellID, cell_region[i].area, cell_region[i].nEdgesOnCell, cell_region[i].meshDensity)
		i = i+1
	end

	-- Close the file
	file_close(ncid)
	
end
regentlib.start(main)
