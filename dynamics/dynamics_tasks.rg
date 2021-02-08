import "regent"
require "data_structures"

local constants = require("constants")
local format = require("std/format")

--------------------------------------------------------------
-------------NOTES ON PORTING OVER FROM MPAS FORTRAN----------
--------------------------------------------------------------

-- Fortran "do" loops are Regent "for" loops. Fortran loops use inclusive stop values, while Regent uses exclusive
-- Arrays are 1-indexed in Fortran, not 0-indexed as in Regent
-- MPAS helper functions such as Edge.verticesOnEdge(1) use 1 indexed values, we need 0-indexed values for our array fields
------ example: MPAS calls cellsOnVertex(j, iVtx), we call vertex_region[iVtx].cellsOnVertex[j-1]
-- Cell, Vertex and Edge ID's are 1-indexed in MPAS, however our regions are 0-indexed, may need to subtract 1 on access
-- It can be helpful to directly copy over lines from MPAS and go through line by line rewriting in Regent

--------------------------------------------------------------

local nCells = constants.nCells
local nEdges = constants.nEdges
local nVertices = constants.nVertices
local maxEdges = constants.maxEdges
local maxEdges2 = constants.maxEdges2
local TWO = constants.TWO
local FIFTEEN = constants.FIFTEEN
local vertexDegree = constants.vertexDegree
local nVertLevels = constants.nVertLevels
local nRelaxZone = constants.nRelaxZone
local seconds_per_day = constants.seconds_per_day


local cio = terralib.includec("stdio.h")

acos = regentlib.acos(double)
asin = regentlib.asin(double)
cos = regentlib.cos(double)
sin = regentlib.sin(double)

copysign = regentlib.copysign(double)
fabs = regentlib.fabs(double)
pow = regentlib.pow(double)
sqrt = regentlib.sqrt(double)

--__demand(__cuda)
task atm_compute_signs(cr : region(ispace(int2d), cell_fs),
                       er : region(ispace(int2d), edge_fs),
                       vr : region(ispace(int2d), vertex_fs))
where
  reads (cr.{edgesOnCell, nEdgesOnCell, verticesOnCell}, er.{cellsOnEdge, verticesOnEdge, zb, zb3}, vr.{cellsOnVertex, edgesOnVertex}),
  writes (cr.{edgesOnCellSign, kiteForCell, zb_cell, zb3_cell}, vr.edgesOnVertexSign)
do

  format.println("Calling atm_compute_signs...")

  var vertex_range = rect2d { int2d {0, 0}, int2d {nVertices - 1, 0} }
  var cell_range_2d = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels} } -- Loop goes to nVertLevels + 1
  var cell_range_1d = rect2d { int2d {0, 0}, int2d {nCells - 1, 0} }

  for iVtx in vertex_range do
    for i = 0, vertexDegree do
      if (vr[iVtx].edgesOnVertex[i] <= nEdges) then
        if (iVtx.x == er[{vr[iVtx].edgesOnVertex[i], 0}].verticesOnEdge[1]) then
          vr[iVtx].edgesOnVertexSign[i] = 1.0
        else
          vr[iVtx].edgesOnVertexSign[i] = -1.0
        end
      else
        vr[iVtx].edgesOnVertexSign[i] = 0.0
      end
    end
  end

  for iCell in cell_range_1d do
    for i = 0, cr[iCell].nEdgesOnCell do
      if (cr[iCell].edgesOnCell[i] <= nEdges) then
        if (iCell.x == er[{cr[iCell].edgesOnCell[i], 0}].cellsOnEdge[0]) then
          cr[iCell].edgesOnCellSign[i] = 1.0
        else
          cr[iCell].edgesOnCellSign[i] = -1.0
        end
      else
        cr[iCell].edgesOnCellSign[i] = 0.0
      end
    end
  end

  for iCell in cell_range_2d do
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      if (cr[{iCell.x, 0}].edgesOnCell[i] <= nEdges) then
        if (iCell.x == er[{cr[{iCell.x, 0}].edgesOnCell[i], 0}].cellsOnEdge[0]) then
          --VERTICAL STRUCTURE ACCESSED--
          cr[iCell].zb_cell[i] = er[{cr[{iCell.x, 0}].edgesOnCell[i], iCell.y}].zb[0]
          cr[iCell].zb3_cell[i] = er[{cr[{iCell.x, 0}].edgesOnCell[i], iCell.y}].zb3[0]

          --cio.printf("zb at cell %d and level %d, index %d, is %f \n", iCell.x, iCell.y, i, cr[iCell].zb_cell[i])
          --cio.printf("zb3 at cell %d and level %d, index %d, is %f \n", iCell.x, iCell.y, i, cr[iCell].zb3_cell[i])

        else
          --VERTICAL--

          cr[iCell].zb_cell[i] = er[{cr[{iCell.x, 0}].edgesOnCell[i], iCell.y}].zb[1]
          cr[iCell].zb3_cell[i] = er[{cr[{iCell.x, 0}].edgesOnCell[i], iCell.y}].zb3[1]
          --cio.printf("zb at cell %d and level %d, index %d, is %f \n", iCell.x, iCell.y, i, cr[iCell].zb_cell[i])
          --cio.printf("zb3 at cell %d and level %d, index %d, is %f \n", iCell.x, iCell.y, i, cr[iCell].zb3_cell[i])

        end
      end
    end
  end


  for iCell in cell_range_1d do
    for i = 0, cr[iCell].nEdgesOnCell do
      var iVtx = cr[iCell].verticesOnCell[i]
      if (iVtx <= nVertices) then
        for j = 1, vertexDegree do
          if (iCell.x == vr[{iVtx, 0}].cellsOnVertex[j]) then
            cr[iCell].kiteForCell[i] = j
            break
          end
        end
      -- trimmed a log statement here
      else
        cr[iCell].kiteForCell[i] = 1
      end
      --cio.printf("cr[{%d, 0}].kiteForCell[%d] is %f \n", iCell, i, cr[{iCell, 0}].kiteForCell[i])
    end
  end
end

--_demand(__cuda)
task atm_adv_coef_compression(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs))
where
  reads (cr.{cellsOnCell, nEdgesOnCell}, er.{cellsOnEdge, dcEdge, deriv_two, dvEdge}),
  writes (er.advCellsForEdge),
  reads writes (er.{adv_coefs, adv_coefs_3rd, nAdvCellsForEdge})
do

  format.println("Calling atm_adv_coef_compression...")

  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, 0} }
  var cell_list : int[maxEdges]

  for iEdge in edge_range do
    er[iEdge].nAdvCellsForEdge = 0
    var cell1 = er[iEdge].cellsOnEdge[0]
    var cell2 = er[iEdge].cellsOnEdge[1]

      --
      -- do only if this edge flux is needed to update owned cells
      --
    if (cell1 <= nCells or cell2 <= nCells) then
      cell_list[0] = cell1
      cell_list[1] = cell2
      var n = 1 -- n is index of cells currently in list

      --  add cells surrounding cell 1.  n is number of cells currently in list
      for i = 0, cr[{cell1, 0}].nEdgesOnCell do
        if (cr[{cell1, 0}].cellsOnCell[i] ~= cell2) then
          n += 1
          cell_list[n] = cr[{cell1, 0}].cellsOnCell[i]
        end
      end

      --  add cells surrounding cell 2 (brute force approach)
      for iCell = 0, cr[{cell2, 0}].nEdgesOnCell do
        var addcell = true
        for i = 0, n do
          if (cell_list[i] == cr[{cell2, 0}].cellsOnCell[iCell]) then
            addcell = false
          end
        end
        if (addcell and n < maxEdges - 1) then -- Temporary solution to avoid segfault: n can go out of bounds
          n += 1
          cell_list[n] = cr[{cell2, 0}].cellsOnCell[iCell]
        end
      end

      er[iEdge].nAdvCellsForEdge = n
      for iCell = 0, er[iEdge].nAdvCellsForEdge do
        er[iEdge].advCellsForEdge[iCell] = cell_list[iCell]
      end

      -- we have the ordered list, now construct coefficients
      for coef = 0, FIFTEEN do
        er[iEdge].adv_coefs[coef] = 0.0
        er[iEdge].adv_coefs_3rd[coef] = 0.0
      end -- initialize list to 0

      -- pull together third and fourth order contributions to the flux
      -- first from cell1

      var j_in = 0
      for j = 0, n do
        if (cell_list[j] == cell1) then
          j_in = j
        end
      end
      er[iEdge].adv_coefs[j_in] += er[iEdge].deriv_two[0]
      er[iEdge].adv_coefs_3rd[j_in] += er[iEdge].deriv_two[0]

      for iCell = 0, cr[{cell1, 0}].nEdgesOnCell do
        j_in = 0
        for j = 0, n do
          if (cell_list[j] == cr[{cell1, 0}].cellsOnCell[iCell]) then
            j_in = j
          end
        end
        er[iEdge].adv_coefs[j_in] += er[iEdge].deriv_two[iCell * FIFTEEN + 0]
        er[iEdge].adv_coefs_3rd[j_in] += er[iEdge].deriv_two[iCell * FIFTEEN + 0]
      end

      -- pull together third and fourth order contributions to the flux
      -- now from cell2

      j_in = 0
      for j = 0, n do
        if (cell_list[j] == cell2) then
          j_in = j
        end
      end
      er[iEdge].adv_coefs[j_in] += er[iEdge].deriv_two[1]
      er[iEdge].adv_coefs_3rd[j_in] += er[iEdge].deriv_two[1]

      for iCell = 0, cr[{cell2, 0}].nEdgesOnCell do
        j_in = 0
        for j = 0, n do
          if (cell_list[j] == cr[{cell2, 0}].cellsOnCell[iCell]) then
            j_in = j
          end
        end
        er[iEdge].adv_coefs[j_in] += er[iEdge].deriv_two[iCell * FIFTEEN + 1]
        er[iEdge].adv_coefs_3rd[j_in] += er[iEdge].deriv_two[iCell * FIFTEEN + 1]
      end

      for j = 0, n do
        er[iEdge].adv_coefs[j] =  -1.0 * pow(er[iEdge].dcEdge, 2) * er[iEdge].adv_coefs[j] / 12 -- this should be a negative number
        er[iEdge].adv_coefs_3rd[j] =  -1.0 * pow(er[iEdge].dcEdge, 2) * er[iEdge].adv_coefs_3rd[j] / 12
      end

      -- 2nd order centered contribution - place this in the main flux weights

      j_in = 0
      for j = 0, n do
        if (cell_list[j] == cell1) then
          j_in = j
        end
      end
      er[iEdge].adv_coefs[j_in] += 0.5

      j_in = 0
      for j = 0, n do
        if (cell_list[j] == cell2) then
          j_in = j
        end
      end
      er[iEdge].adv_coefs[j_in] += 0.5

      --  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply

      for j = 0, n do
        er[iEdge].adv_coefs[j] *= er[iEdge].dvEdge
        er[iEdge].adv_coefs_3rd[j] *= er[iEdge].dvEdge
      end
    end  -- only do for edges of owned-cells
  end -- end loop over edges
end


--config_zd: default 22000.0, config_xnutr: default 0.2. From config
__demand(__cuda)
task atm_compute_damping_coefs(config_zd : double,
                               config_xnutr : double,
                               cr : region(ispace(int2d), cell_fs))
where
  reads (cr.{meshDensity, zgrid}),
  reads writes (cr.dss)
do

  format.println("Calling atm_compute_damping_coefs...")

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }

  var m1 = -1.0
  var pii = acos(m1)

  var dx_scale_power = 1.0

  for iCell in cell_range do
    cr[iCell].dss = 0.0
    var zt = cr[{iCell.x, nVertLevels}].zgrid
    var z = 0.5 * (cr[iCell].zgrid + cr[iCell + {0, 1}].zgrid)
    if (z > config_zd) then
      cr[iCell].dss = config_xnutr * pow(sin(0.5 * pii * (z-config_zd)/(zt-config_zd)), 2.0)
      cr[iCell].dss /= pow(cr[{iCell.x, 0}].meshDensity, (0.25*dx_scale_power))
    end
  end
end

__demand(__cuda)
task atm_couple_coef_3rd_order(config_coef_3rd_order : double,
                               cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs))
where
  reads writes (cr.zb3_cell, er.adv_coefs_3rd)
do

  format.println("Calling atm_couple_coef_3rd_order...")

  var edge_range = rect1d { 0, nEdges - 1 }
  for iEdge in edge_range do
    for i = 0, FIFTEEN do
      er[{iEdge, 0}].adv_coefs_3rd[i] *= config_coef_3rd_order
    end
  end

  var cell_range = rect1d { 0, nCells - 1 }
  for iCell in cell_range do
    for j = 0, maxEdges do
      cr[{iCell, 0}].zb3_cell[j] *= config_coef_3rd_order
    end
  end
end

__demand(__cuda)
task atm_compute_solve_diagnostics(cr : region(ispace(int2d), cell_fs),
                                   er : region(ispace(int2d), edge_fs),
                                   vr : region(ispace(int2d), vertex_fs),
                                   hollingsworth : bool,
                                   rk_step : int)
where
  reads (cr.{edgesOnCell, edgesOnCellSign, h, invAreaCell, kiteForCell, nEdgesOnCell, verticesOnCell}, er.{cellsOnEdge, dcEdge, dvEdge, edgesOnEdge_ECP, nEdgesOnEdge, u, verticesOnEdge, weightsOnEdge}, vr.{edgesOnVertex, edgesOnVertexSign, fVertex, invAreaTriangle, kiteAreasOnVertex}),
  writes (er.{h_edge, pv_edge}),
  reads writes (cr.{divergence, ke}, er.{ke_edge, v}, vr.{ke_vertex, pv_vertex, vorticity})
do

  format.println("Calling atm_compute_solve_diagnostics...")

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }
  var vertex_range = rect2d { int2d {0, 0}, int2d {nVertices - 1, nVertLevels - 1} }

  -- Compute height on cell edges at velocity locations
  for iEdge in edge_range do
    var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]

    er[iEdge].h_edge = 0.5 * (cr[{cell1, iEdge.y}].h + cr[{cell2, iEdge.y}].h)
    var efac = er[{iEdge.x, 0}].dcEdge * er[{iEdge.x, 0}].dvEdge
    er[iEdge].ke_edge = efac * pow(er[iEdge].u, 2)
  end

  -- Compute circulation and relative vorticity at each vertex
  for iVertex in vertex_range do
    vr[iVertex].vorticity = 0.0
    for i = 0, vertexDegree do
      var iEdge = vr[{iVertex.x, 0}].edgesOnVertex[i]
      var s = vr[{iVertex.x, 0}].edgesOnVertexSign[i] * er[{iEdge, 0}].dcEdge

      vr[iVertex].vorticity += s * er[{iEdge, iVertex.y}].u
    end

    vr[iVertex].vorticity *= vr[{iVertex.x, 0}].invAreaTriangle
  end

  -- Compute the divergence at each cell center
  for iCell in cell_range do
    cr[iCell].divergence = 0.0
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
      var s = cr[{iCell.x, 0}].edgesOnCellSign[i] * er[{iEdge, 0}].dvEdge

      cr[iCell].divergence += s + er[{iEdge, iCell.y}].u
    end
    var r = cr[{iCell.x, 0}].invAreaCell
    cr[iCell].divergence *= r
  end

  -- Compute kinetic energy in each cell (Ringler et al JCP 2009)
  for iCell in cell_range do
    cr[iCell].ke = 0.0
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
      cr[iCell].ke += 0.25 * er[{iEdge, iCell.y}].ke_edge
    end

    cr[iCell].ke *= cr[{iCell.x, 0}].invAreaCell
  end

  if (hollingsworth) then
    -- Compute ke at cell vertices - AG's new KE construction, part 1
    -- *** approximation here because we don't have inner triangle areas
    for iVertex in vertex_range do
      var r = 0.25 * vr[{iVertex.x, 0}].invAreaTriangle
      vr[iVertex].ke_vertex = (er[{vr[{iVertex.x, 0}].edgesOnVertex[0], iVertex.y}].ke_edge
                              + er[{vr[{iVertex.x, 0}].edgesOnVertex[1], iVertex.y}].ke_edge
                              + er[{vr[{iVertex.x, 0}].edgesOnVertex[2], iVertex.y}].ke_edge) * r
    end

    -- adjust ke at cell vertices - AG's new KE construction, part 2
    var ke_fact = 1.0 - 0.375

    for iCell in cell_range do
      cr[iCell].ke *= ke_fact
    end

    for iCell in cell_range do
      var r = cr[{iCell.x, 0}].invAreaCell
      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
        var iVertex = cr[{iCell.x, 0}].verticesOnCell[i]
        var j = cr[{iCell.x, 0}].kiteForCell[i]

        cr[iCell].ke += (1.0 - ke_fact)*vr[{iVertex, 0}].kiteAreasOnVertex[j] * vr[{iVertex, iCell.y}].ke_vertex * r
      end
    end
  end

  -- Compute v (tangential) velocities following Thuburn et al JCP 2009
  -- The tangential velocity is only used to compute the Smagorinsky coefficient
  var reconstruct_v = true

  -- The original code checks that rk_step is present and not equal to 3.
  -- -1 is used as a sentinel value, and we are 0-indexing rk_step.
  if (rk_step ~= -1 and rk_step ~= 2) then
    reconstruct_v = false
  end

  if (reconstruct_v) then
    for iEdge in edge_range do
      er[iEdge].v = 0
      for i = 1, er[{iEdge.x, 0}].nEdgesOnEdge do
        var eoe = er[{iEdge.x, 0}].edgesOnEdge_ECP[i]

          er[iEdge].v += er[{iEdge.x, 0}].weightsOnEdge[i] * er[{eoe, iEdge.y}].u
      end
    end
  end

  -- Compute height at vertices, pv at vertices, and average pv to edge locations
  -- ( this computes pv_vertex at all vertices bounding real cells )
  for iVertex in vertex_range do
    vr[iVertex].pv_vertex = vr[{iVertex.x, 0}].fVertex + vr[iVertex].vorticity
  end

  -- Compute pv at the edges
  --  ( this computes pv_edge at all edges bounding real cells )
  for iEdge in edge_range do
    er[iEdge].pv_edge =  0.5 * (vr[{er[{iEdge.x, 0}].verticesOnEdge[0], iEdge.y}].pv_vertex + vr[{er[{iEdge.x, 0}].verticesOnEdge[1], iEdge.y}].pv_vertex)
  end

    --SKIPPED: (config_apvm_upwinding > 0.0) then---
end

-- Comments:
-- This function contains nCellsSolve, moist_start, moist_end, and scalars,
-- which we are currently not sure how to translate
__demand(__cuda)
task atm_compute_moist_coefficients(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs))
where
  reads (er.cellsOnEdge),
  writes (cr.cqw, er.cqu),
  reads writes (cr.qtot)
do

  format.println("Calling atm_compute_moist_coefficients...")

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }

  for iCell in cell_range do
    cr[iCell].qtot = 0.0
    for k = 0, nVertLevels do
      --TODO: What should we use instead of moist_start/moist_end?
      --for iq = moist_start, moist_end do
        --TODO: not sure how to translate: scalars(iq, k, iCell)
        --cr[iCell].qtot += cr[iCell].scalars[iq]
      --end
    end
  end

  for iCell in cell_range do
    if (iCell.y > 0) then -- Original from Fortran: do k = 2, nVertLevels
      var qtotal = 0.5 * (cr[iCell].qtot + cr[iCell - {0, 1}].qtot)
      cr[iCell].cqw = 1.0 / (1.0 + qtotal)
    end
  end

  for iEdge in edge_range do
    var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
    --if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then
      --var qtotal = 0.0
      --for iq = moist_start, moist_end do
        --qtotal += 0.5 * ( cr[cell1, iEdge.y].scalars[iq] + cr[cell2, iEdge.y].scalars[iq] )
      --end
      --er[iEdge].cqu = 1.0 / (1.0 + qtotal)
    --end
  end
end

--Comments for atm_compute_vert_imp_coefs
--Combined this and atm_compute_vert_imp_coefs_work (one was just a helper)
--purpose: "Compute coefficients for vertically implicit gravity-wave/acoustic computations"
--dts is passed in as double now - in code, passed in as rk_sub_timestep(rk_step)
--1.0_RKIND translated as just 1.0
--mpas used cellSolveStart, cellSolveEnd; we believe those deal with how mpas uses parallelization (a start and end for which cell pool is currently being worked on). We will just loop thru the entire cell region instead.
--mpas renames some fields; we will just use the original name --- p : exner, t : theta_m, rb : rho_base, rtb : rtheta_base, pb : exner_base, rt : rtheta_p

__demand(__cuda)
task atm_compute_vert_imp_coefs(cr : region(ispace(int2d), cell_fs),
                                vert_r : region(ispace(int1d), vertical_fs),
                                dts : double)
where
  reads (cr.{cqw, exner, exner_base, qtot, rho_base, rtheta_base, rtheta_p, theta_m, zz},
         vert_r.{rdzu, rdzw, fzm, fzp}),
  reads writes (cr.{a_tri, alpha_tri, b_tri, c_tri, coftz, cofwr, cofwt, cofwz, gamma_tri}, vert_r.cofrz)
do

  format.println("Calling atm_compute_vert_imp_coefs...")

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var cell_range_1d = rect2d { int2d {0, 0}, int2d {nCells - 1, 0} }
  var vert_level_range = rect1d { 0, nVertLevels - 1 }

  --  set coefficients
  var dtseps = .5 * dts * (1.0 + constants.config_epssm)
  var rgas = constants.rgas
  var rcv = rgas / (constants.cp - rgas)
  var c2 = constants.cp * rcv

  var qtotal : double

  -- MGD bad to have all threads setting this variable?
  for k in vert_level_range do
    vert_r[k].cofrz = dtseps * vert_r[k].rdzw
  end

  for iCell in cell_range_1d do
    --cr[iCell].a_tri = 0.0 --a_tri(1,iCell) = 0.  -- note, this value is never used
    --cr[{iCell.x, 0}].b_tri = 1.0    -- note, this value is never used
    --cr[{iCell.x, 0}].c_tri = 0.0    -- note, this value is never used
    cr[iCell].gamma_tri = 0.0
    --cr[iCell].alpha_tri = 0.0  -- note, this value is never used
  end

  -- NB: The following loops must be kept separate, otherwise the parallelism will fail.
  for iCell in cell_range do --  we only need to do cells we are solving for, not halo cells
    if (iCell.y > 0) then
      cr[iCell].cofwr = .5 * dtseps * constants.gravity * (vert_r[iCell.y].fzm * cr[iCell].zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz)
    end

    cr[iCell].coftz = 0.0
    if (iCell.y > 0) then
      cr[iCell].cofwz = dtseps * c2 * (vert_r[iCell.y].fzm * cr[iCell].zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz) * vert_r[iCell.y].rdzu * cr[iCell].cqw * (vert_r[iCell.y].fzm * cr[iCell].exner + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].exner)
      cr[iCell].coftz = dtseps * (vert_r[iCell.y].fzm * cr[iCell].theta_m + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].theta_m)
    end

    var qtotal = cr[iCell].qtot

    cr[iCell].cofwt = .5 * dtseps * rcv * cr[iCell].zz * constants.gravity * cr[iCell].rho_base / ( 1.0 + qtotal) * cr[iCell].exner / ((cr[iCell].rtheta_base + cr[iCell].rtheta_p) * cr[iCell].exner_base)
  end

  for iCell in cell_range do
    if (iCell.y > 0) then --k=2,nVertLevels
      cr[iCell].a_tri = -1.0 * cr[iCell].cofwz * cr[iCell - {0, 1}].coftz * vert_r[iCell.y - 1].rdzw * cr[iCell - {0, 1}].zz
                        + cr[iCell].cofwr * vert_r[iCell.y - 1].cofrz - cr[iCell - {0, 1}].cofwt * cr[iCell - {0, 1}].coftz * vert_r[iCell.y - 1].rdzw

      cr[iCell].b_tri = 1.0 + cr[iCell].cofwz * (cr[iCell].coftz * vert_r[iCell.y].rdzw * cr[iCell].zz + cr[iCell].coftz
                       * vert_r[iCell.y - 1].rdzw * cr[iCell - {0, 1}].zz) - cr[iCell].coftz * (cr[iCell].cofwt * vert_r[iCell.y].rdzw
                       - cr[iCell].cofwt * vert_r[iCell.y - 1].rdzw) + cr[iCell].cofwr * ((vert_r[iCell.y].cofrz - vert_r[iCell.y - 1].cofrz))

      cr[iCell].c_tri = -1.0 * cr[iCell].cofwz * cr[iCell + {0, 1}].coftz * vert_r[iCell.y].rdzw * cr[iCell].zz
                       - cr[iCell].cofwr * vert_r[iCell.y].cofrz + cr[iCell].cofwt * cr[iCell + {0, 1}].coftz * vert_r[iCell.y].rdzw
    end
  end

  for iCell in cell_range do
    if (iCell.y > 0) then --k=2,nVertLevels
      --MGD VECTOR DEPENDENCE
      cr[iCell].alpha_tri = 1.0 / (cr[iCell].b_tri - cr[iCell].a_tri * cr[iCell - {0, 1}].gamma_tri)
    end
  end

  for iCell in cell_range do
    if (iCell.y > 0) then --k=2,nVertLevels
      cr[iCell].gamma_tri = cr[iCell].c_tri * cr[iCell].alpha_tri
    end
  end -- loop over cells
end

__demand(__cuda)
task atm_compute_mesh_scaling(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs),
                              config_h_ScaleWithMesh : bool)
where
  reads (cr.meshDensity, er.cellsOnEdge),
  writes (cr.meshScalingRegionalCell, er.{meshScalingDel2, meshScalingDel4, meshScalingRegionalEdge})
do

  format.println("Calling atm_compute_mesh_scaling...")

  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, 0} }
  -- Compute the scaling factors to be used in the del2 and del4 dissipation
  for iEdge in edge_range do
    er[iEdge].meshScalingDel2 = 1.0
    er[iEdge].meshScalingDel4 = 1.0
  end

  if (config_h_ScaleWithMesh) then
    for iEdge in edge_range do
      var cell1 = er[iEdge].cellsOnEdge[0]
      var cell2 = er[iEdge].cellsOnEdge[1]
      er[iEdge].meshScalingDel2 = 1.0 / pow((cr[{cell1, 0}].meshDensity + cr[{cell2, 0}].meshDensity)/2.0, 0.25)
      er[iEdge].meshScalingDel4 = 1.0 / pow((cr[{cell1, 0}].meshDensity + cr[{cell2, 0}].meshDensity)/2.0, 0.75)
    end
  end

  -- Compute the scaling factors to be used in relaxation zone of regional configuration

  var cell_range = rect2d { int2d{0, 0}, int2d{nCells - 1, 0} }
  for iCell in cell_range do
    cr[iCell].meshScalingRegionalCell = 1.0
  end

  for iEdge in edge_range do
    er[iEdge].meshScalingRegionalEdge = 1.0
  end

  if (config_h_ScaleWithMesh) then
    for iEdge in edge_range do
      var cell1 = er[iEdge].cellsOnEdge[0]
      var cell2 = er[iEdge].cellsOnEdge[1]
      er[iEdge].meshScalingRegionalEdge = 1.0 / pow((cr[{cell1, 0}].meshDensity + cr[{cell2, 0}].meshDensity)/2.0, 0.25)
    end

    for iCell in cell_range do
      cr[iCell].meshScalingRegionalCell = 1.0 / pow(cr[iCell].meshDensity, 0.25)
    end
  end
end

--Not sure how to translate: scalars(index_qv,k,iCell)
--sign(1.0_RKIND,flux) translated as copysign(1.0, flux)
__demand(__cuda)
task atm_init_coupled_diagnostics(cr : region(ispace(int2d), cell_fs),
                                  er : region(ispace(int2d), edge_fs),
                                  vert_r : region(ispace(int1d), vertical_fs))
where
  reads (cr.{edgesOnCell, edgesOnCellSign, nEdgesOnCell, rho_base, theta, theta_base, theta_m, w,
             zb_cell, zb3_cell, zz},
         er.{cellsOnEdge, u},
         vert_r.{fzm, fzp}),
  writes (cr.{pressure_base, pressure_p}),
  reads writes (cr.{exner, exner_base, rho_p, rho_zz, rtheta_base, rtheta_p, rw, theta_m},
                er.ru)
do

  format.println("Calling atm_init_coupled_diagnostics...")
  var rgas = constants.rgas
  var rcv = rgas / (constants.cp - rgas)
  var p0 = 100000

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }

  for iCell in cell_range do
    --cr[iCell].theta_m = cr[iCell].theta * (1.0 + constants.rvord * cr[iCell].scalars[index_qv])
    cr[iCell].rho_zz /= cr[iCell].zz
  end

  for iEdge in edge_range do
    var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
    er[iEdge].ru = 0.5 * er[iEdge].u * (cr[{cell1, iEdge.y}].rho_zz + cr[{cell2, iEdge.y}].rho_zz)
  end

  -- Compute rw (i.e. rho_zz * omega) from rho_zz, w, and ru.
  -- We are reversing the procedure we use in subroutine atm_recover_large_step_variables.
  -- first, the piece that depends on w.
  for iCell in cell_range do
    cr[iCell].rw = 0
    if (iCell.y > 0 and iCell.y < nVertLevels) then -- Original Fortran: do k = 2, nVertLevels
      cr[iCell].rw = cr[iCell].w * (vert_r[iCell.y].fzp * cr[iCell - {0, 1}].rho_zz + vert_r[iCell.y].fzm * cr[iCell].rho_zz)
                     * (vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz + vert_r[iCell.y].fzm * cr[iCell].zz)
    end
  end

  -- next, the piece that depends on ru
  for iCell in cell_range do
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
      if (iCell.y > 0) then -- Original Fortran: do k = 2, nVertLevels
        var flux = vert_r[iCell.y].fzm * er[{iEdge, iCell.y}].ru + vert_r[iCell.y].fzp * er[{iEdge, iCell.y - 1}].ru
        cr[iCell].rw -= cr[{iCell.x, 0}].edgesOnCellSign[i]
                        * (cr[iCell].zb_cell[i] + copysign(1.0, flux) * cr[iCell].zb3_cell[i]) * flux
                        * (vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz + vert_r[iCell.y].fzm * cr[iCell].zz)
      end
    end
  end

  for iCell in cell_range do
    cr[iCell].rho_p = cr[iCell].rho_zz - cr[iCell].rho_base
    cr[iCell].rtheta_base = cr[iCell].theta_base * cr[iCell].rho_base
    cr[iCell].rtheta_p = cr[iCell].theta_m * cr[iCell].rho_p + cr[iCell].rho_base
                         * (cr[iCell].theta_m - cr[iCell].theta_base)
    cr[iCell].exner = pow(cr[iCell].zz * (rgas/p0) * (cr[iCell].rtheta_p + cr[iCell].rtheta_base), rcv)
    cr[iCell].exner_base = pow(cr[iCell].zz * (rgas/p0) * (cr[iCell].rtheta_base), rcv)
    cr[iCell].pressure_p = cr[iCell].zz * rgas * (cr[iCell].exner * cr[iCell].rtheta_p + cr[iCell].rtheta_base
                                                  * (cr[iCell].exner - cr[iCell].exner_base))
    cr[iCell].pressure_base = cr[iCell].zz * rgas * cr[iCell].exner_base * cr[iCell].rtheta_base

    --format.println("zz at cell ({}, {}) is {} \n", iCell.x, iCell.y, cr[iCell].zz)
    --format.println("exner at cell ({}, {}) is {} \n", iCell.x, iCell.y, cr[iCell].exner)
    --format.println("rtheta_p at cell ({}, {}) is {} \n", iCell.x, iCell.y, cr[iCell].rtheta_p)
    --format.println("rtheta_base at cell ({}, {}) is {} \n", iCell.x, iCell.y, cr[iCell].rtheta_base)
    --format.println("exner_base at cell ({}, {}) is {} \n", iCell.x, iCell.y, cr[iCell].exner_base)
    --format.println("Pressure_p at cell ({}, {}) is {} \n", iCell.x, iCell.y, cr[iCell].pressure_p)
  end
end

--Not sure how to translate: scalars(index_qv,k,iCell)
__demand(__cuda)
task atm_compute_output_diagnostics(cr : region(ispace(int2d), cell_fs))
where
  reads (cr.{pressure_base, pressure_p, rho_zz, theta_m, zz}),
  writes (cr.{pressure, rho, theta})
do
  format.println("Calling atm_compute_output_diagnostics...")

  var cell_range = rect2d { int2d{0, 0}, int2d{nCells - 1, nVertLevels - 1} }
  for iCell in cell_range do
    --Original contains scalars(index_qv,k,iCell). Currently translating as follows
    --cr[iCell].theta = cr[iCell].theta_m / (1 + constants.rvord * cr[iCell].scalars[index_qv])
    cr[iCell].rho = cr[iCell].rho_zz * cr[iCell].zz
    cr[iCell].pressure = cr[iCell].pressure_base + cr[iCell].pressure_p
  end

end

__demand(__cuda)
task atm_rk_integration_setup(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs))
where
  reads (cr.{rho_p, rho_zz, rtheta_p, rw, theta_m, w}, er.{ru, u}),
  writes (cr.{rho_p_save, rho_zz_2, rho_zz_old_split, rtheta_p_save, rw_save, theta_m_2, w_2}, er.{ru_save, u_2})
do

  format.println("Calling atm_rk_integration_setup...")
  var edge_range = rect2d { int2d{0, 0}, int2d{nEdges - 1, nVertLevels - 1} }
  var cell_range = rect2d { int2d{0, 0}, int2d{nCells - 1, nVertLevels - 1} }

  for iEdge in edge_range do
    er[iEdge].ru_save = er[iEdge].ru
    er[iEdge].u_2 = er[iEdge].u
  end

  for iCell in cell_range do
    cr[iCell].rw_save = cr[iCell].rw
    cr[iCell].rtheta_p_save = cr[iCell].rtheta_p
    cr[iCell].rho_p_save = cr[iCell].rho_p

    cr[iCell].w_2 = cr[iCell].w
    cr[iCell].theta_m_2 = cr[iCell].theta_m
    cr[iCell].rho_zz_2 = cr[iCell].rho_zz
    cr[iCell].rho_zz_old_split = cr[iCell].rho_zz
    --Not sure how to translate scalars
    --Original: scalars_2(:,:,cellStart:cellEnd) = scalars_1(:,:,cellStart:cellEnd)
    --for i = ???, ??? do -- Scalar bounds
      --cr[iCell].scalars_2[i] = cr[iCell].scalars[i]
    --end
  end
end

__demand(__inline)
task flux4(q_im2 : double, q_im1 : double, q_i : double, q_ip1 : double, ua : double) : double
  return ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0
end

__demand(__inline)
task flux3(q_im2 : double, q_im1 : double, q_i : double, q_ip1 : double, ua : double, coef3 : double) : double
  return flux4(q_im2, q_im1, q_i, q_ip1, ua)
         + coef3*fabs(ua)*((q_ip1-q_im2)-3.*(q_i-q_im1))/12.0
end

__demand(__inline)
task rayleigh_damp_coef(vertLevel : double) : double
  var rayleigh_coef_inverse = 1.0 / ( [double](constants.config_number_rayleigh_damp_u_levels)
                                      * (constants.config_rayleigh_damp_u_timescale_days * seconds_per_day) )
  return [double](vertLevel - (nVertLevels - constants.config_number_rayleigh_damp_u_levels)) * rayleigh_coef_inverse
end

-- Comments for atm_compute_dyn_tend_work
-- A large number of variables appeared to be missing from the registry. I have inferred their types
-- and regions to the best of my ability, but may have gotten some of them wrong.
-- They are as follows:
--    cr: tend_rho_physics, dpdz, rr_save, delsq_divergence, delsq_w,
--    tend_w_euler, theta_m_save, delsq_theta, tend_theta_euler, tend_rtheta_physics
--    er: tend_u_euler, delsq_u, tend_ru_physics
--    vr: delsq_vorticity
-- Variable c_s appears to be a renaming of constants.config_smagorinsky_coef.
-- Some pieces of code were inside an #ifdef CURVATURE. I have ignored the ifdefs and put the code in
-- unconditionally.
-- Some config values should be configurable, rather than constants. These are currently in constants.rg.
-- I have also used previous conventions like removing "_RKIND", looping over all cells instead of
-- cellSolveStart to cellStartEnd, using copysign and fabs for sign and abs, etc.

__demand(__cuda)
task atm_compute_dyn_tend_work(cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs),
                               vr : region(ispace(int2d), vertex_fs),
                               vert_r : region(ispace(int1d), vertical_fs),
                               rk_step : int,
                               dt : double,
                               config_horiz_mixing : regentlib.string,
                               config_mpas_cam_coef : double,
                               config_mix_full : bool,
                               config_rayleigh_damp_u : bool)
where
  reads (cr.{cqw, defc_a, defc_b, divergence, edgesOnCell, edgesOnCell_sign, invAreaCell, ke, lat,
             nEdgesOnCell, pressure_p, qtot, rho_base, rho_zz, rho_p_save, rt_diabatic_tend, rw, rw_save, t_init, tend_rho_physics,
             tend_rtheta_physics, theta_m, theta_m_save, uReconstructZonal, uReconstructMeridional, w, zgrid, zz},
         er.{advCellsForEdge, adv_coefs, adv_coefs_3rd, angleEdge, cellsOnEdge, cqu, dcEdge, dvEdge,
             edgesOnEdge, invDcEdge, invDvEdge, lat, meshScalingDel2, meshScalingDel4, nAdvCellsForEdge,
             nEdgesOnEdge, pv_edge, rho_edge, ru, ru_save, tend_ru_physics, u, v, verticesOnEdge,
             weightsOnEdge, zxu},
         vr.{edgesOnVertex, edgesOnVertex_sign, invAreaTriangle, vorticity},
         vert_r.{fzm, fzp, rdzu, rdzw, u_init, v_init}),
  writes (cr.{rthdynten, tend_rho, tend_rtheta_adv}),
  reads writes (cr.{delsq_divergence, delsq_theta, delsq_w, dpdz, flux_arr, h_divergence, kdiff, ru_edge_w,
                    tend_theta, tend_theta_euler, w, tend_w_euler, wdtz, wdwz},
                er.{delsq_u, q, tend_u, tend_u_euler, u_mix, wduz},
                vr.delsq_vorticity)
do

  format.println("Calling atm_compute_dyn_tend_work...")

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }
  var vertex_range = rect2d { int2d {0, 0}, int2d {nVertices - 1, nVertLevels - 1} }

  var prandtl_inv = 1.0 / constants.prandtl
  -- Can't find dt
  var invDt = 1.0 / dt
  var r_earth = constants.sphere_radius
  var inv_r_earth = 1.0 / r_earth

  var v_mom_eddy_visc2 = constants.config_v_mom_eddy_visc2 -- 0.0
  var v_theta_eddy_visc2 = constants.config_v_theta_eddy_visc2 -- 0.0
  var h_mom_eddy_visc4 = constants.config_h_mom_eddy_visc4 -- 0.0
  var h_theta_eddy_visc4 = constants.config_h_theta_eddy_visc4 --0.0

  if (rk_step == 0) then
    -- Smagorinsky eddy viscosity, based on horizontal deformation (in this case on model coordinate surfaces).
    -- The integration coefficients were precomputed and stored in defc_a and defc_b
    if ([rawstring](config_horiz_mixing) == "2d_smagorinsky") then
      var c_s = constants.config_smagorinsky_coef
      for iCell in cell_range do
        var d_diag : double[nVertLevels]
        var d_off_diag : double[nVertLevels]
        for k = 0, nVertLevels do
          d_diag[k] = 0.0
          d_off_diag[k] = 0.0
        end
        for iEdge = 0, cr[{iCell.x, 0}].nEdgesOnCell do
          for k = 0, nVertLevels do
            var e = cr[{iCell.x, 0}].edgesOnCell[iEdge]
            d_diag[k]     += cr[{iCell.x, 0}].defc_a[iEdge] * er[{e, k}].u
                              - cr[{iCell.x, 0}].defc_b[iEdge] * er[{e, k}].v
            d_off_diag[k] += cr[{iCell.x, 0}].defc_b[iEdge] * er[{e, k}].u
                              + cr[{iCell.x, 0}].defc_a[iEdge] * er[{e, k}].v
          end
        end

        -- here is the Smagorinsky formulation,
        -- followed by imposition of an upper bound on the eddy viscosity

        -- Original: kdiff(k,iCell) = min((c_s * config_len_disp)**2 * sqrt(d_diag(k)**2 + d_off_diag(k)**2),(0.01*config_len_disp**2) * invDt)
        cr[iCell].kdiff = min(pow(c_s * constants.config_len_disp, 2.0)
                              * sqrt(pow(d_diag[iCell.y], 2.0) + pow(d_off_diag[iCell.y], 2.0)),
                              (0.01 * pow(constants.config_len_disp, 2.0)) * invDt)
      end

      h_mom_eddy_visc4 = constants.config_visc4_2dsmag * pow(constants.config_len_disp, 3.0) --0.05 * 120000.0**3
      h_theta_eddy_visc4 = h_mom_eddy_visc4

    elseif ([rawstring](config_horiz_mixing) == "2d_fixed") then
      for iCell in cell_range do
        cr[iCell].kdiff = 0.0
      end
    end

    if (config_mpas_cam_coef > 0.0) then
      for iCell in cell_range do
        -- 2nd-order filter for top absorbing layer as in CAM-SE :  WCS 10 May 2017
        -- From MPAS-CAM V4.0 code, with addition to config-specified coefficient (V4.0_coef = 0.2; SE_coef = 1.0)

        --Original:
        --cr[{iCell.x, nVertLevels-2}].kdiff = max(cr[{iCell.x, nVertLevels-2}].kdiff,
        --                                                2.0833 * constants.config_len_disp * config_mpas_cam_coef)
        --cr[{iCell.x, nVertLevels-1}].kdiff = max(cr[{iCell.x, nVertLevels-1}].kdiff,
        --                                          2.0 * 2.0833 * constants.config_len_disp * config_mpas_cam_coef)
        --cr[{iCell.x, nVertLevels  }].kdiff = max(cr[{iCell.x, nVertLevels  }].kdiff,
        --                                          4.0 * 2.0833 * constants.config_len_disp * config_mpas_cam_coef)

        if (iCell.y >= nVertLevels - 2 and iCell.y <= nVertLevels) then
          var p = iCell.y - (nVertLevels - 2)
          cr[iCell].kdiff = max(cr[iCell].kdiff, pow(2, p) * 2.0833 * constants.config_len_disp * config_mpas_cam_coef)
        end
      end
    end
  end

  -- tendency for density.
  -- accumulate total water here for later use in w tendency calculation.

  -- accumulate horizontal mass-flux

  for iCell in cell_range do
    cr[iCell].h_divergence = 0.0
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
      var edge_sign = cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge
      cr[iCell].h_divergence += edge_sign * er[{iEdge, iCell.y}].ru
    end
  end

  -- compute horiontal mass-flux divergence, add vertical mass flux divergence to complete tend_rho

  for iCell in cell_range do
    var r = cr[{iCell.x, 0}].invAreaCell
    cr[iCell].h_divergence *= r
  end

  -- dp / dz and tend_rho
  -- only needed on first rk_step with pert variables defined a pert from time t
  if (rk_step == 0) then

    var rgas_cprcv = constants.rgas * constants.cp / constants.cv
    for iCell in cell_range do
      --Not sure what tend_rho_physics, dpdz, rb are.
      cr[iCell].tend_rho = -cr[iCell].h_divergence - vert_r[iCell.y].rdzw * (cr[iCell + {0, 1}].rw - cr[iCell].rw + cr[iCell].tend_rho_physics)
      --Original: dpdz(k,iCell) = -gravity*(rb(k,iCell)*(qtot(k,iCell)) + rr_save(k,iCell)*(1.+qtot(k,iCell)))
      cr[iCell].dpdz = -constants.gravity * (cr[iCell].rho_base * (cr[iCell].qtot) + cr[iCell].rho_p_save * (1.0 + cr[iCell].qtot))
    end
  end



-------- BEGIN U SECTION --------

  -- Compute u (normal) velocity tendency for each edge (cell face)
  for iEdge in edge_range do

    var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]

    -- horizontal pressure gradient
    if (rk_step == 0) then
      --Unable to find: tend_u_euler, pp, dpdz
      --It looks like tend_u_euler is from er and the others are from cr
      er[iEdge].tend_u_euler = -er[iEdge].cqu * ( (cr[{cell2, iEdge.y}].pressure_p - cr[{cell1, iEdge.y}].pressure_p) * er[{iEdge.x, 0}].invDcEdge
                                / (0.5 * (cr[{cell2, iEdge.y}].zz + cr[{cell1, iEdge.y}].zz))
                                - 0.5 * er[iEdge].zxu * (cr[{cell1, iEdge.y}].dpdz + cr[{cell2, iEdge.y}].dpdz) )
    end

    -- vertical transport of u
    er[iEdge].wduz = 0.0
    if (iEdge.y == 1 or iEdge.y == nVertLevels - 1) then
      er[iEdge].wduz = 0.5 * (cr[{cell1, iEdge.y}].rw + cr[{cell2, iEdge.y}].rw) * (vert_r[iEdge.y].fzm * er[iEdge].u + vert_r[iEdge.y].fzp * er[iEdge - {0, 1}].u)
    end
    if (iEdge.y > 1 and iEdge.y < nVertLevels - 1) then
      er[iEdge].wduz = flux3( er[iEdge - {0, 2}].u, er[iEdge - {0, 1}].u, er[iEdge].u, er[iEdge + {0, 1}].u,
                                      0.5 * (cr[{cell1, iEdge.y}].rw + cr[{cell2, iEdge.y}].rw), 1.0 )
    end
  end

  for iEdge in edge_range do

    var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
    er[iEdge].tend_u = -vert_r[iEdge.y].rdzw * (er[iEdge + {0, 1}].wduz - er[iEdge].wduz)

    -- Next, nonlinear Coriolis term (q) following Ringler et al JCP 2009

    er[iEdge].q = 0.0

    for j = 0, er[{iEdge.x, 0}].nEdgesOnEdge do
      var eoe = er[{iEdge.x, 0}].edgesOnEdge[j]
      for k = 0, nVertLevels do
        var workpv = 0.5 * (er[iEdge].pv_edge + er[{eoe, iEdge.y}].pv_edge)
        -- the original definition of pv_edge had a factor of 1/density.  We have removed that factor
        -- given that it was not integral to any conservation property of the system
        er[iEdge].q += er[{iEdge.x, 0}].weightsOnEdge[j] * er[{eoe, iEdge.y}].u * workpv
      end
    end

    -- horizontal ke gradient and vorticity terms in the vector invariant formulation
    -- of the horizontal momentum equation
    er[iEdge].tend_u += er[iEdge].rho_edge
                        * ( er[iEdge].q - (cr[{cell2, iEdge.y}].ke - cr[{cell1, iEdge.y}].ke) * er[{iEdge.x, 0}].invDcEdge )
                        - er[iEdge].u * 0.5 * (cr[{cell1, iEdge.y}].h_divergence + cr[{cell2, iEdge.y}].h_divergence)

    -- #ifdef CURVATURE
    -- curvature terms for the sphere
    er[iEdge].tend_u -= ( 2.0 * constants.omega * cos(er[{iEdge.x, 0}].angleEdge)
                        * cos(er[{iEdge.x, 0}].lat) * er[iEdge].rho_edge
                        * 0.25 * (cr[{cell1, iEdge.y}].w + cr[{cell1, iEdge.y + 1}].w
                        + cr[{cell2, iEdge.y}].w + cr[{cell2, iEdge.y + 1}].w) )
                        - ( er[iEdge].u * 0.25 * (cr[{cell1, iEdge.y}].w + cr[{cell1, iEdge.y + 1}].w
                        + cr[{cell2, iEdge.y}].w + cr[{cell2, iEdge.y + 1}].w) * er[iEdge].rho_edge
                        * inv_r_earth )
    -- #endif
  end -- loop over edges

  -- horizontal mixing for u
  -- mixing terms are integrated using forward-Euler, so this tendency is only computed in the
  -- first Runge-Kutta substep and saved for use in later RK substeps 2 and 3.

  if (rk_step == 0) then

    -- del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
    -- First, storage to hold the result from the first del^2 computation.

    for iEdge in edge_range do
      er[iEdge].delsq_u = 0.0
      var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
      var vertex1 = er[{iEdge.x, 0}].verticesOnEdge[0]
      var vertex2 = er[{iEdge.x, 0}].verticesOnEdge[1]
      var r_dc = er[{iEdge.x, 0}].invDcEdge
      var r_dv = min(er[{iEdge.x, 0}].invDvEdge, 4 * r_dc)

      -- Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
      --                    only valid for h_mom_eddy_visc4 == constant
      var u_diffusion = (cr[{cell2, iEdge.y}].divergence - cr[{cell1, iEdge.y}].divergence) * r_dc
                        - (vr[{vertex2, iEdge.y}].vorticity - vr[{vertex1, iEdge.y}].vorticity) * r_dv
      er[iEdge].delsq_u += u_diffusion
      var kdiffu = 0.5 * (cr[{cell1, iEdge.y}].kdiff + cr[{cell2, iEdge.y}].kdiff)
      -- include 2nd-orer diffusion here
      er[iEdge].tend_u_euler += er[iEdge].rho_edge * kdiffu * u_diffusion
                                * er[{iEdge.x, 0}].meshScalingDel2
    end

    if (h_mom_eddy_visc4 > 0.0) then  -- 4th order mixing is active

      for iVertex in vertex_range do
        vr[iVertex].delsq_vorticity = 0.0
        for i = 0, vertexDegree do
          var iEdge = vr[{iVertex.x, 0}].edgesOnVertex[i]
          var edge_sign = vr[{iVertex.x, 0}].invAreaTriangle * er[{iEdge, 0}].dcEdge * vr[{iVertex.x, 0}].edgesOnVertex_sign[i]

          vr[iVertex].delsq_vorticity += edge_sign * er[{iEdge, iVertex.y}].delsq_u
        end
      end

      for iCell in cell_range do
        cr[iCell].delsq_divergence = 0.0
        var r = cr[{iCell.x, 0}].invAreaCell
        for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
          var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
          var edge_sign = r * er[{iEdge, 0}].dvEdge * cr[{iCell.x, 0}].edgesOnCell_sign[i]
          cr[iCell].delsq_divergence += edge_sign * er[{iEdge, iCell.y}].delsq_u
        end
      end

      for iEdge in edge_range do
        var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
        var vertex1 = er[{iEdge.x, 0}].verticesOnEdge[0]
        var vertex2 = er[{iEdge.x, 0}].verticesOnEdge[1]

        var u_mix_scale = er[{iEdge.x, 0}].meshScalingDel4 * h_mom_eddy_visc4
        var r_dc = u_mix_scale * constants.config_del4u_div_factor * er[{iEdge.x, 0}].invDcEdge
        var r_dv = u_mix_scale * min(er[{iEdge.x, 0}].invDvEdge, 4 * er[{iEdge.x, 0}].invDcEdge)

        -- Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
        --                    only valid for h_mom_eddy_visc4 == constant
        -- Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor
        --    relative to the rotational part.  The stability constraint on the divergence component is much less
        --    stringent than the rotational part, and this flexibility may be useful.
        var u_diffusion = er[iEdge].rho_edge * ( (cr[{cell2, iEdge.y}].delsq_divergence - cr[{cell1, iEdge.y}].delsq_divergence) * r_dc
                          - (vr[{vertex2, iEdge.y}].delsq_vorticity - vr[{vertex1, iEdge.y}].delsq_vorticity) * r_dv )
        er[iEdge].tend_u_euler -= u_diffusion
      end
    end -- 4th order mixing is active

    -- vertical mixing for u - 2nd order filter in physical (z) space
    if (v_mom_eddy_visc2 > 0.0) then

      if (config_mix_full) then  -- mix full state

        for iEdge in edge_range do
          var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
          if (iEdge.y > 0 and iEdge.y < nVertLevels - 1) then
            var z1 = 0.5 * (cr[{cell1, iEdge.y - 1}].zgrid + cr[{cell2, iEdge.y - 1}].zgrid)
            var z2 = 0.5 * (cr[{cell1, iEdge.y    }].zgrid + cr[{cell2, iEdge.y    }].zgrid)
            var z3 = 0.5 * (cr[{cell1, iEdge.y + 1}].zgrid + cr[{cell2, iEdge.y + 1}].zgrid)
            var z4 = 0.5 * (cr[{cell1, iEdge.y + 2}].zgrid + cr[{cell2, iEdge.y + 2}].zgrid)

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)

            er[iEdge].tend_u_euler += er[iEdge].rho_edge * v_mom_eddy_visc2
                                      * ( (er[iEdge + {0, 1}].u - er[iEdge].u) / (zp-z0)
                                        - (er[iEdge].u - er[iEdge - {0, 1}].u) / (z0-zm) )
                                      / (0.5 * (zp - zm))
          end
        end

      else  -- idealized cases where we mix on the perturbation from the initial 1-D state

        for iEdge in edge_range do
          er[iEdge].u_mix = er[iEdge].u - vert_r[iEdge.y].u_init * cos(er[{iEdge.x, 0}].angleEdge)
                            - vert_r[iEdge.y].v_init * sin(er[{iEdge.x, 0}].angleEdge)
        end

        for iEdge in edge_range do
          var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]

          if (iEdge.y > 0 and iEdge.y < nVertLevels - 1) then
            var z1 = 0.5 * (cr[{cell1, iEdge.y - 1}].zgrid + cr[{cell2, iEdge.y - 1}].zgrid)
            var z2 = 0.5 * (cr[{cell1, iEdge.y    }].zgrid + cr[{cell2, iEdge.y    }].zgrid)
            var z3 = 0.5 * (cr[{cell1, iEdge.y + 1}].zgrid + cr[{cell2, iEdge.y + 1}].zgrid)
            var z4 = 0.5 * (cr[{cell1, iEdge.y + 2}].zgrid + cr[{cell2, iEdge.y + 2}].zgrid)

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)

            er[iEdge].tend_u_euler += er[iEdge].rho_edge * v_mom_eddy_visc2
                                      * ( (er[iEdge + {0, 1}].u_mix - er[iEdge].u_mix) / (zp-z0)
                                        - (er[iEdge].u_mix - er[iEdge - {0, 1}].u_mix) / (z0-zm) )
                                      / (0.5 * (zp - zm))
          end
        end
      end -- mix perturbation state
    end -- vertical mixing of horizontal momentum
  end -- (rk_step 1 test for computing mixing terms)

  --  add in mixing and physics tendency for u

  --  Rayleigh damping on u
  if (config_rayleigh_damp_u) then

    for iEdge in edge_range do
      if (iEdge.y > nVertLevels - constants.config_number_rayleigh_damp_u_levels + 1) then
        er[iEdge].tend_u -= er[iEdge].rho_edge * er[iEdge].u * rayleigh_damp_coef(iEdge.y)
      end
    end
  end

  for iEdge in edge_range do
    er[iEdge].tend_u += er[iEdge].tend_u_euler + er[iEdge].tend_ru_physics
  end



-------- BEGIN W SECTION ---------

  --  horizontal advection for w
  for iCell in cell_range do
    cr[iCell].w = 0.0
  end

  for iCell in cell_range do
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
      var edge_sign = cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge * 0.5

      if (iCell.y > 0) then
        cr[iCell].ru_edge_w = vert_r[iCell.y].fzm * er[{iEdge, iCell.y}].ru + vert_r[iCell.y].fzp * er[{iEdge, iCell.y - 1}].ru
      end

      for k = 0, nVertLevels do
        cr[iCell].flux_arr = 0.0
      end

      -- flux_arr stores the value of w at the cell edge used in the horizontal transport

      for j = 0, er[{iEdge, 0}].nAdvCellsForEdge do
        var iAdvCell = er[{iEdge, 0}].advCellsForEdge[j]
        if (iCell.y > 0) then
          var scalar_weight = er[{iEdge, 0}].adv_coefs[j] + copysign(1.0, cr[iCell].ru_edge_w) * er[{iEdge, 0}].adv_coefs_3rd[j]
          cr[iCell].flux_arr += scalar_weight * cr[{iAdvCell, iCell.y}].w
        end
      end
    end
  end

  for iCell in cell_range do
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      if (iCell.y > 0) then
        cr[iCell].w -= cr[{iCell.x, 0}].edgesOnCell_sign[i] * cr[iCell].ru_edge_w * cr[iCell].flux_arr
      end
    end
  end

  -- #ifdef CURVATURE
  for iCell in cell_range do
    if (iCell.y > 0) then
      cr[iCell].w += (cr[iCell].rho_zz * vert_r[iCell.y].fzm
                          + cr[iCell - {0, 1}].rho_zz * vert_r[iCell.y].fzp)
                          * ( pow(vert_r[iCell.y].fzm * cr[iCell].uReconstructZonal + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].uReconstructZonal, 2.0)
                            + pow(vert_r[iCell.y].fzm * cr[iCell].uReconstructMeridional + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].uReconstructMeridional, 2.0) ) / r_earth
                          + 2.0 * constants.omega * cos(cr[{iCell.x, 0}].lat)
                          * (vert_r[iCell.y].fzm * cr[iCell].uReconstructZonal + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].uReconstructZonal)
                          * (cr[iCell].rho_zz * vert_r[iCell.y].fzm + cr[iCell - {0, 1}].rho_zz * vert_r[iCell.y].fzp)
    end
  end
  -- #endif

  -- horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),
  -- but here we can also code in hyperdiffusion if we wish (2nd order at present)

  if (rk_step == 0) then

    -- del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).

    -- First, storage to hold the result from the first del^2 computation.
    --  we copied code from the theta mixing, hence the theta* names.

    for iCell in cell_range do
      cr[iCell].delsq_w = 0.0
      cr[iCell].tend_w_euler = 0.0
      var r_areaCell = cr[{iCell.x, 0}].invAreaCell
      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
          var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]

          var edge_sign = 0.5 * r_areaCell * cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge
                          * er[{iEdge, 0}].invDcEdge

          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

        if (iCell.y > 0) then
          var w_turb_flux = edge_sign * (er[{iEdge, iCell.y}].rho_edge + er[{iEdge, iCell.y - 1}].rho_edge)
                            * (cr[{cell2, iCell.y}].w - cr[{cell1, iCell.y}].w)
          cr[iCell].delsq_w += w_turb_flux
          w_turb_flux *= er[{iEdge, 0}].meshScalingDel2 * 0.25
                         * (cr[{cell1, iCell.y}].kdiff + cr[{cell2, iCell.y}].kdiff
                           + cr[{cell1, iCell.y - 1}].kdiff + cr[{cell2, iCell.y - 1}].kdiff)
          cr[iCell].tend_w_euler += w_turb_flux
        end
      end
    end

    if (h_mom_eddy_visc4 > 0.0) then  -- 4th order mixing is active

      for iCell in cell_range do
        var r_areaCell = h_mom_eddy_visc4 * cr[{iCell.x, 0}].invAreaCell
        for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
          var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

          var edge_sign = er[{iEdge, 0}].meshScalingDel4 * r_areaCell * er[{iEdge, 0}].dvEdge
                          * cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].invDcEdge

          if (iCell.y > 0) then
            cr[iCell].tend_w_euler -= edge_sign * (cr[{cell2, iCell.y}].delsq_w - cr[{cell1, iCell.y}].delsq_w)
          end
        end
      end
    end -- 4th order mixing is active
  end -- horizontal mixing for w computed in first rk_step

  --  vertical advection, pressure gradient and buoyancy for w
  for iCell in cell_range do

    cr[iCell].wdwz = 0.0
    if (iCell.y == 1 or iCell.y == nVertLevels - 1) then
      cr[iCell].wdwz = 0.25 * (cr[iCell].rw + cr[iCell - {0, 1}].rw)*(cr[iCell].w + cr[iCell - {0, 1}].w)
    end
    if (iCell.y > 1 and iCell.y < nVertLevels - 1) then
      cr[iCell].wdwz = flux3( cr[iCell - {0, 2}].w, cr[iCell - {0, 1}].w, cr[iCell].w, cr[iCell + {0, 1}].w,
                                     0.5 * (cr[iCell].rw + cr[iCell - {0, 1}].rw), 1.0 )
    end
  end

  for iCell in cell_range do
    -- Note: next we are also dividing through by the cell area after the horizontal flux divergence
    if (iCell.y > 0) then
      cr[iCell].w *= cr[{iCell.x, 0}].invAreaCell - vert_r[iCell.y].rdzu * (cr[iCell + {0, 1}].wdwz - cr[iCell].wdwz)
    end

    if (rk_step == 0) then
      if (iCell.y > 0) then
        cr[iCell].tend_w_euler -= cr[iCell].cqw * ( vert_r[iCell.y].rdzu *
                                  (cr[iCell].pressure_p - cr[iCell - {0, 1}].pressure_p)
                                  - (vert_r[iCell.y].fzm * cr[iCell].dpdz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].dpdz) )  -- dpdz is the buoyancy term here.
      end
    end
  end

  if (rk_step == 0) then
    if (v_mom_eddy_visc2 > 0.0) then
      for iCell in cell_range do
        if (iCell.y > 0) then
          cr[iCell].tend_w_euler += v_mom_eddy_visc2 * (cr[iCell].rho_zz + cr[iCell - {0, 1}].rho_zz)
                                         * 0.5 * ( (cr[iCell + {0, 1}].w - cr[iCell].w) * vert_r[iCell.y].rdzw
                                         - (cr[iCell].w - cr[iCell - {0, 1}].w) * vert_r[iCell.y - 1].rdzw )
                                         * vert_r[iCell.y].rdzu
        end
      end
    end
  end -- mixing term computed first rk_step

  -- add in mixing terms for w
  for iCell in cell_range do
    if (iCell.y > 0) then
      cr[iCell].w += cr[iCell].tend_w_euler
    end
  end



-------- BEGIN THETA SECTION --------
  -- horizontal advection for theta
  for iCell in cell_range do
    cr[iCell].tend_theta = 0.0
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]

      cr[iCell].flux_arr = 0.0

      for j = 0, er[{iEdge, 0}].nAdvCellsForEdge do
        var iAdvCell = er[{iEdge, 0}].advCellsForEdge[j]
        var scalar_weight = er[{iEdge, 0}].adv_coefs[j] + copysign(1.0, er[{iEdge, iCell.y}].ru)
                            * er[{iEdge, 0}].adv_coefs_3rd[j]
        cr[iCell].flux_arr += scalar_weight * cr[{iAdvCell, iCell.y}].theta_m
      end

      cr[iCell].tend_theta -= cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, iCell.y}].ru * cr[iCell].flux_arr
    end
  end

  -- addition to pick up perturbation flux for rtheta_pp equation
  if (rk_step > 0) then
    for iCell in cell_range do
      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

        var flux = cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge
                    * (er[{iEdge, iCell.y}].ru_save - er[{iEdge, iCell.y}].ru) * 0.5
                    * (cr[{cell2, iCell.y}].theta_m_save + cr[{cell1, iCell.y}].theta_m_save)
        cr[iCell].tend_theta -= flux -- division by areaCell picked up down below
      end
    end
  end

  -- horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
  -- but here we can also code in hyperdiffusion if we wish (2nd order at present)
  if (rk_step == 0) then
    for iCell in cell_range do
      cr[iCell].delsq_theta = 0.0
      cr[iCell].tend_theta_euler = 0.0
      var r_areaCell = cr[{iCell.x, 0}].invAreaCell
      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
        var edge_sign = r_areaCell * cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge * er[{iEdge, 0}].invDcEdge
        var pr_scale = prandtl_inv * er[{iEdge, 0}].meshScalingDel2
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

        -- we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below
        var theta_turb_flux = edge_sign * (cr[{cell2, iCell.y}].theta_m - cr[{cell1, iCell.y}].theta_m) * er[{iEdge, iCell.y}].rho_edge
        cr[iCell].delsq_theta += theta_turb_flux
        theta_turb_flux *= 0.5 * (cr[{cell1, iCell.y}].kdiff + cr[{cell2, iCell.y}].kdiff) * pr_scale
        cr[iCell].tend_theta_euler += theta_turb_flux
      end
    end

    if (h_theta_eddy_visc4 > 0.0) then  -- 4th order mixing is active

      for iCell in cell_range do
        var r_areaCell = h_theta_eddy_visc4 * prandtl_inv * cr[{iCell.x, 0}].invAreaCell
        for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do

          var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
          var edge_sign = er[{iEdge, 0}].meshScalingDel4 * r_areaCell * er[{iEdge, 0}].dvEdge
                          * cr[{iCell.x, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].invDcEdge

          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

          cr[iCell].tend_theta_euler -= edge_sign * (cr[{cell2, iCell.y}].delsq_theta - cr[{cell1, iCell.y}].delsq_theta)
        end
      end
    end -- 4th order mixing is active
  end -- theta mixing calculated first rk_step

  -- vertical advection plus diabatic term
  -- Note: we are also dividing through by the cell area after the horizontal flux divergence

  for iCell in cell_range do
    cr[iCell].wdtz = 0.0
    -- Don't change the order of these statements!
    if (iCell.y > 0 and iCell.y < nVertLevels - 1) then
      cr[iCell].wdtz = ( (cr[iCell].rw_save - cr[iCell].rw) * (vert_r[iCell.y].fzm * cr[iCell].theta_m_save
                          + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].theta_m_save) )
    end
    if (iCell.y == 1) then
      cr[iCell].wdtz += cr[iCell].rw * (vert_r[iCell.y].fzm * cr[iCell].theta_m + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].theta_m)
    end
    if (iCell.y == nVertLevels - 1) then
      cr[iCell].wdtz = cr[iCell].rw_save * (vert_r[iCell.y].fzm * cr[iCell].theta_m_save
                       + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].theta_m_save)  -- rtheta_pp redefinition
    end
  end

  for iCell in cell_range do
    cr[iCell].tend_theta *= cr[{iCell.x, 0}].invAreaCell - vert_r[iCell.y].rdzw * (cr[iCell + {0, 1}].wdtz - cr[iCell].wdtz)
    cr[iCell].tend_rtheta_adv = cr[iCell].tend_theta -- this is for the Tiedke scheme -- Note: this value is never used
    cr[iCell].rthdynten = cr[iCell].tend_theta / cr[iCell].rho_zz -- this is for the Grell-Freitas scheme -- Note: this value is never used
    cr[iCell].tend_theta += cr[iCell].rho_zz * cr[iCell].rt_diabatic_tend
  end

  --  vertical mixing for theta - 2nd order
  if (rk_step == 0) then

    if (v_theta_eddy_visc2 > 0.0) then  -- vertical mixing for theta_m

      if (config_mix_full) then

        for iCell in cell_range do
          if (iCell.y > 0 and iCell.y < nVertLevels - 1) then
            var z1 = cr[iCell - {0, 1}].zgrid
            var z2 = cr[iCell         ].zgrid
            var z3 = cr[iCell + {0, 1}].zgrid
            var z4 = cr[iCell + {0, 2}].zgrid

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)

            cr[iCell].tend_theta_euler += v_theta_eddy_visc2 * prandtl_inv * cr[iCell].rho_zz
                                          * ( (cr[iCell + {0, 1}].theta_m - cr[iCell].theta_m) / (zp-z0)
                                            - (cr[iCell].theta_m - cr[iCell - {0, 1}].theta_m) / (z0-zm) )
                                          / (0.5 * (zp - zm))
          end
        end

      else  -- idealized cases where we mix on the perturbation from the initial 1-D state
        for iCell in cell_range do
          if (iCell.y > 0 and iCell.y < nVertLevels - 1) then
            var z1 = cr[iCell - {0, 1}].zgrid
            var z2 = cr[iCell         ].zgrid
            var z3 = cr[iCell + {0, 1}].zgrid
            var z4 = cr[iCell + {0, 2}].zgrid

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)
            cr[iCell].tend_theta_euler += v_theta_eddy_visc2 * prandtl_inv * cr[iCell].rho_zz
                                          * ( ((cr[iCell + {0, 1}].theta_m - cr[iCell + {0, 1}].t_init)
                                                - (cr[iCell].theta_m - cr[iCell].t_init)) / (zp-z0)
                                          - ( (cr[iCell].theta_m - cr[iCell].t_init)
                                              - (cr[iCell - {0, 1}].theta_m - cr[iCell - {0, 1}].t_init)) / (z0-zm) )
                                          / (0.5*(zp-zm))
          end
        end
      end
    end
  end -- compute vertical theta mixing on first rk_step

  for iCell in cell_range do
    cr[iCell].tend_theta += cr[iCell].tend_theta_euler + cr[iCell].tend_rtheta_physics
  end
end

-- There is a special case when rthdynten is not associated with packages. Ignoring for now by
-- assuming that they are associated properly
task atm_compute_dyn_tend(cr : region(ispace(int2d), cell_fs),
                          er : region(ispace(int2d), edge_fs),
                          vr : region(ispace(int2d), vertex_fs),
                          vert_r : region(ispace(int1d), vertical_fs),
                          rk_step : int,
                          dt : double,
                          config_horiz_mixing : regentlib.string,
                          config_mpas_cam_coef : double,
                          config_mix_full : bool,
                          config_rayleigh_damp_u : bool)
where
  reads writes (cr, er, vr, vert_r)
do

  format.println("Calling atm_compute_dyn_tend_work...")
  atm_compute_dyn_tend_work(cr, er, vr, vert_r, rk_step, dt, config_horiz_mixing, config_mpas_cam_coef, config_mix_full, config_rayleigh_damp_u)
end

__demand(__cuda)
task atm_set_smlstep_pert_variables_work(cr : region(ispace(int2d), cell_fs),
                                         er : region(ispace(int2d), edge_fs),
                                         vert_r : region(ispace(int1d), vertical_fs))
where
  reads (cr.{bdyMaskCell, edgesOnCell, edgesOnCell_sign, nEdgesOnCell, zb_cell, zb3_cell, zz},
         er.u_tend, vert_r.{fzm, fzp}),
  reads writes (cr.w)
do

  format.println("Calling atm_set_smlstep_pert_variables_work...")

  var cell_range = rect2d { int2d {0, 1}, int2d {nCells - 1, nVertLevels - 1} } --All vertical loops in the original code go from 2 to nVertLevels

  for iCell in cell_range do
    if (cr[{iCell.x, 0}].bdyMaskCell <= nRelaxZone) then
      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
        var flux = cr[{iCell.x, 0}].edgesOnCell_sign[i] * (vert_r[iCell.y].fzm * er[{iEdge, iCell.y}].u_tend + vert_r[iCell.y].fzp * er[{iEdge, iCell.y - 1}].u_tend)
        --sign function copied as copysign, _RKIND removed as done previously
        cr[iCell].w -= (cr[iCell].zb_cell[i] + copysign(1.0, er[{iEdge, iCell.y}].u_tend) * cr[{iCell, iCell.y}].zb3_cell[i]) * flux
      end

      cr[iCell].w *= ( vert_r[iCell.y].fzm * cr[iCell].zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz )
    end
  end
end

task atm_set_smlstep_pert_variables(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    vert_r : region(ispace(int1d), vertical_fs))
where
  reads writes (cr, er, vert_r)
do

  cio.printf("set small step vars\n")
  atm_set_smlstep_pert_variables_work(cr, er, vert_r)
end

--Comments for atm_advance_acoustic_step_work
--dts is passed in as double now - in code, passed in as rk_sub_timestep(rk_step)
--1.0_RKIND translated as just 1.0
--Tendency variables: tend_rw, tend_rt, tend_rho (added to CR), tend_ru (added to ER). Not in registry
--Other variables not in registry: rs, ts: added to CR as part of vertical grid, ru_Avg (added to ER)
--__demand(__cuda)
task atm_advance_acoustic_step_work(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    vert_r : region(ispace(int1d), vertical_fs),
                                    dts : double,
                                    small_step : int)
                                    -- nCellsSolve : int)
where
  reads (cr.{a_tri, alpha_tri, coftz, cofwr, cofwt, cofwz, dss, edgesOnCell, edgesOnCellSign, exner,
             gamma_tri, invAreaCell, nEdgesOnCell, rho_pp, rho_zz, rtheta_pp, rw, rw_save, specZoneMaskCell,
             tend_rho, theta_m, theta_m, w, zz},
         er.{cellsOnEdge, cqu, dvEdge, invDcEdge, specZoneMaskEdge, tend_ru, zxu},
         vert_r.{cofrz, fzm, fzp, rdzw}),
  writes (cr.rtheta_pp_old),
  reads writes (cr.{rho_pp, rtheta_pp, rw_p, wwAvg},
                er.{ruAvg, ru_p})
do

  format.println("Calling atm_advance_acoustic_step_work...")

  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var cell_range_P1 = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels} }
  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }

  var epssm = constants.config_epssm
  var rgas = constants.rgas
  var rcv = rgas / (constants.cp - rgas)
  var c2 = constants.cp * rcv
  var resm = (1.0 - epssm) / (1.0 + epssm)
  var rdts = 1.0 / dts

  var rs : double[nVertLevels]
  var ts : double[nVertLevels]

  format.println("Calling atm_advance_acoustic_step_work...")

  if (small_step ~= 0) then -- not needed on first small step
    -- forward-backward acoustic step integration.
    -- begin by updating the horizontal velocity u,
    -- and accumulating the contribution from the updated u to the other tendencies.
    for iEdge in edge_range do
      var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]

      -- update edges for block-owned cells
      --if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

        --var pgrad = ((cr[{cell2, iEdge.y}].rtheta_pp - cr[{cell1, iEdge.y}].rtheta_pp) * er[{iEdge.x, 0}].invDcEdge) / (0.5 * (cr[{cell2, iEdge.y}].zz + cr[{cell1, iEdge.y}].zz))
        --pgrad *= er[iEdge].cqu * 0.5 * c2 * (cr[{cell1, iEdge.y}].exner + cr[{cell2, iEdge.y}].exner)
        --pgrad += 0.5 * er[iEdge].zxu * constants.gravity * (cr[{cell1, iEdge.y}].rho_pp + cr[{cell2, iEdge.y}].rho_pp)
        --er[iEdge].ru_p += dts * (er[iEdge].tend_ru - (1.0 - er[{iEdge.x, 0}].specZoneMaskEdge) * pgrad)  --NEEDS FIXING

        -- accumulate ru_p for use later in scalar transport
        --er[iEdge].ruAvg += er[iEdge].ru_p
      --end -- end test for block-owned cells
    end -- end loop over edges

  else -- this is all that us needed for ru_p update for first acoustic step in RK substep
    for iEdge in edge_range do
      var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]

      -- update edges for block-owned cells
      --if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then
        --er[iEdge].ru_p = dts * er[iEdge].tend_ru
        --er[iEdge].ruAvg = er[iEdge].ru_p
      --end -- end test for block-owned cells
    end -- end loop over edges
  end -- test for first acoustic step

  if (small_step == 0) then -- initialize here on first small timestep.
    for iCell in cell_range do
      cr[iCell].rtheta_pp_old = 0
    end
  else
    for iCell in cell_range do
      cr[iCell].rtheta_pp_old = cr[iCell].rtheta_pp
    end
  end

  for iCell in cell_range_P1 do
    if (small_step == 0) then
      cr[iCell].wwAvg = 0
      cr[iCell].rw_p = 0
    end
  end

  for iCell in cell_range do
    if (small_step == 0) then
      cr[iCell].rho_pp = 0
      cr[iCell].rtheta_pp = 0
    end

    if (cr[{iCell, 0}].specZoneMaskCell == 0.0) then -- not specified zone, compute...
      for i = 0, nVertLevels do
        ts[i] = 0
        rs[i] = 0
      end

      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]    --edgesOnCell(i,iCell)
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]    --cell2 = cellsOnEdge(2,iEdge)

        var flux = cr[{iCell.x, 0}].edgesOnCellSign[i] * dts * er[{iEdge, 0}].dvEdge * er[{iEdge, iCell.y}].ru_p * cr[{iCell.x, 0}].invAreaCell
        rs[iCell.y] -= flux
        ts[iCell.y] -= flux * 0.5 * (cr[{cell2, iCell.y}].theta_m + cr[{cell1, iCell.y}].theta_m)
      end

      -- vertically implicit acoustic and gravity wave integration.
      -- this follows Klemp et al MWR 2007, with the addition of an implicit Rayleigh damping of w
      -- serves as a gravity-wave absorbing layer, from Klemp et al 2008.
      rs[iCell.y] = cr[iCell].rho_pp + dts * cr[iCell].tend_rho + rs[iCell.y] - vert_r[iCell.y].cofrz* resm * (cr[iCell + {0, 1}].rw_p - cr[iCell].rw_p)
      ts[iCell.y] = cr[iCell].rtheta_pp + dts * cr[iCell].theta_m + ts[iCell.y] - resm * vert_r[iCell.y].rdzw * (cr[iCell + {0, 1}].coftz * cr[iCell + {0, 1}].rw_p  - cr[iCell].coftz * cr[iCell].rw_p )

      if (iCell.y > 0) then
        cr[iCell].wwAvg += 0.5 * (1.0 - epssm) * cr[iCell].rw_p
        cr[iCell].rw_p += dts * cr[iCell].w - cr[iCell].cofwz
                          * ((cr[iCell].zz * ts[iCell.y] - cr[iCell - {0, 1}].zz * ts[iCell.y - 1]) + resm
                          * (cr[iCell].zz * cr[iCell].rtheta_pp - cr[iCell - {0, 1}].zz * cr[iCell - {0, 1}].rtheta_pp))
                          - cr[iCell].cofwr * ((rs[iCell.y] + rs[iCell.y - 1]) + resm * (cr[iCell].rho_pp + cr[iCell - {0, 1}].rho_pp))
                          + cr[iCell].cofwt * (ts[iCell.y] + resm * cr[iCell].rtheta_pp)
                          + cr[iCell - {0, 1}].cofwt * (ts[iCell.y - 1] +resm * cr[iCell - {0, 1}].rtheta_pp)

        -- tridiagonal solve sweeping up and then down the column
        cr[iCell].rw_p -= cr[iCell].a_tri * cr[iCell - {0, 1}].rw_p
        cr[iCell].rw_p *= cr[iCell].alpha_tri
      end

      --TODO: how to parallelize this???
      --for k = nVertLevels, 1, -1 do
      --  cr[{iCell, k}].rw_p -= cr[{iCell, k}].gamma_tri * cr[{iCell, k+1}].rw_p
      --end


      -- the implicit Rayleigh damping on w (gravity-wave absorbing)
      if (iCell.y > 0) then
        cr[iCell].rw_p += (cr[iCell].rw_save - cr[iCell].rw) - dts * cr[iCell].dss
                          * (vert_r[iCell.y].fzm * cr[iCell].zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz)
                          * (vert_r[iCell.y].fzm * cr[iCell].rho_zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].rho_zz) * cr[iCell].w
        cr[iCell].rw_p /= (1.0 + dts * cr[iCell].dss)
        cr[iCell].rw_p -= (cr[iCell].rw_save - cr[iCell].rw)

        -- accumulate (rho*omega)' for use later in scalar transport
        cr[iCell].wwAvg += 0.5 * (1.0 + epssm) * cr[iCell].rw_p
      end

      -- update rho_pp and theta_pp given updated rw_p

      cr[iCell].rho_pp = rs[iCell.y] - vert_r[iCell.y].cofrz*(cr[iCell + {0, 1}].rw_p - cr[iCell].rw_p)
      cr[iCell].rtheta_pp = ts[iCell.y] - vert_r[iCell.y].rdzw
                            * (cr[iCell + {0, 1}].coftz * cr[iCell + {0, 1}].rw_p - cr[iCell].coftz * cr[iCell].rw_p)

    else -- specifed zone in regional_MPAS
      cr[iCell].rho_pp  = cr[iCell].rho_pp + dts * cr[iCell].tend_rho
      cr[iCell].rtheta_pp = cr[iCell].rtheta_pp + dts * cr[iCell].theta_m
      cr[iCell].rw_p = cr[iCell].rw_p + dts * cr[iCell].w
      cr[iCell].wwAvg = cr[iCell].wwAvg + 0.5 * (1.0+epssm) * cr[iCell].rw_p
    end
  end -- end of loop over cells
end

task atm_advance_acoustic_step(cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs),
                               vert_r : region(ispace(int1d), vertical_fs),
                               dts : double,
                               small_step : int)
                               -- nCellsSolve : int)
where
  reads writes (er, cr, vert_r)
do

  cio.printf("advancing acoustic step\n")
  atm_advance_acoustic_step_work(cr, er, vert_r, dts, small_step)
end

-- Comments:
-- dts is passed in as double - in code, passed in as rk_sub_timestep(rk_step)
-- 1.0_RKIND and 2.0_RKIND translated as 1.0 and 2.0
-- This function also contains nCellsSolve, which has not been resolved yet
__demand(__cuda)
task atm_divergence_damping_3d(cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs),
                               dts : double)
where
  reads (cr.{rtheta_pp, rtheta_pp_old, theta_m}, er.{cellsOnEdge, specZoneMaskEdge}),
  reads writes (er.ru_p)
do

  format.println("Calling atm_divergence_damping_3d...")

  var smdiv = constants.config_smdiv
  var rdts = 1.0 / dts
  var coef_divdamp = 2.0 * smdiv * constants.config_len_disp * rdts

  var edge_range_1d = rect1d { 0, nEdges - 1 }

  for iEdgex in edge_range_1d do

    var cell1 = er[{iEdgex, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdgex, 0}].cellsOnEdge[1]

    -- update edges for block-owned cells
    -- if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

      for k = 0, nVertLevels do

        -- scaled 3d divergence damping
        var divCell1 = -(cr[{cell1, k}].rtheta_pp - cr[{cell1, k}].rtheta_pp_old)
        var divCell2 = -(cr[{cell2, k}].rtheta_pp - cr[{cell2, k}].rtheta_pp_old)
        er[{iEdgex, k}].ru_p += coef_divdamp * (divCell2 - divCell1) *
                              (1.0 - er[{iEdgex, 0}].specZoneMaskEdge)
                              / (cr[{cell1, k}].theta_m + cr[{cell2, k}].theta_m)
      end
    -- end
  end
end

__demand(__cuda)
task atm_recover_large_step_variables_work(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    vert_r : region(ispace(int1d), vertical_fs),
                                    ns : int,
                                    rk_step : int,
                                    dt : double)
where
  reads (cr.{bdyMaskCell, edgesOnCell, edgesOnCell_sign, exner_base, nEdgesOnCell, rho_base, rho_p_save, rho_pp,
             rt_diabatic_tend, rtheta_base, rtheta_p_save, rtheta_pp, rw_p, rw_save, zb_cell, zb3_cell, zz},
         er.{cellsOnEdge, ru_p, ru_save},
         vert_r.{cf1, cf2, cf3, fzm, fzp}),
  writes (cr.{pressure_p, theta_m}, er.u),
  reads writes (cr.{exner, rho_p, rho_zz, rtheta_p, rw, w, wwAvg},
                er.{ruAvg, ru})
do

  var rgas = constants.rgas
  var rcv = rgas / (constants.cp - rgas)
  var p0 = 100000

  var garbage_range = rect2d { int2d {nCells, 0}, int2d {nCells, nVertLevels - 1} }
  var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var edge_range = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }

  -- Avoid FP errors caused by a potential division by zero below by
  -- initializing the "garbage cell" of rho_zz to a non-zero value
  for iCell in garbage_range do
    cr[iCell].rho_zz = 1.0
  end

  -- compute new density everywhere so we can compute u from ru.
  -- we will also need it to compute theta_m below
  var invNs = 1 / [double](ns)

  for iCell in cell_range do
    cr[iCell].rho_p = cr[iCell].rho_p_save + cr[iCell].rho_pp
    cr[iCell].rho_zz = cr[iCell].rho_p + cr[iCell].rho_base

    cr[iCell].w = 0.0

    cr[iCell].wwAvg *= invNs
    cr[iCell].wwAvg += cr[iCell].rw_save
    cr[iCell].rw = cr[iCell].rw_save + cr[iCell].rw_p
    -- pick up part of diagnosed w from omega - divide by density later
    cr[iCell].w = cr[iCell].rw / (vert_r[iCell.y].fzm * cr[iCell].zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].zz)

    if (iCell.y == nVertLevels) then
      cr[iCell].w = 0.0
    end

    if (rk_step == 2) then
      cr[iCell].rtheta_p = cr[iCell].rtheta_p_save + cr[iCell].rtheta_pp - dt * cr[iCell].rho_zz * cr[iCell].rt_diabatic_tend
      cr[iCell].theta_m = (cr[iCell].rtheta_p + cr[iCell].rtheta_base) / cr[iCell].rho_zz
      cr[iCell].exner = cr[iCell].zz * (rgas / p0) * pow((cr[iCell].rtheta_p + cr[iCell].rtheta_base), rcv)
      -- pressure_p is perturbation pressure
      cr[iCell].pressure_p = cr[iCell].zz * rgas * (cr[iCell].exner * cr[iCell].rtheta_p + cr[iCell].rtheta_base * (cr[iCell].exner - cr[iCell].exner_base))
    else
      cr[iCell].rtheta_p = cr[iCell].rtheta_p_save + cr[iCell].rtheta_pp
      cr[iCell].theta_m = (cr[iCell].rtheta_p + cr[iCell].rtheta_base) / cr[iCell].rho_zz
    end
  end

  -- recover time-averaged ruAvg on all edges of owned cells (for upcoming scalar transport).
  -- we solved for these in the acoustic-step loop.
  -- we will compute ru and u here also, given we are here, even though we only need them on nEdgesSolve

  for iEdge in edge_range do
    var cell1 = er[{iEdge.x, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge.x, 0}].cellsOnEdge[1]
    er[iEdge].ruAvg *= invNs
    er[iEdge].ruAvg += er[iEdge].ru_save
    er[iEdge].ru = er[iEdge].ru_save * er[iEdge].ru_p
    er[iEdge].u = 2 * er[iEdge].ru / (cr[{cell1, iEdge.y}].rho_zz + cr[{cell2, iEdge.y}].rho_zz)
  end

  for iCell in cell_range do

    -- finish recovering w from (rho*omega)_p.  as when we formed (rho*omega)_p from u and w, we need
    -- to use the same flux-divergence operator as is used for the horizontal theta transport
    -- (See Klemp et al 2003).

    if (cr[{iCell.x, 0}].bdyMaskCell <= nRelaxZone) then
      for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
        var flux = (vert_r[0].cf1 * er[{iEdge, 0}].ru + vert_r[0].cf2 * er[{iEdge, 1}].ru + vert_r[0].cf3 * er[{iEdge, 2}].ru)
        cr[{iCell.x, 0}].w += cr[{iCell.x, 0}].edgesOnCell_sign[i] * (cr[{iCell.x, 0}].zb_cell[i] + copysign(1.0, flux) * cr[{iCell.x, 0}].zb3_cell[i]) * flux
        var flux2 = vert_r[iCell.y].fzm * er[{iEdge, iCell.y}].ru * (vert_r[iCell.y].fzp * er[{iEdge, iCell.y - 1}].ru)
        cr[iCell].w += cr[{iCell, 0}].edgesOnCell_sign[i] * (cr[iCell].zb_cell[i] + copysign(1.0, flux2) * cr[iCell].zb3_cell[i]) * flux2
      end
    end
  end

  for iCell in cell_range do
    if (cr[{iCell.x, 0}].bdyMaskCell <= nRelaxZone) then
      if (iCell.y == 0) then
        cr[{iCell.x, 0}].w /= (vert_r[0].cf1 * cr[{iCell.x, 0}].rho_zz + vert_r[0].cf2 * cr[{iCell.x, 1}].rho_zz + vert_r[0].cf3 * cr[{iCell.x, 2}].rho_zz)
      end

      if (iCell.y > 0) then
        cr[iCell].w /= (vert_r[iCell.y].fzm * cr[iCell].rho_zz + vert_r[iCell.y].fzp * cr[iCell - {0, 1}].rho_zz)
      end
    end
  end
end


task atm_recover_large_step_variables(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    vert_r : region(ispace(int1d), vertical_fs),
                                    ns : int,
                                    rk_step : int,
                                    dt : double)
where
  reads writes (cr, er, vert_r)
do
  cio.printf("recovering large step vars\n")
  atm_recover_large_step_variables_work(cr, er, vert_r, ns, rk_step, dt)

end

-- Comments on mpas_reconstruct
-- In Fortran, mpas_reconstruct is an interface with 1d and 2d versions.
-- The program decides which to call based on the types of the arguments passed in.
-- In the dynamics loop, only the 2d version is called.
__demand(__cuda)
task mpas_reconstruct_2d(cr : region(ispace(int2d), cell_fs),
                         er : region(ispace(int2d), edge_fs),
                         includeHalos : bool,
                         on_a_sphere : bool)
where
  reads (cr.{coeffs_reconstruct, edgesOnCell, lat, lon, nEdgesOnCell}, er.u),
  writes (cr.{uReconstructMeridional, uReconstructZonal}),
  reads writes (cr.{uReconstructX, uReconstructY, uReconstructZ})
do

  format.println("Calling mpas_reconstruct_2d...")

  -- The Fortran code has includeHalos as an optional argument. For now, we will make it mandatory,
  -- and pass in false whenever it is not passed in in the Fortran code.
  var nCellsReconstruct = nCells
  --if (includeHalos) then
  --  nCellsReconstruct = nCellsSolve
  --end
  var cell_range = rect2d { int2d {0, 0}, int2d {nCellsReconstruct - 1, nVertLevels - 1} }

  -- initialize the reconstructed vectors
  for iCell in cell_range do
    cr[iCell].uReconstructX = 0.0
    cr[iCell].uReconstructY = 0.0
    cr[iCell].uReconstructZ = 0.0
  end

  for iCell in cell_range do
    -- a more efficient reconstruction where rbf_values*matrix_reconstruct
    -- has been precomputed in coeffs_reconstruct
    for i = 0, cr[{iCell.x, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell.x, 0}].edgesOnCell[i]
      cr[iCell].uReconstructX += cr[{iCell.x, 0}].coeffs_reconstruct[i][0] * er[{iEdge, iCell.y}].u
      cr[iCell].uReconstructY += cr[{iCell.x, 0}].coeffs_reconstruct[i][1] * er[{iEdge, iCell.y}].u
      cr[iCell].uReconstructZ += cr[{iCell.x, 0}].coeffs_reconstruct[i][2] * er[{iEdge, iCell.y}].u
    end
  end -- iCell

  if (on_a_sphere) then
    for iCell in cell_range do
      var clat = cos(cr[{iCell.x, 0}].lat)
      var slat = sin(cr[{iCell.x, 0}].lat)
      var clon = cos(cr[{iCell.x, 0}].lon)
      var slon = sin(cr[{iCell.x, 0}].lon)
      cr[iCell].uReconstructZonal = -cr[iCell].uReconstructX * slon + cr[iCell].uReconstructY * clon
      cr[iCell].uReconstructMeridional = -(cr[iCell].uReconstructX * clon + cr[iCell].uReconstructY * slon) * slat
                                         + cr[iCell].uReconstructZ * clat
    end
  else
    for iCell in cell_range do
      cr[iCell].uReconstructZonal = cr[iCell].uReconstructX
      cr[iCell].uReconstructMeridional = cr[iCell].uReconstructY
    end
  end
end

__demand(__cuda)
task atm_rk_dynamics_substep_finish(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    dynamics_substep : int,
                                    dynamics_split : int)
where
  reads (cr.{rho_p, rho_zz_2, rho_zz_old_split, rtheta_p, rw, theta_m_2, w_2}, er.{ru, u_2}),
  writes (cr.{rho_p_save, rho_zz, rtheta_p_save, rw_save, theta_m, w}, er.{ru_save, u}),
  reads writes (cr.{wwAvg, wwAvg_split}, er.{ruAvg, ruAvg_split})
do

  format.println("Calling atm_rk_dynamics_substep_finish...")
  var inv_dynamics_split = 1.0 / [double](dynamics_split)
  var edge_range = rect2d { int2d{0, 0}, int2d{nEdges - 1, nVertLevels - 1} }
  var cell_range = rect2d { int2d{0, 0}, int2d{nCells - 1, nVertLevels - 1} }

  if (dynamics_substep < dynamics_split) then
    for iEdge in edge_range do
      er[iEdge].ru_save = er[iEdge].ru
      er[iEdge].u = er[iEdge].u_2
    end
    for iCell in cell_range do
      cr[iCell].rw_save = cr[iCell].rw
      cr[iCell].rtheta_p_save = cr[iCell].rtheta_p
      cr[iCell].rho_p_save = cr[iCell].rho_p

      cr[iCell].w = cr[iCell].w_2
      cr[iCell].theta_m = cr[iCell].theta_m_2
      cr[iCell].rho_zz = cr[iCell].rho_zz_2
    end
  end

  if (dynamics_substep == 1) then
    for iEdge in edge_range do
      er[iEdge].ruAvg_split = er[iEdge].ruAvg
    end
    for iCell in cell_range do
      cr[iCell].wwAvg_split = cr[iCell].wwAvg
    end
  else
    for iEdge in edge_range do
      er[iEdge].ruAvg_split = er[iEdge].ruAvg + er[iEdge].ruAvg_split
    end
    for iCell in cell_range do
      cr[iCell].wwAvg_split = cr[iCell].wwAvg + cr[iCell].wwAvg_split
    end
  end

  if (dynamics_substep == dynamics_split) then
    for iEdge in edge_range do
      er[iEdge].ruAvg = er[iEdge].ruAvg_split * inv_dynamics_split
    end
    for iCell in cell_range do
      cr[iCell].wwAvg = cr[iCell].wwAvg_split * inv_dynamics_split
      cr[iCell].rho_zz = cr[iCell].rho_zz_old_split
    end
  end
end
