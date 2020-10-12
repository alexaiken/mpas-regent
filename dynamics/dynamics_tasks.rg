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

local nCells = 2562
local nEdges = 7680
local nVertices = 5120
local maxEdges = 10
local maxEdges2 = 20
local TWO = 2
local FIFTEEN = 15
local vertexDegree = 3
local nVertLevels = 1
local nRelaxZone = 5
local seconds_per_day = 86400.0


local cio = terralib.includec("stdio.h")
local cmath = terralib.includec("math.h")

__demand(__cuda)
task atm_compute_signs_pt1(cr : region(ispace(int2d), cell_fs),
                        er : region(ispace(int2d), edge_fs),
                        vr : region(ispace(int2d), vertex_fs))
where reads writes(vr, cr), reads (er) do

    var range = rect2d { int2d {0, 0}, int2d {nVertices - 1, 0} }
    for iVtx in range do
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
end


task atm_compute_signs_pt2(cr : region(ispace(int2d), cell_fs),
                          er : region(ispace(int2d), edge_fs),
                          vr : region(ispace(int2d), vertex_fs))
where reads writes(vr, cr), reads (er) do


      for iCell=0,nCells do
         for i=0, cr[{iCell, 0}].nEdgesOnCell do
            if (cr[{iCell, 0}].edgesOnCell[i] <= nEdges) then
              if (iCell == er[{cr[{iCell, 0}].edgesOnCell[i], 0}].cellsOnEdge[0]) then
                  cr[{iCell, 0}].edgesOnCellSign[i] = 1.0
                  --VERTICAL STRUCTURE ACCESSED--
                  for vLevel = 0, nVertLevels + 1 do

                    cr[{iCell, vLevel}].zb_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb[0]
                    cr[{iCell, vLevel}].zb3_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb3[0]

                    --cio.printf("zb at cell %d and level %d, index %d, is %f \n", iCell, vLevel, i, cr[{iCell, vLevel}].zb_cell[i])
                    --cio.printf("zb3 at cell %d and level %d, index %d, is %f \n", iCell, vLevel, i, cr[{iCell, vLevel}].zb3_cell[i])

                  end

               else
                  cr[{iCell, 0}].edgesOnCellSign[i] = -1.0
                  --VERTICAL--
                  for vLevel = 0, nVertLevels + 1 do

                    cr[{iCell, vLevel}].zb_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb[1]
                    cr[{iCell, vLevel}].zb3_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb3[1]
                    --cio.printf("zb at cell %d and level %d, index %d, is %f \n", iCell, vLevel, i, cr[{iCell, vLevel}].zb_cell[i])
                    --cio.printf("zb3 at cell %d and level %d, index %d, is %f \n", iCell, vLevel, i, cr[{iCell, vLevel}].zb3_cell[i])

                  end
               end
            else
               cr[{iCell, 0}].edgesOnCellSign[i] = 0.0
            end
         end
      end


    for iCell=0, nCells do
        for i=0, cr[{iCell, 0}].nEdgesOnCell do
            var iVtx = cr[{iCell, 0}].verticesOnCell[i]
            if (iVtx <= nVertices) then
                for j=1,vertexDegree do
                    if (iCell == vr[{iVtx, 0}].cellsOnVertex[j]) then
                        cr[{iCell, 0}].kiteForCell[i] = j
                        break
                    end
                end
                -- trimmed a log statement here
            else
                cr[{iCell, 0}].kiteForCell[i] = 1
            end
            --cio.printf("cr[{%d, 0}].kiteForCell[%d] is %f \n", iCell, i, cr[{iCell, 0}].kiteForCell[i])
        end
    end
end


task atm_adv_coef_compression(cr : region(ispace(int2d), cell_fs),
                              er : region(ispace(int2d), edge_fs))
where reads writes(er), reads(cr) do

    var cell_list : int[maxEdges]

    for iEdge = 0, nEdges do
        er[{iEdge, 0}].nAdvCellsForEdge = 0
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

         --
         -- do only if this edge flux is needed to update owned cells
         --
        if (cell1 <= nCells or cell2 <= nCells) then
            cell_list[0] = cell1
            cell_list[1] = cell2
            var n = 2 -- n is number of cells currently in list

          --  add cells surrounding cell 1.  n is number of cells currently in list
            for i = 0, cr[{cell1, 0}].nEdgesOnCell do
               if (cr[{cell1, 0}].cellsOnCell[i] ~= cell2) then
                  n = n + 1
                  cell_list[n] = cr[{cell1, 0}].cellsOnCell[i]
               end
            end

          --  add cells surrounding cell 2 (brute force approach)
            for iCell = 0, cr[{cell2, 0}].nEdgesOnCell do
               var addcell = true
               for i=0, n do
                  if (cell_list[i] == cr[{cell2, 0}].cellsOnCell[iCell]) then addcell = false end
               end
               if (addcell) then
                  n = n+1
                  cell_list[n] = cr[{cell2, 0}].cellsOnCell[iCell]
               end
            end


            er[{iEdge, 0}].nAdvCellsForEdge = n
            for iCell = 0, er[{iEdge, 0}].nAdvCellsForEdge do
               er[{iEdge, 0}].advCellsForEdge[iCell] = cell_list[iCell]
            end

          -- we have the ordered list, now construct coefficients
            for coef = 0, FIFTEEN do
                er[{iEdge, 0}].adv_coefs[coef] = 0.0
                er[{iEdge, 0}].adv_coefs_3rd[coef] = 0.0
            end -- initialize list to 0

          -- pull together third and fourth order contributions to the flux
          -- first from cell1

            var j_in = 0
            for j=0, n do
               if (cell_list[j] == cell1 ) then j_in = j end
            end
            er[{iEdge, 0}].adv_coefs[j_in]= er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[0]
            er[{iEdge, 0}].adv_coefs_3rd[j_in]= er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[0]

            for iCell = 0, cr[{cell1, 0}].nEdgesOnCell do
               j_in = 0
               for j=0, n do
                 if( cell_list[j] == cr[{cell1, 0}].cellsOnCell[iCell]) then j_in = j end
               end
               er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[iCell * FIFTEEN + 0]
               er[{iEdge, 0}].adv_coefs_3rd[j_in] = er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[iCell * FIFTEEN + 0]
            end

          -- pull together third and fourth order contributions to the flux
          -- now from cell2

            j_in = 0
            for j=0, n do
               if( cell_list[j] == cell2 ) then j_in = j end
            end
              er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[1]
              er[{iEdge, 0}].adv_coefs_3rd[j_in] = er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[1]

            for iCell = 0, cr[{cell2, 0}].nEdgesOnCell do
               j_in = 0
               for j=0, n do
                  if( cell_list[j] == cr[{cell2, 0}].cellsOnCell[iCell] ) then j_in = j end
               end
               er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[iCell * FIFTEEN + 1]
               er[{iEdge, 0}].adv_coefs_3rd[j_in] = er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[iCell * FIFTEEN + 1]
            end

            for j = 0,n do
              -- TODO: how to calculate exponent?
              -- er[{iEdge, 0}].adv_coefs[j] =  -1.0 * er[{iEdge, 0}].dcEdge**2 * er[{iEdge, 0}].adv_coefs[j] / 12 -- this should be a negative number
              -- er[{iEdge, 0}].adv_coefs_3rd[j] =  -1.0 * er[{iEdge, 0}].dcEdge**2 * er[{iEdge, 0}].adv_coefs_3rd[j] / 12
--               adv_coefs    (j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs    (j,iEdge) / 12.
--               adv_coefs_3rd(j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
            end

          -- 2nd order centered contribution - place this in the main flux weights

            j_in = 0
            for j=0, n do
               if( cell_list[j] == cell1 ) then j_in = j end
            end
            er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + 0.5

            j_in = 0
            for j=0, n do
               if( cell_list[j] == cell2 ) then j_in = j end
            end
            er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + 0.5

          --  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply

            for j=0,n do
               er[{iEdge, 0}].adv_coefs[j] = er[{iEdge, 0}].dvEdge * er[{iEdge, 0}].adv_coefs[j]
               er[{iEdge, 0}].adv_coefs_3rd[j] = er[{iEdge, 0}].dvEdge * er[{iEdge, 0}].adv_coefs_3rd[j]
            end

         end  -- only do for edges of owned-cells

    end -- end loop over edges


end


--config_zd: default 22000.0, config_xnutr: default 0.2. From config
task atm_compute_damping_coefs(config_zd : double, config_xnutr : double, cr : region(ispace(int2d), cell_fs))
where reads writes (cr) do
  var m1 = -1.0
  var pii = cmath.acos(m1) -- find equivelelt transformation in Regent for acos()
  --cio.printf("pii = %f\n", pii)

  var dx_scale_power = 1.0
  fill(cr.dss, 0)
  --cio.printf("cr[{%d, %d}].dss is %f\n", 10, 3, cr[{10, 3}].dss)

  for iCell = 0, nCells do
    var zt = cr[{iCell, nVertLevels}].zgrid
    for k = 0, nVertLevels do
      var z = 0.5 * (cr[{iCell, k}].zgrid + cr[{iCell, k+1}].zgrid)
      if (z > config_zd) then
        cr[{iCell, k}].dss = config_xnutr * cmath.pow(cmath.sin(0.5 * pii * (z-config_zd)/(zt-config_zd)), 2.0)
        cr[{iCell, k}].dss = cr[{iCell, k}].dss / cmath.pow(cr[{iCell, 0}].meshDensity, (0.25*dx_scale_power))
      end
    end
  end
end

task atm_couple_coef_3rd_order(config_coef_3rd_order : double,
                               cr : region(ispace(int2d), cell_fs),
                               er : region(ispace(int2d), edge_fs))
where reads writes (er, cr) do
  for iEdge = 0, nEdges do
    for i = 0, FIFTEEN do
      er[{iEdge, 0}].adv_coefs_3rd[i] = config_coef_3rd_order * er[{iEdge, 0}].adv_coefs_3rd[i]
    end
  end


  for iCell = 0, nCells do
    for vLevel = 0, nVertLevels + 1 do
      for i = 0, maxEdges do
        cr[{iCell, vLevel}].zb3_cell[i] = config_coef_3rd_order * cr[{iCell, vLevel}].zb3_cell[i]
      end
    end
  end
end

task atm_compute_solve_diagnostics(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), vr : region(ispace(int2d), vertex_fs), hollingsworth : bool)
where reads writes(vr, cr, er) do

  cio.printf("computing solve diagnostics\n")

  for iEdge = 0, nEdges do
    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

    for k = 0, nVertLevels do
      er[{iEdge, k}].h_edge = 0.5 * (cr[{cell1, k}].h + cr[{cell2, k}].h)
    end

    var efac = er[{iEdge, 0}].dcEdge * er[{iEdge, 0}].dvEdge

    for k = 0, nVertLevels do
      er[{iEdge, k}].ke_edge = efac * cmath.pow(er[{k, iEdge}].u, 2)
    end
  end

  for iVertex = 0, nVertices do
    for j = 0, nVertLevels do
      vr[{iVertex, j}].vorticity
    end
    for i = 0, vertexDegree do
      var iEdge = vr[{iVertex, 0}].edgesOnVertex[i]
      var s = vr[{iVertex, 0}].edgesOnVertexSign[i] * er[{iEdge, 0}].dcEdge

      for k = 0, nVertLevels do
        vr[{iVertex, k}].vorticity = vr[{iVertex, k}].vorticity + s * er[{iEdge, k}].u
      end
    end

    for k = 0, nVertLevels do
      vr[{iVertex, k}].vorticity = vr[{iVertex, k}].vorticity * vr[{iVertex, 0}].invAreaTriangle
    end
  end

  for iCell = 0, nCells do
    for j = 0, nVertLevels do
      cr[{iCell, j}].divergence = 0
    end
    for i = 1, cr[{iCell, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell, 0}].edgesOnCell[i]
      var s = cr[{iCell, 0}].edgesOnCellSign[i] * er[{iEdge, 0}].dvEdge

      for k = 0, nVertLevels do
        cr[{iCell, k}].divergence = cr[{iCell, k}].divergence + s + er[{iEdge, k}].u
      end
    end
    var r = cr[{iCell, 0}].invAreaCell
    for k = 0, nVertLevels do
      cr[{iCell, k}].divergence = cr[{iCell, k}].divergence * r
    end
  end

  for iCell = 0, nCells do
    for j = 1, nVertLevels do
      cr[{iCell, j}].ke = 0
    end
    for i = 0, cr[{iCell, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell, 0}].edgesOnCell[i]
      for k = 0, nVertLevels do
        cr[{iCell, k}].ke = cr[{iCell, k}].ke + 0.25 * er[{iEdge, k}].ke_edge
      end
    end

    for k = 0, nVertLevels do
      cr[{iCell, k}].ke = cr[{iCell, k}].ke * cr[{iCell,0}].invAreaCell
    end
  end

  if (hollingsworth) then
    for iVertex = 0, nVertices do
      var r = 0.25 * vr[{iVertex, 0}].invAreaTriangle
      for k = 0, nVertLevels do
        vr[{iVertex, k}].ke_vertex = (er[{vr[{iVertex, 0}].edgesOnVertex[0], k}].ke_edge + er[{vr[{iVertex, 0}].edgesOnVertex[1], k}].ke_edge + er[{vr[{iVertex, 0}].edgesOnVertex[2], k}].ke_edge)*r
      end
    end

    var ke_fact = 1.0 - 0.375

    for iCell = 0, nCells do
      for k = 0, nVertLevels do
        cr[{iCell, k}].ke = ke_fact * cr[{iCell, k}].ke
      end
    end

    for iCell = 0, nCells do
      var r = cr[{iCell, 0}].invAreaCell
      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
       var iVertex = cr[{iCell, 0}].verticesOnCell[i]
       var j = cr[{iCell, 0}].kiteForCell[i]

       for k = 0, nVertLevels do
        cr[{iCell, k}].ke = cr[{iCell, k}].ke + (1.0 - ke_fact)*vr[{iVertex, 0}].kiteAreasOnVertex[j] * vr[{iVertex, k}].ke_vertex * r
       end
      end
    end
  end

  var reconstruct_v = true

  ----Comment (Arjun): What to do with present(rk_step)?
  --if(present(rk_step)) then
  --  if(rk_step /= 3) reconstruct_v = .false.
  --end if

  if(reconstruct_v) then
    for iEdge = 0, nEdges do
      for j = 0, nVertLevels do
        er[{iEdge, j}].v = 0
      end
      for i = 1, er[{iEdge, 0}].nEdgesOnEdge do
        var eoe = er[{iEdge, 0}].edgesOnEdge_ECP[i]

        for k = 0, nVertLevels do
          er[{iEdge, k}].v = er[{iEdge, k}].v + er[{iEdge, 0}].weightsOnEdge[i] * er[{eoe, k}].u
        end
      end
    end
  end

  for iVertex = 0, nVertices do
    for k = 0, nVertLevels do
      vr[{iVertex, k}].pv_vertex = vr[{iVertex, 0}].fVertex + vr[{iVertex, k}].vorticity
    end
  end

  for iEdge = 0, nEdges do
    for k =0, nVertLevels do
      er[{iEdge, k}].pv_edge =  0.5 * (vr[{er[{iEdge, 0}].verticesOnEdge[0],k}].pv_vertex + vr[{er[{iEdge, 0}].verticesOnEdge[1],k}].pv_vertex)
    end
  end

    --SKIPPED: (config_apvm_upwinding > 0.0) then---
end

-- Comments:
-- This function contains nCellsSolve, moist_start, moist_end, and scalars,
-- which we are currently not sure how to translate 
task atm_compute_moist_coefficients(cr : region(ispace(int2d), cell_fs), 
                                    er : region(ispace(int2d), edge_fs))
where reads writes(cr, er) do 

  cio.printf("computing moist coefficients\n")

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].qtot = 0.0
    end
    for k = 0, nVertLevels do
      --TODO: What should we use instead of moist_start/moist_end?
      --for iq = moist_start, moist_end do
        --TODO: not sure how to translate: scalars(iq, k, iCell)
        --cr[{iCell, k}].qtot = cr[{iCell, k}].qtot + scalars(iq, k, iCell)
      --end
    end
  end

  for iCell = 0, nCells do
    for k = 1, nVertLevels do
      var qtotal = 0.5 * (cr[{iCell, k}].qtot + cr[{iCell, k - 1}].qtot)
      cr[{iCell, k}].cqw = 1.0 / (1.0 + qtotal)
      cio.printf("cr[{%d, %d}].cqw = %f", iCell, k, cr[{iCell, k}].cqw)
    end
  end

  for iEdge = 0, nEdges do
    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
    --if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then
      --for k = 0, nVertLevels do
        --var qtotal = 0.0
        --for iq = moist_start, moist_end do
          --qtotal = qtotal + 0.5 * ( scalars(iq, k, cell1) + scalars(iq, k, cell2) )
        --end
        --er[{iEdge, k}].cqu = 1.0 / (1.0 + qtotal)
      --end
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

task atm_compute_vert_imp_coefs(cr : region(ispace(int2d), cell_fs),
                                vert_r : region(ispace(int1d), vertical_fs),
                                dts : double)
where reads writes (cr, vert_r) do
      --  set coefficients
      var dtseps = .5 * dts * (1.0 + constants.config_epssm)
      var rgas = constants.rgas
      var rcv = rgas / (constants.cp - rgas)
      var c2 = constants.cp * rcv

      var qtotal : double
      var b_tri : double[nVertLevels]
      var c_tri : double[nVertLevels]


-- MGD bad to have all threads setting this variable?
      for k = 0, nVertLevels do
         vert_r[k].cofrz= dtseps * vert_r[k].rdzw
      end

      for iCell = 0, nCells do --  we only need to do cells we are solving for, not halo cells

--DIR$ IVDEP
         for k=1, nVertLevels do
            cr[{iCell, k}].cofwr = .5 * dtseps * constants.gravity * (vert_r[k].fzm * cr[{iCell, k}].zz + vert_r[k].fzp * cr[{iCell, k-1}].zz)
         end
         cr[{iCell, 0}].coftz = 0.0 --coftz(1,iCell) = 0.0
--DIR$ IVDEP
         for k=1, nVertLevels do
            cr[{iCell, k}].cofwz = dtseps * c2 * (vert_r[k].fzm * cr[{iCell, k}].zz + vert_r[k].fzp * cr[{iCell, k-1}].zz) * vert_r[k].rdzu * cr[{iCell, k}].cqw * (vert_r[k].fzm * cr[{iCell, k}].exner + vert_r[k].fzp * cr[{iCell, k-1}].exner)
            cr[{iCell, k}].coftz = dtseps * (vert_r[k].fzm * cr[{iCell, k}].theta_m + vert_r[k].fzp * cr[{iCell, k-1}].theta_m)
         end
         cr[{iCell, nVertLevels}].coftz = 0.0 -- coftz(nVertLevels+1,iCell)
--DIR$ IVDEP
         for k=0, nVertLevels do

            qtotal = cr[{iCell, k}].qtot

            cr[{iCell, k}].cofwt = .5 * dtseps * rcv * cr[{iCell, k}].zz * constants.gravity * cr[{iCell, k}].rho_base / ( 1.0 + qtotal) * cr[{iCell, k}].exner / ((cr[{iCell, k}].rtheta_base + cr[{iCell, k}].rtheta_p) * cr[{iCell, k}].exner_base)
         end

         cr[{iCell, 0}].a_tri = 0.0 --a_tri(1,iCell) = 0.  -- note, this value is never used
         b_tri[0] = 1.0    -- note, this value is never used
         c_tri[0] = 0.0    -- note, this value is never used
         cr[{iCell, 0}].gamma_tri = 0.0
         cr[{iCell, 0}].alpha_tri = 0.0  -- note, this value is never used

--DIR$ IVDEP
         for k=1, nVertLevels do --k=2,nVertLevels
            cr[{iCell,k}].a_tri = -1.0 * cr[{iCell, k}].cofwz * cr[{iCell, k-1}].coftz * vert_r[k-1].rdzw * cr[{iCell, k-1}].zz + cr[{iCell, k}].cofwr * vert_r[k-1].cofrz - cr[{iCell, k-1}].cofwt * cr[{iCell, k-1}].coftz * vert_r[k-1].rdzw

            b_tri[k] = 1.0 + cr[{iCell, k}].cofwz * (cr[{iCell, k}].coftz * vert_r[k].rdzw * cr[{iCell, k}].zz +  cr[{iCell, k}].coftz * vert_r[k-1].rdzw * cr[{iCell, k-1}].zz) -  cr[{iCell, k}].coftz * (cr[{iCell, k}].cofwt * vert_r[k].rdzw - cr[{iCell, k-1}].cofwt * vert_r[k-1].rdzw) + cr[{iCell, k}].cofwr * ((vert_r[k].cofrz- vert_r[k-1].cofrz))

            c_tri[k] =   -1.0 * cr[{iCell, k}].cofwz * cr[{iCell, k+1}].coftz * vert_r[k].rdzw * cr[{iCell, k}].zz - cr[{iCell, k}].cofwr * vert_r[k].cofrz+ cr[{iCell, k}].cofwt * cr[{iCell, k+1 }].coftz * vert_r[k].rdzw
         end
--MGD VECTOR DEPENDENCE
         for k=1, nVertLevels do -- k=2, nVertLevels
            cr[{iCell, k}].alpha_tri = 1.0/ (b_tri[k]-cr[{iCell, k}].a_tri * cr[{iCell, k-1}].gamma_tri)
            cr[{iCell, k}].gamma_tri = c_tri[k] * cr[{iCell, k}].alpha_tri
         end

      end -- loop over cells
end


task atm_compute_mesh_scaling(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), config_h_ScaleWithMesh : bool)
where reads writes (cr, er) do

  for iEdge = 0, nEdges do
    er[{iEdge, 0}].meshScalingDel2 = 1.0
    er[{iEdge, 0}].meshScalingDel4 = 1.0
  end

  if (config_h_ScaleWithMesh) then
    for iEdge = 0, nEdges do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
      er[{iEdge, 0}].meshScalingDel2 = 1.0 / cmath.pow((cr[{cell1, 0}].meshDensity + cr[{cell2, 0}].meshDensity)/2.0, 0.25)
      er[{iEdge, 0}].meshScalingDel4 = 1.0 / cmath.pow((cr[{cell1, 0}].meshDensity + cr[{cell2, 0}].meshDensity)/2.0, 0.75)
    end
  end

  for iCell = 0, nCells do
    cr[{iCell, 0}].meshScalingRegionalCell = 1.0
  end

  for iEdge = 0, nEdges do
    er[{iEdge, 0}].meshScalingRegionalEdge = 1.0
  end

  if (config_h_ScaleWithMesh) then
    for iEdge = 0, nEdges do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
      er[{iEdge, 0}].meshScalingRegionalEdge = 1.0 / cmath.pow((cr[{cell1, 0}].meshDensity + cr[{cell2, 0}].meshDensity)/2.0, 0.25)
    end

    for iCell = 0, nCells do
      cr[{iCell, 0}].meshScalingRegionalCell = 1.0/ cmath.pow(cr[{iCell, 0}].meshDensity, 0.25)
    end
  end
end

--Not sure how to translate: scalars(index_qv,k,iCell)
--sign(1.0_RKIND,flux) translated as cmath.copysign(1.0, flux)

task atm_init_coupled_diagnostics(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), vert_r : region(ispace(int1d), vertical_fs))
where reads writes (cr, er, vert_r) do
  cio.printf("initializing coupled diagnostics\n")
  var rgas = constants.rgas
  var rcv = rgas / (constants.cp - rgas)
  var p0 = 100000

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      --cr[{iCell, k}].theta_m = cr[{iCell, k}].theta * (1.0 + constants.rvord * scalars(index_qv,k,iCell))
      cr[{iCell, k}].rho_zz = cr[{iCell, k}].rho_zz / cr[{iCell, k}].zz
    end
  end

  for iEdge = 0, nEdges do
    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
    for k = 0, nVertLevels do
      er[{iEdge, k}].ru = 0.5 * er[{iEdge, k}].u * (cr[{cell1, k}].rho_zz + cr[{cell2, k}].rho_zz)
    end
  end

  for iCell = 0, nCells do
    cr[{iCell, 0}].rw = 0
    cr[{iCell, nVertLevels}].rw = 0
    for k = 1, nVertLevels do
      cr[{iCell, k}].rw = cr[{iCell, k}].w * (vert_r[k].fzp * cr[{iCell, k-1}].rho_zz + vert_r[k].fzm * cr[{iCell, k}].rho_zz) * (vert_r[k].fzp * cr[{iCell, k-1}].zz + vert_r[k].fzm * cr[{iCell, k}].zz)
    end
  end

  for iCell = 0, nCells do
    for i = 0, cr[{iCell, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell, 0}].edgesOnCell[i]
      for k = 1, nVertLevels do
        var flux = vert_r[k].fzm * er[{iEdge, k}].ru + vert_r[k].fzp * er[{iEdge, k-1}].ru
        cr[{iCell, k}].rw = cr[{iCell, k}].rw - cr[{iCell, 0}].edgesOnCellSign[i] * (cr[{iCell, k}].zb_cell[i] + cmath.copysign(1.0, flux) * cr[{iCell, k}].zb3_cell[i]) * flux * (vert_r[k].fzp * cr[{iCell, k-1}].zz + vert_r[k].fzm * cr[{iCell, k}].zz)
      end
    end
  end

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].rho_p = cr[{iCell, k}].rho_zz - cr[{iCell, k}].rho_base
    end
  end

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].rtheta_base = cr[{iCell, k}].theta_base * cr[{iCell, k}].rho_base
    end
  end

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].rtheta_p = cr[{iCell, k}].theta_m * cr[{iCell, k}].rho_p + cr[{iCell, k}].rho_base * (cr[{iCell, k}].theta_m - cr[{iCell, k}].theta_base)
    end
  end



  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].exner = cmath.pow(cr[{iCell, k}].zz * (rgas/p0) * (cr[{iCell, k}].rtheta_p + cr[{iCell, k}].rtheta_base), rcv)
      cr[{iCell, k}].exner_base = cmath.pow(cr[{iCell, k}].zz * (rgas/p0) * (cr[{iCell, k}].rtheta_base), rcv)
    end
  end

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].pressure_p = cr[{iCell, k}].zz * rgas * (cr[{iCell, k}].exner * cr[{iCell, k}].rtheta_p + cr[{iCell, k}].rtheta_base * (cr[{iCell, k}].exner - cr[{iCell, k}].exner_base))
      cr[{iCell, k}].pressure_base = cr[{iCell, k}].zz * rgas * cr[{iCell, k}].exner_base * cr[{iCell, k}].rtheta_base
      --cio.printf("zz at cell %d and %d is %f \n", iCell, k, cr[{iCell, k}].zz)
      --cio.printf("exner at cell %d and %d is %f \n", iCell, k, cr[{iCell, k}].exner)
      --cio.printf("rtheta_p at cell %d and %d is %f \n", iCell, k, cr[{iCell, k}].rtheta_p)
      --cio.printf("rtheta_base at cell %d and %d is %f \n", iCell, k, cr[{iCell, k}].rtheta_base)
      --cio.printf("exner_base at cell %d and %d is %f \n", iCell, k, cr[{iCell, k}].exner_base)
      --cio.printf("Pressure_p at cell %d and %d is %f \n", iCell, k, cr[{iCell, k}].pressure_p)
    end
  end

end

--Not sure how to translate: scalars(index_qv,k,iCell)
--__demand(__cuda)
task atm_compute_output_diagnostics(cr : region(ispace(int2d), cell_fs))
where reads writes (cr) do

  for iCell = 1, nCells do
    for k = 0, nVertLevels do
      --cr[{iCell, k}].theta = cr[{iCell, k}].theta_m / (1 + constants.rvord * scalars(index_qv,k,iCell))
      cr[{iCell, k}].rho = cr[{iCell, k}].rho_zz * cr[{iCell, k}].zz
      cr[{iCell, k}].pressure = cr[{iCell, k}].pressure_base + cr[{iCell, k}].pressure_p
    end
  end

end

task atm_rk_integration_setup(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs))
where reads writes (cr, er) do
  format.println("setting up rk integration")
  var edge_range = rect2d { int2d{0, 0}, int2d{nEdges - 1, nVertLevels - 1} }
  var cell_range = rect2d { int2d{0, 0}, int2d{nCells - 1, nVertLevels - 1} }

  for i in edge_range do
    er[i].ru_save = er[i].ru
    er[i].u_2 = er[i].u
  end

  for i in cell_range do
    cr[i].rw_save = cr[i].rw
    cr[i].rtheta_p_save = cr[i].rtheta_p
    cr[i].rho_p_save = cr[i].rho_p

    cr[i].w_2 = cr[i].w
    cr[i].theta_m_2 = cr[i].theta_m
    cr[i].rho_zz_2 = cr[i].rho_zz
    cr[i].rho_zz_old_split = cr[i].rho_zz
    --Not sure how to translate scalars
    --scalars_2(:,:,cellStart:cellEnd) = scalars_1(:,:,cellStart:cellEnd)
  end
end

__demand(__inline)
task flux4(q_im2 : double, q_im1 : double, q_i : double, q_ip1 : double, ua : double) : double
  return ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0
end

__demand(__inline)
task flux3(q_im2 : double, q_im1 : double, q_i : double, q_ip1 : double, ua : double, coef3 : double) : double
  return flux4(q_im2, q_im1, q_i, q_ip1, ua) 
         + coef3*cmath.fabs(ua)*((q_ip1-q_im2)-3.*(q_i-q_im1))/12.0
end

-- Comments for atm_compute_dyn_tend_work
-- A large number of variables appeared to be missing from the registry. I have inferred their types
-- and regions to the best of my ability, but may have gotten some of them wrong.
-- They are as follows:
--    cr: tend_rho_physics, dpdz, rb, rr_save, pp, delsq_divergence, ur_cell, vr_cell, delsq_w, 
--    tend_w_euler, theta_m_save, delsq_theta, tend_theta_euler, tend_rtheta_physics
--    er: tend_u_euler, delsq_u, tend_ru_physics
--    vr: delsq_vorticity
-- Variable c_s appears to be a renaming of constants.config_smagorinsky_coef.
-- Some pieces of code were inside an #ifdef CURVATURE. I have ignored the ifdefs and put the code in 
-- unconditionally.
-- Some config values should be configurable, rather than constants. These are currently in constants.rg.
-- I have also used previous conventions like removing "_RKIND", looping over all cells instead of 
-- cellSolveStart to cellStartEnd, using cmath.copysign and cmath.fabs for sign and abs, etc.

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
where reads writes (cr, er, vr, vert_r) do
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
      for iCell = 0, nCells do
        var d_diag : double[nVertLevels]
        var d_off_diag : double[nVertLevels]
        for k = 0, nVertLevels do
          d_diag[k] = 0.0
        end
        for k = 0, nVertLevels do
          d_off_diag[k] = 0.0
        end
        for iEdge = 0, cr[{iCell, 0}].nEdgesOnCell do
          for k = 0, nVertLevels do
            var e = cr[{iCell, 0}].edgesOnCell[iEdge]
            d_diag[k]     += cr[{iCell, 0}].defc_a[iEdge] * er[{e, k}].u
                              - cr[{iCell, 0}].defc_b[iEdge] * er[{e, k}].v
            d_off_diag[k] += cr[{iCell, 0}].defc_b[iEdge] * er[{e, k}].u 
                              + cr[{iCell, 0}].defc_a[iEdge] * er[{e, k}].v
          end
        end

        for k = 0, nVertLevels do
          -- here is the Smagorinsky formulation, 
          -- followed by imposition of an upper bound on the eddy viscosity

          -- Original: kdiff(k,iCell) = min((c_s * config_len_disp)**2 * sqrt(d_diag(k)**2 + d_off_diag(k)**2),(0.01*config_len_disp**2) * invDt)
          cr[{iCell, k}].kdiff = min(cmath.pow(c_s * constants.config_len_disp, 2.0) 
                                     * cmath.pow(cmath.pow(d_diag[k], 2.0) + cmath.pow(d_off_diag[k], 2.0), 0.5),
                                     (0.01 * cmath.pow(constants.config_len_disp, 2.0)) * invDt)
        end
      end

      h_mom_eddy_visc4 = constants.config_visc4_2dsmag * cmath.pow(constants.config_len_disp, 3.0) --0.05 * 120000.0**3
      h_theta_eddy_visc4 = h_mom_eddy_visc4

    elseif ([rawstring](config_horiz_mixing) == "2d_fixed") then
      fill(cr.kdiff, 0.0)
    end

    if (config_mpas_cam_coef > 0.0) then

      for iCell = 0, nCells do
        -- 2nd-order filter for top absorbing layer as in CAM-SE :  WCS 10 May 2017
        -- From MPAS-CAM V4.0 code, with addition to config-specified coefficient (V4.0_coef = 0.2; SE_coef = 1.0)
        cr[{iCell, nVertLevels-2}].kdiff = max(cr[{iCell, nVertLevels-2}].kdiff,
                                                        2.0833 * constants.config_len_disp * config_mpas_cam_coef)
        cr[{iCell, nVertLevels-1}].kdiff = max(cr[{iCell, nVertLevels-1}].kdiff,
                                                  2.0 * 2.0833 * constants.config_len_disp * config_mpas_cam_coef)
        cr[{iCell, nVertLevels  }].kdiff = max(cr[{iCell, nVertLevels  }].kdiff,
                                                  4.0 * 2.0833 * constants.config_len_disp * config_mpas_cam_coef)
      end
    end
  end

  -- tendency for density.
  -- accumulate total water here for later use in w tendency calculation.

  -- accumulate horizontal mass-flux

  fill(cr.h_divergence, 0.0)
  for iCell = 0, nCells do
    for i = 0, cr[{iCell, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell, 0}].edgesOnCell[i]
      var edge_sign = cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge
      for k = 0, nVertLevels do
        cr[{iCell, k}].h_divergence += edge_sign * er[{iEdge, k}].ru
      end
    end
  end

  -- compute horiontal mass-flux divergence, add vertical mass flux divergence to complete tend_rho

  for iCell = 0, nCells do
    var r = cr[{iCell, 0}].invAreaCell
    for k = 0, nVertLevels do
      cr[{iCell, k}].h_divergence *= r
    end
  end

  -- dp / dz and tend_rho
  -- only needed on first rk_step with pert variables defined a pert from time t
  if (rk_step == 0) then

    var rgas_cprcv = constants.rgas * constants.cp / constants.cv
    for iCell = 0, nCells do

      for k = 0, nVertLevels do
        --Not sure what tend_rho_physics, dpdz, rb, rr_save are.
        cr[{iCell, k}].tend_rho = -cr[{iCell, k}].h_divergence - vert_r[k].rdzw * (cr[{iCell, k + 1}].rw - cr[{iCell, k}].rw + cr[{iCell, k}].tend_rho_physics)
        --Original: dpdz(k,iCell) = -gravity*(rb(k,iCell)*(qtot(k,iCell)) + rr_save(k,iCell)*(1.+qtot(k,iCell)))
        cr[{iCell, k}].dpdz = -constants.gravity * (cr[{iCell, k}].rb * (cr[{iCell, k}].qtot) + cr[{iCell, k}].rr_save * (1.0 + cr[{iCell, k}].qtot))
      end
    end
  end



-------- BEGIN U SECTION --------

  -- Compute u (normal) velocity tendency for each edge (cell face)
  for iEdge = 0, nEdges do

    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

    -- horizontal pressure gradient 
    if (rk_step == 0) then
      for k = 0, nVertLevels do
        --Unable to find: tend_u_euler, pp, dpdz
        --It looks like tend_u_euler is from er and the others are from cr
        er[{iEdge, k}].tend_u_euler = -er[{iEdge, k}].cqu * ( (cr[{cell2, k}].pp - cr[{cell1, k}].pp) * er[{iEdge, 0}].invDcEdge 
                                      / (0.5 * (cr[{cell2, k}].zz + cr[{cell1, k}].zz)) 
                                      - 0.5 * er[{iEdge, k}].zxu * (cr[{cell1, k}].dpdz + cr[{cell2, k}].dpdz) )
      end

    end

    -- vertical transport of u
    var wduz : double[nVertLevels + 1]
    wduz[0] = 0.0

    var k = 1
    wduz[k] = 0.5 * (cr[{cell1, k}].rw + cr[{cell2, k}].rw) * (vert_r[k].fzm * er[{iEdge, k}].u + vert_r[k].fzp * er[{iEdge, k - 1}].u)
    for k = 2, nVertLevels - 1 do
      wduz[k] = flux3( er[{iEdge,k-2}].u, er[{iEdge,k-1}].u, er[{iEdge,k}].u, er[{iEdge,k+1}].u, 0.5*(cr[{cell1,k}].rw+cr[{cell2,k}].rw), 1.0 )
    end
    k = nVertLevels - 1
    wduz[k] = 0.5 * (cr[{cell1, k}].rw + cr[{cell2, k}].rw) * (vert_r[k].fzm * er[{iEdge, k}].u + vert_r[k].fzp * er[{iEdge, k - 1}].u)

    wduz[nVertLevels] = 0.0

    for k = 0, nVertLevels do
      er[{iEdge, k}].tend_u = -vert_r[k].rdzw * (wduz[k+1] - wduz[k])
    end

    -- Next, nonlinear Coriolis term (q) following Ringler et al JCP 2009

    var q : double[nVertLevels]
    for k = 0, nVertLevels do
      q[k] = 0.0
    end

    for j = 0, er[{iEdge, 0}].nEdgesOnEdge do
      var eoe = er[{iEdge, 0}].edgesOnEdge[j]
      for k = 0, nVertLevels do
        var workpv = 0.5 * (er[{iEdge, k}].pv_edge + er[{eoe, k}].pv_edge)
        -- the original definition of pv_edge had a factor of 1/density.  We have removed that factor
        -- given that it was not integral to any conservation property of the system
        q[k] += er[{iEdge, 0}].weightsOnEdge[j] * er[{eoe, k}].u * workpv
      end
    end

    for k = 0, nVertLevels do

      -- horizontal ke gradient and vorticity terms in the vector invariant formulation
      -- of the horizontal momentum equation
      er[{iEdge, k}].tend_u += er[{iEdge, k}].rho_edge 
                               * ( q[k] - (cr[{cell2, k}].ke - cr[{cell1, k}].ke) * er[{iEdge, 0}].invDcEdge )
                               - er[{iEdge, k}].u * 0.5 * (cr[{cell1, k}].h_divergence + cr[{cell2, k}].h_divergence)

      -- #ifdef CURVATURE
      -- curvature terms for the sphere
      er[{iEdge, k}].tend_u -= ( 2.0 * constants.omega * cmath.cos(er[{iEdge, 0}].angleEdge) 
                               * cmath.cos(er[{iEdge, 0}].latEdge) * er[{iEdge, k}].rho_edge 
                               * 0.25 * (cr[{cell1, k}].w + cr[{cell1, k + 1}].w 
                               + cr[{cell2, k}].w + cr[{cell2, k + 1}].w) )
                               - ( er[{iEdge, k}].u * 0.25 * (cr[{cell1, k}].w + cr[{cell1, k + 1}].w 
                               + cr[{cell2, k}].w + cr[{cell2, k + 1}].w) * er[{iEdge, k}].rho_edge 
                               * inv_r_earth )
      -- #endif
    end
  end

  -- horizontal mixing for u
  -- mixing terms are integrated using forward-Euler, so this tendency is only computed in the
  -- first Runge-Kutta substep and saved for use in later RK substeps 2 and 3.

  if (rk_step == 0) then

    -- del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).
    -- First, storage to hold the result from the first del^2 computation.
    fill(er.delsq_u, 0.0)

    for iEdge = 0, nEdges do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
      var vertex1 = er[{iEdge, 0}].verticesOnEdge[0]
      var vertex2 = er[{iEdge, 0}].verticesOnEdge[1]
      var r_dc = er[{iEdge, 0}].invDcEdge
      var r_dv = min(er[{iEdge, 0}].invDvEdge, 4 * r_dc)

      for k = 0, nVertLevels do

          -- Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
          --                    only valid for h_mom_eddy_visc4 == constant
        var u_diffusion = (cr[{cell2, k}].divergence - cr[{cell1, k}].divergence) * r_dc 
                          - (vr[{vertex2, k}].vorticity - vr[{vertex1, k}].vorticity) * r_dv
        er[{iEdge, k}].delsq_u += u_diffusion
        var kdiffu = 0.5 * (cr[{cell1, k}].kdiff + cr[{cell2, k}].kdiff)
        -- include 2nd-orer diffusion here 
        er[{iEdge, k}].tend_u_euler += er[{iEdge, k}].rho_edge * kdiffu * u_diffusion 
                                       * er[{iEdge, 0}].meshScalingDel2

      end
    end

    if (h_mom_eddy_visc4 > 0.0) then  -- 4th order mixing is active

      fill(vr.delsq_vorticity, 0.0)
      for iVertex = 0, nVertices do
        for i = 0, vertexDegree do
          var iEdge = vr[{iVertex, 0}].edgesOnVertex[i]
          var edge_sign = vr[{iVertex, 0}].invAreaTriangle * er[{iEdge, 0}].dcEdge * vr[{iVertex, 0}].edgesOnVertex_sign[i]
          for k = 0, nVertLevels do
              vr[{iVertex, k}].delsq_vorticity += edge_sign * er[{iEdge, k}].delsq_u
          end
        end
      end

      fill(cr.delsq_divergence, 0.0)
      for iCell = 0, nCells do
        var r = cr[{iCell, 0}].invAreaCell
        for i = 0, cr[{iCell, 0}].nEdgesOnCell do
          var iEdge = cr[{iCell, 0}].edgesOnCell[i]
          var edge_sign = r * er[{iEdge, 0}].dvEdge * cr[{iCell, 0}].edgesOnCell_sign[i]
          for k = 0, nVertLevels do
            cr[{iCell, k}].delsq_divergence += edge_sign * er[{iEdge, k}].delsq_u
          end
        end
      end

      for iEdge = 0, nEdges do
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
        var vertex1 = er[{iEdge, 0}].verticesOnEdge[0]
        var vertex2 = er[{iEdge, 0}].verticesOnEdge[1]

        var u_mix_scale = er[{iEdge, 0}].meshScalingDel4 * h_mom_eddy_visc4
        var r_dc = u_mix_scale * constants.config_del4u_div_factor * er[{iEdge, 0}].invDcEdge
        var r_dv = u_mix_scale * min(er[{iEdge, 0}].invDvEdge, 4 * er[{iEdge, 0}].invDcEdge)

        for k = 0, nVertLevels do

          -- Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
          --                    only valid for h_mom_eddy_visc4 == constant
          -- Here, we scale the diffusion on the divergence part a factor of config_del4u_div_factor 
          --    relative to the rotational part.  The stability constraint on the divergence component is much less
          --    stringent than the rotational part, and this flexibility may be useful.
          var u_diffusion =  er[{iEdge, k}].rho_edge * ( (cr[{cell2, k}].delsq_divergence - cr[{cell1, k}].delsq_divergence) * r_dc
                            - (vr[{vertex2, k}].delsq_vorticity - vr[{vertex1, k}].delsq_vorticity) * r_dv )
          er[{iEdge, k}].tend_u_euler -= u_diffusion
          
        end
      end
    
    end -- 4th order mixing is active 

    -- vertical mixing for u - 2nd order filter in physical (z) space
    if (v_mom_eddy_visc2 > 0.0) then

      if (config_mix_full) then  -- mix full state

        for iEdge = 0, nEdges do
          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
          for k = 1, nVertLevels - 1 do
            var z1 = 0.5 * (cr[{cell1, k - 1}].zgrid + cr[{cell2, k - 1}].zgrid)
            var z2 = 0.5 * (cr[{cell1, k    }].zgrid + cr[{cell2, k    }].zgrid)
            var z3 = 0.5 * (cr[{cell1, k + 1}].zgrid + cr[{cell2, k + 1}].zgrid)
            var z4 = 0.5 * (cr[{cell1, k + 2}].zgrid + cr[{cell2, k + 2}].zgrid)

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)

            er[{iEdge, k}].tend_u_euler += er[{iEdge, k}].rho_edge * v_mom_eddy_visc2
                                           * ( (er[{iEdge, k+1}].u - er[{iEdge, k  }].u) / (zp-z0)
                                             - (er[{iEdge, k  }].u - er[{iEdge, k-1}].u) / (z0-zm) )
                                           / (0.5 * (zp - zm))
          end
        end

      else  -- idealized cases where we mix on the perturbation from the initial 1-D state

        var u_mix : double[nVertLevels]

        for iEdge = 0, nEdges do
          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

          for k = 0, nVertLevels do
            u_mix[k] = er[{iEdge, k}].u - vert_r[k].u_init * cmath.cos(er[{iEdge, 0}].angleEdge)
                       - vert_r[k].v_init * cmath.sin(er[{iEdge, 0}].angleEdge)
          end

          for k = 1, nVertLevels - 1 do
            var z1 = 0.5 * (cr[{cell1, k - 1}].zgrid + cr[{cell2, k - 1}].zgrid)
            var z2 = 0.5 * (cr[{cell1, k    }].zgrid + cr[{cell2, k    }].zgrid)
            var z3 = 0.5 * (cr[{cell1, k + 1}].zgrid + cr[{cell2, k + 1}].zgrid)
            var z4 = 0.5 * (cr[{cell1, k + 2}].zgrid + cr[{cell2, k + 2}].zgrid)

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)

            er[{iEdge, k}].tend_u_euler += er[{iEdge, k}].rho_edge * v_mom_eddy_visc2
                                           * ( (u_mix[k+1] - u_mix[k  ]) / (zp-z0)
                                             - (u_mix[k  ] - u_mix[k-1]) / (z0-zm) )
                                           / (0.5 * (zp - zm))
          end
        end
      end -- mix perturbation state
    end -- vertical mixing of horizontal momentum
  end -- (rk_step 1 test for computing mixing terms)

  --  add in mixing and physics tendency for u

  --  Rayleigh damping on u
  if (config_rayleigh_damp_u) then
    var rayleigh_damp_coef : double[nVertLevels]
    var rayleigh_coef_inverse = 1.0 / ( [double](constants.config_number_rayleigh_damp_u_levels)
                                        * (constants.config_rayleigh_damp_u_timescale_days * seconds_per_day) )

    for k = nVertLevels - constants.config_number_rayleigh_damp_u_levels + 1, nVertLevels do
      rayleigh_damp_coef[k] = [double](k - (nVertLevels - constants.config_number_rayleigh_damp_u_levels))
                              * rayleigh_coef_inverse
    end

    for iEdge = 0, nEdges do
      for k = nVertLevels - constants.config_number_rayleigh_damp_u_levels + 1, nVertLevels do
        er[{iEdge, k}].tend_u -= er[{iEdge, k}].rho_edge * er[{iEdge, k}].u * rayleigh_damp_coef[k]
      end
    end
  end

  for iEdge = 0, nEdges do
    for k = 0, nVertLevels do
      er[{iEdge, k}].tend_u += er[{iEdge, k}].tend_u_euler + er[{iEdge, k}].tend_ru_physics
    end
  end



-------- BEGIN W SECTION ---------

  --  horizontal advection for w
  var ru_edge_w : double[nVertLevels]

  for iCell = 0, nCells do
    for k = 0, nVertLevels + 1 do
      cr[{iCell, k}].tend_w = 0.0
    end
    for i = 0, cr[{iCell, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell, 0}].edgesOnCell[i]
      var edge_sign = cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge * 0.5

      for k = 1, nVertLevels do
        ru_edge_w[k] = vert_r[k].fzm * er[{iEdge, k}].ru + vert_r[k].fzp * er[{iEdge, k - 1}].ru
      end

      var flux_arr : double[nVertLevels]
      for k = 0, nVertLevels do
        flux_arr[k] = 0.0
      end

      -- flux_arr stores the value of w at the cell edge used in the horizontal transport

      for j = 0, er[{iEdge, 0}].nAdvCellsForEdge do
        var iAdvCell = er[{iEdge, 0}].advCellsForEdge[j]
        for k = 1, nVertLevels do
          var scalar_weight = er[{iEdge, 0}].adv_coefs[j] + cmath.copysign(1.0, ru_edge_w[k]) * er[{iEdge, 0}].adv_coefs_3rd[j]
          flux_arr[k] += scalar_weight * cr[{iAdvCell, k}].w
        end
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].tend_w -= cr[{iCell, 0}].edgesOnCell_sign[i] * ru_edge_w[k] * flux_arr[k]
      end
    end
  end

  -- #ifdef CURVATURE
  for iCell = 0, nCells do
    for k = 1, nVertLevels do
      cr[{iCell, k}].tend_w += (cr[{iCell, k}].rho_zz * vert_r[k].fzm 
                                + cr[{iCell, k-1}].rho_zz * vert_r[k].fzp)
                               * ( cmath.pow(vert_r[k].fzm * cr[{iCell, k}].ur_cell 
                                             + vert_r[k].fzp * cr[{iCell, k-1}].ur_cell, 2.0)
                                  + cmath.pow(vert_r[k].fzm * cr[{iCell, k}].vr_cell
                                             + vert_r[k].fzp * cr[{iCell, k-1}].vr_cell, 2.0) ) / r_earth
                               + 2.0 * constants.omega * cmath.cos(cr[{iCell, 0}].latCell)
                               * (vert_r[k].fzm * cr[{iCell, k}].ur_cell 
                                  + vert_r[k].fzp * cr[{iCell, k-1}].ur_cell)
                               * (cr[{iCell, k}].rho_zz * vert_r[k].fzm 
                                  + cr[{iCell, k-1}].rho_zz * vert_r[k].fzp)
    end
  end
  -- #endif

  -- horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),
  -- but here we can also code in hyperdiffusion if we wish (2nd order at present)

  if (rk_step == 0) then

    -- del^4 horizontal filter.  We compute this as del^2 ( del^2 (u) ).

    -- First, storage to hold the result from the first del^2 computation.
    --  we copied code from the theta mixing, hence the theta* names.

    fill(cr.delsq_w, 0.0)

    for iCell = 0, nCells do
      for k = 0, nVertLevels + 1 do
        cr[{iCell, k}].tend_w_euler = 0.0
      end
      var r_areaCell = cr[{iCell, 0}].invAreaCell
      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
          var iEdge = cr[{iCell, 0}].edgesOnCell[i]

          var edge_sign = 0.5 * r_areaCell * cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge 
                          * er[{iEdge, 0}].invDcEdge

          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

        for k = 1, nVertLevels do
          var w_turb_flux = edge_sign * (er[{iEdge, k}].rho_edge + er[{iEdge, k-1}].rho_edge)
                            * (cr[{cell2, k}].w - cr[{cell1, k}].w)
          cr[{iCell, k}].delsq_w += w_turb_flux
          w_turb_flux *= er[{iEdge, 0}].meshScalingDel2 * 0.25
                          * (cr[{cell1, k}].kdiff + cr[{cell2, k}].kdiff 
                            + cr[{cell1, k-1}].kdiff + cr[{cell2, k-1}].kdiff)
          cr[{iCell, k}].tend_w_euler += w_turb_flux
        end
      end
    end

    if (h_mom_eddy_visc4 > 0.0) then  -- 4th order mixing is active

      for iCell = 0, nCells do
        var r_areaCell = h_mom_eddy_visc4 * cr[{iCell, 0}].invAreaCell
        for i = 0, cr[{iCell, 0}].nEdgesOnCell do
          var iEdge = cr[{iCell, 0}].edgesOnCell[i]
          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

          var edge_sign = er[{iEdge, 0}].meshScalingDel4 * r_areaCell * er[{iEdge, 0}].dvEdge 
                          * cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].invDcEdge

          for k = 1, nVertLevels do
            cr[{iCell, k}].tend_w_euler -= edge_sign * (cr[{cell2, k}].delsq_w - cr[{cell1, k}].delsq_w)
          end
        end
      end
    end -- 4th order mixing is active 
  end -- horizontal mixing for w computed in first rk_step

  --  vertical advection, pressure gradient and buoyancy for w
  for iCell = 0, nCells do

    var wdwz : double[nVertLevels + 1]
    wdwz[0] = 0.0
    var k = 1
    wdwz[k] = 0.25 * (cr[{iCell, k}].rw + cr[{iCell, k-1}].rw)*(cr[{iCell, k}].w + cr[{iCell, k-1}].w)
    for k = 2, nVertLevels - 1 do
      wdwz[k] = flux3( cr[{iCell, k-2}].w, cr[{iCell, k-1}].w, cr[{iCell, k}].w, cr[{iCell, k+1}].w, 
                       0.5 * (cr[{iCell, k}].rw + cr[{iCell, k-1}].rw), 1.0 )
    end
    k = nVertLevels - 1
    wdwz[k] = 0.25 * (cr[{iCell, k}].rw + cr[{iCell, k-1}].rw) * (cr[{iCell, k}].w + cr[{iCell, k-1}].rw)
    wdwz[nVertLevels] = 0.0

    -- Note: next we are also dividing through by the cell area after the horizontal flux divergence
    for k = 1, nVertLevels do
      cr[{iCell, k}].tend_w *= cr[{iCell, 0}].invAreaCell - vert_r[k].rdzu * (wdwz[k+1] - wdwz[k])
    end

    if (rk_step == 0) then
      for k = 1, nVertLevels do
        cr[{iCell, k}].tend_w_euler -= cr[{iCell, k}].cqw * ( vert_r[k].rdzu * 
                                       (cr[{iCell, k}].pp - cr[{iCell, k-1}].pp)
                                       - (vert_r[k].fzm * cr[{iCell, k}].dpdz + vert_r[k].fzp * cr[{iCell, k-1}].dpdz) )  -- dpdz is the buoyancy term here.
      end
    end
  end

  if (rk_step == 0) then
    if (v_mom_eddy_visc2 > 0.0) then
      for iCell = 0, nCells do
        for k = 1, nVertLevels do
          cr[{iCell, k}].tend_w_euler += v_mom_eddy_visc2 * (cr[{iCell, k}].rho_zz + cr[{iCell, k-1}].rho_zz)
                                         * 0.5 * ( (cr[{iCell, k+1}].w - cr[{iCell, k}].w) * vert_r[k].rdzw
                                         - (cr[{iCell, k}].w - cr[{iCell, k-1}].w) * vert_r[k-1].rdzw ) 
                                         * vert_r[k].rdzu
        end
      end
    end
  end -- mixing term computed first rk_step

  -- add in mixing terms for w
  for iCell = 0, nCells do
    for k = 1, nVertLevels do
      cr[{iCell, k}].tend_w += cr[{iCell, k}].tend_w_euler
    end
  end



-------- BEGIN THETA SECTION --------

  -- horizontal advection for theta
  fill(cr.tend_theta, 0.0)
  for iCell = 0, nCells do
    for i = 0, cr[{iCell, 0}].nEdgesOnCell do
      var iEdge = cr[{iCell, 0}].edgesOnCell[i]

      var flux_arr : double[nVertLevels]
      for k = 0, nVertLevels do
        flux_arr[k] = 0.0
      end

      for j = 0, er[{iEdge, 0}].nAdvCellsForEdge do
        var iAdvCell = er[{iEdge, 0}].advCellsForEdge[j]
        for k = 0, nVertLevels do
          var scalar_weight = er[{iEdge, 0}].adv_coefs[j] + cmath.copysign(1.0, er[{iEdge, k}].ru) 
                              * er[{iEdge, 0}].adv_coefs_3rd[j]
          flux_arr[k] += scalar_weight * cr[{iAdvCell, k}].theta_m
        end
      end

      for k = 0, nVertLevels do
        cr[{iCell, k}].tend_theta -= cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, k}].ru * flux_arr[k]
      end
    end
  end

  -- addition to pick up perturbation flux for rtheta_pp equation
  if (rk_step > 0) then
    for iCell = 0, nCells do
      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell, 0}].edgesOnCell[i]
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

        for k = 0, nVertLevels do
          var flux = cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge
                     * (er[{iEdge, k}].ru_save - er[{iEdge, k}].ru) * 0.5
                     * (cr[{cell2, k}].theta_m_save + cr[{cell1, k}].theta_m_save)
          cr[{iCell, k}].tend_theta -= flux -- division by areaCell picked up down below
        end
      end
    end
  end

  -- horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
  -- but here we can also code in hyperdiffusion if we wish (2nd order at present)
  if (rk_step == 0) then
    fill(cr.delsq_theta, 0.0)
    fill(cr.tend_theta_euler, 0.0)
    for iCell = 0, nCells do
      var r_areaCell = cr[{iCell, 0}].invAreaCell
      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell, 0}].edgesOnCell[i]
        var edge_sign = r_areaCell * cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].dvEdge * er[{iEdge, 0}].invDcEdge
        var pr_scale = prandtl_inv * er[{iEdge, 0}].meshScalingDel2
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

        for k = 0, nVertLevels do
          -- we are computing the Smagorinsky filter at more points than needed here so as to pick up the delsq_theta for 4th order filter below
          var theta_turb_flux = edge_sign * (cr[{cell2, k}].theta_m - cr[{cell1, k}].theta_m) * er[{iEdge, k}].rho_edge
          cr[{iCell, k}].delsq_theta += theta_turb_flux
          theta_turb_flux *= 0.5 * (cr[{cell1, k}].kdiff + cr[{cell2, k}].kdiff) * pr_scale
          cr[{iCell, k}].tend_theta_euler += theta_turb_flux
        end
      end
    end

    if (h_theta_eddy_visc4 > 0.0) then  -- 4th order mixing is active

      for iCell = 0, nCells do
        var r_areaCell = h_theta_eddy_visc4 * prandtl_inv * cr[{iCell, 0}].invAreaCell
        for i = 0, cr[{iCell, 0}].nEdgesOnCell do

          var iEdge = cr[{iCell, 0}].edgesOnCell[i]
          var edge_sign = er[{iEdge, 0}].meshScalingDel4 * r_areaCell * er[{iEdge, 0}].dvEdge 
                          * cr[{iCell, 0}].edgesOnCell_sign[i] * er[{iEdge, 0}].invDcEdge

          var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
          var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

          for k = 0, nVertLevels do
            cr[{iCell, k}].tend_theta_euler -= edge_sign * (cr[{cell2, k}].delsq_theta - cr[{cell1, k}].delsq_theta)
          end
        end
      end
    end -- 4th order mixing is active 
  end -- theta mixing calculated first rk_step

  -- vertical advection plus diabatic term
  -- Note: we are also dividing through by the cell area after the horizontal flux divergence

  var wdtz : double[nVertLevels + 1]
  for iCell = 0, nCells do

    wdtz[0] = 0.0
    var k = 1
    wdtz[k] = cr[{iCell, k}].rw * (vert_r[k].fzm * cr[{iCell, k}].theta_m + vert_r[k].fzp * cr[{iCell, k-1}].theta_m) 
              + ( (cr[{iCell, k}].rw_save - cr[{iCell, k}].rw) 
              * (vert_r[k].fzm * cr[{iCell, k}].theta_m_save + vert_r[k].fzp * cr[{iCell, k-1}].theta_m_save) )
    for k = 2, nVertLevels - 1 do
      wdtz[k] = flux3( cr[{iCell, k-2}].theta_m, cr[{iCell, k-1}].theta_m, cr[{iCell, k}].theta_m, 
                       cr[{iCell, k+1}].theta_m, cr[{iCell, k}].rw, constants.config_coef_3rd_order )
      wdtz[k] += (cr[{iCell, k}].rw_save - cr[{iCell, k}].rw)
                 * (vert_r[k].fzm * cr[{iCell, k}].theta_m_save + vert_r[k].fzp * cr[{iCell, k-1}].theta_m_save)  -- rtheta_pp redefinition
    end
    k = nVertLevels - 1
    wdtz[k] =  cr[{iCell, k}].rw_save * (vert_r[k].fzm * cr[{iCell, k}].theta_m_save 
                             + vert_r[k].fzp * cr[{iCell, k-1}].theta_m_save)  -- rtheta_pp redefinition
    wdtz[nVertLevels] = 0.0

    for k = 0, nVertLevels do
      cr[{iCell, k}].tend_theta *= cr[{iCell, 0}].invAreaCell - vert_r[k].rdzw * (wdtz[k+1] - wdtz[k])
      cr[{iCell, k}].tend_rtheta_adv = cr[{iCell, k}].tend_theta -- this is for the Tiedke scheme
      cr[{iCell, k}].rthdynten = cr[{iCell, k}].tend_theta / cr[{iCell, k}].rho_zz -- this is for the Grell-Freitas scheme
      cr[{iCell, k}].tend_theta += cr[{iCell, k}].rho_zz * cr[{iCell, k}].rt_diabatic_tend
    end
  end

  --  vertical mixing for theta - 2nd order
  if (rk_step == 0) then

    if (v_theta_eddy_visc2 > 0.0) then  -- vertical mixing for theta_m

      if (config_mix_full) then

        for iCell = 0, nCells do
          for k = 1, nVertLevels - 1 do
            var z1 = cr[{iCell, k - 1}].zgrid
            var z2 = cr[{iCell, k    }].zgrid
            var z3 = cr[{iCell, k + 1}].zgrid
            var z4 = cr[{iCell, k + 2}].zgrid

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)

            cr[{iCell, k}].tend_theta_euler += v_theta_eddy_visc2 * prandtl_inv * cr[{iCell, k}].rho_zz
                                               * ( (cr[{iCell, k+1}].theta_m - cr[{iCell, k}].theta_m) / (zp-z0)
                                               - (cr[{iCell, k}].theta_m-cr[{iCell, k-1}].theta_m) / (z0-zm) )
                                               / (0.5*(zp-zm))
          end
        end

      else  -- idealized cases where we mix on the perturbation from the initial 1-D state
        for iCell = 0, nCells do
          for k = 1, nVertLevels - 1 do
            var z1 = cr[{iCell, k - 1}].zgrid
            var z2 = cr[{iCell, k    }].zgrid
            var z3 = cr[{iCell, k + 1}].zgrid
            var z4 = cr[{iCell, k + 2}].zgrid

            var zm = 0.5 * (z1 + z2)
            var z0 = 0.5 * (z2 + z3)
            var zp = 0.5 * (z3 + z4)
            cr[{iCell, k}].tend_theta_euler += v_theta_eddy_visc2 * prandtl_inv * cr[{iCell, k}].rho_zz
                                              * ( ((cr[{iCell, k+1}].theta_m - cr[{iCell, k+1}].t_init) 
                                                    - (cr[{iCell, k}].theta_m - cr[{iCell, k}].t_init)) / (zp-z0)
                                              - ( (cr[{iCell, k}].theta_m - cr[{iCell, k}].t_init)
                                                  - (cr[{iCell, k-1}].theta_m - cr[{iCell, k-1}].t_init)) / (z0-zm) )
                                              / (0.5*(zp-zm))
          end
        end
      end
    end
  end -- compute vertical theta mixing on first rk_step

  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].tend_theta += cr[{iCell, k}].tend_theta_euler + cr[{iCell, k}].tend_rtheta_physics
    end
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
where reads writes (cr, er, vr, vert_r) do
  cio.printf("computing dynamic tendencies\n")
  atm_compute_dyn_tend_work(cr, er, vr, vert_r, rk_step, dt, config_horiz_mixing, config_mpas_cam_coef, config_mix_full, config_rayleigh_damp_u)
end

task atm_set_smlstep_pert_variables_work(cr : region(ispace(int2d), cell_fs),
                                         er : region(ispace(int2d), edge_fs),
                                         vert_r : region(ispace(int1d), vertical_fs))
where reads writes (cr, er, vert_r) do

  for iCell = 0, nCells do
    if (cr[{iCell, 0}].bdyMaskCell <= nRelaxZone) then
      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell, 0}].edgesOnCell[i]
        for k = 1, nVertLevels do
            var flux = cr[{iCell, 0}].edgesOnCell_sign[i] * (vert_r[k].fzm * er[{iEdge, k}].u_tend + vert_r[k].fzp * er[{iEdge, k - 1}].u_tend)
            --sign function copied as cmath.copysign, _RKIND removed as done previously
            cr[{iCell, k}].w_tend = cr[{iCell, k}].w_tend - (cr[{iCell, k}].zb_cell[i] + cmath.copysign(1.0, er[{iEdge, k}].u_tend) * cr[{iCell, k}].zb3_cell[i]) * flux
        end
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].w_tend = ( vert_r[k].fzm * cr[{iCell, k}].zz + vert_r[k].fzp * cr[{iCell, k - 1}].zz ) * cr[{iCell, k}].w_tend
      end
    end
  end
end

task atm_set_smlstep_pert_variables(cr : region(ispace(int2d), cell_fs),
                                    er : region(ispace(int2d), edge_fs),
                                    vert_r : region(ispace(int1d), vertical_fs))
where reads writes (cr, er, vert_r) do
  cio.printf("set small step vars\n")
  atm_set_smlstep_pert_variables_work(cr, er, vert_r)
end

--Comments for atm_advance_acoustic_step_work
--dts is passed in as double now - in code, passed in as rk_sub_timestep(rk_step)
--1.0_RKIND translated as just 1.0
--Tendency variables: tend_rw, tend_rt, tend_rho (added to CR), tend_ru (added to ER). Not in registry
--Other variables not in registry: rs, ts: added to CR as part of vertical grid, ru_Avg (added to ER)

task atm_advance_acoustic_step_work(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), vert_r : region(ispace(int1d), vertical_fs), dts : double, small_step : int) -- nCellsSolve : int)
where reads writes (er, cr, vert_r) do
  var epssm = constants.config_epssm
  var rgas = constants.rgas
  var rcv = rgas / (constants.cp - rgas)
  var c2 = constants.cp * rcv
  var resm = (1.0 - epssm) / (1.0 + epssm)
  var rdts = 1.0 / dts

  var rs : double[nVertLevels]
  var ts : double[nVertLevels]

  cio.printf("advancing acoustic step\n")

  if (small_step == 1) then
    for iEdge = 0, nEdges do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

      --if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

        --for k = 0, nVertLevels do
          --var pgrad = ((cr[{cell2, k}].rtheta_pp - cr[{cell1, k}].rtheta_pp) * er[{iEdge, 0}].invDcEdge) / (0.5 * (cr[{cell2, k}].zz +cr[{cell1, k}].zz))
          --pgrad = er[{iEdge, k}].cqu * 0.5 * c2 * (cr[{cell1, k}].exner + cr[{cell2, k}].exner) * pgrad
          --pgrad = pgrad + 0.5 * er[{iEdge, k}].zxu * constants.gravity * (cr[{cell1, k}].rho_pp + cr[{cell2, k}].rho_pp)
          --er[{iEdge, k}].ru_p = er[{iEdge, k}].ru_p + dts * (er[{iEdge, k}].tend_ru - (1.0 - er[{iEdge, 0}].specZoneMaskEdge) * pgrad)  --NEEDS FIXING
        --end
        --for k = 0, nVertLevels do
          --er[{iEdge, k}].ruAvg = er[{iEdge, k}].ruAvg + er[{iEdge, k}].ru_p
        --end
      --end
    end

  else
    for iEdge = 0, nEdges do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

      --if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

        --for k = 0, nVertLevels do
          --er[{iEdge, k}].ru_p = dts * er[{iEdge, k}].tend_ru
        --end

        --for k = 0, nVertLevels do
          --er[{iEdge, k}].ruAvg = er[{iEdge, k}].ru_p
        --end

      --end
    end
  end

  if (small_step == 1) then
    for iCell = 0, nCells do
      for j = 0, nVertLevels do
        cr[{iCell, j}].rtheta_pp_old = 0
      end
    end
  else
    for iCell = 0, nCells do
      for j = 0, nVertLevels do
        cr[{iCell, j}].rtheta_pp_old = cr[{iCell, j}].rtheta_pp
      end
    end
  end

  for iCell = 0, nCells do
    if(small_step == 1) then
      for j = 0, nVertLevels do
        cr[{iCell, j}].wwAvg = 0
        cr[{iCell, j}].rho_pp = 0
        cr[{iCell, j}].rtheta_pp = 0
        cr[{iCell, j}].rw_p = 0
      end
      cr[{iCell, nVertLevels}].wwAvg = 0
      cr[{iCell, nVertLevels}].rw_p = 0
    end

    if(cr[{iCell, 0}].specZoneMaskCell == 0.0) then
      for i = 0, nVertLevels do
        ts[i] = 0
        rs[i] = 0
      end

      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell, 0}].edgesOnCell[i]    --edgesOnCell(i,iCell)
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]    --cell2 = cellsOnEdge(2,iEdge)

        for k = 0, nVertLevels do
          var flux = cr[{iCell, 0}].edgesOnCellSign[i] * dts * er[{iEdge, 0}].dvEdge * er[{iEdge, k}].ru_p * cr[{iCell, 0}].invAreaCell
          rs[k] = rs[k] - flux
          ts[k] = ts[k] - flux * 0.5 * (cr[{cell2, k}].theta_m + cr[{cell1, k}].theta_m)
        end
      end

      for k = 0, nVertLevels do
         rs[k] = cr[{iCell, k}].rho_pp + dts * cr[{iCell, k}].tend_rho + rs[k] - vert_r[k].cofrz* resm * (cr[{iCell, k+1}].rw_p - cr[{iCell, k}].rw_p)
         ts[k] = cr[{iCell, k}].rtheta_pp + dts * cr[{iCell, k}].tend_rt + ts[k] - resm * vert_r[k].rdzw * (cr[{iCell, k+1}].coftz * cr[{iCell, k+1}].rw_p  - cr[{iCell, k}].coftz * cr[{iCell, k}].rw_p )
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].wwAvg = cr[{iCell, k}].wwAvg + 0.5 * (1.0 - epssm) * cr[{iCell, k}].rw_p
      end

      for k = 1, nVertLevels do
         cr[{iCell, k}].rw_p = cr[{iCell, k}].rw_p +  dts * cr[{iCell, k}].tend_rw - cr[{iCell, k}].cofwz * ((cr[{iCell, k}].zz * ts[k] - cr[{iCell, k-1}].zz * ts[k-1]) + resm * (cr[{iCell, k}].zz * cr[{iCell, k}].rtheta_pp - cr[{iCell, k-1}].zz * cr[{iCell, k-1}].rtheta_pp)) - cr[{iCell, k}].cofwr * ((rs[k] + rs[k-1]) + resm * (cr[{iCell, k}].rho_pp + cr[{iCell, k-1}].rho_pp))  + cr[{iCell, k}].cofwt * (ts[k] + resm * cr[{iCell, k}].rtheta_pp) + cr[{iCell, k-1}].cofwt * (
         ts[k-1] +resm * cr[{iCell, k-1}].rtheta_pp)
      end

      for k = 1, nVertLevels do
         cr[{iCell, k}].rw_p = (cr[{iCell, k}].rw_p - cr[{iCell, k}].a_tri * cr[{iCell, k-1}].rw_p) * cr[{iCell, k}].alpha_tri
      end

      for k = nVertLevels, 1, -1 do
        cr[{iCell, k}].rw_p = cr[{iCell, k}].rw_p - cr[{iCell, k}].gamma_tri * cr[{iCell, k+1}].rw_p
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].rw_p = (cr[{iCell, k}].rw_p + (cr[{iCell, k}].rw_save - cr[{iCell, k}].rw) - dts * cr[{iCell, k}].dss * (vert_r[k].fzm * cr[{iCell, k}].zz + vert_r[k].fzp * cr[{iCell, k-1}].zz)*(vert_r[k].fzm * cr[{iCell, k}].rho_zz + vert_r[k].fzp * cr[{iCell, k-1}].rho_zz) * cr[{iCell, k}].w)/(1.0 + dts * cr[{iCell, k}].dss)  - (cr[{iCell, k}].rw_save - cr[{iCell, k}].rw)
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].wwAvg = cr[{iCell, k}].wwAvg + 0.5 * (1.0 + epssm) * cr[{iCell, k}].rw_p
      end

      for k=0, nVertLevels do
         cr[{iCell, k}].rho_pp  = rs[k] - vert_r[k].cofrz*(cr[{iCell, k+1}].rw_p - cr[{iCell, k}].rw_p)
         cr[{iCell, k}].rtheta_pp = ts[k]  - vert_r[k].rdzw * (cr[{iCell, k+1}].coftz * cr[{iCell, k+1}].rw_p - cr[{iCell, k}].coftz * cr[{iCell, k}].rw_p)
      end

    else
      for k=0, nVertLevels do
         cr[{iCell, k}].rho_pp  = cr[{iCell, k}].rho_pp + dts * cr[{iCell, k}].tend_rho
         cr[{iCell, k}].rtheta_pp = cr[{iCell, k}].rtheta_pp + dts * cr[{iCell, k}].tend_rt
         cr[{iCell, k}].rw_p = cr[{iCell, k}].rw_p + dts * cr[{iCell, k}].tend_rw
         cr[{iCell, k}].wwAvg = cr[{iCell, k}].wwAvg + 0.5 * (1.0+epssm) * cr[{iCell, k}].rw_p
      end
    end
  end
end

task atm_advance_acoustic_step(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), vert_r : region(ispace(int1d), vertical_fs), dts : double, small_step : int) -- nCellsSolve : int)
where reads writes (er, cr, vert_r) do
  cio.printf("advancing acoustic step\n")
  atm_advance_acoustic_step_work(cr, er, vert_r, dts, small_step)
end

-- Comments:
-- dts is passed in as double - in code, passed in as rk_sub_timestep(rk_step)
-- 1.0_RKIND and 2.0_RKIND translated as 1.0 and 2.0
-- This function also contains nCellsSolve, which has not been resolved yet
task atm_divergence_damping_3d(cr : region(ispace(int2d), cell_fs),
                               er: region(ispace(int2d), edge_fs),
                               dts : double)
where reads writes (cr, er) do
  cio.printf("update horizontal momentum\n")

  var smdiv = constants.config_smdiv
  var rdts = 1.0 / dts
  var coef_divdamp = 2.0 * smdiv * constants.config_len_disp * rdts

  for iEdge = 0, nEdges do

    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

    -- update edges for block-owned cells
    -- if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

      for k = 0, nVertLevels do

        -- scaled 3d divergence damping
        var divCell1 = -(cr[{cell1, k}].rtheta_pp - cr[{cell1, k}].rtheta_pp_old)
        var divCell2 = -(cr[{cell2, k}].rtheta_pp - cr[{cell2, k}].rtheta_pp_old)
        er[{iEdge, k}].ru_p = er[{iEdge, k}].ru_p + coef_divdamp * (divCell2 - divCell1) * 
                              (1.0 - er[{iEdge, 0}].specZoneMaskEdge) 
                              / (cr[{cell1, k}].theta_m + cr[{cell2, k}].theta_m)
      end
    -- end
  end
end

task atm_recover_large_step_variables()
  cio.printf("recovering large step vars\n")
end

task atm_rk_dynamics_substep_finish()
  cio.printf("finishing substep\n")
end

--__demand(__cuda)
task atm_core_init(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), vr : region(ispace(int2d), vertex_fs), vert_r : region(ispace(int1d), vertical_fs))
where reads writes (cr, er, vr, vert_r) do

  atm_compute_signs_pt1(cr, er, vr)

  atm_compute_signs_pt2(cr, er, vr)

  --atm_adv_coef_compression(cr, er)

  --config_coef_3rd_order = 0.25 in namelist
  atm_couple_coef_3rd_order(0.25, cr, er)

  atm_init_coupled_diagnostics(cr, er, vert_r)

  atm_compute_solve_diagnostics(cr, er, vr, false) --last param is hollingsworth

  --mpas_reconstruct()

  atm_compute_mesh_scaling(cr, er, true)

  --config_zd: default 22000.0, config_xnutr: default 0.2. From config
  atm_compute_damping_coefs(22000, 0.2, cr)

end
