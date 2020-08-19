import "regent"
require "data_structures"

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


local cio = terralib.includec("stdio.h")
local cmath = terralib.includec("math.h")

task atm_compute_signs(vr : region(ispace(int1d), vertex_fs),
                       er : region(ispace(int2d), edge_fs),
                       cr : region(ispace(int2d), cell_fs))
where reads writes(vr, cr), reads (er) do 

    for iVtx = 1, nVertices do -- TODO: change bounds once you know whether vOnEdge contains ID's or indices
        for i = 0, vertexDegree do
            if (vr[iVtx].edgesOnVertex[i] <= nEdges) then
                --note: when MPAS calls vOnEdge(2, ...), we access vOnEdge index 1!
                if (iVtx == er[{vr[iVtx].edgesOnVertex[i], 0}].verticesOnEdge[1]) then
                    vr[iVtx].edgesOnVertexSign[i] = 1.0
                else
                    vr[iVtx].edgesOnVertexSign[i] = -1.0
                end
            else
                vr[iVtx].edgesOnVertexSign[i] = 0.0
            end
        end
    end


      for iCell=1,nCells do
         for i=1, cr[{iCell, 0}].nEdgesOnCell do
            if (cr[{iCell, 0}].edgesOnCell[i] <= nEdges) then
              if (iCell == er[{cr[{iCell, 0}].edgesOnCell[i], 0}].cellsOnEdge[0]) then
                  cr[{iCell, 0}].edgesOnCellSign[i] = 1.0
                  --VERTICAL STRUCTURE ACCESSED--
                  for vLevel = 0, nVertLevels + 1 do
                  
                    cr[{iCell, vLevel}].zb_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb[0]
                    cr[{iCell, vLevel}].zb3_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb3[0]

                  end 

               else
                  cr[{iCell, 0}].edgesOnCellSign[i] = -1.0
                  --VERTICAL--
                  for vLevel = 0, nVertLevels + 1 do
                  
                    cr[{iCell, vLevel}].zb_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb[1]
                    cr[{iCell, vLevel}].zb3_cell[i] = er[{cr[{iCell, 0}].edgesOnCell[i], vLevel}].zb3[1]

                  end 
               end 
            else
               cr[{iCell, 0}].edgesOnCellSign[i] = 0.0
            end 
         end 
      end 


    for iCell=1, nCells do
        for i=1, cr[{iCell, 0}].nEdgesOnCell do
            var iVtx = cr[{iCell, 0}].verticesOnCell[i]
            if (iVtx <= nVertices) then
                for j=1,vertexDegree do
                    if (iCell == vr[iVtx].cellsOnVertex[j]) then
                        cr[{iCell, 0}].kiteForCell[i] = j
                        break
                    end
                end 
                -- trimmed a log statement here
            else
                cr[{iCell, 0}].kiteForCell[i] = 1
            end 
        end 
    end 
end


task atm_adv_coef_compression(er : region(ispace(int2d), edge_fs),
                              cr : region(ispace(int2d), cell_fs))
where reads writes(er), reads(cr) do
    
    var cell_list : int[maxEdges]

    for iEdge = 1, nEdges do
        er[{iEdge, 0}].nAdvCellsForEdge = 0
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0] 
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1] 
        
         --
         -- do only if this edge flux is needed to update owned cells
         --
        if (cell1 <= nCells or cell2 <= nCells) then
            cell_list[1] = cell1
            cell_list[2] = cell2
            var n = 2 -- n is number of cells currently in list 
  
          --  add cells surrounding cell 1.  n is number of cells currently in list
            for i = 1, cr[{cell1, 0}].nEdgesOnCell do
               if (cr[{cell1, 0}].cellsOnCell[i] ~= cell2) then
                  n = n + 1
                  cell_list[n] = cr[{cell1, 0}].cellsOnCell[i]
               end 
            end 
  
          --  add cells surrounding cell 2 (brute force approach)
            for iCell = 1, cr[{cell2, 0}].nEdgesOnCell do
               var addcell = true
               for i=1, n do
                  if (cell_list[i] == cr[{cell2, 0}].cellsOnCell[iCell]) then addcell = false end
               end 
               if (addcell) then
                  n = n+1
                  cell_list[n] = cr[{cell2, 0}].cellsOnCell[iCell]
               end 
            end 
  
  
            er[{iEdge, 0}].nAdvCellsForEdge = n
            for iCell = 1, er[{iEdge, 0}].nAdvCellsForEdge do
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
            for j=1, n do
               if (cell_list[j] == cell1 ) then j_in = j end
            end 
            er[{iEdge, 0}].adv_coefs[j_in]= er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[0][0]
            er[{iEdge, 0}].adv_coefs_3rd[j_in]= er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[0][0]
  
            for iCell = 1, cr[{cell1, 0}].nEdgesOnCell do
               j_in = 0
               for j=1, n do
                 if( cell_list[j] == cr[{cell1, 0}].cellsOnCell[iCell]) then j_in = j end
               end
               er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[iCell][0]
               er[{iEdge, 0}].adv_coefs_3rd[j_in] = er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[iCell][0]
            end 
  
          -- pull together third and fourth order contributions to the flux
          -- now from cell2
  
            j_in = 0
            for j=1, n do
               if( cell_list[j] == cell2 ) then j_in = j end
            end 
              er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[0][1]
              er[{iEdge, 0}].adv_coefs_3rd[j_in] = er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[0][1]
  
            for iCell = 1, cr[{cell2, 0}].nEdgesOnCell do
               j_in = 0
               for j=1, n do
                  if( cell_list[j] == cr[{cell2, 0}].cellsOnCell[iCell] ) then j_in = j end
               end 
               er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + er[{iEdge, 0}].deriv_two[iCell][1]
               er[{iEdge, 0}].adv_coefs_3rd[j_in] = er[{iEdge, 0}].adv_coefs_3rd[j_in] + er[{iEdge, 0}].deriv_two[iCell][1]
            end 
  
            for j = 1,n do
              -- TODO: how to calculate exponent?
              -- er[{iEdge, 0}].adv_coefs[j] =  -1.0 * er[{iEdge, 0}].dcEdge**2 * er[{iEdge, 0}].adv_coefs[j] / 12 -- this should be a negative number
              -- er[{iEdge, 0}].adv_coefs_3rd[j] =  -1.0 * er[{iEdge, 0}].dcEdge**2 * er[{iEdge, 0}].adv_coefs_3rd[j] / 12 
--               adv_coefs    (j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs    (j,iEdge) / 12.
--               adv_coefs_3rd(j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
            end 
  
          -- 2nd order centered contribution - place this in the main flux weights
  
            j_in = 0
            for j=1, n do
               if( cell_list[j] == cell1 ) then j_in = j end
            end
            er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + 0.5
  
            j_in = 0
            for j=1, n do
               if( cell_list[j] == cell2 ) then j_in = j end
            end
            er[{iEdge, 0}].adv_coefs[j_in] = er[{iEdge, 0}].adv_coefs[j_in] + 0.5
  
          --  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply
  
            for j=1,n do
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
                               er : region(ispace(int2d), edge_fs),
                               cr : region(ispace(int2d), cell_fs))
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

--Comments for atm_advance_acoustic_step_work
--Constants declared in task: rgas, cp, epssm, dts
--Constants from mpas_constants.F : rgas, cp, gravity. These are passed in to the task as doubles at the moment
--Contants from config: config_epssm: default = 0.1. Passed in as double at the moment.
--dts is passed in as double now - in code, passed in as rk_sub_timestep(rk_step)
--1.0_RKIND translated as just 1.0
--Tendency variables: tend_rw, tend_rt, tend_rho (added to CR), tend_ru (added to ER). Not in registry
--Other variables not in registry: rs, ts: added to CR as part of vertical grid, ru_Avg (added to ER)

task atm_advance_acoustic_step_work(cr : region(ispace(int2d), cell_fs), er : region(ispace(int2d), edge_fs), rgas : double, cp : double, gravity : double, dts : double, config_epssm: double, edgeStart : int, edgeEnd : int, cellStart : int, cellEnd : int, cellSolveStart : int, cellSolveEnd : int, small_step : int, nCellsSolve : int)
where reads writes (er, cr) do
  var epssm = config_epssm
  var rcv = rgas / (cp - rgas)
  var c2 = cp * rcv
  var resm = (1.0 - epssm) / (1.0 + epssm)
  var rdts = 1./dts
  cio.printf("advancing acoustic step\n")

  if (small_step == 1) then
    for iEdge = edgeStart-1, edgeEnd do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

      if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

        for k = 0, nVertLevels do
          var pgrad = ((cr[{cell2, k}].rtheta_pp - cr[{cell1, k}].rtheta_pp) * er[{iEdge, 0}].invDcEdge) / (0.5 * (cr[{cell2, k}].zz +cr[{cell1, k}].zz))
          pgrad = er[{iEdge, k}].cqu * 0.5 * c2 * (cr[{cell1, k}].exner + cr[{cell2, k}].exner) * pgrad
          pgrad = pgrad + 0.5 * er[{iEdge, k}].zxu * gravity * (cr[{cell1, k}].rho_pp + cr[{cell2, k}].rho_pp)
          er[{iEdge, k}].ru_p = er[{iEdge, k}].ru_p + dts * (er[{iEdge, k}].tend_ru - (1.0 - er[{iEdge, 0}].specZoneMaskEdge) * pgrad)  --NEEDS FIXING
        end
        for k = 0, nVertLevels do
          er[{iEdge, k}].ruAvg = er[{iEdge, k}].ruAvg + er[{iEdge, k}].ru_p
        end
      end
    end

  else
    for iEdge = edgeStart-1, edgeEnd do
      var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
      var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

      if (cell1 <= nCellsSolve or cell2 <= nCellsSolve) then

        for k = 0, nVertLevels do
          er[{iEdge, k}].ru_p = dts * er[{iEdge, k}].tend_ru
        end

        for k = 0, nVertLevels do
          er[{iEdge, k}].ruAvg = er[{iEdge, k}].ru_p
        end

      end
    end
  end

  if (small_step == 1) then
    for iCell = cellStart - 1, cellEnd do
      for j = 0, nVertLevels do
        cr[{iCell, j}].rtheta_pp_old = 0
      end
    end
  else
    for iCell = cellStart, cellEnd do
      for j = 0, nVertLevels do
        cr[{iCell, j}].rtheta_pp_old = cr[{iCell, j}].rtheta_pp
      end
    end
  end

  for iCell = cellSolveStart, cellSolveEnd do
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
        cr[{0, i}].ts = 0
        cr[{0, i}].rs = 0
      end

      for i = 0, cr[{iCell, 0}].nEdgesOnCell do
        var iEdge = cr[{iCell, 0}].edgesOnCell[i]    --edgesOnCell(i,iCell)
        var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
        var cell2 = er[{iEdge, 0}].cellsOnEdge[1]    --cell2 = cellsOnEdge(2,iEdge)

        for k = 0, nVertLevels do
          var flux = cr[{iCell, 0}].edgesOnCell_sign[i] * dts * er[{iEdge, 0}].dvEdge * er[{iEdge, k}].ru_p * cr[{iCell, 0}].invAreaCell
          cr[{0, k}].rs = cr[{0, k}].rs - flux
          cr[{0, k}].ts = cr[{0, k}].ts - flux * 0.5 * (cr[{cell2, k}].theta_m + cr[{cell1, k}].theta_m)
        end
      end

      for k = 0, nVertLevels do
         cr[{0, k}].rs = cr[{iCell, k}].rho_pp + dts * cr[{iCell, k}].tend_rho + cr[{0, k}].rs - cr[{0, k}].cofrz * resm * (cr[{iCell, k+1}].rw_p - cr[{iCell, k}].rw_p)
         cr[{0, k}].ts = cr[{iCell, k}].rtheta_pp + dts * cr[{iCell, k}].tend_rt + cr[{0, k}].ts - resm * cr[{0, k}].rdzw * (cr[{iCell, k+1}].coftz * cr[{iCell, k+1}].rw_p  - cr[{iCell, k}].coftz * cr[{iCell, k}].rw_p )
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].wwAvg = cr[{iCell, k}].wwAvg + 0.5 * (1.0 - epssm) * cr[{iCell, k}].rw_p
      end

      for k = 1, nVertLevels do
         cr[{iCell, k}].rw_p = cr[{iCell, k}].rw_p +  dts * cr[{iCell, k}].tend_rw - cr[{iCell, k}].cofwz * ((cr[{iCell, k}].zz * cr[{0, k}].ts - cr[{iCell, k-1}].zz * cr[{0, k-1}].ts) + resm * (cr[{iCell, k}].zz * cr[{iCell, k}].rtheta_pp - cr[{iCell, k-1}].zz * cr[{iCell, k-1}].rtheta_pp)) - cr[{iCell, k}].cofwr * ((cr[{0, k}].rs + cr[{0, k-1}].rs) + resm * (cr[{iCell, k}].rho_pp + cr[{iCell, k-1}].rho_pp))  + cr[{iCell, k}].cofwt * (cr[{0, k}].ts + resm * cr[{iCell, k}].rtheta_pp) + cr[{iCell, k-1}].cofwt * (cr[{0, k-1}].ts +resm * cr[{iCell, k-1}].rtheta_pp)
      end

      for k = 1, nVertLevels do
         cr[{iCell, k}].rw_p = (cr[{iCell, k}].rw_p - cr[{iCell, k}].a_tri * cr[{iCell, k-1}].rw_p) * cr[{iCell, k}].alpha_tri
      end

      for k = nVertLevels, 1, -1 do
        cr[{iCell, k}].rw_p = cr[{iCell, k}].rw_p - cr[{iCell, k}].gamma_tri * cr[{iCell, k+1}].rw_p
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].rw_p = (cr[{iCell, k}].rw_p + (cr[{iCell, k}].rw_save - cr[{iCell, k}].rw) - dts * cr[{iCell, k}].dss * (cr[{0, k}].fzm * cr[{iCell, k}].zz + cr[{0, k}].fzp * cr[{iCell, k-1}].zz)*(cr[{0, k}].fzm * cr[{iCell, k}].rho_zz + cr[{0, k}].fzp * cr[{iCell, k-1}].rho_zz) * cr[{iCell, k}].w)/(1.0 + dts * cr[{iCell, k}].dss)  - (cr[{iCell, k}].rw_save - cr[{iCell, k}].rw)
      end

      for k = 1, nVertLevels do
        cr[{iCell, k}].wwAvg = cr[{iCell, k}].wwAvg + 0.5 * (1.0 + epssm) * cr[{iCell, k}].rw_p
      end

      for k=0, nVertLevels do
         cr[{iCell, k}].rho_pp  = cr[{0, k}].rs - cr[{0, k}].cofrz *(cr[{iCell, k+1}].rw_p - cr[{iCell, k}].rw_p)
         cr[{iCell, k}].rtheta_pp = cr[{0, k}].ts  - cr[{0, k}].rdzw * (cr[{iCell, k+1}].coftz * cr[{iCell, k+1}].rw_p - cr[{iCell, k}].coftz * cr[{iCell, k}].rw_p)
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



task atm_compute_moist_coefficients()
  cio.printf("computing moist coefficients\n")
end  

--Comments for atm_compute_vert_imp_coefs
--Combined this and atm_compute_vert_imp_coefs_work (one was just a helper)
--purpose: "Compute coefficients for vertically implicit gravity-wave/acoustic computations"
--Constants declared in task: 
--Constants from mpas_constants.F : . These are passed in to the task as doubles at the moment
--Contants from config: config_epssm: default = 0.1. Passed in as double at the moment.
--dts is passed in as double now - in code, passed in as rk_sub_timestep(rk_step)
--1.0_RKIND translated as just 1.0
--mpas used cellSolveStart, cellSolveEnd; we believe those deal with how mpas uses parallelization (a start and end for which cell pool is currently being worked on). We will just loop thru the entire cell region instead.
--mpas renames some fields; we will just use the original name --- p : exner, t : theta_m, rb : rho_base, rtb : rtheta_base, pb : exner_base, rt : rtheta_p

task atm_compute_vert_imp_coefs(cr : region(ispace(int2d), cell_fs),
                                dts : double,
                                epssm : double,
                                rgas : double,
                                cp : double,
                                gravity : double)  
where reads writes (cr) do
      --  set coefficients
      var dtseps = .5*dts*(1.+epssm)
      var rcv = rgas/(cp-rgas)
      var c2 = cp*rcv
      
      var qtotal : double
      var b_tri : double[nVertLevels]
      var c_tri : double[nVertLevels]


-- MGD bad to have all threads setting this variable?
      for k = 0, nVertLevels do
         cr[{0, k}].cofrz = dtseps * cr[{0, k}].rdzw
      end 

      for iCell = 0, nCells do --  we only need to do cells we are solving for, not halo cells

--DIR$ IVDEP
         for k=1, nVertLevels do
            cr[{iCell, k}].cofwr = .5 * dtseps * gravity * (cr[{0, k}].fzm * cr[{iCell, k}].zz + cr[{0, k}].fzp * cr[{iCell, k-1}].zz)
         end 
         cr[{iCell, 0}].coftz = 0.0 --coftz(1,iCell) = 0.0
--DIR$ IVDEP
         for k=1, nVertLevels do
            cr[{iCell, k}].cofwz = dtseps * c2 * (cr[{0, k}].fzm * cr[{iCell, k}].zz + cr[{0, k}].fzp * cr[{iCell, k-1}].zz) * cr[{0, k}].rdzu * cr[{iCell, k}].cqw * (cr[{0, k}].fzm * cr[{iCell, k}].exner + cr[{0, k}].fzp * cr[{iCell, k-1}].exner)
            cr[{iCell, k}].coftz = dtseps * (cr[{0, k}].fzm * cr[{iCell, k}].theta_m + cr[{0, k}].fzp * cr[{iCell, k-1}].theta_m)
         end 
         cr[{iCell, nVertLevels}].coftz = 0.0 -- coftz(nVertLevels+1,iCell) 
--DIR$ IVDEP
         for k=0, nVertLevels do

            qtotal = cr[{iCell, k}].qtot

            cr[{iCell, k}].cofwt = .5 * dtseps * rcv * cr[{iCell, k}].zz * gravity * cr[{iCell, k}].rho_base / ( 1.0 + qtotal) * cr[{iCell, k}].exner / ((cr[{iCell, k}].rtheta_base + cr[{iCell, k}].rtheta_p) * cr[{iCell, k}].exner_base)
         end 

         cr[{iCell, 0}].a_tri = 0.0 --a_tri(1,iCell) = 0.  -- note, this value is never used
         b_tri[0] = 1.0    -- note, this value is never used
         c_tri[0] = 0.0    -- note, this value is never used
         cr[{iCell, 0}].gamma_tri = 0.0
         cr[{iCell, 0}].alpha_tri = 0.0  -- note, this value is never used

--DIR$ IVDEP
         for k=1, nVertLevels do --k=2,nVertLevels
            cr[{iCell,k}].a_tri = -1.0 * cr[{iCell, k}].cofwz * cr[{iCell, k-1}].coftz * cr[{0, k-1}].rdzw * cr[{iCell, k-1}].zz + cr[{iCell, k}].cofwr * cr[{0, k-1}].cofrz - cr[{iCell, k-1}].cofwt * cr[{iCell, k-1}].coftz * cr[{0, k-1}].rdzw
         
            b_tri[k] = 1.0 + cr[{iCell, k}].cofwz * (cr[{iCell, k}].coftz * cr[{0, k}].rdzw * cr[{iCell, k}].zz +  cr[{iCell, k}].coftz * cr[{0, k-1}].rdzw * cr[{iCell, k-1}].zz) -  cr[{iCell, k}].coftz * (cr[{iCell, k}].cofwt * cr[{0, k}].rdzw - cr[{iCell, k-1}].cofwt * cr[{0, k-1}].rdzw) + cr[{iCell, k}].cofwr * ((cr[{0, k}].cofrz - cr[{0, k-1}].cofrz))
            
            c_tri[k] =   -1.0 * cr[{iCell, k}].cofwz * cr[{iCell, k+1}].coftz * cr[{0, k}].rdzw * cr[{iCell, k}].zz - cr[{iCell, k}].cofwr * cr[{0, k}].cofrz + cr[{iCell, k}].cofwt * cr[{iCell, k+1 }].coftz * cr[{iCell, k}].rdzw
         end 
--MGD VECTOR DEPENDENCE
         for k=1, nVertLevels do -- k=2, nVertLevels
            cr[{iCell, k}].alpha_tri = 1.0/ (b_tri[k]-cr[{iCell, k}].a_tri * cr[{iCell, k-1}].gamma_tri)
            cr[{iCell, k}].gamma_tri = c_tri[k] * cr[{iCell, k}].alpha_tri
         end 

      end -- loop over cells
end

task atm_compute_dyn_tend()
  cio.printf("computing dynamic tendencies\n")
end

task atm_set_smlstep_pert_variables()
  cio.printf("set small step vars\n")
end

task atm_advance_acoustic_step()
  cio.printf("advancing acoustic step\n")
end

task atm_divergence_damping_3d()
  cio.printf("update horizontal momentum\n")
end

task atm_recover_large_step_variables()
  cio.printf("recovering large step vars\n")
end

task atm_rk_dynamics_substep_finish()
  cio.printf("finishing substep\n")
end
