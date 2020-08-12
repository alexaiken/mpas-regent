import "regent"
require "data_structures"


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

task atm_compute_signs(vr : region(ispace(int1d), vertex_fs),
                       er : region(ispace(int1d), edge_fs),
                       cr : region(ispace(int1d), cell_fs))
where reads writes(vr, cr), reads (er) do 

    for iVtx=1, nVertices do
        for i=1, vertexDegree do
            if (vr[iVtx].edgesOnVertex[i] <= nEdges) then
                --note: when MPAS calls vOnEdge(2, ...), we access vOnEdge index 1!
                if (iVtx == er[vr[iVtx].edgesOnVertex[i]].verticesOnEdge[1]) then
                    vr[iVtx].edgesOnVertexSign[i] = 1.0
                else
                    vr[iVtx].edgesOnVertexSign[i] = -1.0
                end
            else
                vr[iVtx].edgesOnVertexSign[i] = 0.0
            end
        end
    end


-- TODO: port over to Regent once vertical structure representation is figured out (need zb)
--      for iCell=1,nCells do
--         for i=1, cr[iCell].nEdgesOnCell do
--            if (cr[iCell].edgesOnCell[i] <= nEdges) then
--              if (iCell == er[cr[iCell].edgesOnCell[i]].cellsOnEdge[0]) then
--                  cr[iCell].edgesOnCellSign[i] = 1.0
--                  --VERTICAL STRUCTURE NEEDED--
--                  zb_cell(:,i,iCell) = zb(:,1,edgesOnCell(i,iCell))
--                  zb3_cell(:,i,iCell) = zb3(:,1,edgesOnCell(i,iCell))
--               else
--                  cr[iCell].edgesOnCellSign[i] = -1.0
--                  --VERTICAL--
--                  zb_cell(:,i,iCell) = zb(:,2,edgesOnCell(i,iCell))
--                  zb3_cell(:,i,iCell) = zb3(:,2,edgesOnCell(i,iCell))
--               end 
--            else
--               edgesOnCell_sign(i,iCell) = 0.0
--            end 
--         end 
--      end 


    for iCell=1, nCells do
        for i=1, cr[iCell].nEdgesOnCell do
            var iVtx = cr[iCell].verticesOnCell[i]
            if (iVtx <= nVertices) then
                for j=1,vertexDegree do
                    if (iCell == vr[iVtx].cellsOnVertex[j]) then
                        cr[iCell].kiteForCell[i] = j
                        break
                    end
                end 
                -- trimmed a log statement here
            else
                cr[iCell].kiteForCell[i] = 1
            end 
        end 
    end 
end


task atm_adv_coef_compression(er : region(edge_fs),
                              cr : region(cell_fs))
where reads writes(er), reads(cr) do
    
    var cell_list : int[maxEdges]

    for iEdge = 1, nEdges do
        er[iEdge].nAdvCellsForEdge = 0
        var cell1 = er[iEdge].cellsOnEdge[0] 
        var cell2 = er[iEdge].cellsOnEdge[1] 
        
         --
         -- do only if this edge flux is needed to update owned cells
         --
        if (cell1 <= nCells or cell2 <= nCells) then
            cell_list[1] = cell1
            cell_list[2] = cell2
            var n = 2 -- n is number of cells currently in list 
  
          --  add cells surrounding cell 1.  n is number of cells currently in list
            for i = 1, cr[cell1].nEdgesOnCell do
               if (cr[cell1].cellsOnCell[i] ~= cell2) then
                  n = n + 1
                  cell_list[n] = cr[cell1].cellsOnCell[i]
               end 
            end 
  
          --  add cells surrounding cell 2 (brute force approach)
            for iCell = 1, cr[cell2].nEdgesOnCell do
               var addcell = true
               for i=1, n do
                  if (cell_list[i] == cr[cell2].cellsOnCell[iCell]) then addcell = false end
               end 
               if (addcell) then
                  n = n+1
                  cell_list[n] = cr[cell2].cellsOnCell[iCell]
               end 
            end 
  
  
            er[iEdge].nAdvCellsForEdge = n
            for iCell = 1, er[iEdge].nAdvCellsForEdge do
               er[iEdge].advCellsForEdge[iCell] = cell_list[iCell]
            end 
  
          -- we have the ordered list, now construct coefficients
          --TODO: add adv_coefs to fspaces
            for coef = 0, FIFTEEN do
                er[iEdge].adv_coefs[coef] = 0.0
                er[iEdge].adv_coefs_3rd[coef] = 0.0
            end -- initialize list to 0
          
          -- pull together third and fourth order contributions to the flux
          -- first from cell1
  
            var j_in = 0
            for j=1, n do
               if (cell_list[j] == cell1 ) then j_in = j end
            end 
            er[iEdge].adv_coefs[j_in]= er[iEdge].adv_coefs[j_in] + er[iEdge].deriv_two[0][0]
            er[iEdge].adv_coefs_3rd[j_in]= er[iEdge].adv_coefs_3rd[j_in] + er[iEdge].deriv_two[0][0]
  
            for iCell = 1, cr[cell1].nEdgesOnCell do
               j_in = 0
               for j=1, n do
                 if( cell_list[j] == cr[cell1].cellsOnCell[iCell]) then j_in = j end
               end
               er[iEdge].adv_coefs[j_in] = er[iEdge].adv_coefs[j_in] + er[iEdge].deriv_two[iCell][0]
               er[iEdge].adv_coefs_3rd[j_in] = er[iEdge].adv_coefs_3rd[j_in] + er[iEdge].deriv_two[iCell][0]
            end 
  
          -- pull together third and fourth order contributions to the flux
          -- now from cell2
  
            j_in = 0
            for j=1, n do
               if( cell_list[j] == cell2 ) then j_in = j end
            end 
              er[iEdge].adv_coefs[j_in] = er[iEdge].adv_coefs[j_in] + er[iEdge].deriv_two[0][1]
              er[iEdge].adv_coefs_3rd[j_in] = er[iEdge].adv_coefs_3rd[j_in] + er[iEdge].deriv_two[0][1]
  
            for iCell = 1, cr[cell2].nEdgesOnCell do
               j_in = 0
               for j=1, n do
                  if( cell_list[j] == cr[cell2].cellsOnCell[iCell] ) then j_in = j end
               end 
               er[iEdge].adv_coefs[j_in] = er[iEdge].adv_coefs[j_in] + er[iEdge].deriv_two[iCell][1]
               er[iEdge].adv_coefs_3rd[j_in] = er[iEdge].adv_coefs_3rd[j_in] + er[iEdge].deriv_two[iCell][1]
            end 
  
            for j = 1,n do
              -- er[iEdge].adv_coefs[j] =  - er[iEdge].dcEdge**2 * er[iEdge].adv_coefs[j] / 12 -- this should be a negative number
              -- er[iEdge].adv_coefs_3rd[j] =  - er[iEdge].dcEdge**2 * er[iEdge].adv_coefs_3rd[j] / 12 
--               adv_coefs    (j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs    (j,iEdge) / 12.
--               adv_coefs_3rd(j,iEdge) = - (dcEdge(iEdge) **2) * adv_coefs_3rd(j,iEdge) / 12.
            end 
  
          -- 2nd order centered contribution - place this in the main flux weights
  
            j_in = 0
            for j=1, n do
               if( cell_list[j] == cell1 ) then j_in = j end
            end
            er[iEdge].adv_coefs[j_in] = er[iEdge].adv_coefs[j_in] + 0.5
  
            j_in = 0
            for j=1, n do
               if( cell_list[j] == cell2 ) then j_in = j end
            end
            er[iEdge].adv_coefs[j_in] = er[iEdge].adv_coefs[j_in] + 0.5
  
          --  multiply by edge length - thus the flux is just dt*ru times the results of the vector-vector multiply
  
            for j=1,n do
               er[iEdge].adv_coefs[j] = er[iEdge].dvEdge * er[iEdge].adv_coefs[j]
               er[iEdge].adv_coefs_3rd[j] = er[iEdge].dvEdge * er[iEdge].adv_coefs_3rd[j]
            end 
 
         end  -- only do for edges of owned-cells
         
    end -- end loop over edges 


end





--task atm_compute_damping_coefs()
--    var m1 = -1.0
--    var pii = acos(m1) -- find equivelent transformation in Regent for acos()
--    var dx_scale_power = 1.0
--    dss(:,:) = 0.0 --dss should come from the mesh (dimensions vertlevels x ncells)
--    for iCell=1, nCells do
--        zt = zgrid(nVertLevels+1,iCell)
--        for k=1, nVertLevels do
--            z = 0.5*(zgrid(k,iCell) + zgrid(k+1,iCell))
--            if (z > config_zd) then
--                dss(k,iCell) = config_xnutr*sin(0.5*pii*(z-config_zd)/(zt-config_zd))**2.0
 --               dss(k,iCell) = dss(k,iCell) / meshDensity(iCell)**(0.25*dx_scale_power)
 --           end 
 --       end 
--    end 
--end

task atm_couple_coef_3rd_order()
  --  adv_coefs_3rd(:,:) = config_coef_3rd_order * adv_coefs_3rd(:,:)
   -- zb3_cell(:,:,:) = config_coef_3rd_order * zb3_cell(:,:,:)
end


task atm_compute_moist_coefficients()
  cio.printf("computing moist coefficients\n")
end  

task atm_compute_vert_imp_coefs()  
  cio.printf("computing vertical coefficients\n")
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

task atm_compute_solve_diagnostics()
  cio.printf("computing solve diagnostics\n")
end

task atm_rk_dynamics_substep_finish()
  cio.printf("finishing substep\n")
end
