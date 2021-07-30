import "regent"
require "data_structures"

local constants = require("constants")
local format = require("std/format")
local nCells = constants.nCells
local nEdges = constants.nEdges
local nVertices = constants.nVertices
local maxEdges = constants.maxEdges
local maxEdges2 = constants.maxEdges2
local TWO = constants.TWO
local FIFTEEN = constants.FIFTEEN
local vertexDegree = constants.vertexDegree
local nVertLevels = constants.nVertLevels
local sphere_radius = constants.sphere_radius
local nlat = constants.nlat

local cio = terralib.includec("stdio.h")
-- local cmath = terralib.includec("math.h")
-- Can't use terralib.includec("math.h") with Cuda.
local pow = regentlib.pow(double)
local cos = regentlib.cos(double)
local sin = regentlib.sin(double)
local exp = regentlib.exp(double)
local sqrt = regentlib.sqrt(double)


__demand(__cuda)
task init_atm_case_jw(cr : region(ispace(int2d), cell_fs),
                      cpr : region(ispace(int2d), cell_fs),
                      csr : region(ispace(int2d), cell_fs),
                      cgr : region(ispace(int2d), cell_fs),
                      er : region(ispace(int2d), edge_fs),
                      vr : region(ispace(int2d), vertex_fs),
                      vertr : region(ispace(int1d), vertical_fs))
where
  reads (cr.lat,
         er.{cellsOnEdge, edgesOnEdge_ECP, lat, lon, nEdgesOnEdge, verticesOnEdge, weightsOnEdge}, 
         vr.{lat, lon}),
  writes (cr.{qsat, relhum, rho, rtheta_base, rtheta_p, surface_pressure, theta, w},
          er.{fEdge, zxu},
          vr.fVertex,
          vertr.{rdzu, rdzw}),
  reads writes (cr.{areaCell, dss, exner, hx, pressure_base, pressure_p, qv, rho_base, rho_p, rho_zz, 
                    rw, surface_pressure, theta_base, theta_m, x, y, z, zgrid, zz},
                er.{dvEdge, dcEdge, ru, u, v, x, y, z, zb, zb3},
                vr.{areaTriangle, kiteAreasOnVertex, x, y, z},
                vertr.{dzu, fzm, fzp, cf1, cf2, cf3})
do

-- Ranges for Cuda code generation
  var cell_range_1d = rect1d {0, nCells - 1 }
  var cell_range_1d_0 = rect2d { int2d {0, 0}, int2d {nCells - 1, 0} }

  var edge_range_1d = rect1d { 0, nEdges - 1 }
  var edge_range_1d_0 = rect2d { int2d {0, 0}, int2d {nEdges - 1, 0} }
  var edge_range_2d = rect2d { int2d {0, 0}, int2d {nEdges - 1, nVertLevels - 1} }

  var vertex_range_1d = rect2d { int2d {0, 0}, int2d {nVertices - 1, 0} } 

  var vertical_range_1d = rect1d {0, nVertLevels - 1 } 
  var vertical_range_s1_1d = rect1d { 1, nVertLevels - 1 } 
  var vertical_range_extra_1d = rect1d { 0, nVertLevels } -- Loop goes to nVertLevels + 1. Should be replaced with rect1d.
  var vertical_nlat_range = rect2d { int2d {0, 0}, int2d {nVertLevels - 1, nlat - 1} } 

  var nlat_range_1d = rect1d { 0, nlat - 1 }

  var cell_range_2d = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
  var cell_range_extra_2d = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels} }
  
-- local vars/constants defined at beginning of subroutine
  var cp = constants.cp
  var rgas = constants.rgas
  var gravity = constants.gravity
  var pii = constants.pii
  var u0 = 35.0
  var alpha_grid = 0.0 -- no grid rotation
  var omega_e : double
  -- Moved because of Cuda. 
  -- var flux : double 
  -- var u_pert : double
  -- var phi : double 
  -- var eoe : int
  -- var ztemp : double 
  var t0b = 250.0
  var t0 = 288.0
  var delta_t = 4.8E+05 -- TODO: scientific notation in regent?
  var dtdz = 0.005
  var eta_t = 0.2
  var u_perturbation = 1.0
  var pert_radius = 0.1
  var latitude_pert = 40.0
  var longitude_pert = 20.0
  var theta_c = pii/4.0
  var lambda_c = 3.0*pii/2.0
  var k_x = 9.0           -- Normal mode wave number

--initialization of moisture:
  var rh_max = 0.40 -- Maximum relative humidity
  -- scalars, index_qv, qsat, relhum
--end initialization

  var rebalance = true
  --var nlat = 721
  var moisture = false
  var nz1 = nVertLevels
  var nz = nz1 + 1
  var arr : double[2][3]

  var qv_2d : double[nVertLevels * nlat]
---- Scale all distances and areas from a unit sphere to one with radius sphere_radius

  for iCell in cell_range_1d_0 do
    -- iCell is the same as { iCell.x, 0 }
    cr[iCell].x = cr[iCell].x * sphere_radius
    cr[iCell].y = cr[iCell].y * sphere_radius
    cr[iCell].z = cr[iCell].z * sphere_radius
    --cr[iCell].areaCell = cr[iCell].areaCell * cmath.pow(sphere_radius, 2.0)
    cr[iCell].areaCell = cr[iCell].areaCell * sphere_radius * sphere_radius
  end

  for iVert in vertex_range_1d do
    -- iVert is the same as { iVert.x, 0 }
    vr[iVert].x = vr[iVert].x * sphere_radius
    vr[iVert].y = vr[iVert].y * sphere_radius
    vr[iVert].z = vr[iVert].z * sphere_radius
    --vr[iVert].areaTriangle = vr[iVert].areaTriangle * cmath.pow(sphere_radius, 2.0)
    vr[iVert].areaTriangle = vr[iVert].areaTriangle * sphere_radius * sphere_radius
    for vDeg = 0, vertexDegree do
      --vr[iVert].kiteAreasOnVertex[vertexDegree] = vr[iVert].kiteAreasOnVertex[vertexDegree] * cmath.pow(sphere_radius, 2.0)
      vr[iVert].kiteAreasOnVertex[vertexDegree] = vr[iVert].kiteAreasOnVertex[vertexDegree] * sphere_radius * sphere_radius
    end
  end

  for iEdge in edge_range_1d_0 do
    -- iEdge is the same as { iEdge.x, 0 }
    er[iEdge].x = er[iEdge].x * sphere_radius
    er[iEdge].y = er[iEdge].y * sphere_radius
    er[iEdge].z = er[iEdge].z * sphere_radius
    er[iEdge].dvEdge = er[iEdge].dvEdge * sphere_radius
    er[iEdge].dcEdge = er[iEdge].dcEdge * sphere_radius
  end




-- initialization of moisture:
  for iCell in cell_range_2d do
    -- iCell = { nCells, nVertLevels }
    cr[iCell].qv = 0.0
    cr[iCell].qsat = 0.0
    cr[iCell].relhum = 0.0
  end

  for iVertLat in vertical_nlat_range do
    qv_2d[iVertLat.x * nlat + iVertLat.y] = 0.0
  end

-- end initialization of moisture.


  for iCell in cell_range_1d_0 do
    -- iCell is the same as { iCell.x, 0 }
    cr[iCell].surface_pressure = 0.0
  end

--      call atm_initialize_advection_rk(mesh, nCells, nEdges, maxEdges, on_a_sphere, sphere_radius )
--      call atm_initialize_deformation_weights(mesh, nCells, on_a_sphere, sphere_radius)

--      call mpas_pool_get_dimension(state, 'index_qv', index_qv)


  var xnutr = 0.0
  var zd = 12000.0
  var znut = eta_t
  var etavs = (1.0 - 0.252) * pii/2.
  var r_earth = sphere_radius
  omega_e = constants.omega
  var p0 = 1.0E+05

--! We may pass in an hx(:,:) that has been precomputed elsewhere.
--! For now it is independent of k
  for iCell in cell_range_extra_2d do
    -- Declared on each iteration to avoid compiler warning for loop-carried dependecy.
    var phi = cr[{iCell.x, 0}].lat
    cr[iCell].hx = u0 / gravity * pow(cos(etavs), 1.5) * ((-2.0 * pow(sin(phi), 6) * (pow(cos(phi), 2.0) + 1.0/3.0) + 10.0/63.0) * u0 * pow(cos(etavs), 1.5) + (1.6 * pow(cos(phi), 3) * (pow(sin(phi), 2) + 2.0/3.0) - pii/4.0) * r_earth * omega_e)
  end


--      !     Metrics for hybrid coordinate and vertical stretching

  var str = 1.5
  var zt = 45000.0
  var dz = zt / nz1
  var sh : double[nVertLevels + 1]
  var zw : double[nVertLevels + 1]
  var ah : double[nVertLevels + 1]
  var dzw : double[nVertLevels]
  var zu : double[nVertLevels]
  var rdzwp : double[nVertLevels]
  var rdzwm : double[nVertLevels]

  for iVert in vertical_range_extra_1d do
    var k = int(iVert)
    if k == 0 then
      sh[k] = -1
    else
  --!  sh[k] is the stretching specified for height surfaces
      sh[k] = pow(([double](k - 1.0) * dz / zt), str) --- this was (real(k-1)) in mpas; do we need to cast it here too?
    end

    --cio.printf("dz is %f \n", dz)
    --cio.printf("zt is %f \n", zt)
    --cio.printf("str is %f \n", str)
    --cio.printf("cmath.pow(([float](%d - 1.0) * dz=%f / zt=%f), str) is %f \n", k, dz, zt, cmath.pow(([float](k - 1.0) * dz / zt), str))


--!           to specify specific heights zc(k) for coordinate surfaces,
--!           input zc(k) and define sh[k] = zc(k)/zt
--!           zw(k) is the hieght of zeta surfaces
--!                zw[k] = (k-1)*dz yields constant dzeta
--!                        and nonconstant dzeta/dz
--!                zw[k] = sh[k]*zt yields nonconstant dzeta
--!                        and nearly constant dzeta/dz

    zw[k] = (k - 1) * dz -- in mpas they cast to float

--!            zw[k] = sh[k]*zt --- see above comments for which version you want
--!
--!           ah[k] governs the transition between terrain-following
--!           and pureheight coordinates
--!                ah[k] = 0 is a terrain-following coordinate
--!                ah[k] = 1 is a height coordinate
--
    ah[k] = 1.0 - pow(cos(.5 * pii * (k - 1) * dz / zt), 6.0)
--!            ah[k] = 0.

    --cio.printf("sh[%d] is %f \n", k, sh[k])
    --cio.printf("ah[%d] is %f \n", k, ah[k])
    --cio.printf("zw[%d] is %f \n", k, zw[k])
  end

  for iVert in vertical_range_1d do -- Old comment: nz1 is just nVertLevels, idk why mpas renamed it
    var k = int(iVert)
    dzw[k] = zw[k + 1] - zw[k]
    -- Can't mix assignments to stack-allocated arrays and regions.
    -- vertr[k].rdzw = 1.0 / dzw[k] 
    zu[k] = .5 * (zw[k] + zw[k + 1])
  end
  -- Split for loop into two to allow for Cuda code generation.
  for iVert in vertical_range_1d do -- Old comment: nz1 is just nVertLevels, idk why mpas renamed it
    var k = int(iVert)
    vertr[k].rdzw = 1.0 / dzw[k]
  end

  for kIndex in vertical_range_s1_1d do -- k=2,nz1 in mpas
    var k = int(kIndex)
    vertr[k].dzu = .5 * (dzw[k] + dzw[k - 1])
    vertr[k].rdzu = 1.0 / vertr[k].dzu
    vertr[k].fzp = .5 * dzw[k] / vertr[k].dzu
    vertr[k].fzm = .5 * dzw[k - 1] / vertr[k].dzu
  end
  -- Split for loop into two to allow for Cuda code generation.
  for kIndex in vertical_range_s1_1d do -- k=2,nz1 in mpas
    var k = int(kIndex)
    rdzwp[k] = dzw[k - 1] / (dzw[k] * (dzw[k] + dzw[k - 1]))
    rdzwm[k] = dzw[k] / (dzw[k - 1] * (dzw[k] + dzw[k - 1]))
  end



--!**********  how are we storing cf1, cf2 and cf3?

-- TODO: What to do with this?
--  var COF1 = (2. * vertr[1].dzu + vertr[2].dzu) / (vertr[1].dzu + vertr[2].dzu) * dzw[0] / vertr[1].dzu
--  var COF2 = vertr[1].dzu / (vertr[1].dzu + vertr[2].dzu) * dzw[0] / vertr[2].dzu
--  vertr[0].cf1 = vertr[1].fzp + COF1
--  vertr[0].cf2 = vertr[1].fzm - COF1 - COF2
--  vertr[0].cf3 = COF2


  for iCell in cell_range_extra_2d do
    cr[iCell].zgrid = (1.0 - ah[iCell.y]) * (sh[iCell.y] * (zt - cr[iCell].hx) + cr[iCell].hx) + ah[iCell.y] * sh[iCell.y] * zt
    --cio.printf("ah[%d] is %f \n", k, ah[k]) --cio.printf("sh[%d] is %f \n", k, sh[k])
    --cio.printf("cr[{%d, %d}].hx is %f\n", iCell, k, cr[{iCell, k}].hx)
  end
  -- Split for loop into two to allow for Cuda code generation.
  for iCell in cell_range_2d do
    cr[iCell].zz = (zw[iCell.y + 1] - zw[iCell.y]) / (cr[{iCell.x, iCell.y + 1}].zgrid - cr[iCell].zgrid)
    --cio.printf("cr[{iCell = %d, k = %d}].zz is %f \n", iCell, k, cr[{iCell, k}].zz)
    --cio.printf("zw[%d] is %f \n", k, zw[k]) : zw is set
    --cio.printf("cr[{%d, %d}].zgrid is %f", iCell, k, cr[{iCell,k}].zgrid) : zgrid is not set
  end

  for iCell in cell_range_2d do
    var iCell1 = er[{iCell.x, 0}].cellsOnEdge[0] --cellsOnEdge(1,i)
    var iCell2 = er[{iCell.x, 0}].cellsOnEdge[1] --cellsOnEdge(2,i)
    er[iCell].zxu = 0.5 * (cr[{iCell2, iCell.y}].zgrid - cr[{iCell1, iCell.y}].zgrid + cr[{iCell2, iCell.y + 1}].zgrid - cr[{iCell1, iCell.y + 1}].zgrid) / er[{iCell.x, 0}].dcEdge
  end
--  for index in cell_range_1d do
--    var i = int(index)
--    var iCell1 = er[{i, 0}].cellsOnEdge[0] --cellsOnEdge(1,i)
--    var iCell2 = er[{i, 0}].cellsOnEdge[1] --cellsOnEdge(2,i)
--    for k=1, nz1 do -- TODO: should this be 0?
--      er[{i, k}].zxu = 0.5 * (cr[{iCell2, k}].zgrid - cr[{iCell1, k}].zgrid + cr[{iCell2, k + 1}].zgrid - cr[{iCell1, k + 1}].zgrid) / er[{i, 0}].dcEdge
--    end
--  end
  for iCell in cell_range_2d do
    var ztemp = .5 * (cr[{iCell.x, iCell.y + 1}].zgrid + cr[iCell].zgrid)
    cr[iCell].dss = 0.0
    ztemp = cr[iCell].zgrid -- TODO: This looks wrong. But, it's identical to the Fortran code.
    if (ztemp > zd + .1) then
       cr[iCell].dss = cr[iCell].dss + xnutr * pow(sin( .5 * pii * (ztemp - zd) / (zt - zd)), 2)
    end
  end





--!**************  section for 2d (z,lat) calc for zonal velocity

  var dlat = 0.5 * pii / float(nlat - 1)
  var lat_2d : double[nlat]
  var zgrid_2d : double[(nVertLevels + 1) * nlat]
  var zz_2d : double[nVertLevels * nlat]
  var ppb_2d : double[nVertLevels * nlat]
  var pb_2d : double[nVertLevels * nlat]
  var rb_2d : double[nVertLevels * nlat]
  var tb_2d : double[nVertLevels * nlat]
  var etavs_2d : double[nVertLevels * nlat]
  var u_2d : double[nVertLevels * nlat]
  var rho_2d : double[nVertLevels * nlat]
  var p_2d : double[nVertLevels * nlat]
  var pp_2d : double[nVertLevels * nlat]
  var rr_2d : double[nVertLevels * nlat]
  var rtb_2d : double[nVertLevels * nlat]


  -- TODO: This loops seems to accomplish nothing.
  for iNlat in nlat_range_1d do
    var i = int(iNlat)
    -- Moved declarations of helper variables inside of for loop to avoid loop-carried dependencies.
    var eta : double[nVertLevels]
    var teta : double[nVertLevels]
    var tt : double[nVertLevels]
    var etav : double[nVertLevels]
    var temperature_1d : double[nVertLevels]
    var ptemp : double
    var ztemp : double
    var phi : double

    lat_2d[i] = float(i-1) * dlat
    phi = lat_2d[i]
    var hx_1d = u0 / gravity * pow(cos(etavs), 1.5) * ((-2.0 * pow(sin(phi), 6) * (pow(cos(phi), 2) + 1.0 / 3.0) + 10.0 / 63.0) 
                * (u0) * pow(cos(etavs), 1.5) + (1.6 * pow(cos(phi), 3) * (pow(sin(phi), 2) + 2.0 / 3.0) - pii / 4.0) * r_earth * omega_e)
    for k=0, nz do
      zgrid_2d[k * nlat + i] = (1. - ah[k]) * (sh[k] * (zt - hx_1d) + hx_1d) + ah[k] * sh[k]* zt
    end
    for k=0, nz1 do
      zz_2d[k * nlat + i] = (zw[k + 1] - zw[k]) / (zgrid_2d[(k + 1) * nlat + i] - zgrid_2d[k * nlat + i])
    end

    for k=1,nz1 do
      ztemp = .5 * (zgrid_2d[(k + 1) * nlat + i] + zgrid_2d[k * nlat + i])
      ppb_2d[k * nlat + i] = p0 * exp(-gravity * ztemp / (rgas * t0b))
      pb_2d[k * nlat + i] = pow((ppb_2d[k * nlat + i] / p0), (rgas / cp))
      rb_2d[k * nlat + i] = ppb_2d[k * nlat + i] / (rgas * t0b * zz_2d[k * nlat + i])
      tb_2d[k * nlat + i] = t0b / pb_2d[k * nlat + i]
      rtb_2d[k * nlat + i] = rb_2d[k * nlat + i] * tb_2d[k * nlat + i]
      p_2d[k * nlat + i] = pb_2d[k * nlat + i]
      pp_2d[k * nlat + i] = 0.0
      rr_2d[k * nlat + i] = 0.0
    end


    for itr = 0, 10 do

      for k=0, nz1 do
        eta[k] = (ppb_2d[k * nlat + i] + pp_2d[k * nlat + i]) / p0
        etav[k] = (eta[k] - .252) * pii / 2.0
        if(eta[k] >= znut)  then
          teta[k] = t0 * pow(eta[k], (rgas * dtdz / gravity))
        else
          teta[k] = t0 * pow(eta[k], (rgas * dtdz / gravity)) + delta_t * pow((znut - eta[k]), 5)
        end
      end

      phi = lat_2d[i]
      for k=1, nz1 do
        temperature_1d[k] = teta[k] + .75 * eta[k] * pii * u0 / rgas * sin(etav[k]) * sqrt(cos(etav[k])) * ((-2. * pow(sin(phi), 6) * (pow(cos(phi), 2) + 1.0 / 3.0) + 10.0 / 63.0) 
                            * 2.0 * u0 * pow(cos(etav[k]), 1.5) + (1.6 * pow(cos(phi), 3) * (pow(sin(phi), 2) + 2.0 / 3.0) - pii / 4.0) * r_earth * omega_e) / (1.0 + 0.61 * qv_2d[nlat * k + i])

        -- Only used in skipped conditional below. Would change to var ... = ... for cuda generation.
        -- ztemp = .5 * (zgrid_2d[k * nlat + i] + zgrid_2d[(k + 1) * nlat + i])
        -- ptemp = ppb_2d[k * nlat + i] + pp_2d[k * nlat + i]

        --get moisture
        ----SKIPPING THIS CONDITIONAL FOR NOW----
        --if (moisture) then
          --qv_2d[k * nlat + i] = env_qv( ztemp, temperature_1d[k], ptemp, rh_max )
        --end

        tt[k] = temperature_1d[k] * (1.0 + 1.61 * qv_2d[k * nlat + i])
      end

      for itrp = 0,25 do
        for k=0,nz1 do
          rr_2d[k * nlat + i]  = (pp_2d[k * nlat + i] / (rgas * zz_2d[k * nlat + i]) - rb_2d[k * nlat + i] * (tt[k] - t0b)) / tt[k]
        end

        var ppi : double[nVertLevels]
        -- TODO: I believe this should be ppi[0].
        ppi[1] = p0 - .5 * dzw[1] * gravity * (1.25 * (rr_2d[1 * nlat + i] + rb_2d[1 * nlat + i]) * (1.0 + qv_2d[1 * nlat + i]) - .25* (rr_2d[2 * nlat + i] + rb_2d[2 * nlat + i]) * (1.0 + qv_2d[2 * nlat + i]))

        ppi[1] = ppi[1] - ppb_2d[1 * nlat + i]

        for k=0, nz1-1 do
          -- TODO: Can't access vertr here.
--          ppi[k + 1] = ppi[k] - vertr[k+1].dzu * gravity * ((rr_2d[k * nlat + i] + (rr_2d[k * nlat + i] + rb_2d[k * nlat + i]) * qv_2d[k * nlat + i]) * vertr[k+1].fzp + (rr_2d[(k + 1) * nlat + i] + (rr_2d[(k + 1) * nlat + i] + rb_2d[(k + 1) * nlat + i]) * qv_2d[(k + 1) * nlat + i]) * vertr[k + 1].fzm)
        end

        for k=0, nz1 do
          pp_2d[k * nlat + i] = .2 * ppi[k] + .8 * pp_2d[k * nlat + i]
        end
--
      end   -- end inner iteration loop itrp

    end   -- end outer iteration loop itr
--
--
    for k = 0, nz1 do
      rho_2d[k * nlat + i] = rr_2d[k * nlat + i] + rb_2d[k * nlat + i]
      etavs_2d[k * nlat + i] = ((ppb_2d[k * nlat + i] + pp_2d[k * nlat + i]) / p0 - 0.252) * pii / 2.0
      u_2d[k * nlat + i] = u0 * (pow(sin(2. * lat_2d[i]), 2)) * (pow(cos(etavs_2d[k * nlat + i]), 1.5))
    end

  end   -- end loop over latitudes for 2D zonal wind field calc


------------TODO: skipping this conditional for now ---------------


--
--      !SHP-balance:: in case of rebalacing for geostrophic wind component
--      if (rebalance) then
--
--        do i=1,nlat-1
--          do k=1,nz1
--            zx_2d(k,i) = (zgrid_2d(k,i+1)-zgrid_2d(k,i))/(dlat*r_earth)
--          end do
--        end do
--
--        call init_atm_recompute_geostrophic_wind(u_2d, rho_2d, pp_2d, qv_2d, lat_2d, zz_2d, zx_2d,     &
--                                        cf1, cf2, cf3, fzm, fzp, rdzw, nz1, nlat, dlat, sphere_radius)
--
--      end if
--
--!******************************************************************
--




--!******************************************************************

--!
--!---- baroclinc wave initialization ---------------------------------
--!
--!     reference sounding based on dry isothermal atmosphere
--!
  for iCell in cell_range_1d do
    var i = int(iCell)
    -- Moved declarations of helper variables inside for loop to avoid loop-carried dependencies.
    var eta : double[nVertLevels]
    var teta : double[nVertLevels]
    var tt : double[nVertLevels]
    var etav : double[nVertLevels]
    var temperature_1d : double[nVertLevels]
    var ptemp : double
    var phi : double

    for k=0, nz1 do
      -- TODO: Seems incorrect! 
      -- Fortran: ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
      var ztemp = .5 * (cr[{i, k + 1}].zgrid + cr[{i, k}].zgrid)
      cr[{i, k}].pressure_base = p0 * exp(-gravity * ztemp / (rgas * t0b))
      cr[{i, k}].pressure_p = pow((cr[{i, k}].pressure_base / p0), (rgas / cp))
      cr[{i, k}].rho_base = cr[{i, k}].pressure_base / (rgas * t0b * cr[{i, k}].zz)
      cr[{i, k}].theta_base = t0b / cr[{i, k}].pressure_p
      cr[{i, k}].rtheta_base = cr[{i, k}].rho_base * cr[{i, k}].theta_base
      cr[{i, k}].exner = cr[{i, k}].pressure_p
      cr[{i, k}].pressure_p = 0.0
      cr[{i, k}].rho_p = 0.0
    end

--!     iterations to converge temperature as a function of pressure
--!
    for itr = 0, 10 do

      for k=0, nz1 do
        eta[k] = (cr[{i, k}].pressure_base + cr[{i, k}].pressure_p) / p0
        etav[k] = (eta[k] - .252) * pii / 2.
        if (eta[k] >= znut) then
          teta[k] = t0 * pow(eta[k], (rgas * dtdz / gravity))
        else
          teta[k] = t0 * pow(eta[k], (rgas * dtdz / gravity)) + delta_t * pow((znut - eta[k]), 5)
        end
      end
      phi = cr[{i, 0}].lat
      for k=0,nz1 do
        temperature_1d[k] = teta[k] + .75 * eta[k] * pii * u0 / rgas * sin(etav[k]) * sqrt(cos(etav[k])) * ((-2.0 * pow(sin(phi), 6) * (pow(cos(phi), 2) + 1.0 / 3.0) + 10.0 / 63.0) * 2. * u0 * pow(cos(etav[k]), 1.5) + (1.6 * pow(cos(phi), 3) * (pow(sin(phi), 2) + 2.0 / 3.0) - pii / 4.0) * r_earth * omega_e) / (1. + 0.61 * cr[{i, k}].qv)
        
        -- Only needed for conditional below.
        --ztemp = .5 * (cr[{k,i}].zgrid + cr[{i, k + 1}].zgrid)
        --ptemp = cr[{i, k}].pressure_base + cr[{i, k}].pressure_p

------SKIPPING BECAUSE CONDITIONAL ----------
--  --!get moisture
--        if (moisture) then

--            !scalars(index_qv,k,i) = env_qv( ztemp, temperature_1d[k], ptemp, rh_max )
--
--          if(ptemp < 50000.) then
--              relhum(k,i) = 0.0
--           elseif(ptemp > p0) then
--              relhum(k,i) = 1.0
--           else
--              relhum(k,i) = (1.-((p0-ptemp)/50000.)**1.25)
--           end
--           relhum(k,i) = min(rh_max,relhum(k,i))
--
--           !.. calculation of water vapor mixing ratio:
--           if (temperature_1d[k] > 273.15) then
--               es  = 1000.*0.6112*cmath.exp(17.67*(temperature_1d[k]-273.15)/(temperature_1d[k]-29.65))
--           else
--               es  = 1000.*0.6112*cmath.exp(21.8745584*(temperature_1d[k]-273.15)/(temperature_1d[k]-7.66))
--           end
--           qsat(k,i) = (287.04/461.6)*es/(ptemp-es)
--           if(relhum(k,i) .eq. 0.0) then qsat(k,i) = 0.0
--           end
--           scalars(index_qv,k,i) = relhum(k,i)*qsat(k,i)
--        end

        tt[k] = temperature_1d[k] * (1. + 1.61 * cr[{i, k}].qv)

      end




      for itrp = 0,25 do
        for k=0, nz1 do
          cr[{i, k}].rho_p = (cr[{i, k}].pressure_p / (rgas * cr[{i, k}].zz) - cr[{i, k}].rho_base * (tt[k] - t0b)) / tt[k]
        end

        var ppi : double[nVertLevels]
        -- TODO: should be 0.
        ppi[1] = p0 - .5 * dzw[1] * gravity * (1.25 * (cr[{i, 1}].rho_p + cr[{i, 1}].rho_base) * (1. + cr[{i, 0}].qv ) - .25 * (cr[{i, 2}].rho_p + cr[{i, 2}].rho_base) * (1. + cr[{i, 1}].qv))

        ppi[1] = ppi[1] - cr[{i, 1}].pressure_base
        for k=0, nz1-1 do
           ppi[k + 1] = ppi[k] - vertr[k + 1].dzu * gravity * ( (cr[{i, k}].rho_p + (cr[{i, k}].rho_p + cr[{i, k}].rho_base) 
                        * cr[{i, k}].qv) * vertr[k + 1].fzp + (cr[{i, k + 1}].rho_p + (cr[{i, k + 1}].rho_p + cr[{i, k + 1}].rho_base) * cr[{i, k + 1}].qv) * vertr[k + 1].fzm)
        end

        for k=0, nz1 do
          cr[{i, k}].pressure_p = .2 * ppi[k] + .8 * cr[{i, k}].pressure_p
        end

      end   -- end inner iteration loop itrp

    end   -- end outer iteration loop itr



    for k=0, nz1 do
      cr[{i, k}].exner = pow(((cr[{i, k}].pressure_base + cr[{i, k}].pressure_p) / p0), (rgas / cp))
      cr[{i, k}].theta_m = tt[k] / cr[{i, k}].exner
      cr[{i, k}].rtheta_p = cr[{i, k}].theta_m * cr[{i,k}].rho_p + cr[{i, k}].rho_base * (cr[{i, k}].theta_m - cr[{i, k}].theta_base)
      cr[{i, k}].rho_zz = cr[{i, k}].rho_base + cr[{i,k}].rho_p
    end

    --calculation of surface pressure:
    cr[{i, 0}].surface_pressure = 0.5 * dzw[1] * gravity * (1.25 * (cr[{i, 1}].rho_p + cr[{i, 1}].rho_base) * (1.0 + cr[{i, 0}].qv) - 0.25 * (cr[{i, 2}].rho_p + cr[{i, 2}].rho_base) * (1.0 + cr[{i, 1}].qv))
    cr[{i, 0}].surface_pressure = cr[{i, 0}].surface_pressure + cr[{i, 1}].pressure_p + cr[{i, 1}].pressure_base

  end   -- end loop over cells

  var lat_pert = latitude_pert * pii / 180.0
  var lon_pert = longitude_pert * pii / 180.0



  for i in edge_range_1d do
    var iEdge = int(i)
    var u_pert : double

    var vtx1 = er[{iEdge, 0}].verticesOnEdge[0] --verticesOnEdge(1,iEdge)
    var vtx2 = er[{iEdge, 0}].verticesOnEdge[1]
    var lat1 = vr[{vtx1, 0}].lat
    var lat2 = vr[{vtx2, 0}].lat
    var iCell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var iCell2 = er[{iEdge, 0}].cellsOnEdge[1]
    var flux = (0.5 * (lat2 - lat1) - 0.125 * (sin(4. * lat2) - sin(4. * lat1))) * sphere_radius / er[{iEdge, 0}].dvEdge

--     if (config_init_case == 2) then
--        r_pert = sphere_distance( er[{iEdge, 0}].lat, er[{iEdge, 0}].lon, &
--                                  lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
--        u_pert = u_perturbation*cmath.exp(-r_pert**2)*(lat2-lat1) * sphere_radius / dvEdge(iEdge)
--
--     elseif (config_init_case == 3) then
--        lon_Edge = er[{iEdge, 0}].lon
--        u_pert = u_perturbation*cmath.cos(k_x*(lon_Edge - lon_pert)) &
--                     *(0.5*(lat2-lat1) - 0.125*(cmath.sin(4.*lat2) - cmath.sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)
--     else
--        u_pert = 0.0
--     end

----SKIPPING THE CONDITIONAL above, defaulting to the ELSE case
       u_pert = 0.0
-------end conditional else----------

----SKIPPING THE CONDITIONAL, ONLY PORTING THE ELSE STATEMENT FOR NOW--------
--     if (rebalance) then
--
--       call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,dvEdge(iEdge),sphere_radius,u0,nz1,nlat)
--       do k=1,nVertLevels
--         fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
--         u(k,iEdge) = fluxk + u_pert
--       end do

--     else

--       do k=1,nVertLevels
--         etavs = (0.5*(pressure_base(k,iCell1)+pressure_base(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
--         fluxk = u0*flux*(cmath.cos(etavs)**1.5)
--         u(k,iEdge) = fluxk + u_pert
--       end do

--     end if


    for k=0, nVertLevels do
      var etavs = (0.5 * (cr[{iCell1, k}].pressure_base + cr[{iCell2, k}].pressure_base + cr[{iCell1, k}].pressure_p + cr[{iCell2, k}].pressure_p) / p0 - 0.252) * pii / 2.0
      var fluxk = u0 * flux * (pow(cos(etavs), 1.5))
      er[{iEdge, k}].u = fluxk + u_pert
    end

------end conditional else---------

    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
    for k=0, nz1 do
       er[{iEdge, k}].ru = 0.5 * (cr[{cell1, k}].rho_zz + cr[{cell2, k}].rho_zz) * er[{iEdge, k}].u
    end

--  !
--  ! Generate rotated Coriolis field
--  !

    er[{iEdge, 0}].fEdge = 2.0 * omega_e * ( -1.0 * cos(er[{iEdge, 0}].lon) * cos(er[{iEdge, 0}].lat) * sin(alpha_grid) + sin(er[{iEdge, 0}].lat) * cos(alpha_grid) )
  end


  for iVtx=0, nVertices do
    vr[{iVtx, 0}].fVertex = 2.0 * omega_e * (-1.0 * cos(vr[{iVtx, 0}].lon) * cos(vr[{iVtx, 0}].lat) * sin(alpha_grid) + sin(vr[{iVtx, 0}].lat) * cos(alpha_grid))
  end



--  !
--  !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
--  !
--
--  !
--  !     pre-calculation z-metric terms in omega eqn.
--  !



  for i in edge_range_1d do
    var iEdge = int(i)
    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
    -- Moved z_edges declaration inside because of Cuda.
    var z_edge : double
    var z_edge3 : double

    --Below conditional is unnecessary because we use regions, not pools/threads
     --! Avoid a potential divide by zero below if areaCell(nCells+1) is used in the denominator
    -- if (cell1 <= nCellsSolve or cell2 <= nCellsSolve ) then
        for k = 0, nVertLevels do
          z_edge = (cr[{cell1, k}].zgrid + cr[{cell2, k}].zgrid) / 2.0
          -- if (config_theta_adv_order == 2) then

            --  z_edge = (cr[{cell1, k}].zgrid+cr[{cell2, k}].zgrid)/2.

          -- elseif (config_theta_adv_order == 3 or config_theta_adv_order ==4) then --theta_adv_order == 3 or 4

          --    d2fdx2_cell1 = deriv_two(1,1,iEdge) * cr[{cell1, k}].zgrid
          --    d2fdx2_cell2 = deriv_two(1,2,iEdge) * cr[{cell2, k}].zgrid

--  WCS fix 20120711

          --    for i=0, nEdgesOnCell(cell1) do
          --       if ( cellsOnCell(i,cell1) > 0) then d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
          --       end
          --    end
          --    for i=0, nEdgesOnCell(cell2) do
          --       if ( cellsOnCell(i,cell2) > 0) then  d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
          --       end
          --    end

          --    z_edge =  0.5*(cr[{cell1, k}].zgrid + cr[{cell2, k}].zgrid)         &
          --                  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

          --    if (config_theta_adv_order == 3) then
          --       z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
          --    else
          --       z_edge3 = 0.
          --    end
              -- only using the ELSE default for now
               z_edge3 = 0.0
     --      end

           er[{iEdge, k}].zb[0] = (z_edge - cr[{cell1, k}].zgrid) * er[{iEdge, 0}].dvEdge / cr[{cell1, 0}].areaCell
           er[{iEdge, k}].zb[1] = (z_edge - cr[{cell2, k}].zgrid) * er[{iEdge, 0}].dvEdge / cr[{cell2, 0}].areaCell
           er[{iEdge, k}].zb3[0] = z_edge3 * er[{iEdge, 0}].dvEdge / cr[{cell1, 0}].areaCell
           er[{iEdge, k}].zb3[1] = z_edge3 * er[{iEdge, 0}].dvEdge / cr[{cell2, 0}].areaCell

        end
    -- end

  end

  --! for including terrain
  for i in edge_range_1d do
    var iEdge = int(i)

    var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
    var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

-------------------------------------------------------------------------------
----------TODO: past this point, conversion from Fortran is incomplete---------
-------------------------------------------------------------------------------

    for k = 1, nVertLevels do
       var flux = (vertr[k].fzm * er[{iEdge, k}].ru + vertr[k].fzp * er[{iEdge, k - 1}].ru)
       cr[{cell2, k}].rw = cr[{cell2, k}].rw + (vertr[k].fzm * cr[{cell2, k}].zz + vertr[k].fzp * cr[{cell2, k-1}].zz) * er[{iEdge, k}].zb[1] * flux
       cr[{cell1, k}].rw = cr[{cell1, k}].rw - (vertr[k].fzm * cr[{cell1, k}].zz + vertr[k].fzp * cr[{cell1, k-1}].zz) * er[{iEdge, k}].zb[0] * flux

--        if (config_theta_adv_order ==3) then
--           cr[{cell2, k}].rw = cr[{cell2, k}].rw    &
--                                        - sign(1.0_RKIND,er[{iEdge, k}].ru)*config_coef_3rd_order* &
--                                          (vertr[k].fzm*zz(k,cell2)+vertr[k].fzp*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
--           cr[{cell1, k}].rw = cr[{cell1, k}].rw    &
--                                        + sign(1.0_RKIND,er[{iEdge, k}].ru)*config_coef_3rd_order* &
--                                          (vertr[k].fzm*zz(k,cell1)+vertr[k].fzp*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
--        end

    end
  end

  --! Compute w from rho_zz and rw
  var cell_range_2d_special = rect2d { int2d {0, 1}, int2d {nCells - 1, nVertLevels - 1} }
  for iCell in cell_range_2d_special do
    cr[iCell].w = cr[iCell].rw / (vertr[iCell.y].fzp * cr[{iCell.x, iCell.y - 1}].rho_zz + vertr[iCell.y].fzm * cr[iCell].rho_zz)
  end


  --!
  --! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
  --!
  for iEdge in edge_range_2d do
    er[iEdge].v = 0.0
  end

  for i in edge_range_1d do
    var iEdge = int(i)
    for i = 0, er[{iEdge, 0}].nEdgesOnEdge do
      var eoe = er[{iEdge, 0}].edgesOnEdge_ECP[i]
      for k = 0, nVertLevels do
        er[{iEdge, k}].v = er[{iEdge, k}].v + er[{iEdge, 0}].weightsOnEdge[i] * er[{eoe, k}].u
      end
    end
  end

--------PSURF is never used, in mpas it only exists for purposes of logging. excluding because of this
--  for i=0,10 do
--    psurf = (vertr[0].cf1*(pressure_base(1,i)+pp(1,i)) + vertr[0].cf2*(pressure_base(2,i)+pp(2,i)) + vertr[0].cf3*(pressure_base(3,i)+pp(3,i)))/100.
--
--        psurf = (pressure_base(1,i)+pp(1,i)) + .5*dzw[1]*gravity        &
--                      *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
--                       -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))
--
--  end

  --! Compute rho and theta from rho_zz and theta_m
  for iCell in cell_range_2d do
      cr[iCell].rho = cr[iCell].rho_zz * cr[iCell].zz
      cr[iCell].theta = cr[iCell].theta_m / (1.0 + 1.61 * cr[iCell].qv)
  end
end
