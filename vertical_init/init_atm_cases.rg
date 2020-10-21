import "regent"
require "data_structures"

local constants = require("constants")

local nCells = 2562
local nEdges = 7680
local nVertices = 5120
local maxEdges = 10
local maxEdges2 = 20
local TWO = 2
local FIFTEEN = 15
local vertexDegree = 3
local nVertLevels = 1
local sphere_radius = 6371229.0 --from ncdump of the grid, default is 1.0, but set to "a" from constants.F in init_atm_core.F
local nlat = 721

local cio = terralib.includec("stdio.h")
local cmath = terralib.includec("math.h")


--__demand(__cuda)
task init_atm_case_jw(cr : region(ispace(int2d), cell_fs),
                      er : region(ispace(int2d), edge_fs),
                      vr : region(ispace(int2d), vertex_fs),
                      vertr : region(ispace(int1d), vertical_fs))
where reads (cr.lat, er.cellsOnEdge, er.edgesOnEdge_ECP, er.lat, er.lon, er.nEdgesOnEdge, er.verticesOnEdge, er.weightsOnEdge, vr.lat, vr.lon),
writes (cr.qsat, cr.relhum, cr.rho, cr.rtheta_base, cr.rtheta_p, cr.surface_pressure, cr.theta, cr.w, er.fEdge, er.zxu, vr.fVertex, vertr.rdzu, vertr.rdzw),
reads writes (cr.areaCell, cr.dss, cr.exner, cr.hx, cr.pressure_base, cr.pressure_p, cr.qv, cr.rho_base, cr.rho_p, cr.rho_zz, cr.rw, cr.surface_pressure, cr.theta_base, cr.theta_m, cr.x, cr.y, cr.z, cr.zgrid, cr.zz, er.dvEdge, er.dcEdge, er.ru, er.u, er.v, er.x, er.y, er.z, er.zb, er.zb3, vr.areaTriangle, vr.kiteAreasOnVertex, vr.x, vr.y, vr.z, vertr.dzu, vertr.fzm, vertr.fzp) do

-- local vars/constants defined at beginning of subroutine
  var cp = constants.cp
  var rgas = constants.rgas
  var gravity = constants.gravity
  var pii = constants.pii
  var u0 = 35.0
  var alpha_grid = 0.0 -- no grid rotation
  var omega_e : double
  var flux : double
  var u_pert : double
  var phi : double
  var eoe : int
  var ztemp : double
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

  for iCell = 0, nCells do
    cr[{iCell, 0}].x = cr[{iCell, 0}].x * sphere_radius
    cr[{iCell, 0}].y = cr[{iCell, 0}].y * sphere_radius
    cr[{iCell, 0}].z = cr[{iCell, 0}].z * sphere_radius
    cr[{iCell, 0}].areaCell = cr[{iCell, 0}].areaCell * cmath.pow(sphere_radius, 2.0)
  end

  for iVert = 0, nVertices do
    vr[{iVert, 0}].x = vr[{iVert, 0}].x * sphere_radius
    vr[{iVert, 0}].y = vr[{iVert, 0}].y * sphere_radius
    vr[{iVert, 0}].z = vr[{iVert, 0}].z * sphere_radius
    vr[{iVert, 0}].areaTriangle = vr[{iVert, 0}].areaTriangle * cmath.pow(sphere_radius, 2.0)
    for vDeg = 0, vertexDegree do
      vr[{iVert, 0}].kiteAreasOnVertex[vertexDegree] = vr[{iVert, 0}].kiteAreasOnVertex[vertexDegree] * cmath.pow(sphere_radius, 2.0)
    end
  end

  for iEdge = 0, nEdges do
    er[{iEdge, 0}].x = er[{iEdge, 0}].x * sphere_radius
    er[{iEdge, 0}].y = er[{iEdge, 0}].y * sphere_radius
    er[{iEdge, 0}].z = er[{iEdge, 0}].z * sphere_radius
    er[{iEdge, 0}].dvEdge = er[{iEdge, 0}].dvEdge * sphere_radius
    er[{iEdge, 0}].dcEdge = er[{iEdge, 0}].dcEdge * sphere_radius

  end




-- initialization of moisture:
  for iCell = 0, nCells do
    for k = 0, nVertLevels do
      cr[{iCell, k}].qv = 0.0
      cr[{iCell, k}].qsat = 0.0
      cr[{iCell, k}].relhum = 0.0
    end
  end

  for i = 0, nVertLevels do
    for j = 0, nlat do
      qv_2d[i*nlat + j] = 0.0
    end
  end

-- end initialization of moisture.


  for iCell = 0, nCells do
    cr[{iCell, 0}].surface_pressure = 0.0
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

  for iCell=0, nCells do
    for k = 0, nz do
      phi = cr[{iCell, 0}].lat
      cr[{iCell, k}].hx = u0 / gravity * cmath.pow(cmath.cos(etavs), 1.5) * ((-2.0 * cmath.pow(cmath.sin(phi), 6) * (cmath.pow(cmath.cos(phi), 2.0) + 1.0/3.0) + 10.0/63.0) * u0*cmath.pow(cmath.cos(etavs), 1.5) + (1.6 * cmath.pow(cmath.cos(phi), 3) * (cmath.pow(cmath.sin(phi), 2) + 2.0/3.0) - pii/4.0)*r_earth*omega_e)
    end
  end


--      !     Metrics for hybrid coordinate and vertical stretching

  var str = 1.5
  var zt = 45000.0
  var dz = zt/nz1
  var sh : double[nVertLevels + 1]
  var zw : double[nVertLevels + 1]
  var ah : double[nVertLevels + 1]
  var dzw : double[nVertLevels]
  var zu : double[nVertLevels]
  var rdzwp : double[nVertLevels]
  var rdzwm : double[nVertLevels]

  for k=0, nz do
    if k == 0 then
      sh[k] = -1
    else
  --!  sh[k] is the stretching specified for height surfaces
      sh[k] = cmath.pow(([double](k - 1.0) * dz / zt), str) --- this was (real(k-1)) in mpas; do we need to cast it here too?
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

    zw[k] = (k-1)*dz -- in mpas they cast to float

--!            zw[k] = sh[k]*zt --- see above comments for which version you want
--!
--!           ah[k] governs the transition between terrain-following
--!           and pureheight coordinates
--!                ah[k] = 0 is a terrain-following coordinate
--!                ah[k] = 1 is a height coordinate
--
    ah[k] = 1.0 - cmath.pow(cmath.cos(.5*pii*(k-1)*dz/zt), 6.0)
--!            ah[k] = 0.

    --cio.printf("sh[%d] is %f \n", k, sh[k])
    --cio.printf("ah[%d] is %f \n", k, ah[k])
    --cio.printf("zw[%d] is %f \n", k, zw[k])
  end
  for k=0, nz1 do -- nz1 is just nVertLevels, idk why mpas renamed it
    dzw[k] = zw[k+1]-zw[k]
    vertr[k].rdzw = 1.0/dzw[k]
    zu[k] = .5*(zw[k]+zw[k+1])
  end

  for k=1, nz1 do -- k=2,nz1 in mpas
    vertr[k].dzu = .5*(dzw[k]+dzw[k-1])
    vertr[k].rdzu  =  1.0/vertr[k].dzu
    vertr[k].fzp = .5* dzw[k]/vertr[k].dzu
    vertr[k].fzm =.5* dzw[k-1]/vertr[k].dzu
    rdzwp[k] = dzw[k-1]/(dzw[k]*(dzw[k]+dzw[k-1]))
    rdzwm[k] = dzw[k]/(dzw[k-1]*(dzw[k]+dzw[k-1]))
  end


--!**********  how are we storing cf1, cf2 and cf3?

  var COF1 = (2.*vertr[1].dzu+vertr[2].dzu)/(vertr[1].dzu + vertr[2].dzu) * dzw[0]/ vertr[1].dzu
  var COF2 = vertr[1].dzu / (vertr[1].dzu + vertr[2].dzu)*dzw[0]/ vertr[2].dzu
  var CF1  = vertr[1].fzp + COF1
  var CF2  = vertr[1].fzm - COF1 - COF2
  var CF3  = COF2




  for iCell=0, nCells do
    for k=0, nz do
        cr[{iCell,k}].zgrid = (1.0-ah[k])*(sh[k]*(zt-cr[{iCell, k}].hx)+cr[{iCell, k}].hx) + ah[k] * sh[k]* zt
        --cio.printf("ah[%d] is %f \n", k, ah[k])
        --cio.printf("sh[%d] is %f \n", k, sh[k])
        --cio.printf("cr[{%d, %d}].hx is %f\n", iCell, k, cr[{iCell, k}].hx)
    end
    for k=0, nz1 do
      cr[{iCell, k}].zz = (zw[k+1]-zw[k])/(cr[{iCell, k+1}].zgrid-cr[{iCell,k}].zgrid)
      --cio.printf("cr[{iCell = %d, k = %d}].zz is %f \n", iCell, k, cr[{iCell, k}].zz)
      --cio.printf("zw[%d] is %f \n", k, zw[k]) : zw is set
      --cio.printf("cr[{%d, %d}].zgrid is %f", iCell, k, cr[{iCell,k}].zgrid) : zgrid is not set
    end
  end

  for i=0, nEdges do
    var iCell1 = er[{i, 0}].cellsOnEdge[0] --cellsOnEdge(1,i)
    var iCell2 = er[{i, 0}].cellsOnEdge[1] --cellsOnEdge(2,i)
    for k=1,nz1 do
      er[{i, k}].zxu = 0.5 * (cr[{iCell2, k}].zgrid-cr[{iCell1, k}].zgrid + cr[{iCell2, k+1}].zgrid-cr[{iCell1, k+1}].zgrid) / er[{i, 0}].dcEdge
    end
  end
  for i=0, nCells do
    for k=0, nz1 do
      ztemp = .5*(cr[{i, k+1}].zgrid+cr[{k,i}].zgrid)
      cr[{i, k}].dss = 0.0
      ztemp = cr[{k,i}].zgrid
      if(ztemp > zd+.1)  then
         cr[{i, k}].dss = cr[{i, k}].dss+xnutr*cmath.pow(cmath.sin(.5*pii*(ztemp-zd)/(zt-zd)), 2)
      end
    end
  end




--!**************  section for 2d (z,lat) calc for zonal velocity

  var dlat = 0.5*pii / float(nlat-1)
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
  var eta : double[nVertLevels]
  var ppi : double[nVertLevels]
  var teta : double[nVertLevels]
  var tt : double[nVertLevels]
  var etav : double[nVertLevels]
  var temperature_1d : double[nVertLevels]
  var ptemp : double
  for i = 0,  nlat do

    lat_2d[i] = float(i-1)*dlat
    phi = lat_2d[i]
    var hx_1d = u0 / gravity * cmath.pow(cmath.cos(etavs),1.5) * ((-2.0 * cmath.pow(cmath.sin(phi), 6) * (cmath.pow(cmath.cos(phi),2)+1.0/3.0)+10.0/63.0) *(u0)*cmath.pow(cmath.cos(etavs),1.5) +(1.6*cmath.pow(cmath.cos(phi),3) *(cmath.pow(cmath.sin(phi),2)+2.0/3.0)-pii/4.0)*r_earth*omega_e)
    for k=0, nz do
      zgrid_2d[k * nlat + i] = (1.-ah[k])*(sh[k]*(zt-hx_1d)+hx_1d) + ah[k] * sh[k]* zt
    end
    for k=0, nz1 do
      zz_2d[k* nlat + i] = (zw[k+1]-zw[k])/(zgrid_2d[(k+1) * nlat + i]-zgrid_2d[k*nlat + i])
    end

    for k=1,nz1 do
      ztemp = .5*(zgrid_2d[(k+1)*nlat + i]+zgrid_2d[k*nlat + i])
      ppb_2d[k * nlat + i] = p0*cmath.exp(-gravity*ztemp/(rgas*t0b))
      pb_2d[k * nlat + i] = cmath.pow((ppb_2d[k * nlat + i]/p0),(rgas/cp))
      rb_2d[k * nlat + i] = ppb_2d[k * nlat + i]/(rgas*t0b*zz_2d[k * nlat + i])
      tb_2d[k * nlat + i] = t0b/pb_2d[k * nlat + i]
      rtb_2d[k * nlat + i] = rb_2d[k * nlat + i]*tb_2d[k * nlat + i]
      p_2d[k * nlat + i] = pb_2d[k * nlat + i]
      pp_2d[k * nlat + i] = 0.0
      rr_2d[k * nlat + i] = 0.0
    end


    for itr = 0,10 do

      for k=0,nz1 do
        eta[k] = (ppb_2d[k * nlat + i]+pp_2d[k * nlat + i])/p0
        etav[k] = (eta[k]-.252)*pii/2.0
        if(eta[k] >= znut)  then
          teta[k] = t0*cmath.pow(eta[k],(rgas*dtdz/gravity))
        else
          teta[k] = t0*cmath.pow(eta[k],(rgas*dtdz/gravity)) + delta_t*cmath.pow((znut-eta[k]),5)
        end
      end

      phi = lat_2d[i]
      for k=1,nz1 do
        temperature_1d[k] = teta[k]+.75*eta[k]*pii*u0/rgas*cmath.sin(etav[k])  *cmath.sqrt(cmath.cos(etav[k]))* ((-2.*cmath.pow(cmath.sin(phi),6)  *(cmath.pow(cmath.cos(phi),2)+1.0/3.0)+10.0/63.0)  *2.0*u0*cmath.pow(cmath.cos(etav[k]),1.5) +(1.6*cmath.pow(cmath.cos(phi),3) *(cmath.pow(cmath.sin(phi),2)+2.0/3.0)-pii/4.0)*r_earth*omega_e)/(1.0+0.61*qv_2d[nlat*k + i])

        ztemp   = .5*(zgrid_2d[k * nlat + i]+zgrid_2d[(k+1) * nlat + i])
        ptemp   = ppb_2d[k * nlat + i] + pp_2d[k * nlat + i]

        --get moisture
        ----SKIPPING THIS CONDITIONAL FOR NOW----
        --if (moisture) then
          --qv_2d[k * nlat + i] = env_qv( ztemp, temperature_1d[k], ptemp, rh_max )
        --end

        tt[k] = temperature_1d[k]*(1.0+1.61*qv_2d[k * nlat + i])
      end

      for itrp = 0,25 do
        for k=0,nz1 do
          rr_2d[k * nlat + i]  = (pp_2d[k * nlat + i]/(rgas*zz_2d[k * nlat + i]) - rb_2d[k * nlat + i]*(tt[k]-t0b))/tt[k]
        end

        ppi[1] = p0-.5*dzw[1]*gravity *(1.25*(rr_2d[1 * nlat + i]+rb_2d[1 * nlat + i])*(1.0+qv_2d[1*nlat + i])  -.25*(rr_2d[2 * nlat + i]+rb_2d[2 * nlat + i])*(1.0+qv_2d[2 * nlat + i]))

        ppi[1] = ppi[1]-ppb_2d[1 * nlat + i]
        for k=0, nz1-1 do

          ppi[k+1] = ppi[k]-vertr[k+1].dzu*gravity*  ( (rr_2d[k * nlat + i]+(rr_2d[k * nlat + i] +rb_2d[k * nlat + i])*qv_2d[k*nlat + i])*vertr[k+1].fzp  + (rr_2d[(k+1) * nlat + i]+(rr_2d[(k+1) * nlat + i]+rb_2d[(k+1) * nlat + i])*qv_2d[(k+1) * nlat + i])*vertr[k+1].fzm )
        end

        for k=0, nz1 do
          pp_2d[k * nlat + i] = .2*ppi[k]+.8*pp_2d[k * nlat + i]
        end

      end   -- end inner iteration loop itrp

    end   -- end outer iteration loop itr


    for k = 0, nz1 do
      rho_2d[k * nlat + i] = rr_2d[k * nlat + i]+rb_2d[k * nlat + i]
      etavs_2d[k * nlat + i] = ((ppb_2d[k * nlat + i]+pp_2d[k * nlat + i])/p0 - 0.252)*pii/2.0
      u_2d[k * nlat + i] = u0*(cmath.pow(cmath.sin(2.*lat_2d[i]),2)) * (cmath.pow(cmath.cos(etavs_2d[k * nlat + i]),1.5))
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
  for i=0, nCells do
    for k=0,nz1 do
      ztemp    = .5*(cr[{i, k+1}].zgrid+cr[{k,i}].zgrid)
      cr[{i, k}].pressure_base = p0*cmath.exp(-gravity*ztemp/(rgas*t0b))
      cr[{i, k}].pressure_p = cmath.pow((cr[{i, k}].pressure_base/p0),(rgas/cp))
      cr[{i, k}].rho_base = cr[{i, k}].pressure_base/(rgas*t0b*cr[{i, k}].zz)
      cr[{i, k}].theta_base = t0b/cr[{i, k}].pressure_p
      cr[{i, k}].rtheta_base = cr[{i, k}].rho_base*cr[{i, k}].theta_base
      cr[{i, k}].exner = cr[{i, k}].pressure_p
      cr[{i, k}].pressure_p = 0.0
      cr[{i, k}].rho_p = 0.0
    end

--!     iterations to converge temperature as a function of pressure
--!
    for itr = 0, 10 do

      for k=0, nz1 do
        eta[k] = (cr[{i, k}].pressure_base+cr[{i, k}].pressure_p)/p0
        etav[k] = (eta[k]-.252)*pii/2.
        if(eta[k] >= znut)  then
          teta[k] = t0*cmath.pow(eta[k],(rgas*dtdz/gravity))
        else
          teta[k] = t0*cmath.pow(eta[k],(rgas*dtdz/gravity)) + delta_t*cmath.pow((znut-eta[k]),5)
        end
      end
      phi = cr[{i, 0}].lat
      for k=0,nz1 do
        temperature_1d[k] = teta[k]+.75*eta[k]*pii*u0/rgas*cmath.sin(etav[k])  *cmath.sqrt(cmath.cos(etav[k]))* ((-2.0*cmath.pow(cmath.sin(phi),6)  *(cmath.pow(cmath.cos(phi),2)+1.0/3.0)+10.0/63.0) *2.*u0*cmath.pow(cmath.cos(etav[k]),1.5)   +(1.6*cmath.pow(cmath.cos(phi),3)   *(cmath.pow(cmath.sin(phi),2)+2.0/3.0)-pii/4.0)*r_earth*omega_e)/(1.+0.61*cr[{i, k}].qv)

        ztemp   = .5*(cr[{k,i}].zgrid+cr[{i, k+1}].zgrid)
        ptemp   = cr[{i, k}].pressure_base + cr[{i, k}].pressure_p

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

        tt[k] = temperature_1d[k]*(1.+1.61*cr[{i, k}].qv)

      end




      for itrp = 0,25 do
        for k=0,nz1 do
          cr[{i,k}].rho_p  = (cr[{i, k}].pressure_p/(rgas*cr[{i, k}].zz) - cr[{i, k}].rho_base*(tt[k]-t0b))/tt[k]
        end

        ppi[1] = p0-.5*dzw[1]*gravity *(1.25*(cr[{i, 1}].rho_p+cr[{i, 1}].rho_base)*(1.+cr[{i, 0}].qv ) -.25*(cr[{i, 2}].rho_p+cr[{i, 2}].rho_base)*(1.+cr[{i, 1}].qv))

        ppi[1] = ppi[1]-cr[{i, 1}].pressure_base
        for k=0,nz1-1 do


           ppi[k+1] = ppi[k]-vertr[k+1].dzu*gravity* ( (cr[{i,k}].rho_p+(cr[{i,k}].rho_p+cr[{i, k}].rho_base)*cr[{i, k}].qv)*vertr[k+1].fzp   + (cr[{i, k+1}].rho_p+(cr[{i, k+1}].rho_p+cr[{i, k+1}].rho_base)*cr[{i, k+1}].qv)*vertr[k+1].fzm)

        end

        for k=0 ,nz1 do
          cr[{i, k}].pressure_p = .2*ppi[k]+.8*cr[{i, k}].pressure_p
        end

      end   -- end inner iteration loop itrp

    end   -- end outer iteration loop itr





    for k=0,nz1 do
      cr[{i, k}].exner = cmath.pow(((cr[{i, k}].pressure_base+cr[{i, k}].pressure_p)/p0),(rgas/cp))
      cr[{i, k}].theta_m = tt[k]/cr[{i, k}].exner
      cr[{i, k}].rtheta_p = cr[{i, k}].theta_m * cr[{i,k}].rho_p+cr[{i, k}].rho_base*(cr[{i, k}].theta_m-cr[{i, k}].theta_base)
      cr[{i, k}].rho_zz = cr[{i, k}].rho_base + cr[{i,k}].rho_p
    end

    --calculation of surface pressure:
    cr[{i, 0}].surface_pressure = 0.5*dzw[1]*gravity * (1.25*(cr[{i, 1}].rho_p + cr[{i, 1}].rho_base) * (1.0 + cr[{i, 0}].qv) -  0.25*(cr[{i, 2}].rho_p + cr[{i, 2}].rho_base) * (1.0 + cr[{i, 1}].qv))
    cr[{i, 0}].surface_pressure = cr[{i, 0}].surface_pressure + cr[{i, 1}].pressure_p + cr[{i, 1}].pressure_base

  end   -- end loop over cells

  var lat_pert = latitude_pert*pii/180.0
  var lon_pert = longitude_pert*pii/180.0




  for iEdge=0,nEdges do

     var vtx1 = er[{iEdge, 0}].verticesOnEdge[0] --verticesOnEdge(1,iEdge)
     var vtx2 = er[{iEdge, 0}].verticesOnEdge[1]
     var lat1 = vr[{vtx1, 0}].lat
     var lat2 = vr[{vtx2, 0}].lat
     var iCell1 = er[{iEdge, 0}].cellsOnEdge[0]
     var iCell2 = er[{iEdge, 0}].cellsOnEdge[1]
     var flux = (0.5*(lat2-lat1) - 0.125*(cmath.sin(4.*lat2) - cmath.sin(4.*lat1))) * sphere_radius / er[{iEdge, 0}].dvEdge

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
       etavs = (0.5*(cr[{iCell1, k}].pressure_base+cr[{iCell2, k}].pressure_base+cr[{iCell1, k}].pressure_p+cr[{iCell2, k}].pressure_p)/p0 - 0.252)*pii/2.0
       var fluxk = u0*flux*(cmath.pow(cmath.cos(etavs),1.5))
       er[{iEdge, k}].u = fluxk + u_pert
     end

------end conditional else---------

     var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
     var cell2 = er[{iEdge, 0}].cellsOnEdge[1]
     for k=0,nz1 do
        er[{iEdge, k}].ru  = 0.5*(cr[{cell1, k}].rho_zz+cr[{cell2, k}].rho_zz)*er[{iEdge, k}].u
     end

--  !
--  ! Generate rotated Coriolis field
--  !

     er[{iEdge, 0}].fEdge = 2.0 * omega_e *  ( -1.0*cmath.cos(er[{iEdge, 0}].lon) * cmath.cos(er[{iEdge, 0}].lat) * cmath.sin(alpha_grid) +  cmath.sin(er[{iEdge, 0}].lat) * cmath.cos(alpha_grid) )
  end


  for iVtx=0,nVertices do
     vr[{iVtx, 0}].fVertex = 2.0 * omega_e *  (-1.0*cmath.cos(vr[{iVtx, 0}].lon) * cmath.cos(vr[{iVtx, 0}].lat) * cmath.sin(alpha_grid) +  cmath.sin(vr[{iVtx, 0}].lat) * cmath.cos(alpha_grid))
  end



--  !
--  !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
--  !
--
--  !
--  !     pre-calculation z-metric terms in omega eqn.
--  !


  var z_edge : double
  var z_edge3 : double
  for iEdge = 0,nEdges do
     var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
     var cell2 = er[{iEdge, 0}].cellsOnEdge[1]

    --Below conditional is unnecessary because we use regions, not pools/threads
     --! Avoid a potential divide by zero below if areaCell(nCells+1) is used in the denominator
    -- if (cell1 <= nCellsSolve or cell2 <= nCellsSolve ) then
        for k = 0, nVertLevels do
          z_edge = (cr[{cell1, k}].zgrid+cr[{cell2, k}].zgrid)/2.0
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

           er[{iEdge, k}].zb[0] = (z_edge-cr[{cell1, k}].zgrid)*er[{iEdge, 0}].dvEdge/cr[{cell1, 0}].areaCell
           er[{iEdge, k}].zb[1] = (z_edge-cr[{cell2, k}].zgrid)*er[{iEdge, 0}].dvEdge/cr[{cell2, 0}].areaCell
           er[{iEdge, k}].zb3[0]=  z_edge3*er[{iEdge, 0}].dvEdge/cr[{cell1, 0}].areaCell
           er[{iEdge, k}].zb3[1] =  z_edge3*er[{iEdge, 0}].dvEdge/cr[{cell2, 0}].areaCell

        end
    -- end

  end

  --! for including terrain
  for iEdge = 0,nEdges do

     var cell1 = er[{iEdge, 0}].cellsOnEdge[0]
     var cell2 = er[{iEdge, 0}].cellsOnEdge[1]



-------------------------------------------------------------------------------
----------TODO: past this point, conversion from Fortran is incomplete---------
-------------------------------------------------------------------------------



     for k = 1, nVertLevels do
        flux =  (vertr[k].fzm*er[{iEdge, k}].ru+vertr[k].fzp*er[{iEdge, k-1}].ru)
        cr[{cell2, k}].rw = cr[{cell2, k}].rw + (vertr[k].fzm*cr[{cell2, k}].zz+vertr[k].fzp*cr[{cell2, k-1}].zz)*er[{iEdge, k}].zb[1]*flux
        cr[{cell1, k}].rw = cr[{cell1, k}].rw - (vertr[k].fzm*cr[{cell1, k}].zz+vertr[k].fzp*cr[{cell1, k-1}].zz)*er[{iEdge, k}].zb[0]*flux

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
  for iCell=0,nCells do
     for k=1,nVertLevels do
        cr[{iCell, k}].w = cr[{iCell, k}].rw / (vertr[k].fzp * cr[{iCell, k-1}].rho_zz + vertr[k].fzm * cr[{iCell, k}].rho_zz)
     end
  end


  --!
  --! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
  --!
  for iEdge = 0, nEdges do
    for k =0, nVertLevels do
      er[{iEdge, k}].v = 0.0
    end
  end

  for iEdge = 0, nEdges do
     for i=0, er[{iEdge, 0}].nEdgesOnEdge do
        eoe = er[{iEdge, 0}].edgesOnEdge_ECP[i]
        for k = 0, nVertLevels do
           er[{iEdge, k}].v = er[{iEdge, k}].v + er[{iEdge, 0}].weightsOnEdge[i] * er[{eoe, k}].u
        end
     end
  end

--------PSURF is never used, in mpas it only exists for purposes of logging. excluding because of this
--  for i=0,10 do
--    psurf = (cf1*(pressure_base(1,i)+pp(1,i)) + cf2*(pressure_base(2,i)+pp(2,i)) + cf3*(pressure_base(3,i)+pp(3,i)))/100.
--
--        psurf = (pressure_base(1,i)+pp(1,i)) + .5*dzw[1]*gravity        &
--                      *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
--                       -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))
--
--  end

  --! Compute rho and theta from rho_zz and theta_m
  for iCell=0,nCells do
     for k=0,nVertLevels do
        cr[{iCell, k}].rho = cr[{iCell, k}].rho_zz * cr[{iCell, k}].zz
        cr[{iCell, k}].theta = cr[{iCell, k}].theta_m / (1.0 + 1.61 * cr[{iCell, k}].qv)
     end
  end

end
