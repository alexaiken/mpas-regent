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
local sphere_radius = 6371229.0 --from ncdump of the grid, default is 1.0, but set to "a" from constants.F in init_atm_core.F


local cio = terralib.includec("stdio.h")
local cmath = terralib.includec("math.h")



task init_atm_case_jw(vr : region(ispace(int2d), vertex_fs),
                      er : region(ispace(int2d), edge_fs),
                      cr : region(ispace(int2d), cell_fs))
where reads writes(vr, er, cr) do 
-- constants from constants.F --- TODO: combine constants in one place so you don't redefine every time you need

  var pii = 3.141592653589793
  var omega = 7.29212E-5 

-- local vars/constants defined at beginning of subroutine
  var u0 = 35.0
  var alpha_grid = 0.0 -- no grid rotation
  var omega_e : double
 
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
  var nlat = 721
  var moisture = false
  var nz1 = nVertLevels 
  var nz = nz1 + 1


---- Scale all distances and areas from a unit sphere to one with radius sphere_radius

  for iCell = 0, nCells do
    cr[{iCell, 0}].x = cr[{iCell, 0}].x * sphere_radius 
    cr[{iCell, 0}].y = cr[{iCell, 0}].y * sphere_radius 
    cr[{iCell, 0}].z = cr[{iCell, 0}].z * sphere_radius 
    cr[{iCell, 0}].areaCell = cr[{iCell, 0}].areaCell * cmath.pow(sphere_radius, 2.0) 
  end 
  for iVert = 0, nVertices do
    vr[iVert].x = vr[iVert].x * sphere_radius
    vr[iVert].y = vr[iVert].y * sphere_radius
    vr[iVert].z = vr[iVert].z * sphere_radius
    vr[iVert].areaTriangle = vr[iVert].areaTriangle * cmath.pow(sphere_radius, 2.0)
    vr[iVert].kiteAreasOnVertex = vr[iVert].kiteAreasOnVertex * cmath.pow(sphere_radius, 2.0)
  end
  for iEdge = 0, nEdges do

    er[{iEdge, 0}].x = er[{iEdge, 0}].x * sphere_radius 
    er[{iEdge, 0}].y = er[{iEdge, 0}].y * sphere_radius 
    er[{iEdge, 0}].z = er[{iEdge, 0}].z * sphere_radius 
    er[{iEdge, 0}].dvEdge = er[{iEdge, 0}].dvEdge * sphere_radius      
    er[{iEdge, 0}].dcEdge = er[{iEdge, 0}].dcEdge * sphere_radius      

  end



-- initialization of moisture:
      scalars(:,:,:) = 0.0
      qsat(:,:)      = 0.0
      relhum(:,:)    = 0.0
      qv_2d(:,:)     = 0.0 -- dimensions = "nVertLevels, nlat"
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
  var omega_e = omega
  var p0 = 1.0E+05

--! We may pass in an hx(:,:) that has been precomputed elsewhere.
--! For now it is independent of k

  for iCell=0, nCells do
    for k = 0, nz do
      phi = cr[{iCell, 0}].latCell
      cr[{iCell, k}].hx = u0 / gravity * cmath.pow(cmath.cos(etavs), 1.5) * ((-2.0 * cmath.pow(cmath.sin(phi), 6) * (cmath.pow(cmath.cos(phi), 2.0) + 1.0/3.0) + 10.0/63.0) * u0*cmath.pow(cmath.cos(etavs), 1.5) + (1.6 * cmath.pow(cmath.cos(phi), 3) * (cmath.pow(cmath.sin(phi), 2) + 2.0/3.0) - pii/4.0)*r_earth*omega_e)
    end
  end


--
--      !     Metrics for hybrid coordinate and vertical stretching
--
--      str = 1.5
--      zt = 45000.
--      dz = zt/float(nz1)
--
--      call mpas_log_write(' hx computation complete ')
--
--      do k=1,nz
--
--!           sh(k) is the stretching specified for height surfaces
--
--            sh(k) = (real(k-1)*dz/zt)**str 
--
--!           to specify specific heights zc(k) for coordinate surfaces,
--!           input zc(k) and define sh(k) = zc(k)/zt
--!           zw(k) is the hieght of zeta surfaces
--!                zw(k) = (k-1)*dz yields constant dzeta
--!                        and nonconstant dzeta/dz
--!                zw(k) = sh(k)*zt yields nonconstant dzeta
--!                        and nearly constant dzeta/dz 
--
--            zw(k) = float(k-1)*dz
--!            zw(k) = sh(k)*zt
--!
--!           ah(k) governs the transition between terrain-following 
--!           and pureheight coordinates
--!                ah(k) = 0 is a terrain-following coordinate
--!                ah(k) = 1 is a height coordinate
-- 
--            ah(k) = 1.-cos(.5*pii*(k-1)*dz/zt)**6
--!            ah(k) = 0.
--            call mpas_log_write(' k, sh, zw, ah = $i $r $r $r', intArgs=(/k/), realArgs=(/sh(k),zw(k),ah(k)/))
--      end do
--      do k=1,nz1
--         dzw (k) = zw(k+1)-zw(k)
--         rdzw(k) = 1./dzw(k)
--         zu(k  ) = .5*(zw(k)+zw(k+1))
--      end do
--      do k=2,nz1
--         dzu (k)  = .5*(dzw(k)+dzw(k-1))
--         rdzu(k)  =  1./dzu(k)
--         fzp (k)  = .5* dzw(k  )/dzu(k)
--         fzm (k)  = .5* dzw(k-1)/dzu(k)
--         rdzwp(k) = dzw(k-1)/(dzw(k  )*(dzw(k)+dzw(k-1)))
--         rdzwm(k) = dzw(k  )/(dzw(k-1)*(dzw(k)+dzw(k-1)))
--      end do
--
--!**********  how are we storing cf1, cf2 and cf3?
--
--      COF1 = (2.*DZU(2)+DZU(3))/(DZU(2)+DZU(3))*DZW(1)/DZU(2) 
--      COF2 =     DZU(2)        /(DZU(2)+DZU(3))*DZW(1)/DZU(3) 
--      CF1  = FZP(2) + COF1
--      CF2  = FZM(2) - COF1 - COF2
--      CF3  = COF2       
--
--!      d1  = .5*dzw(1)
--!      d2  = dzw(1)+.5*dzw(2)
--!      d3  = dzw(1)+dzw(2)+.5*dzw(3)
--!      cf1 = d2*d3*(d3-d2)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
--!      cf2 = d1*d3*(d1-d3)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
--!      cf3 = d1*d2*(d2-d1)/(d2*d3*(d3-d2)+d1*d3*(d1-d3)+d1*d2*(d2-d1))
--
--      call mpas_log_write(' cf1, cf2, cf3 = $r $r $r', realArgs=(/cf1,cf2,cf3/))
--
--      do iCell=1,nCells
--        do k=1,nz
--          zgrid(k,iCell) = (1.-ah(k))*(sh(k)*(zt-hx(k,iCell))+hx(k,iCell))  &
--                         + ah(k) * sh(k)* zt
--        end do
--        do k=1,nz1
--          zz (k,iCell) = (zw(k+1)-zw(k))/(zgrid(k+1,iCell)-zgrid(k,iCell))
--        end do
--      end do
--
--      do i=1, nEdges
--        iCell1 = cellsOnEdge(1,i)
--        iCell2 = cellsOnEdge(2,i)
--        do k=1,nz1
--          zxu (k,i) = 0.5 * (zgrid(k,iCell2)-zgrid(k,iCell1) + zgrid(k+1,iCell2)-zgrid(k+1,iCell1)) / dcEdge(i)
--        end do
--      end do
--      do i=1, nCells
--        do k=1,nz1
--          ztemp = .5*(zgrid(k+1,i)+zgrid(k,i))
--          dss(k,i) = 0.
--          ztemp = zgrid(k,i)
--          if(ztemp.gt.zd+.1)  then
--             dss(k,i) = dss(k,i)+xnutr*sin(.5*pii*(ztemp-zd)/(zt-zd))**2
--          end if
--        end do
--      end do
--
--      !do k=1,nz1
--      !  call mpas_log_write(' k, zgrid(k,1),hx(k,1) = $i $r $r', intArgs=(/k/), realArgs=(/zgrid(k,1),hx(k,1)/))
--      !end do
--
--      !do k=1,nz1
--      !  call mpas_log_write(' k, zxu(k,1) = $i $r', intArgs=(/k/), realArgs=(/zxu(k,1)/))
--      !end do
--
--      call mpas_log_write(' grid metrics setup complete ')
--
--!**************  section for 2d (z,lat) calc for zonal velocity
--
--      dlat = 0.5*pii/float(nlat-1)
--      do i = 1,nlat
--
--        lat_2d(i) = float(i-1)*dlat
--        phi = lat_2d(i)
--        hx_1d    = u0/gravity*cos(etavs)**1.5                           &
--                   *((-2.*sin(phi)**6                                   &
--                         *(cos(phi)**2+1./3.)+10./63.)                  &
--                         *(u0)*cos(etavs)**1.5                          &
--                    +(1.6*cos(phi)**3                                   &
--                         *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)
--
--        do k=1,nz
--          zgrid_2d(k,i) = (1.-ah(k))*(sh(k)*(zt-hx_1d)+hx_1d)  &
--                         + ah(k) * sh(k)* zt
--        end do
--        do k=1,nz1
--          zz_2d(k,i) = (zw(k+1)-zw(k))/(zgrid_2d(k+1,i)-zgrid_2d(k,i))
--        end do
--
--        do k=1,nz1
--          ztemp    = .5*(zgrid_2d(k+1,i)+zgrid_2d(k,i))
--          ppb_2d(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
--          pb_2d(k,i) = (ppb_2d(k,i)/p0)**(rgas/cp)
--          rb_2d(k,i) = ppb_2d(k,i)/(rgas*t0b*zz_2d(k,i))
--          tb_2d(k,i) = t0b/pb_2d(k,i)
--          rtb_2d(k,i) = rb_2d(k,i)*tb_2d(k,i)
--          p_2d(k,i) = pb_2d(k,i)
--          pp_2d(k,i) = 0.0
--          rr_2d(k,i) = 0.0
--        end do
--
--
--        do itr = 1,10
--
--          do k=1,nz1
--            eta (k) = (ppb_2d(k,i)+pp_2d(k,i))/p0
--            etav(k) = (eta(k)-.252)*pii/2.
--            if(eta(k) >= znut)  then
--              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
--            else
--              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
--            end if
--          end do
--
--          phi = lat_2d(i)
--          do k=1,nz1
--            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
--                            *sqrt(cos(etav(k)))*                   &
--                              ((-2.*sin(phi)**6                    &
--                                   *(cos(phi)**2+1./3.)+10./63.)   &
--                                   *2.*u0*cos(etav(k))**1.5        &
--                              +(1.6*cos(phi)**3                    &
--                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*qv_2d(k,i))
--
--
--            ztemp   = .5*(zgrid_2d(k,i)+zgrid_2d(k+1,i))
--            ptemp   = ppb_2d(k,i) + pp_2d(k,i)
--
--            !get moisture 
--            if (moisture) then
--              qv_2d(k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
--            end if
--
--            tt(k) = temperature_1d(k)*(1.+1.61*qv_2d(k,i))
--          end do
--
--          do itrp = 1,25
--            do k=1,nz1
--              rr_2d(k,i)  = (pp_2d(k,i)/(rgas*zz_2d(k,i)) - rb_2d(k,i)*(tt(k)-t0b))/tt(k)
--            end do
--
--            ppi(1) = p0-.5*dzw(1)*gravity                            &
--                          *(1.25*(rr_2d(1,i)+rb_2d(1,i))*(1.+qv_2d(1,i))   &
--                            -.25*(rr_2d(2,i)+rb_2d(2,i))*(1.+qv_2d(2,i)))
--
--            ppi(1) = ppi(1)-ppb_2d(1,i)
--            do k=1,nz1-1
--
--              ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                       &
--                            ( (rr_2d(k  ,i)+(rr_2d(k  ,i)+rb_2d(k  ,i))*qv_2d(k  ,i))*fzp(k+1)   &
--                            + (rr_2d(k+1,i)+(rr_2d(k+1,i)+rb_2d(k+1,i))*qv_2d(k+1,i))*fzm(k+1) )
--            end do
--
--            do k=1,nz1
--              pp_2d(k,i) = .2*ppi(k)+.8*pp_2d(k,i)
--            end do
--
--          end do  ! end inner iteration loop itrp
--
--        end do  ! end outer iteration loop itr
--
--        do k=1,nz1
--          rho_2d(k,i) = rr_2d(k,i)+rb_2d(k,i)
--          etavs_2d(k,i) = ((ppb_2d(k,i)+pp_2d(k,i))/p0 - 0.252)*pii/2.
--          u_2d(k,i) = u0*(sin(2.*lat_2d(i))**2) *(cos(etavs_2d(k,i))**1.5)
--        end do
--
--      end do  ! end loop over latitudes for 2D zonal wind field calc
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
--!
--!---- baroclinc wave initialization ---------------------------------
--!
--!     reference sounding based on dry isothermal atmosphere
--!
--      do i=1, nCells
--        do k=1,nz1
--          ztemp    = .5*(zgrid(k+1,i)+zgrid(k,i))
--          ppb(k,i) = p0*exp(-gravity*ztemp/(rgas*t0b)) 
--          pb (k,i) = (ppb(k,i)/p0)**(rgas/cp)
--          rb (k,i) = ppb(k,i)/(rgas*t0b*zz(k,i))
--          tb (k,i) = t0b/pb(k,i)
--          rtb(k,i) = rb(k,i)*tb(k,i)
--          p  (k,i) = pb(k,i)
--          pp (k,i) = 0.
--          rr (k,i) = 0.
--        end do
--
--!       if(i == 1) then
--!         do k=1,nz1
--!           call mpas_log_write(' k, ppb, pb, rb, tb (k,1) = $i $r $r $r $r', intArgs=(/k/), realArgs=(/ppb(k,1),pb(k,1),rb(k,1)*zz(k,1),tb(k,1)/))
--!         end do
--!       end if
--
--!     iterations to converge temperature as a function of pressure
--!
--        do itr = 1,10
--
--          do k=1,nz1
--            eta (k) = (ppb(k,i)+pp(k,i))/p0
--            etav(k) = (eta(k)-.252)*pii/2.
--            if(eta(k).ge.znut)  then
--              teta(k) = t0*eta(k)**(rgas*dtdz/gravity)
--            else
--              teta(k) = t0*eta(k)**(rgas*dtdz/gravity) + delta_t*(znut-eta(k))**5
--            end if
--          end do
--          phi = latCell(i)
--          do k=1,nz1
--            temperature_1d(k) = teta(k)+.75*eta(k)*pii*u0/rgas*sin(etav(k))      &
--                            *sqrt(cos(etav(k)))*                   &
--                              ((-2.*sin(phi)**6                    &
--                                   *(cos(phi)**2+1./3.)+10./63.)   &
--                                   *2.*u0*cos(etav(k))**1.5        &
--                              +(1.6*cos(phi)**3                    &
--                                *(sin(phi)**2+2./3.)-pii/4.)*r_earth*omega_e)/(1.+0.61*scalars(index_qv,k,i))
--
--            ztemp   = .5*(zgrid(k,i)+zgrid(k+1,i))
--            ptemp   = ppb(k,i) + pp(k,i)
--
--            !get moisture 
--            if (moisture) then
-- 
--                !scalars(index_qv,k,i) = env_qv( ztemp, temperature_1d(k), ptemp, rh_max )
--
--               if(ptemp < 50000.) then
--                  relhum(k,i) = 0.0
--               elseif(ptemp > p0) then
--                  relhum(k,i) = 1.0
--               else
--                  relhum(k,i) = (1.-((p0-ptemp)/50000.)**1.25)
--               end if
--               relhum(k,i) = min(rh_max,relhum(k,i))
--
--               !.. calculation of water vapor mixing ratio:
--               if (temperature_1d(k) > 273.15) then
--                   es  = 1000.*0.6112*exp(17.67*(temperature_1d(k)-273.15)/(temperature_1d(k)-29.65))
--               else
--                   es  = 1000.*0.6112*exp(21.8745584*(temperature_1d(k)-273.15)/(temperature_1d(k)-7.66))
--               end if
--               qsat(k,i) = (287.04/461.6)*es/(ptemp-es)
--               if(relhum(k,i) .eq. 0.0) qsat(k,i) = 0.0
--               scalars(index_qv,k,i) = relhum(k,i)*qsat(k,i)
--            end if
--
--            tt(k) = temperature_1d(k)*(1.+1.61*scalars(index_qv,k,i))
--
--          end do
--
--          do itrp = 1,25
--            do k=1,nz1
--              rr(k,i)  = (pp(k,i)/(rgas*zz(k,i)) - rb(k,i)*(tt(k)-t0b))/tt(k)
--            end do
--
--            ppi(1) = p0-.5*dzw(1)*gravity                         &
--                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
--                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))
--
--            ppi(1) = ppi(1)-ppb(1,i)
--            do k=1,nz1-1
--
--!              ppi(k+1) = ppi(k)-.5*dzu(k+1)*gravity*                     &
--!                            (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i)   &
--!                            +rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))
--
--               ppi(k+1) = ppi(k)-dzu(k+1)*gravity*                                                  &
--                             ( (rr(k  ,i)+(rr(k  ,i)+rb(k  ,i))*scalars(index_qv,k  ,i))*fzp(k+1)   &
--                             + (rr(k+1,i)+(rr(k+1,i)+rb(k+1,i))*scalars(index_qv,k+1,i))*fzm(k+1) )
--
--            end do
--
--            do k=1,nz1
--              pp(k,i) = .2*ppi(k)+.8*pp(k,i)
--            end do
--
--          end do  ! end inner iteration loop itrp
--
--        end do  ! end outer iteration loop itr
--
--        do k=1,nz1
--          p (k,i) = ((ppb(k,i)+pp(k,i))/p0)**(rgas/cp)
--          t (k,i) = tt(k)/p(k,i)
--          rt (k,i) = t(k,i)*rr(k,i)+rb(k,i)*(t(k,i)-tb(k,i))
--          rho_zz (k,i) = rb(k,i) + rr(k,i)
--        end do
--
--        !calculation of surface pressure:
--        surface_pressure(i) = 0.5*dzw(1)*gravity                                    &
--                        * (1.25*(rr(1,i) + rb(1,i)) * (1. + scalars(index_qv,1,i))  &
--                        -  0.25*(rr(2,i) + rb(2,i)) * (1. + scalars(index_qv,2,i)))
--        surface_pressure(i) = surface_pressure(i) + pp(1,i) + ppb(1,i)
--
--      end do  ! end loop over cells
--
--      !call mpas_log_write('')
--      !call mpas_log_write('--- initialization of water vapor:')
--      !do iCell = 1, nCells
--      !   if(iCell == 1 .or. iCell == nCells) then
--      !      do k = nz1, 1, -1
--      !         call mpas_log_write('$i $i $r $r $r $r', intArgs=(/iCell,k/), realArgs=(/t(k,iCell),relhum(k,iCell),qsat(k,iCell),scalars(index_qv,k,iCell)/))
--      !      end do
--       !     call mpas_log_write('')
--      !   end if
--      !end do
--
--      lat_pert = latitude_pert*pii/180.
--      lon_pert = longitude_pert*pii/180.
--
--      do iEdge=1,nEdges
--
--         vtx1 = verticesOnEdge(1,iEdge)
--         vtx2 = verticesOnEdge(2,iEdge)
--         lat1 = latVertex(vtx1)
--         lat2 = latVertex(vtx2)
--         iCell1 = cellsOnEdge(1,iEdge)
--         iCell2 = cellsOnEdge(2,iEdge)
--         flux = (0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)
--
--         if (config_init_case == 2) then
--            r_pert = sphere_distance( latEdge(iEdge), lonEdge(iEdge), &
--                                      lat_pert, lon_pert, 1.0_RKIND)/(pert_radius)
--            u_pert = u_perturbation*exp(-r_pert**2)*(lat2-lat1) * sphere_radius / dvEdge(iEdge)
--
--         else if (config_init_case == 3) then
--            lon_Edge = lonEdge(iEdge)
--            u_pert = u_perturbation*cos(k_x*(lon_Edge - lon_pert)) &
--                         *(0.5*(lat2-lat1) - 0.125*(sin(4.*lat2) - sin(4.*lat1))) * sphere_radius / dvEdge(iEdge)
--         else
--            u_pert = 0.0
--         end if
--
--         if (rebalance) then
--
--           call init_atm_calc_flux_zonal(u_2d,etavs_2d,lat_2d,flux_zonal,lat1,lat2,dvEdge(iEdge),sphere_radius,u0,nz1,nlat)
--           do k=1,nVertLevels
--             fluxk = u0*flux_zonal(k)/(0.5*(rb(k,iCell1)+rb(k,iCell2)+rr(k,iCell1)+rr(k,iCell2)))
--             u(k,iEdge) = fluxk + u_pert
--           end do
--
--         else 
--
--           do k=1,nVertLevels
--             etavs = (0.5*(ppb(k,iCell1)+ppb(k,iCell2)+pp(k,iCell1)+pp(k,iCell2))/p0 - 0.252)*pii/2.
--             fluxk = u0*flux*(cos(etavs)**1.5)
--             u(k,iEdge) = fluxk + u_pert
--           end do
--
--         end if
--
--         cell1 = cellsOnEdge(1,iEdge)
--         cell2 = cellsOnEdge(2,iEdge)
--         do k=1,nz1
--            ru(k,iEdge)  = 0.5*(rho_zz(k,cell1)+rho_zz(k,cell2))*u(k,iEdge)
--         end do
--
--      !
--      ! Generate rotated Coriolis field
--      !
--
--         fEdge(iEdge) = 2.0 * omega_e * &
--                                       ( -cos(lonEdge(iEdge)) * cos(latEdge(iEdge)) * sin(alpha_grid) + &
--                                         sin(latEdge(iEdge)) * cos(alpha_grid) &
--                                       )
--      end do
--
--      do iVtx=1,nVertices
--         fVertex(iVtx) = 2.0 * omega_e * &
--                                         (-cos(lonVertex(iVtx)) * cos(latVertex(iVtx)) * sin(alpha_grid) + &
--                                          sin(latVertex(iVtx)) * cos(alpha_grid) &
--                                         )
--      end do
--
--      !
--      !  CALCULATION OF OMEGA, RW = ZX * RU + ZZ * RW
--      !
--
--      !
--      !     pre-calculation z-metric terms in omega eqn.
--      !
--      do iEdge = 1,nEdges
--         cell1 = cellsOnEdge(1,iEdge)
--         cell2 = cellsOnEdge(2,iEdge)
--
--         
--         ! Avoid a potential divide by zero below if areaCell(nCells+1) is used in the denominator
--         if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then
--            do k = 1, nVertLevels
--
--               if (config_theta_adv_order == 2) then
--
--                  z_edge = (zgrid(k,cell1)+zgrid(k,cell2))/2.
--
--               else if (config_theta_adv_order == 3 .or. config_theta_adv_order ==4) then !theta_adv_order == 3 or 4 
--
--                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * zgrid(k,cell1)
--                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * zgrid(k,cell2)
--
--!  WCS fix 20120711
--
--                  do i=1, nEdgesOnCell(cell1)
--                     if ( cellsOnCell(i,cell1) > 0)       &
--                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * zgrid(k,cellsOnCell(i,cell1))
--                  end do
--                  do i=1, nEdgesOnCell(cell2)
--                     if ( cellsOnCell(i,cell2) > 0)       &
--                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * zgrid(k,cellsOnCell(i,cell2))
--                  end do             
--
--                  z_edge =  0.5*(zgrid(k,cell1) + zgrid(k,cell2))         &
--                                - (dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.
--
--                  if (config_theta_adv_order == 3) then
--                     z_edge3 =  - (dcEdge(iEdge) **2) * (d2fdx2_cell1 - d2fdx2_cell2) / 12.
--                  else
--                     z_edge3 = 0.
--                  end if
--
--               end if
--
--               zb(k,1,iEdge) = (z_edge-zgrid(k,cell1))*dvEdge(iEdge)/areaCell(cell1)
--               zb(k,2,iEdge) = (z_edge-zgrid(k,cell2))*dvEdge(iEdge)/areaCell(cell2)
--               zb3(k,1,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell1)
--               zb3(k,2,iEdge)=  z_edge3*dvEdge(iEdge)/areaCell(cell2)
--
--            end do
--         end if
--
--      end do
--
--      ! for including terrain
--      rw = 0.0
--      w = 0.0
--      do iEdge = 1,nEdges
--
--         cell1 = cellsOnEdge(1,iEdge)
--         cell2 = cellsOnEdge(2,iEdge)
--
--         do k = 2, nVertLevels
--            flux =  (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
--            rw(k,cell2) = rw(k,cell2) + (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb(k,2,iEdge)*flux
--            rw(k,cell1) = rw(k,cell1) - (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb(k,1,iEdge)*flux
--
--            if (config_theta_adv_order ==3) then 
--               rw(k,cell2) = rw(k,cell2)    &
--                                            - sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
--                                              (fzm(k)*zz(k,cell2)+fzp(k)*zz(k-1,cell2))*zb3(k,2,iEdge)*flux
--               rw(k,cell1) = rw(k,cell1)    &
--                                            + sign(1.0_RKIND,ru(k,iEdge))*config_coef_3rd_order* &
--                                              (fzm(k)*zz(k,cell1)+fzp(k)*zz(k-1,cell1))*zb3(k,1,iEdge)*flux
--            end if
--
--         end do
--
--      end do
--
--      ! Compute w from rho_zz and rw
--      do iCell=1,nCells
--         do k=2,nVertLevels
--            w(k,iCell) = rw(k,iCell) / (fzp(k) * rho_zz(k-1,iCell) + fzm(k) * rho_zz(k,iCell))
--         end do
--      end do
--
--
--      !
--      ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
--      !
--      v(:,:) = 0.0
--      do iEdge = 1, nEdges
--         do i=1,nEdgesOnEdge(iEdge)
--            eoe = edgesOnEdge(i,iEdge)
--            do k = 1, nVertLevels
--               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
--           end do
--         end do
--      end do
--
--      do i=1,10
--        psurf = (cf1*(ppb(1,i)+pp(1,i)) + cf2*(ppb(2,i)+pp(2,i)) + cf3*(ppb(3,i)+pp(3,i)))/100.
--
--            psurf = (ppb(1,i)+pp(1,i)) + .5*dzw(1)*gravity        &
--                          *(1.25*(rr(1,i)+rb(1,i))*(1.+scalars(index_qv,1,i))   &
--                            -.25*(rr(2,i)+rb(2,i))*(1.+scalars(index_qv,2,i)))
--
--        call mpas_log_write(' i, psurf, lat = $i $r $r', intArgs=(/i/), realArgs=(/psurf,latCell(i)*180./3.1415828/))
--      end do
--
--      ! Compute rho and theta from rho_zz and theta_m
--      do iCell=1,nCells
--         do k=1,nVertLevels
--            rho(k,iCell) = rho_zz(k,iCell) * zz(k,iCell)
--            theta(k,iCell) = t(k,iCell) / (1.0 + 1.61 * scalars(index_qv,k,iCell))
--         end do
--      end do
--
--   end subroutine init_atm_case_jw
