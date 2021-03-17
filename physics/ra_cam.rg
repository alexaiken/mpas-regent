import "regent"
require "data_structures"
require "physics_data_structures"
require "physics/ra_cam_cld_support"
require "physics/ra_cam_radctl_support"

local constants = require("constants")
local c = regentlib.c
local format = require("std/format")

task param_cldoptics_calc()
  cldefr()
  cldems()
  cldovrlap()
end

------------------------------------------------------------------------- 
-- MPAS-Model/src/core_atmosphere/physics/physics_wrf/module_ra_cam.F
------------------------------------------------------------------------- 
-- 
-- Purpose: 
-- Driver for radiation computation.
-- 
-- Method: 
-- Radiation uses cgs units, so conversions must be done from
-- model fields to radiation fields.
--
-- Author: CCM1,  CMS Contact: J. Truesdale
-- 
-------------------------------------------------------------------------
task radctl(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
            camrad_1d_r : region(ispace(int1d), camrad_1d_fs),
            camrad_2d_r : region(ispace(int2d), camrad_2d_fs),
            lchnk : int, ncol : int, pcols : int, pver : int, pverp : int, pverr : int, pverrp : int, ppcnst : int, pcnst : int,
            qm1 : region(ispace(int3d), double),
            julian : double,
            ozmixmj : region(ispace(int3d), double), ozmix : region(ispace(int2d), double),
            levsiz : double,
            pin : region(ispace(int1d), double),
            ozncyc : bool,
            m_psp : region(ispace(int1d), double),
            m_psn : region(ispace(int1d), double),
            aerosoljp : region(ispace(int3d), double),
            aerosoljn : region(ispace(int3d), double),
            m_hybi : region(ispace(int1d), double),
            paerlev : int,
            dolw : bool, dosw : bool,
            swcf : region(ispace(int1d), double),
            lwcf : region(ispace(int1d), double),
            flut : region(ispace(int1d), double))
where reads (
        camrad_2d_fs.{pmid, pint, t}, 
        phys_tbls, 
        ozmixmj, 
        ozmix, 
        pin,
        m_psp, m_psn,
        aerosoljp, aerosoljn,
        m_hybi
      ),
      writes (
        camrad_2d_fs.pmid,
        ozmix,
        swcf
      )
do

  -----------------------------Input variables-----------------------------

   var nspint : int = 19           -- Num of spctrl intervals across solar spectrum
   var naer_groups : int = 7       -- Num of aerosol groups for optical diagnostics
          -- current groupings are sul, sslt, all carbons, all dust, background, and all aerosols

  -----------------------------Local variables-----------------------------

  var i : int
  var k : int

  -- indexes of gases in constituent array
  var in2o : int
  var ich4 : int
  var if11 : int
  var if12 : int

  var eccf : double       -- Earth/sun distance factor

  var radctl_1d_r = region(ispace(int1d, pcols), radctl_1d_fs)
  var radctl_2d_pver_r = region(ispace(int2d, {pcols, pver}), radctl_2d_pver_fs)
  var radctl_2d_pverr_r = region(ispace(int2d, {pcols, pverr}), radctl_2d_pverr_fs)
  var radctl_3d_aero_r = region(ispace(int3d, {pcols, nspint, naer_groups}), radctl_3d_aero_fs)
  
  var pnm = region(ispace(int2d, {pcols, pverrp}), double)      -- Model interface pressures (dynes/cm2)
  
  var aerosol = region(ispace(int3d, {pcols, pver, constants.naer_all}), double) -- aerosol mass mixing ratios
  fill(aerosol, 0.0) -- TODO: TEMP
  
  var scales = region(ispace(int1d, constants.naer_all), double)  -- scaling factors for aerosols
  -------------------------------------------------------------------------

  --
  -- Interpolate ozone volume mixing ratio to model levels
  --
  -- WRF: added pin, levsiz, ozmix here
  oznint(julian, ozmixmj, ozmix, levsiz, pcols, ozncyc)

  radozn(cr, camrad_2d_fs, radctl_2d_pverr_r, ncol, pcols, pver, pin, levsiz, ozmix)

  --
  -- Set chunk dependent radiation input
  --
  radinp(cr, radctl_2d_pverr_r, ncol, pver, pverp, pnm)

  --
  -- Solar radiation computation
  --
  if (dosw) then
    --
    -- calculate heating with aerosols
    --
    aqsat(cr, phys_tbls, radctl_2d_pverr_r, ncol, pver)

    for i = 0, ncol do
      for j = 0, pver do
        radctl_2d_pverr_r[{i, j}].rh = qm1[{i, j, 0}] / radctl_2d_pverr_r[{i, j}].qsat[{i, j}] *
            ((1.0 - phys_tbls[0].epsilo) * radctl_2d_pverr_r[{i, j}].qsat + phys_tbls[0].epsilo) /
            ((1.0 - phys_tbls[0].epsilo) * qm1[{i, j, 0}] + phys_tbls[0].epsilo)
      end
    end

    -- REGENT NOTE: Skipping lines 1662-1715 in module_ra_cam.F because radforce always == false

    get_int_scales() -- TODO

    get_aerosol(cr, phys_tbls, camrad_1d_r, lchnk, julian, aerosoljp, aerosoljn, m_hybi, paerlev, 
                pcols, pver, pverp, pverr, pverrp, aerosol, scales)

    aerosol_indirect() -- TODO
  
    radcswmx() -- TODO

    -- -- tls ---------------------------------------------------------------2
    --
    -- Convert units of shortwave fields needed by rest of model from CGS to MKS
    --
    for i = 0, ncol do
      radctl_1d_r[i].solin = radctl_1d_r[i].solin * 1.e-3
      camrad_1d_r[i].fsds = camrad_1d_r[i].fsds * 1.e-3
      radctl_1d_r[i].fsnirt = radctl_1d_r[i].fsnirt * 1.e-3
      radctl_1d_r[i].fsnrtc = radctl_1d_r[i].fsnrtc * 1.e-3
      radctl_1d_r[i].fsnirtsq = radctl_1d_r[i].fsnirtsq * 1.e-3
      camrad_1d_r[i].fsnt = camrad_1d_r[i].fsnt * 1.e-3
      camrad_1d_r[i].fsns = camrad_1d_r[i].fsns * 1.e-3
      radctl_1d_r[i].fsntc = radctl_1d_r[i].fsntc * 1.e-3
      radctl_1d_r[i].fsnsc = radctl_1d_r[i].fsnsc * 1.e-3
      radctl_1d_r[i].fsdsc = radctl_1d_r[i].fsdsc * 1.e-3
      radctl_1d_r[i].fsntoa = radctl_1d_r[i].fsntoa * 1.e-3
      radctl_1d_r[i].fsntoac = radctl_1d_r[i].fsntoac * 1.e-3
      swcf[i] = radctl_1d_r[i].fsntoa - radctl_1d_r[i].fsntoac
    end
    for i = 0, ncol do
      for j = 0, pver do
        radctl_2d_pver_r[{i, j}].ftem = camrad_2d_r[{i, j}].qrs / phys_tbls[0].cpair
      end
    end

    -- Added upward/downward total and clear sky fluxes
    for k = 0, pverp do
      for i = 0, ncol do
        camrad_2d_r[{i, k}].fsup  = camrad_2d_r[{i, k}].fsup * 1.e-3
        camrad_2d_r[{i, k}].fsupc = camrad_2d_r[{i, k}].fsupc * 1.e-3
        camrad_2d_r[{i, k}].fsdn  = camrad_2d_r[{i, k}].fsdn * 1.e-3
        camrad_2d_r[{i, k}].fsdnc = camrad_2d_r[{i, k}].fsdnc * 1.e-3
      end
    end
  end

  --
  -- Longwave radiation computation
  --
  if (dolw) then
    get_int_scales() -- TODO

    get_aerosol(cr, phys_tbls, lchnk, julian, m_psp, m_psn, aerosoljp, aerosoljn, m_hybi, paerlev, 
                pcols, pver, pverp, pverr, pverrp, aerosol, scales)

    --
    -- Convert upward longwave flux units to CGS
    --
    for i = 0, ncol do
      radctl_1d_r[i].lwupcgs = camrad_1d_r[i].lwups
    end

    --
    -- Do longwave computation. If not implementing greenhouse gas code then
    -- first specify trace gas mixing ratios. If greenhouse gas code then:
    --  o ixtrcg   => indx of advected n2o tracer
    --  o ixtrcg+1 => indx of advected ch4 tracer
    --  o ixtrcg+2 => indx of advected cfc11 tracer
    --  o ixtrcg+3 => indx of advected cfc12 tracer
    --
    trcmix() -- TODO

    radclwmx() -- TODO

    --
    -- Convert units of longwave fields needed by rest of model from CGS to MKS
    --
    for i = 0, ncol do
      camrad_1d_r[i].flnt  = camrad_1d_r[i].flnt * 1.e-3
      flut[i]  = flut[i] * 1.e-3
      radctl_1d_r[i].flutc = radctl_1d_r[i].flutc * 1.e-3
      camrad_1d_r[i].flns = camrad_1d_r[i].flns * 1.e-3
      radctl_1d_r[i].flntc = radctl_1d_r[i].flntc * 1.e-3
      radctl_1d_r[i].flnsc = radctl_1d_r[i].flnsc * 1.e-3
      camrad_1d_r[i].flwds = camrad_1d_r[i].flwds * 1.e-3
      lwcf[i]  = radctl_1d_r[i].flutc - flut[i]
    end

    -- Added upward/downward total and clear sky fluxes
    for k = 0, pverp do
      for i = 0, ncol do
        camrad_2d_r[i].flup  = camrad_2d_r[i].flup * 1.e-3
        camrad_2d_r[i].flupc = camrad_2d_r[i].flupc * 1.e-3
        camrad_2d_r[i].fldn  = camrad_2d_r[i].fldn * 1.e-3
        camrad_2d_r[i].fldnc = camrad_2d_r[i].fldnc * 1.e-3
      end
    end
  end
end

task camrad(cr : region(ispace(int2d), cell_fs),
            phys_tbls : region(ispace(int1d), phys_tbls_fs),
            levsiz : int,
            m_psp : region(ispace(int2d), double),
            m_psn : region(ispace(int2d), double),
            julian : double,
            ozncyc : bool,
            paerlev : int,
            naer_c : int
            --abstot_3d : region(isapce(int3d), double)
          )
where 
  reads (phys_tbls)
do
  -----------------------------Local variables-----------------------------

  var lchnk : int
  var ncol : int
  var pcols : int
  var pver : int
  var pverp : int = 10
  var pverr : int = 10
  var pverrp : int = 10

  var pcnst : int
  var pnats : int
  var ppcnst : int
  var i : int
  var j : int
  var k : int
  var ii : int
  var kk : int
  var kk1 : int
  var m : int
  var n : int

  var begchunk : int
  var endchunk : int

  var nyrm : int
  var nyrp : int

  var doymodel : double
  var doydatam : double
  var doydatap : double
  var deltat : double
  var fact1 : double
  var fact2 : double

  var xt24 : double
  var tloctm : double
  var hrang : double
  var xxlat : double
  var oldxt24 : double

  var camrad_1d_r = region(ispace(int1d, constants.nCells), camrad_1d_fs)
  var m_psjp = region(ispace(int1d, constants.nCells), double) -- MATCH surface pressure
  var m_psjn = region(ispace(int1d, constants.nCells), double) -- MATCH surface pressure
  var swcftoa = region(ispace(int1d, constants.nCells), double)
  var lwcftoa = region(ispace(int1d, constants.nCells), double)
  var olrtoa = region(ispace(int1d, constants.nCells), double)

  var m_hybi = region(ispace(int1d, paerlev), double)
  var pin = region(ispace(int1d, levsiz), double)

  var camrad_2d_r = region(ispace(int2d, {constants.nCells, constants.nVertLevels}), camrad_2d_fs)

  var ozmix = region(ispace(int2d, {constants.nCells, levsiz}), double)

  var q = region(isapce(int3d, {constants.nCells, constants.nVertLevels, n_cldadv}), double)

  var ozmixmj = region(ispace(int3d, {constants.nCells, levsiz, constants.nMonths}), double)

  var aerosoljp = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)
  var aerosoljn = region(ispace(int3d, {constants.nCells, paerlev, naer_c}), double)

  var abstot = region(isapce(int3d, {constants.nCells, constants.nVertLevels, constants.nVertLevels}), double) -- Total absorptivity
  var absnxt = region(isapce(int3d, {constants.nCells, constants.nVertLevels, 4}), double) -- Total nearest layer absorptivity
  var emstot = region(isapce(int3d, {constants.nCells, constants.nVertLevels + 1}), double) -- Total emissivity

  -----------------------------

  format.println("Calling camrad...")

  lchnk    = 1
  begchunk = ims
  endchunk = ime
  ncol     = ite - its + 1
  pcols    = ite - its + 1
  pver     = kte - kts + 1
  pverp    = pver + 1
  pverr    = kte - kts + 1
  pverrp   = pverr + 1
  -- number of advected constituents and non-advected constituents (including water vapor)
  ppcnst   = n_cldadv
  -- number of non-advected constituents
  pnats    = 0
  pcnst    = ppcnst - pnats

  -- check the # species defined for the input climatology and naer

  if (naer_c ~= naer_all) then
    format.println("Physics Fatal Error: naer_c-1 ~= naer_all, {}, {}", naer_c, naer_all)
  end

  -- update CO2 volume mixing ratio (co2vmr)
  
  -- determine time interpolation factors, check sanity
  -- of interpolation factors to within 32-bit roundoff
  -- assume that day of year is 1 for all input data
  --
  nyrm     = yr - yrdata(1) + 1
  nyrp     = nyrm + 1
  doymodel = yr * 365.0 + julian
  doydatam = yrdata(nyrm) * 365. + 1.
  doydatap = yrdata(nyrp) * 365. + 1.
  deltat   = doydatap - doydatam
  fact1    = (doydatap - doymodel) / deltat
  fact2    = (doymodel - doydatam) / deltat
  co2vmr   = (co2(nyrm) * fact1 + co2(nyrp) * fact2) * 1.e-06

  co2mmr   = co2vmr * mwco2 / mwdry

  --===================================================
  -- Radiation computations
  --===================================================

  for k=0, levsiz do
    pin(k) = pin0(k)
  end

  for k=0, paerlev do
    m_hybi(k) = m_hybi0(k)
  end

  -- check for uninitialized arrays
  if (abstot_3d(its, kts, kts, jts) == 0.0 and doabsems == false and dolw) then
    format.println("Physics Message: camrad lw: CAUTION: re-calculating abstot,absnxt, on restart")
    doabsems = true
  end

  for j = jts, jte do

    --
    -- Cosine solar zenith angle for current time step
    --

    for i = its, ite do
      ii = i - its + 1
      -- XT24 is the fractional part of simulation days plus half of RADT expressed in 
      -- units of minutes
      -- JULIAN is in days
      -- RADT is in minutes
      XT24 = MOD(XTIME + RADT * 0.5, 1440.0)
      TLOCTM = GMT + XT24 / 60.0 + XLONG(I, J) / 15.
      HRANG = 15.0 * (TLOCTM - 12.0) * DEGRAD
      XXLAT = XLAT(I, J) * DEGRAD
      clat(ii) = xxlat
      coszrs(II) = SIN(XXLAT) * SIN(DECLIN) + COS(XXLAT) * COS(DECLIN) * COS(HRANG)
    end

    -- moist variables

    for k = kts, kte do
      kk = kte - k + kts

      for i = its, ite do
        ii = i - its + 1

        -- convert to specific humidity
        q(ii, kk, 1) = max(1.e-10, 
                           qv3d(i, k, j) / (1.0 + qv3d(i, k, j)))
        
        if F_QI and F_QC and F_QS then
          q(ii, kk, ixcldliq) = max(0.0, 
                                    qc3d(i, k, j) / (1.0 + qv3d(i, k, j)))
          q(ii, kk, ixcldice) = max(0.0, 
                                    (qi3d(i, k, j) + qs3d(i, k, j)) / (1.0 + qv3d(i, k, j)))
        elseif F_QC and F_QR  then
          -- Warm rain or simple ice
          q(ii, kk, ixcldliq) = 0.
          q(ii, kk, ixcldice) = 0.
          if t_phy(i,k,j) > 273.15 then
            q(ii,kk,ixcldliq) = max(0.0,
                                    qc3d(i, k, j) / (1.0 + qv3d(i, k, j)))
          else
            q(ii,kk,ixcldice) = max(0.0,
                                    qc3d(i, k, j) / (1.0 + qv3d(i, k, j)))
          end
        elseif F_QC and F_QS then
          -- For Ferrier (note that currently Ferrier has QI, so this section will not be used)
          q(ii,kk,ixcldice) = max(0.0,
                                  qc3d(i, k, j) / (1.0 + qv3d(i, k, j)) * f_ice_phy(i, k, j))
          q(ii,kk,ixcldliq) = max(0.0,
                                  qc3d(i, k, j) / (1.0 + qv3d(i, k, j)) 
                                  * (1.0 - f_ice_phy(i, k, j)) * (1.0 - f_rain_phy(i, k, j)))
        else
          q(ii, kk, ixcldliq) = 0.0
          q(ii, kk, ixcldice) = 0.0
        end

        cld(ii, kk) = CLDFRA(I, K, J)
      end
    end

    for i = its, ite do
      ii = i - its + 1
      landfrac(ii) = 2.0 - XLAND(I, J)
      landm(ii) = landfrac(ii)
      snowh(ii) = 0.001 * SNOW(I, J)
      icefrac(ii) = XICE(I, J)
    end

    -- ldf (05-15-2011): In MPAS num_months ranges from 1 to 12 (instead of 2 to 13 in WRF):
    -- REGENT NOTE: num_months now ranges from 0 to 11
    for m = 0, num_months do
      for k = 0, levsiz do
        for i = its, ite do
          ii = i - its + 1
          ozmixmj(ii,k,m) = ozmixm(i,k,j,m)
        end
      end
    end

    for i = its, ite do
      ii = i - its + 1
      camrad_1d_r[ii].m_psjp = m_psp[{i,j}]
      camrad_1d_r[ii].m_psjn = m_psn[{i,j}]
    end

    for n = 1, naer_c do
      for k = 1, paerlev do
        for i = its, ite do
          ii = i - its + 1
          aerosoljp(ii, k, n) = aerosolcp(i, k, j, n)
          aerosoljn(ii, k, n) = aerosolcn(i, k, j, n)
        end
      end
    end

    --
    -- Complete radiation calculations
    --
    for i = its, ite do
      ii = i - its + 1
      lwups(ii) = stebol * EMISS(I, J) * pow(TSK(I, J), 4)
    end

    for k = kts, kte + 1 do
      kk = kte - k + kts + 1
      for i = its, ite do
        ii = i - its + 1
        pint(ii, kk) = p8w(i, k, j)
        if k == kts then
          ps(ii) = pint(ii, kk)
        end
        lnpint(ii, kk) = log(pint(ii, kk))
      end
    end

    if doabsems == false and dolw then
      for kk = 0, cam_abs_dim2 do
        for kk1 = kts, kte+1 do
          for i = its, ite do
            abstot(i, kk1, kk) = abstot_3d(i, kk1, kk, j)
          end
        end
      end
      for kk = 0, cam_abs_dim1 do
        for kk1 = kts, kte do
          for i = its, ite do
            absnxt(i, kk1, kk) = absnxt_3d(i, kk1, kk, j)
          end
        end
      end
      for kk = kts, kte + 1 do
        for i = its, ite do
          emstot(i, kk) = emstot_3d(i, kk, j)
        end
      end
    end

    for k = kts, kte do
      kk = kte - k + kts 
      for i = its, ite do
        ii = i - its + 1
        pmid(ii, kk) = p_phy(i, k, j)
        lnpmid(ii, kk) = log(pmid(ii, kk))
        lnpint(ii, kk) = log(pint(ii, kk))
        pdel(ii, kk) = pint(ii, kk + 1) - pint(ii, kk)
        t(ii, kk) = t_phy(i, k, j)
        zm(ii, kk) = z(i, k, j)
      end
    end

    -- Compute cloud water/ice paths and optical properties for input to radiation
    param_cldoptics_calc() -- TODO

    for i = its, ite do
      ii = i - its + 1
      -- use same albedo for direct and diffuse
      -- change this when separate values are provided
      asdir(ii) = albedo(i,j)
      asdif(ii) = albedo(i,j)
      aldir(ii) = albedo(i,j)
      aldif(ii) = albedo(i,j)
    end

    radctl(
      cr, phys_tbls, camrad_1d_r, camrad_2d_r,
      lchnk, ncol, pcols, pver, pverp, pverr, pverrp, ppcnst, pcnst,
      q,
      julian,
      ozmixmj, ozmix, 
      levsiz, 
      pin, 
      ozncyc,
      m_psjp, m_psjn,
      aerosoljp, aerosoljn, 
      m_hybi, 
      paerlev,
      dolw, dosw,
      swcftoa, lwcftoa, olrtoa
    )

    for k = kts, kte do
      kk = kte - k + kts 
      for i = its, ite do
        ii = i - its + 1
        if dolw then
          RTHRATENLW(I,K,J) = 1.e4 * qrl(ii, kk) / (phys_tbls[0].cpair * pi_phy(i, k, j))
        end
        if dosw then
          RTHRATENSW(I,K,J) = 1.e4 * qrs(ii, kk) / (phys_tbls[0].cpair * pi_phy(i, k, j))
        end
        cemiss(i,k,j)     = emis(ii,kk)
        taucldc(i,k,j)    = tauxcl(ii,kk)
        taucldi(i,k,j)    = tauxci(ii,kk)
      end
    end

    if doabsems and dolw then
      for kk = 0, cam_abs_dim2 do
        for kk1 = kts, kte + 1 do
          for i = its, ite do
            abstot_3d(i, kk1, kk, j) = abstot(i, kk1, kk)
          end
        end
      end
      for kk = 0, cam_abs_dim1 do
        for kk1 = kts, kte do
          for i = its, ite do
            absnxt_3d(i, kk1, kk, j) = absnxt(i, kk1, kk)
          end
        end
      end
      for kk = kts, kte + 1 do
        for i = its, ite do
          emstot_3d(i, kk, j) = emstot(i, kk)
        end
      end
    end

    if PRESENT(SWUPT) then
      if dosw then
        -- Added shortwave and longwave upward/downward total and clear sky fluxes
        for k = kts, kte + 1 do
          kk = kte + 1 - k + kts
          for i = its, ite do
            ii = i - its + 1
            if k == kte + 1 then
              swupt(i, j)     = fsup(ii, kk)
              swuptc(i, j)    = fsupc(ii, kk)
              swdnt(i, j)     = fsdn(ii, kk)
              swdntc(i, j)    = fsdnc(ii, kk)
            end
            if k == kts then
              swupb(i, j)     = fsup(ii, kk)
              swupbc(i, j)    = fsupc(ii, kk)
              swdnb(i, j)     = fsdn(ii, kk)
              swdnbc(i, j)    = fsdnc(ii, kk)
            end
          end
        end
      end
      if dolw then
        -- Added shortwave and longwave upward/downward total and clear sky fluxes
        for k = kts, kte + 1 do
          kk = kte + 1 - k + kts
          for i = its, ite do
            ii = i - its + 1
            if k == kte + 1 then
              lwupt(i, j)     = flup(ii, kk)
              lwuptc(i, j)    = flupc(ii, kk)
              lwdnt(i, j)     = fldn(ii, kk)
              lwdntc(i, j)    = fldnc(ii, kk)
            end
            if k == kts then
              lwupb(i, j)     = flup(ii, kk)
              lwupbc(i, j)    = flupc(ii, kk)
              lwdnb(i, j)     = fldn(ii, kk)
              lwdnbc(i, j)    = fldnc(ii, kk)
            end
          end
        end
      end
    end

    for i = its, ite do
      ii = i - its + 1
      -- Added shortwave and longwave cloud forcing at TOA and surface
      if dolw then
        glw(i, j) = flwds(ii)
        lwcf(i, j) = lwcftoa(ii)
        olr(i, j)  = olrtoa(ii)
      end
      if dosw then
        gsw(i, j) = fsns(ii)
        swcf(i, j) = swcftoa(ii)
        coszr(i, j) = coszrs(ii)
      end
    end
  end    -- j-loop
  format.println("Camrad done")
end