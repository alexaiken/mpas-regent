import "regent"
require "data_structures"
require "physics/ra_cam_cld_support"
require "physics/ra_cam_radctl_support"

fspace camrad_1d_fs {
  coszrs : double,
  landfrac : double,
  landm : double,
  snowh : double,
  icefrac : double,
  lwups : double,
  asdir : double,
  asdif : double,
  aldir : double,
  aldif : double,
  ps : double,
  nmxrgn : int,       -- Number of maximally overlapped regions

  fsns : double,      -- Surface absorbed solar flux
  fsnt : double,      -- Net column abs solar flux at model top
  flns : double,      -- Srf longwave cooling (up-down) flux
  flnt : double,      -- Net outgoing lw flux at model top

  swcftoa : double,   -- Top of the atmosphere solar cloud forcing
  lwcftoa : double,   -- Top of the atmosphere longwave cloud forcing
  olrtoa : double,    -- Top of the atmosphere outgoing longwave

  sols : double,      -- Downward solar rad onto surface (sw direct)
  soll : double,      -- Downward solar rad onto surface (lw direct)
  solsd : double,     -- Downward solar rad onto surface (sw diffuse)
  solld : double,     -- Downward solar rad onto surface (lw diffuse)
  fsds : double,      -- Flux Shortwave Downwelling Surface
  flwds : double,     -- Surface down longwave flux
  m_psjp : double,    -- MATCH surface pressure
  m_psjn : double,
  clat : double,      -- latitude in radians for columns
}

fspace camrad_2d_fs {
  cld : double, 
  pmid : double, 
  lnpmid : double, 
  pdel : double, 
  zm : double, 
  t : double,
  cicewp : double,    -- in-cloud cloud ice water path
  cliqwp : double,    -- in-cloud cloud liquid water path
  emis : double,      -- cloud emissivity
  rel : double,       -- effective drop radius (microns)
  rei : double,       -- ice effective drop size (microns)

  qrs : double,       -- Solar heating rate
  qrl : double,       -- Longwave cooling rate
}

fspace camrad_2d_extended_fs {
  ---------------------------------
  -- extended (1 to kte-kts + 2) --
  pint : double, 
  lnpint : double,
  pmxrgn : double,    -- Maximum values of pressure for each
  
  -- Added outputs of total and clearsky fluxes etc
  fsup : double,      -- Upward total sky solar
  fsupc : double,     -- Upward clear sky solar
  fsdn : double,      -- Downward total sky solar
  fsdnc : double,     -- Downward clear sky solar
  flup : double,      -- Upward total sky longwave
  flupc : double,     -- Upward clear sky longwave
  fldn : double,      -- Downward total sky longwave  
  fldnc : double,     -- Downward clear sky longwave

  -----------------------------------
  -- left shifted (0 to kte-kts+1) --
  tauxcl : double,    -- cloud water optical depth
  tauxci : double,    -- cloud ice optical depth
}

fspace ozone_mix_fs {
  ozmixmj : double[constants.num_months],     -- monthly ozone mixing ratio
  ozmix : double,                             -- ozone mixing ratio (time interpolated)
}

fspace aerosol_fs {
  aerosoljp : double, -- monthly aerosol concentrations
  aerosoljn : double, -- monthly aerosol concentrations
}

task param_cldoptics_calc()
  cldefr()
  cldems()
  cldovrlap()
end

task radctl()

  -- passing temp variables for now
  var julian : double = 1.0
  var levsiz : int = 1
  var num_months : int = 1
  var pcols : int = 1
  var ozncyc : bool = true
  var ozmixmj = region(ispace(int3d, {pcols, levsiz, num_months}), double)
  for i=0, pcols do
    for j=0, levsiz do
      for k=0, num_months do
        ozmixmj[{i, j, k}] = .5
      end
    end
  end
  var ozmix = region(ispace(int2d, {pcols, levsiz}), double) -- ozone mixing ratio
  for i=0, pcols do
    for j=0, levsiz do
      ozmix[{i, j}] = .5
    end
  end
  oznint(julian, ozmixmj, ozmix, levsiz, num_months, pcols, ozncyc)

  -- passing temp variables for now
  var ncol : int = 1     -- number of atmospheric columns
  var pver : int = 1
  var pmid = region(ispace(int2d, {pcols, pver}), double)    -- level pressures (mks)
  for i=0, pcols do
    for j=0, levsiz do
      pmid[{i, j}] = .5
    end
  end
  var pin = region(ispace(int1d, levsiz), double)            -- ozone data level pressures (mks)
  for i=0, pcols do
    pin[i] = .5
  end
  var o3vmr = region(ispace(int2d, {pcols, pver}), double)
  for i=0, pcols do
    for j=0, levsiz do
      o3vmr[{i, j}] = .5
    end
  end
  radozn(ncol, pcols, pver, pmid, pin, levsiz, ozmix, o3vmr)

  radinp()
  aqsat()
  get_rf_scales()
  get_aerosol()
  aerosol_indirect()
  radcswmx()
  get_int_scales()
  radclwmx()
  trcmix()
end

task camrad(rthratenlw,
            rthatensw,
            dolw,
            dosw,
            swupt,
            swuptc,
            swdnt,
            swdntc,
            lwupt,
            lwuptc,
            lwdnt,
            lwdntc,
            swupb,
            swupbc,
            swdnb,
            swdnbc,
            lwupb,
            lwupbc,
            lwdnb,
            lwdnbc,
            olr,
            cemiss,
            taucldc,
            taucldi,
            xlat,
            xlong,
            albedo,
            t_phy,
            tsk,
            emiss,
            qv3d,
            qc3d,
            qr3d,
            qi3d,
            qs3d,
            qg3d,
            f_qv,
            f_qc,
            f_qr,
            f_qi,
            f_qs,
            f_qg,
            f_ice_phy,
            f_rain_phy,
            p_phy,
            p8w,
            z,
            pi_phy,
            rho_phy,
            dz8w,
            cldfra,
            xland,
            xice,
            snow,
            ozmixm,
            pin0,
            levsiz,
            num_months,
            m_psp,
            m_psn,
            aerosolcp,
            aerosolcn,
            m_hybi0,
            cam_abs_dim1,
            cam_abs_dim2,
            paerlev,
            naer_c,
            gmt,
            julday,
            julian,
            yr,
            dt,
            xtime,
            declin,
            solcon,
            radt,
            degrad,
            n_cldadv,
            abstot_3d,
            absnxt_3d,
            emstot_3d,
            doabsems,
            ids,
            ide, 
            jds,
            jde, 
            kds,
            kde,
            ims,
            ime, 
            jms,
            jme, 
            kms,
            kme,
            its,
            ite,
            jts,
            jte,
            kts,
            kte)

  ------------------
  -- START LOCALS --

  lchnk : int
  ncol : int
  pcols : int
  pver : int
  pverp : int
  pverr : int
  pverrp : int

  pcnst : int
  pnats : int
  ppcnst : int
  i : int
  j : int
  k : int
  ii : int
  kk : int
  kk1 : int
  m : int
  n : int

  begchunk : int
  endchunk : int

  nyrm : int
  nyrp : int

  doymodel : double
  doydatam : double
  doydatap : double
  deltat : double
  fact1 : double
  fact2 : double

  xt24 : double
  tloctm : double
  hrang : double
  xxlat : double
  oldxt24 : double

  q = region(ispace(int3d, {ite-its, kte-kts, n_cldadv}), double)
  ozone_mix_locals = region(ispace(int2d, {ite-its, levsiz}), ozone_mix_fs)
  pin = region(ispace(int1d, levsiz), double)
  aerosol_locals = region(ispace(int3d, {ite-its, paerlev, naer_c}), aerosol_fs)
  m_hybi = region(ispace(int1d, paerlev), double)
  abstot = region(ispace(int3d, {its:its, kts:kte+1, kts:kte+1}), double)   -- Total absorptivity
  absnxt = region(ispace(int3d, {its:ite, kts:kte, 4}), double)             -- Total nearest layer absorptivity
  emstot = region(ispace(int2d, {its:ite, kts:kte+1}), double)              -- Total emissivity
  camrad_1d_locals = region(ispace(int1d, ite-its), camrad_1d_fs)
  camrad_2d_locals = region(ispace(int2d, {ite-its, kte-kts}), camrad_2d_fs)
  camrad_2dextended_locals = region(ispace(int2d, {ite-its, kte-kts+1}), camrad_2d_extended_fs)

  -- END LOCALS --
  ----------------

  param_cldoptics_calc()
  radctl()
end