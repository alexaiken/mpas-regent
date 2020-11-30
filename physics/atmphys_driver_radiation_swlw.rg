import "regent"
require "data_structures"
require "physics/ra_cam"

local constants = require("constants")

local levsiz = constants.noznLevels

-- Math library imports
min = regentlib.fmin(double, double)

task radconst()
end

----------------
-- SHORT WAVE --
----------------

task allocate_radiation_sw()
end

task deallocate_radiation_sw()
end

task radiation_sw_from_MPAS()
end

task radiation_sw_to_MPAS()
end

task driver_radiation_sw()
  radiation_sw_from_MPAS()
  radconst()
  camrad()
  radiation_sw_to_MPAS()
end

---------------
-- LONG WAVE --
---------------

task allocate_radiation_lw()
end

task deallocate_radiation_lw()
end

task vinterp_ozn()
-- reads: pin_in, ozmix_in
-- writes: o3vmr_out
-- reads writes: o3vmr, ozmix, pin, pmid
  -- Initialize variables
  var kkstart : int -- level index
  var kount : int -- counter
  var dpu : double -- upper level pressure difference
  var dpl : double -- lower level pressure difference
  var kupper : int[pcols] -- level indices for interpolation

  -- ldf begin:
  for k = 0, levsiz do
    pin[k] = pin_in[k]
  end
  for i = 0, pcols do
    for k = 0, levsiz do
      ozmix[{i, k}] = ozmix_in[{i, k}]
    end
  end
  for i = 0, pcols do
    for k = 0, pver do
      var kk = pver - k + 1
      pmid[{i, kk}] = pmid_in[{i, k}]
    end
  end
  -- ldf end.

  -- Initialize index array
  for i = 0, ncol do
    kupper[i] = 1 -- TODO: 0 or 1?
  end

  for k = 0, pver do
    -- Top level we need to start looking is the top level for the previous k
    -- for all longitude points
    kkstart = levsiz
    for i = 0, ncol do
      kkstart = fmin(kkstart, kupper[i])
    end
    kount = 0

    -- Store level indices for interpolation
    for kk = kkstart, levsiz - 1 do -- Original: kkstart to levsiz - 1
      for i = 0, ncol do
        if (pin[kk] < pmid[{i, k}] and pmid[{i, k}] <= pin[kk + 1]) then
          kupper[i] = kk
          kount += 1
        end
      end

      -- If all indices for this level have been found, do the interpolation and
      -- go to the next level
      if (kount == ncol) then
        for i = 0, ncol do
          dpu = pmid[{i, k}] - pin[kupper[i]]
          dpl = pin[kupper[i] + 1] - pmid[{i, k}]
          o3vmr[{i, k}] = (ozmix[{i, kupper[i]}] * dpl +
                            ozmix[{i, kupper[i]+1}] * dpu) / (dpl + dpu)
        end
        --goto 35
      end
    end

    -- If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
    -- must extrapolate from the bottom or top ozone data level for at least some
    -- of the longitude points.
    for i = 0, ncol do
      if (pmid[{i, k}] < pin[0]) then -- Originally: .lt. pin(1)
        o3vmr[{i, k}] = ozmix[{i,0}] * pmid[{i, k}] / pin[0]
      else if (pmid[{i, k}] > pin[levsiz]) then
        o3vmr[{i,k}] = ozmix[{i, levsiz}]
      else
        dpu = pmid[{i, k}] - pin[kupper[i]]
        dpl = pin[kupper[i]+1] - pmid[{i, k}]
        o3vmr[{i, k}] = (ozmix[{i, kupper[i]}] * dpl +
                          ozmix[{i, kupper[i] + 1}] * dpu) / (dpl + dpu)
      end
    end

    if (kount > ncol) then
      -- call endrun ('VINTERP_OZN: Bad ozone data: non-monotonicity suspected')
    end
--35    continue
  end

  --ldf begin:
  for i = 0, pcols do
    for k = 0, pver do
      var kk = pver - k + 1
      o3vmr_out[{i, kk}] = o3vmr[{i, k}]
    end
  end
  --ldf end.
end

task radiation_lw_from_MPAS(cr : region (ispace(int2d), cell_fs),
                            radt_lw_scheme : regentlib.string,
                            microp_scheme : regentlib.string,
                            config_microp_re : bool,
                            config_o3climatology : bool,
                            xtime_s : double)
where reads (cr.cldfrac, cr.lat, cr.lon, cr.m_ps, cr.re_cloud, cr.re_ice, cr.re_snow, cr.sfc_emiss, cr.skintemp, cr.snow, cr.xice, cr.xland)
do
  for j = jts, jte do
    for i = its, ite do
      sfc_emiss_p[{i, j}] = cr[{i, 0}].sfc_emiss
      tsk_p[{i, j}]       = cr[{i, 0}].skintemp
      snow_p[{i, j}]      = cr[{i, 0}].snow
      xice_p[{i, j}]      = cr[{i, 0}].xice
      xland_p[{i, j}]     = cr[{i, 0}].xland
    end
  end
  for j = jts, jte do
    for k = kts, kte do
      for i = its, ite do
        cldfrac_p[{i,k,j}] = cr[{i, k}].cldfrac
      end
    end
  end

  -- initialization:
  for j = jts, jte do
    for k = kts, kte do
      for i = its, ite do
        f_ice[{i,k,j}]  = 0.0
        f_rain[{i,k,j}] = 0.0
      end
    end
  end

  for j = jts, jte do
    for i = its, ite do
      glw_p[{i, j}]      = 0.0
      lwcf_p[{i, j}]     = 0.0
      lwdnb_p[{i, j}]    = 0.0
      lwdnbc_p[{i, j}]   = 0.0
      lwdnt_p[{i, j}]    = 0.0
      lwdntc_p[{i, j}]   = 0.0
      lwupb_p[{i, j}]    = 0.0
      lwupbc_p[{i, j}]   = 0.0
      lwupt_p[{i, j}]    = 0.0
      lwuptc_p[{i, j}]   = 0.0
      olrtoa_p[{i, j}]   = 0.0
    end
  
    for k = kts, kte do
      for i = its, ite do
        rthratenlw_p[{i,k,j}] = 0.0
      end
    end
  end

  if ([rawstring](radt_lw_scheme) == "rrtmg_lw") then -- string should be trimmed
    if ([rawstring](microp_scheme) == "mp_thompson" or [rawstring](microp_scheme) == "mp_wsm6") then
      if (config_microp_re) then

        for j = jts, jte do
          for k = kts, kte do
            for i = its, ite do
              recloud_p[{i, k, j}] = cr[{i, k}].re_cloud
              reice_p[{i, k, j}]   = cr[{i, k}].re_ice
              resnow_p[{i, k, j}]  = cr[{i, k}].re_snow
            end
          end
        end
      else
        for j = jts, jte do
          for k = kts, kte do
            for i = its, ite do
              recloud_p[{i, k, j}] = 0.0
              reice_p[{i, k, j}]   = 0.0
              resnow_p[{i, k, j}]  = 0.0
            end
          end
        end
      end
      for j = jts, jte do
        for k = kts, kte do
          for i = its, ite do
            rrecloud_p[{i, k, j}] = 0.0
            rreice_p[{i, k, j}]   = 0.0
            rresnow_p[{i, k, j}]  = 0.0
          end
        end
      end
    end

    for j = jts, jte do
      for k = kts, kte + 2 do
        for i = its, ite do
          lwdnflx_p[{i, k, j}]  = 0.0
          lwdnflxc_p[{i, k, j}] = 0.0
          lwupflx_p[{i, k, j}]  = 0.0
          lwupflxc_p[{i, k, j}] = 0.0
        end
      end
    end

    if (config_o3climatology) then
      -- ozone mixing ratio:
      for k = 0, num_oznLevels do
        pin_p[k] = pin[k]
      end
      for j = jts, jte do
        for k = 0, num_oznLevels do
          for i = its, ite do
            o3clim_p[{i, k, j}] = cr[{i, 0}].o3clim[k]
          end
        end
      end

      var nlevs = kte - kts + 1
      var ncols = ite - its + 1
      --if(.not.allocated(p2d) ) allocate(p2d(its:ite,kts:kte) )
      --if(.not.allocated(o32d)) allocate(o32d(its:ite,kts:kte))
      for j = jts, jte do
        for i = its, ite do
          for k = kts, kte do
            o32d[{i, k}] = 0.0
            p2d[{i, k}]  = pres_hyd_p[{i, k, j}] / 100.0
          end
        end
        vinterp_ozn(1, ncols, ncols, nlevs, p2d, pin_p, num_oznlevels, o3clim_p[{1, 1, j}], o32d)
        for i = its, ite do
          for k = kts, kte do
            o3vmr[{k,i}] = o32d[{i,k}]
          end
        end
      end
      --if(allocated(p2d))  deallocate(p2d)
      --if(allocated(o32d)) deallocate(o32d)
    else
      for k = 0, num_oznLevels do
        pin_p[k] = 0.0
      end
      for j = jts, jte do
        for k = 0, num_oznLevels do
          for i = its, ite do
            o3clim_p[{i,k,j}] = 0.0
          end
        end
      end
    end

  else if ([rawstring](radt_lw_scheme) == "cam_lw") then
    for j = jts, jte do
      for i = its, ite do
        xlat_p[{i,j}] = cr[{i, 0}].lat / degrad
        xlon_p[{i,j}] = cr[{i, 0}].lon / degrad
        sfc_albedo_p[{i, j}] = cr[{i, 0}].sfc_albedo

        coszr_p[{i, j}]      = 0.0
        gsw_p[{i, j}]        = 0.0
        swcf_p[{i, j}]       = 0.0
        swdnb_p[{i, j}]      = 0.0
        swdnbc_p[{i, j}]     = 0.0
        swdnt_p[{i, j}]      = 0.0
        swdntc_p[{i, j}]     = 0.0
        swupb_p[{i, j}]      = 0.0
        swupbc_p[{i, j}]     = 0.0
        swupt_p[{i, j}]      = 0.0
        swuptc_p[{i, j}]     = 0.0
      end
      for k = kts, kte do
        for i = its, ite do
          rthratensw_p[{i, k, j}] = 0.0
          cemiss_p[{i, k, j}]     = 0.0
          taucldc_p[{i, k, j}]    = 0.0
          taucldi_p[{i, k, j}]    = 0.0
        end
      end
    end

    -- On the first time-step of each model run, the local arrays absnxt_p, absnst_p,
    -- and emstot_p are filled with the MPAS arrays abstot, absnxt, and emstot. If it
    -- is a new run, these three arrays will be initialized to zero;If a restart run,
    -- these three arrays will be filled with the restart values.
    if (xtime_s < 1.0e-12) then
      for j = jts, jte do
        for n = 0, cam_abs_dim1 do -- Originally from 1 to cam_abs_dim1
          for k = kts, kte do
            for i = its, ite do
              absnxt_p[{i,k,n,j}] = 0.0
            end
          end
        end
        for n = 0, cam_abs_dim2 do -- Originally from 1 to cam_abs_dim2
          for k = kts, kte + 1 do
            for i = its, ite do
              abstot_p[{i,k,n,j}] = 0.0
            end
          end
        end
        for k = kts, kte + 1 do
          for i = its, ite do
            emstot_p[{i, k, j}] = 0.0
          end
        end
      end
    end

    -- ozone mixing ratio:
    for k = 0, num_oznlevels do
      pin_p[k] = pin[k]
    end
    for n = 0, num_months do
      for j = jts, jte do
        for k = 0, num_oznlevels do
          for i = its, ite do
            ozmixm_p(i,k,j,n) = ozmixm(n,k,i)
          end
        end
      end
    end
    --aerosol mixing ratio:
    for k = 0, num_aerlevels do
      m_hybi_p[k] = cr[{0, k}].m_hybi
    end
    for i = its, ite do
      for j = jts, jte do
        m_psp_p[{i, j}] = cr[{i, 0}].m_ps
        m_psn_p[{i, j}] = cr[{i, 0}].m_ps
      end
    end
    for n = 0, num_aerosols do
      for j = jts, jte do
        for k = 0, num_aerlevels do
          for i = its, ite do
            aerosolcp_p[{i,k,j,n}] = cr[{i, k}].aerosols[n]
            aerosolcn_p[{i,k,j,n}] = cr[{i, k}].aerosols[n]
          end
        end
      end
    end
  end
end

task radiation_lw_to_MPAS()
end

task driver_radiation_lw()
  radiation_lw_from_MPAS()
  radconst()
  camrad()
  radiation_lw_to_MPAS()
end