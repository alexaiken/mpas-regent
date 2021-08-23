import "regent"
require "data_structures"

local constants = require("constants")
local c = regentlib.c
local inf = 1.0/0.0
local format = require("std/format")

fabs = regentlib.fabs(double)
floor = regentlib.floor(double)
fmod = regentlib.fmod(double, double)

-- Struct of two factors
struct two_factors {
  fact1 : double,
  fact2 : double,
}

-- check sanity of time interpolation factors to within 32-bit roundoff
task validfactors(fact1 : double,
                  fact2 : double)
  var delta : double = 1.e-6
  if (fabs(fact1 + fact2 - 1.0) > delta or
      fact1 > 1.0 + delta or fact1 < -1 * delta or
      fact2 > 1.0 + delta or fact2 < -1 * delta) then
    return false
  end
  return true
end

-- Purpose: Determine time interpolation factors (normally for a boundary dataset)
--          for linear interpolation.
--
-- Method:  Assume 365 days per year.  Output variable fact1 will be the weight to
--          apply to data at calendar time "cdayminus", and fact2 the weight to apply
--          to data at time "cdayplus".  Combining these values will produce a result
--          valid at time "cday".  Output arguments fact1 and fact2 will be between
--          0 and 1, and fact1 + fact2 = 1 to roundoff.
task getfactors(cycflag : bool,
                np1 : int,
                cdayminus : double,
                cdayplus : double,
                cday : double)
  
  var tf : two_factors

  -- Determine time interpolation factors.  Account for December-January
  -- interpolation if dataset is being cycled yearly.
  if cycflag and np1 == 1 then
    var deltat : double = cdayplus + constants.daysperyear - cdayminus
    if cday > cdayplus then
      tf.fact1 = (cdayplus + constants.daysperyear - cday) / deltat
      tf.fact2 = (cday - cdayminus) / deltat
    else
      tf.fact1 = (cdayplus - cday) / deltat
      tf.fact2 = (cday + constants.daysperyear - cdayminus) / deltat
    end
  else
    var deltat : double = cdayplus - cdayminus
    tf.fact1 = (cdayplus - cday) / deltat
    tf.fact2 = (cday - cdayminus) / deltat
  end

  if validfactors(tf.fact1, tf.fact2) ~= true then
    c.printf("Bad fact1 and/or fact2=%.3f,%.3f\n", tf.fact1, tf.fact2)
  end

  return tf
end

task oznint(julian : double,
            ozmixmj : region(ispace(int3d), double),
            ozmix : region(ispace(int2d), double),
            levsiz : int,
            pcols : int,
            ozncyc : bool)
where
  reads (ozmixmj),
  writes (ozmix)
do

  -- julian starts from 0.0 at 0Z on 1 Jan.
  var intjulian : double = julian + 1.0    -- offset by one day
  -- jan 1st 00z is julian = 1.0 here
  var ijul : int = floor(intjulian)
  -- Note that following will drift. 
  -- Need to use actual month/day info to compute julian.
  intjulian = intjulian - ijul
  ijul = fmod(ijul, constants.daysperyear)
  if (ijul == 0) then 
    ijul = constants.daysperyear 
  end
  intjulian = intjulian + ijul

  var date_oz = array(16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350)
  var np1 : int = 0
  var finddate : bool = false
  -- do m = 1, constants.nMonths
  for m = 0, constants.nMonths do
    if (date_oz[m] > intjulian and finddate == false) then
      np1 = m
      finddate = true
    end
  end

  var cdayozp : double = date_oz[np1]
  var cdayozm : double
  var np : int
  var nm : int
  if np1 > 0 then
    cdayozm = date_oz[np1 - 1]
    np = np1
    nm = np-1
  else
    cdayozm = date_oz[11]
    np = np1
    nm = 11
  end

  var tf : two_factors = getfactors(ozncyc, np1, cdayozm, cdayozp, intjulian)
  var fact1 : double = tf.fact1
  var fact2 : double = tf.fact2

  -- Time interpolation
  for k = 0, levsiz do 
    for i = 0, pcols do
      ozmix[{i, k}] = ozmixmj[{i, k, nm}] * fact1 + ozmixmj[{i, k, np}] * fact2
    end
  end

end

-- Purpose: Interpolate ozone from current time-interpolated values to model levels
--
-- Method: Use pressure values to determine interpolation levels
task radozn(cr : region(ispace(int2d), cell_fs),
            camrad_2d_r : region(ispace(int2d), camrad_2d_fs),
            radctl_2d_pverr_r : region(ispace(int2d), radctl_2d_pverr_fs),
            ncol : int,     -- number of atmospheric columns
            pcols : int,
            pver : int,
            pin : region(ispace(int1d), double),     -- ozone data level pressures (mks)
            levsiz : int,                            -- number of ozone layers
            ozmix : region(ispace(int2d), double))   -- ozone mixing ratio
where
  reads (camrad_2d_r.pmid, pin, ozmix),
  writes (radctl_2d_pverr_r.o3vmr)
do
  --
  -- Initialize index array
  --
  var kupper = region(ispace(int1d, pcols), int)     -- Level indices for interpolation
  for i=0, ncol do
    kupper[i] = 1
  end

  for k=0, pver do
    --
    -- Top level we need to start looking is the top level for the previous k
    -- for all longitude points
    --
    var kkstart : int = levsiz
    for i=0, ncol do
      kkstart = min(kkstart, kupper[i])
    end

    --
    -- Store level indices for interpolation
    --
    var kount : int = 0
    var iter_done : bool = false
    for kk=kkstart, levsiz - 1 do
      for i=0, ncol do
        if ((pin[kk] < camrad_2d_r[{i, k}].pmid) and (camrad_2d_r[{i, k}].pmid < pin[kk + 1])) then
          kupper[i] = kk
          kount = kount + 1
        end
      end

      --
      -- If all indices for this level have been found, do the interpolation and
      -- go to the next level
      --
      if (kount == ncol) then
        iter_done = true
        for i=0, ncol do
          var dpu : double = camrad_2d_r[{i, k}].pmid - pin[kupper[i]]
          var dpl : double = pin[kupper[i] + 1] - camrad_2d_r[{i, k}].pmid
          radctl_2d_pverr_r[{i, k}].o3vmr = (ozmix[{i, kupper[i]}] * dpl + ozmix[{i, kupper[i] + 1}] * dpu) / (dpl + dpu)
        end
        break
      end
    end

    if not iter_done then
      --
      -- If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
      -- must extrapolate from the bottom or top ozone data level for at least some
      -- of the longitude points.
      --
      for i=0, ncol do
        if (camrad_2d_r[{i, k}].pmid < pin[0]) then
          radctl_2d_pverr_r[{i, k}].o3vmr = ozmix[{i, 0}] * camrad_2d_r[{i, k}].pmid / pin[0]
        elseif (camrad_2d_r[{i, k}].pmid > pin[levsiz]) then
          radctl_2d_pverr_r[{i, k}].o3vmr = ozmix[{i, levsiz}]
        else
          var dpu : double = camrad_2d_r[{i, k}].pmid - pin[kupper[i]]
          var dpl : double = pin[kupper[i] + 1] - camrad_2d_r[{i, k}].pmid
          radctl_2d_pverr_r[{i, k}].o3vmr = (ozmix[{i, kupper[i]}] * dpl + ozmix[{i, kupper[i] + 1}] * dpu) / (dpl + dpu)
        end
      end
    end

  end
end


-- From MPAS-Model/src/core_atmosphere/physics/physics_wrf/module_ra_cam.F
--
-- Purpose: 
-- Set latitude and time dependent arrays for input to solar
-- and longwave radiation.
-- Convert model pressures to cgs, and compute ozone mixing ratio, needed for
-- the solar radiation.
--
-- NOTE: variable eccf is unused everywhere, so not included in regent version
__demand(__cuda)
task radinp(cr : region(ispace(int2d), cell_fs),
            camrad_2d_r : region(ispace(int2d), camrad_2d_fs),
            radctl_2d_pverr_r : region(ispace(int2d), radctl_2d_pverr_fs),
            ncol : int,         -- number of atmospheric columns
            pver : int,
            pverp : int,
            pnm : region(ispace(int2d), double))
where
  reads (
    camrad_2d_r.{pmid, pint}, 
    radctl_2d_pverr_r.o3vmr
  ),
  writes (
    camrad_2d_r.pmid,
    radctl_2d_pverr_r.{pbr, o3mmr}, 
    pnm
  )
do
  -- Mainly unused currently.
  ---------------------------Local variables-----------------------------
  --var i : int           -- Longitude loop index
  --var k : int           -- Vertical loop index

  --var calday : double   -- current calendar day
  var vmmr : double     -- Ozone volume mixing ratio
  --var delta : double    -- Solar declination angle
  -----------------------------------------------------------------------

  -- Variables for CUDA
  var pver_ncol_2d = rect2d{ {0, 0}, {pver - 1, ncol - 1} }
  var ncol_1d = rect2d{ {0, pverp}, {ncol - 1, pverp} }

  -- Convert pressure from pascals to dynes/cm2
  for i in pver_ncol_2d do
    radctl_2d_pverr_r[i].pbr = camrad_2d_r[i].pmid * 10.0
    pnm[i] = camrad_2d_r[i].pint * 10.0
  end
  for iNcol in ncol_1d do
    --pnm[{i, pverp}] = camrad_2d_r[{i, pverp}].pint * 10.0
    pnm[iNcol] = camrad_2d_r[iNcol].pint * 10.0
  end

  -- Convert ozone volume mixing ratio to mass mixing ratio:
  vmmr = constants.amo / constants.amd
  for cell in radctl_2d_pverr_r do
    radctl_2d_pverr_r[cell].o3mmr = vmmr * radctl_2d_pverr_r[cell].o3vmr
  end

  return
end

--
-- Saturation vapor pressure table lookup
--
__demand(__inline)
task estblf(td : double,            -- Temperature for saturation lookup  
            phys_tbls : region(ispace(int1d), phys_tbls_fs))
where
  reads (phys_tbls.{tmin, tmax, estbl})
do
  var e : double = max(min(td, phys_tbls[0].tmax), phys_tbls[0].tmin)   -- partial pressure
  var i : int = int(e - phys_tbls[0].tmin) + 1
  var ai : double = int(e - phys_tbls[0].tmin)

  return (phys_tbls[0].tmin + ai - e + 1.0) * 
         phys_tbls[0].estbl[i] - (phys_tbls[0].tmin + ai - e) * 
         phys_tbls[0].estbl[i + 1]
end

-- From MPAS-Model/src/core_atmosphere/physics/physics_wrf/module_ra_cam_support.F
--
-- Purpose: 
-- Utility procedure to look up and return saturation vapor pressure from
-- precomputed table, calculate and return saturation specific humidity
-- (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
-- This routine is useful for evaluating only a selected region in the
-- vertical.
__demand(__cuda)
task aqsat(cr : region(ispace(int2d), cell_fs),
           phys_tbls : region(ispace(int1d), phys_tbls_fs),
           camrad_2d_r : region(ispace(int2d), camrad_2d_fs),
           radctl_2d_pverr_r : region(ispace(int2d), radctl_2d_pverr_fs),
           ilen : int,          -- Length of vectors in I direction which
           klen : int)          -- Length of K direction
where
  reads (camrad_2d_r.t, phys_tbls.{tmin, tmax, estbl}),
  reads writes (
    camrad_2d_r.pmid,
    radctl_2d_pverr_r.{esat, qsat}
  )
do
  format.println("{}, {}", klen, ilen)
  var omeps = 1.0 - constants.ep_2
  var ik_range = rect2d{ {0, 0}, {ilen - 1, klen - 1} }
  for ik in ik_range do
    radctl_2d_pverr_r[ik].esat = estblf(camrad_2d_r[ik].t, phys_tbls)

    --
    -- Saturation specific humidity
    --
    radctl_2d_pverr_r[ik].qsat = 
      constants.ep_2 * radctl_2d_pverr_r[ik].esat 
      / (camrad_2d_r[ik].pmid - omeps * radctl_2d_pverr_r[ik].esat)

    --
    -- The following check is to avoid the generation of negative values
    -- that can occur in the upper stratosphere and mesosphere
    --
    radctl_2d_pverr_r[ik].qsat = min(1.0, radctl_2d_pverr_r[ik].qsat)

    if (radctl_2d_pverr_r[ik].qsat < 0.0) then
      radctl_2d_pverr_r[ik].qsat = 1.0
      radctl_2d_pverr_r[ik].esat = camrad_2d_r[ik].pmid
    end
  end
end

-- From MPAS-Model/src/core_atmosphere/physics/physics_wrf/module_ra_cam_support.F
--
-- Input: match surface pressure, cam interface pressure,
--        month index, number of columns, chunk index
--
-- Output: Aerosol mass mixing ratio (AEROSOL_mmr)
--
-- Method:
--         interpolate column mass (cumulative) from match onto
--           cam's vertical grid (pressure coordinate)
--         convert back to mass mixing ratio
--
task vert_interpolate(cr : region(ispace(int2d), cell_fs),              -- surface pressure at a particular month
                      phys_tbls : region(ispace(int1d), phys_tbls_fs),
                      camrad_2d_r : region(ispace(int2d), camrad_2d_fs),
                      Match_ps : region(ispace(int1d), double),
                      aerosolc : region(ispace(int3d), double),
                      m_hybi : region(ispace(int1d), double),
                      paerlev : int,
                      AEROSOL_mmr : region(ispace(int3d), double),      -- aerosol mmr from MATCH
                      pcols : int,
                      pver : int,
                      pverp : int,
                      ncol : int)       -- chunk index and number of columns
where
  reads (
    camrad_2d_r.pint,                            -- interface pressure from CAM
    phys_tbls.idxVOLC,
    Match_ps,
    aerosolc,
    m_hybi
  ),
  writes (
    AEROSOL_mmr
  )
do

  format.println("Calling vert_interpolate...")

  -----------------------------Local variables-----------------------------

  var m : int                           -- index to aerosol species
  var kupper = region(ispace(int1d, pcols), int)   -- last upper bound for interpolation
  var i : int                           -- loop vars for interpolation
  var k : int 
  var kk : int
  var kkstart : int 
  var kount : int
  var isv : int                         -- loop indices to save
  var ksv : int 
  var msv : int            

  var bad : bool                        -- indicates a bad point found
  var lev_interp_comp : bool            -- interpolation completed for a level

  var AEROSOL = region(ispace(int3d, {pcols,pverp,constants.naer}), double)  
                                        -- cumulative mass of aerosol in column beneath upper
                                        -- interface of level in column at particular month
  var dpl : double                      -- lower and upper intepolation factors
  var dpu : double
  var v_coord : double                  -- vertical coordinate
  var m_to_mmr : double                 -- mass to mass mixing ratio conversion factor
  var AER_diff : double                 -- temp var for difference between aerosol masses

  -------------------------------------------------------------------------

  -- Initialize index array
  for i = 0, ncol do
    kupper[i] = 1
  end

  format.println("Calling vert_interpolate... 1")

  -- assign total mass to topmost level   
  for i = 0, ncol do
    for m = 0, constants.naer do
      format.println("{} {} : {}", i, m, pcols)
      AEROSOL[{i,0,m}] = aerosolc[{i,0,m}]
    end
  end

  format.println("Calling vert_interpolate... 2")

  -- At every pressure level, interpolate onto that pressure level
  for k = 1, pver do

    -- Top level we need to start looking is the top level for the previous k
    -- for all longitude points

    kkstart = paerlev
    for i = 0, ncol do
      kkstart = min(kkstart, kupper[i])
    end
    kount = 0

    -- Store level indices for interpolation

    -- for the pressure interpolation should be comparing
    -- pint(column,lev) with m_hybi(lev)*M_ps_cam_col(month,column,chunk)

    lev_interp_comp = false
    for kk = kkstart, paerlev - 1 do
      if (not lev_interp_comp) then
        for i = 0, ncol do
          v_coord = camrad_2d_r[{i,k}].pint
          if (m_hybi[kk] * Match_ps[i] < v_coord and 
              v_coord <= m_hybi[kk + 1] * Match_ps[i]) then
            kupper[i] = kk
            kount = kount + 1
          end
        end

        -- If all indices for this level have been found, do the interpolation and
        -- go to the next level

        -- Interpolate in pressure.
        if (kount == ncol) then
          for i = 0, ncol do
            for m = 0, constants.naer do
              dpu = camrad_2d_r[{i,k}].pint - m_hybi[kupper[i]] * Match_ps[i]
              dpl = m_hybi[kupper[i]+1] * Match_ps[i] - camrad_2d_r[{i,k}].pint
              AEROSOL[{i,k,m}] =
                (aerosolc[{i, kupper[i], m}] * dpl +
                aerosolc[{i, kupper[i] + 1, m}] * dpu) / (dpl + dpu)
            end
          end
          lev_interp_comp = true
        end
      end
    end

    -- If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
    -- must extrapolate from the bottom or top pressure level for at least some
    -- of the longitude points.

    if (not lev_interp_comp) then
      for i = 0, ncol do
        for m = 1, constants.naer do
          if (camrad_2d_r[{i,k}].pint < m_hybi[0] * Match_ps[i]) then
            AEROSOL[{i,k,m}] =  aerosolc[{i,0,m}]
          elseif (camrad_2d_r[{i,k}].pint > m_hybi[paerlev] * Match_ps[i]) then
            AEROSOL[{i,k,m}] = 0.0
          else
            dpu = camrad_2d_r[{i,k}].pint - m_hybi[kupper[i]] * Match_ps[i]
            dpl = m_hybi[kupper[i] + 1] * Match_ps[i] - camrad_2d_r[{i,k}].pint
            AEROSOL[{i,k,m}] =
              (aerosolc[{i, kupper[i], m}] * dpl +
              aerosolc[{i, kupper[i] + 1, m}] * dpu) / (dpl + dpu)
          end
        end
      end
    end
  end

  format.println("Calling vert_interpolate... 3")
  
  -- aerosol mass beneath lowest interface (pverp) must be 0
  for m = 0, constants.naer do
    for i = 0, ncol do
      AEROSOL[{i, pverp, m}] = 0.
    end
  end

  format.println("Calling vert_interpolate... 4")

  -- Set mass in layer to zero whenever it is less than
  -- 1.e-40 kg/m^2 in the layer
  for cell in AEROSOL do
    if (AEROSOL[{i,k,m}] < 1.e-40) then
      AEROSOL[{i,k,m}] = 0.
    end
  end

  -- Set mass in layer to zero whenever it is less than
  --   10^-15 relative to column total mass
  -- convert back to mass mixing ratios.
  -- exit if mmr is negative
  for m = 0, constants.naer do
    for k = 0, pver do
      for i = 0, ncol do
        AER_diff = AEROSOL[{i,k,m}] - AEROSOL[{i,k+1,m}]
        if (fabs(AER_diff) < 1e-15 * AEROSOL[{i,0,m}]) then
          AER_diff = 0.
        end
        m_to_mmr = constants.gravity / (camrad_2d_r[{i,k+1}].pint - camrad_2d_r[{i,k}].pint)
        AEROSOL_mmr[{i,k,m}] = AER_diff * m_to_mmr
      end
    end
  end
  
  format.println("Ending vert_interpolate...")
end

task background()
end

task scale_aerosols()
end

task get_aerosol()
end

task aerosol_indirect()
end

task radcswmx()
end

task get_int_scales()
end

task radclwmx()
end

task trcmix()
end
