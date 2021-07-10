import "regent"
require "data_structures"

local constants = require("constants")

-- Math library imports
round = regentlib.round(double)
exp = regentlib.exp(double)
log = regentlib.log(double)
log10 = regentlib.log10(double)
pow = regentlib.pow(double, double)
fabs = regentlib.fabs(double)

fspace doublefield {
  x         : double;
}

-- Purpose: 
-- Computes saturation vapor pressure over water and/or over ice using
-- Goff & Gratch (1946) relationships. 
-- <Say what the routine does> 
-- 
-- Method: 
-- T (temperature), and itype are input parameters, while es (saturation
-- vapor pressure) is an output parameter.  The input parameter itype
-- serves two purposes: a value of zero indicates that saturation vapor
-- pressures over water are to be returned (regardless of temperature),
-- while a value of one indicates that saturation vapor pressures over
-- ice should be returned when t is less than freezing degrees.  If itype
-- is negative, its absolute value is interpreted to define a temperature
-- transition region below freezing in which the returned
-- saturation vapor pressure is a weighted average of the respective ice
-- and water value.  That is, in the temperature range 0 => -itype
-- degrees c, the saturation vapor pressures are assumed to be a weighted
-- average of the vapor pressure over supercooled water and ice (all
-- water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
task gffgch(phys_tbls : region(ispace(int1d), phys_tbls_fs),
            t : double,
            es_i : int,
            itype : double,
            set_estblh2o : bool)
where
  reads writes (phys_tbls.{estbl, estblh2o})
do

  -- Local variables --
  var eswtr : double    -- Saturation vapor pressure over water
  var tr : double       -- Transition range for es over water to es over ice
  var itypo : int       -- Intermediate scratch variable for holding itype

  --
  -- Check on whether there is to be a transition region for es
  --
  if (itype < 0) then
    tr    = fabs(itype)
    itypo = itype
    itype = 1
  else
    tr    = 0.0
    itypo = itype
  end
  --if (tr > 40.0) then
  --  write(6,900) tr         -- TODO print instead WRITES OUT TO A FILE
  --end if

  if (t >= (constants.tmelt - tr) or itype ~= 1) then 

    -- Water -- 

    var ps = 1013.246               -- Reference pressure (mb)
    var ts = 373.16                 -- Reference temperature (boiling point of water)
    var e1 = 11.344 * (1.0 - t / ts)
    var e2 = -3.49149 * (ts/t - 1.0)
    var f1 = -7.90298 * (ts/t - 1.0)
    var f2 = 5.02808 * log10(ts / t)
    var f3 = -1.3816 * (pow(10.0, e1) - 1.0) / 10000000.0
    var f4 = 8.1328 * (pow(10.0, e2) - 1.0) / 1000.0
    var f5 = log10(ps)
    var f  = f1 + f2 + f3 + f4 + f5
    if set_estblh2o then
      phys_tbls[0].estblh2o[es_i] = pow(10.0, f) * 100.0
      eswtr = phys_tbls[0].estblh2o[es_i]
    else -- set estbl
      phys_tbls[0].estbl[es_i] = pow(10.0, f) * 100.0
      eswtr = phys_tbls[0].estbl[es_i]
    end

    if (t >= constants.tmelt or itype == 0) then
      itype = itypo
      return itype
    end

  end

  -- Ice -- 

  var t0             = constants.tmelt        -- Reference temperature (freezing point of water)
  var term1 : double = 2.01889049 / (t0 / t)
  var term2 : double = 3.56654 * log(t0 / t)
  var term3 : double = 20.947031 * (t0 / t)
  if set_estblh2o then
    phys_tbls[0].estblh2o[es_i] = 575.185606e10 * exp(-(term1 + term2 + term3))
  else -- set estbl
    phys_tbls[0].estbl[es_i] = 575.185606e10 * exp(-(term1 + term2 + term3))
  end

  if (t < (constants.tmelt - tr)) then
    itype = itypo
    return itype
  end

  -- Weighted transition between water and ice --

  var weight : double = min((constants.tmelt - t) / tr, 1.0)
  if set_estblh2o then
    phys_tbls[0].estblh2o[es_i] = weight * phys_tbls[0].estblh2o[es_i] + (1.0 - weight) * eswtr
  else
    phys_tbls[0].estbl[es_i] = weight * phys_tbls[0].estbl[es_i] + (1.0 - weight) * eswtr
  end
  itype = itypo
  return itype

end

task radaeini(phys_tbls : region(ispace(int1d), phys_tbls_fs),
              pstdx : double,     -- Standard pressure (dynes/cm^2)
              mwdryx : double,    -- Molecular weight of dry air 
              mwco2x : double)    -- Molecular weight of carbon dioxide
where
  reads writes (phys_tbls)
do

  -- Reads in CAM_ABS_DATA.DBL file? 

  -- Set up table of H2O saturation vapor pressures for use in calculation effective path RH.  
  -- Need separate table from table in wv_saturation because:
  -- (1. Path temperatures can fall below minimum of that table; and
  -- (2. Abs/Emissivity tables are derived with RH for water only.

  var tmin : int = round(constants.min_tp_h2o)
  var tmax : int = round(constants.max_tp_h2o) + 1
  var itype : double = 0
  -- Fortran loop runs from tmin to tmax inclusive and is 1-indexed
  for t = tmin - 1, tmax, 1 do
    var tdbl = t
    itype = gffgch(phys_tbls, tdbl, t-tmin, itype, true)
  end
end

-- Initialize various constants for radiation scheme; note that
-- the radiation scheme uses cgs units.
task radini(phys_tbls : region(ispace(int1d), phys_tbls_fs),
            gravx : double, 
            cpairx : double, 
            epsilox : double, 
            stebolx : double, 
            pstdx : double)
where
  reads writes (phys_tbls)
do
  radaeini(phys_tbls, pstdx, constants.mwdry, constants.mwco2)
end

-- Purpose:
-- Builds saturation vapor pressure table for later lookup procedure.
--
-- Method:
-- Uses Goff & Gratch (1946) relationships to generate the table
-- according to a set of free parameters defined below.  Auxiliary
-- routines are also included for making rapid estimates (well with 1%)
-- of both es and d(es)/dt for the particular table configuration.
--
-- Author: J. Hack
task gestbl(phys_tbls : region(ispace(int1d), phys_tbls_fs))
where 
  reads writes (phys_tbls)
do
  -- Set es table parameters
  phys_tbls[0].tmin = 173.16            -- Minimum temperature (K) in table
  phys_tbls[0].tmax = 375.16            -- Maximum temperature (K) in table
  phys_tbls[0].ttrice = 20.00           -- Trans. range from es over h2o to es over ice
  phys_tbls[0].icephs = true            -- Ice phase (true or false)

  -- Set physical constants required for es calculation
  phys_tbls[0].epsqs = constants.ep_2
  phys_tbls[0].hlatv = 2.501e6          -- latent heat of evaporation ~ J/kg
  phys_tbls[0].hlatf = 3.336e5          -- latent heat of fusion ~ J/kg
  phys_tbls[0].rgasv = constants.R_v
  phys_tbls[0].cp = constants.cp
  phys_tbls[0].tmelt = 273.16           -- freezing T of fresh water ~ K

  phys_tbls[0].lentbl = [int](phys_tbls[0].tmax - phys_tbls[0].tmin + 2.000001)
  --if (phys_tbls[0].lentbl > constants.plenest) then
    --write(6,9000) phys_tbls[0].tmax, phys_tbls[0].tmin, constants.plenest
    --call endrun ('GESTBL')    -- Abnormal termination
  --end

  -- Begin building es table.
  -- Check whether ice phase requested.
  -- If so, set appropriate transition range for temperature
  if (phys_tbls[0].icephs) then
    if (phys_tbls[0].ttrice ~= 0.0) then
      phys_tbls[0].itype = -1.0 * phys_tbls[0].ttrice
    else
      phys_tbls[0].itype = 1
    end
  else
    phys_tbls[0].itype = 0
  end
  --
  var t = phys_tbls[0].tmin - 1.0
  for n = 0, phys_tbls[0].lentbl do
    t += 1.0
    gffgch(phys_tbls, t, n, phys_tbls[0].itype, false)
  end
  --
  for n = phys_tbls[0].lentbl, constants.plenest do
    phys_tbls[0].estbl[n] = -99999.0
  end

  -- Table complete -- Set coefficients for polynomial approximation of
  -- difference between saturation vapor press over water and saturation
  -- pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
  -- is valid in the range -40 < t < 0 (degrees C).

  --                  --- Degree 5 approximation ---
  phys_tbls[0].pcf[0] =  5.04469588506e-01
  phys_tbls[0].pcf[1] = -5.47288442819e+00
  phys_tbls[0].pcf[2] = -3.67471858735e-01
  phys_tbls[0].pcf[3] = -8.95963532403e-03
  phys_tbls[0].pcf[4] = -7.78053686625e-05

  --#if !defined(mpas)
    --if (masterproc) then
      --write(6,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
    --end
  --#endif

  --return

  --9000 format('GESTBL: FATAL ERROR *********************************',/, &
  --            ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
  --            ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
  --            ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
end

-- initialization of saturation vapor pressures:
task esinti(phys_tbls : region(ispace(int1d), phys_tbls_fs))
where 
  reads writes (phys_tbls)
do

  -- Call gestbl to build saturation vapor pressure table.
  gestbl(phys_tbls)
end

-- initialization of climatological monthly-mean ozone profiles.
-- Fortran equivalent can be found in src/core_atmosphere/mpas_atmphys_camrad_init.F
-- TODO: Check input parameters in other files.
task oznini(cr : region(ispace(int2d), cell_fs),
            ozn_region : region(ispace(int2d), ozn_fs))
where
  reads (cr.{lat, lon}),
  reads writes (ozn_region.{pin, ozmixm})
do
  -- Potentially outdated Fortran comment:
  -- This subroutine assumes uniform distribution of ozone concentration.
  -- It should be replaced by monthly climatology that varies latitudinally and vertically

  -- nCells = constants.nCells
  -- num_Months = constants.nMonths
  -- levsiz = constants.nOznLevels
  -- latsiz = constants.latsiz
  -- lonsiz = constants.lonsiz

  -- Assumes TBL files have one value per line.

  -- read in ozone pressure data
  var file = constants.cio.fopen("OZONE_PLEV.TBL", "r")
  regentlib.assert(not isnull(file), "failed to open OZONE_PLEV.TBL file")
  var buff : int8[constants.MAXCHAR] 
  for k = 0, constants.nOznLevels do
    if (not isnull(constants.cio.fgets(buff, constants.MAXCHAR, file))) then 
      var pressure : double = constants.clib.atof(buff)
      -- TODO: How to set up pin?
      -- ozn_region is indexed by {constants.nCells, constants.nOznLevels + 1}
      ozn_region[{0, k}].pin = pressure * 100
    end
  end
  constants.cio.fclose(file)

  -- read in ozone lat data
  var lat_ozone : double[constants.latsiz]

  file = constants.cio.fopen("OZONE_LAT.TBL", "r")
  regentlib.assert(not isnull(file), "failed to open OZONE_LAT.TBL file")
  for i = 0, constants.latsiz do
    if (not isnull(constants.cio.fgets(buff, constants.MAXCHAR, file))) then
      lat_ozone[i] = constants.clib.atof(buff)
    end
  end
  constants.cio.fclose(file)

  -- read in ozone data

  -- TODO: delete this comment.
  -- Can't use region instead of 4d array for simplicity.
  --var i_ozmixin = ispace(int4d, {constants.lonsiz, constants.nOznLevels, constants.latsiz, constants.nMonths })
  --var ozmixin = region(ispace(int4d, {constants.lonsiz, constants.nOznLevels, constants.latsiz, constants.nMonths }), double)
  --var ozmixin = region(i_ozmixin, double)

  -- Using 1d array since indexing into 4d array is complicated and int4d regions need special flags at compile and runtime.
  -- Array represents double[lonsiz][nOznLevels][latsiz][nMonths]
  var ozmixin : double[constants.lonsiz * constants.nOznLevels * constants.latsiz * constants.nMonths]

  file = constants.cio.fopen("OZONE_DAT.TBL", "r")
  regentlib.assert(not isnull(file), "failed to open OZONE_DAT.TBL file")
  for m = 0, constants.nMonths do
    for j = 0, constants.latsiz do
      for k = 0, constants.nOznLevels do
        for i = 0, constants.lonsiz do
          if (not isnull(constants.cio.fgets(buff, constants.MAXCHAR, file))) then
            -- ozmixin[{i, k, j, m}] = constants.clib.atof(buff)
            ozmixin[i + (constants.lonsiz * k) + (constants.lonsiz * constants.nOznLevels * j) 
                    + (constants.lonsiz * constants.nOznLevels * constants.latsiz * m)] = constants.clib.atof(buff)
          end
        end
      end
    end
  end
  constants.cio.fclose(file)

  -- Interpolation of input ozone data to mpas grid.
  var i1 : int
  var i2 : int

  for iCell = 0, constants.nCells do
    -- cell_id_space = {constants.nCells, constants.nVertlevels + 1}
    var lat : double = cr[{iCell, 0}].lat / constants.degrad
    var lon : double = cr[{iCell, 0}].lon / constants.degrad
    if (lat > lat_ozone[constants.latsiz - 1]) then
      i1 = constants.latsiz - 1
      i2 = constants.latsiz - 1
    elseif (lat < lat_ozone[0]) then
      i1 = 0
      i2 = 0
    else
      for i = 0, constants.latsiz do
        i1 = i
        i2 = i + 1
        if (lat >= lat_ozone[i] and lat < lat_ozone[i + 1]) then
          break
        end
      end
    end

    for m = 0, constants.nMonths do
      for k = 0, constants.nOznLevels do
        for j = 0, constants.lonsiz do
          var dlat = lat_ozone[i2] - lat_ozone[i1]
          var dlatCell = lat - lat_ozone[i1]
          -- Future improvement: Also check whether dflat is close enough to 0 to get underflow.
          if (dlat == 0.) then
            --ozn_region[{iCell, k}].ozmixm[m] = ozmixin[{j, k, i1, m}]
            ozn_region[{iCell, k}].ozmixm[m] = ozmixin[j + (constants.lonsiz * k) + (constants.lonsiz * constants.nOznLevels * i1) 
                                                       + (constants.lonsiz * constants.nOznLevels * constants.latsiz * m)]
          else
            --ozn_region[{iCell, k}].ozmixm[m] = ozmixin[{j, k, i1, m}] + (ozmixin[{j, k, i2, m}] 
            --                                 - ozmixin[{j, k, i1, m}]) * dlatCell / dlat
            ozn_region[{iCell, k}].ozmixm[m] = ozmixin[j + (constants.lonsiz * k) + (constants.lonsiz * constants.nOznLevels * i1) 
                                                       + (constants.lonsiz * constants.nOznLevels * constants.latsiz * m)]
                                               + (ozmixin[j + (constants.lonsiz * k) + (constants.lonsiz * constants.nOznLevels * i2) 
                                                            + (constants.lonsiz * constants.nOznLevels * constants.latsiz * m)]
                                               - ozmixin[j + (constants.lonsiz * k) + (constants.lonsiz * constants.nOznLevels * i1) 
                                                           + (constants.lonsiz * constants.nOznLevels * constants.latsiz * m)])
                                               * dlatCell / dlat



          end 
        end
      end
    end

  end
end

task aer_optics_initialize()
end

-- initialization of aerosol concentrations:
task aerosol_init(cr : region(ispace(int2d), cell_fs),
                  er : region(ispace(int2d), edge_fs))
where
  reads writes (cr)
do

  aer_optics_initialize()
end

-- Counterpart of mpas_atmphys_camrad_init.F:camradinit(...)
-- ------
-- Initialization of CAM radiation codes using MPAS MPI decomposition.
-- Laura D. Fowler (send comments to laura@ucar.edu).
-- 2013-05-01.
-- 
--  subroutine camradinit calls the main subroutines needed to initialize the long- and short-wave
--  CAM radiation codes, and read input data from auxillary files.
-- 
--  subroutines called in mpas_atmphys_camrad_init:
--  -----------------------------------------------
--  radini      :initialization of radiation constants.
--  esinti      :initialization of saturation vapor pressures.
--  oznini      :initialization of climatological monthly-mean ozone profiles.
--  aerosol_init:initialization of aerosol optical properties.
-- 
--  add-ons and modifications to sourcecode:
--  ----------------------------------------
--  * added initialization of variable mxaerl which is the number of layers below 900 hPa in
--    which background aerosols are present. mxaerl is computed using the pressure-base array.
--    -> added diag in the argument list of subroutines camradinit and aerosol_init.
--    -> in subroutine aerosol_init, added initialization of variable mxaerl.
--    Laura D. Fowler (birch.ucar.edu) / 2013-07-01.
--  * moved the arrays pin and ozmixm from the mesh structure to the atm_input structure in
--    subroutine oznini.
--    Laura D. Fowler (birch.ucar.edu) / 2013-07-08.
--  * Replaced the variable g (that originally pointed to gravity) with gravity, for simplicity.
--    Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
--  * Modified sourcecode to use pools.
--    Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
-- 
task camradinit(cr : region(ispace(int2d), cell_fs),
                er : region(ispace(int2d), edge_fs),
                phys_tbls : region(ispace(int1d), phys_tbls_fs),
                ozn_region : region(ispace(int2d), ozn_fs))
where reads writes (cr, phys_tbls, ozn_region)
do

  var pstd : double = 101325.0

  radini(phys_tbls, constants.gravity, constants.cp, constants.ep_2, constants.stbolt, pstd*10.0)
  esinti(phys_tbls)
  oznini(cr, ozn_region)
  aerosol_init(cr, er)
end
