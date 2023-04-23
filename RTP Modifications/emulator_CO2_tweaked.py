# python3 source code for CO2/climate emulator
# author: Knut von Salzen, email: knut.vonsalzen@canada.ca
# Edited by Victoria Spada, 2023
# RTP coefficients were modified from the emulator default settings for a 
# set increase to the Arctic RTP coefficients, and corresponding decreases 
# to the other regional coefficients.

import numpy as np
import os 
import sys
import xarray as xr

### read scenario information
# possible choices are:
#   scenario = SSP126 SSP245 SSP370 SSP585 CLE MFR MFR_SDS CFM
#   realm = mean low_ecs high_ecs low_aci high_aci low_ari_bc high_ari_bc low_ari_s high_ari_s low_alb high_alb

lat_ranges = [[0,30],[30,62],[62,118],[118,180]]
scales = np.zeros(4)
for r in range(4):
    scales[r] = np.cos(-lat_ranges[r][1]+lat_ranges[r][0])/2

data_in = np.genfromtxt('emulator.csv', dtype="<U100,<U100",names=True, delimiter=',')
# scenario = str(data_in['scenario'])
# realm = str(data_in['realm'])

scenario = 'SSP585'
realm = 'high_ecs'

# check if SSP or AMAP scenario

cmip6_scenario = [ 'SSP119', 'SSP245', 'SSP585', 'SSP126', 'SSP370' ]
cmip6 = False
for nsc in range(len(cmip6_scenario)):
  # print(scenario, cmpi6_scenario[nsc])
  if scenario == cmip6_scenario[nsc]:
    cmip6 = True

# emulator input: files for emissions, name of the directory that contains the files

if not cmip6:
  emidir = './emissions/AMAP/'
else:
  emidir = './emissions/CMIP6/'

### emulator configuration

scaling = True                     # scale RTPs to match equilibrium climate sensitivity
fires  = False
reference = 'zero'  # reference year

# natural background emissions of CO2 (Tg/year)

emis_co2_nat = 0.

# initial CO2 mixing ratio at the beginning of the time series (ppm)

co2_start_ppm = 284.32              # Meinshausen (2017) estimate for 1850

# diagnostic parameters

thorizon = [20, 100]      # time horizon for ARTP and AGTP diagnostic calculation (years)
year_ref = [2015, 2015]   # reference year for ARTP and AGTP diagnostic calculation (years)

# factors for unit conversion

ppm2Tg = 7.76e+3                              # convert ppm to tera gram (Tg)
kt2Tg = 1.e-03                                # convert ktons to tera gram (Tg)
Tg2Gg = 1.e+03                                # convert tera gram (Tg) to giga gram (Gg)

# specified CO2 emissions (True) or concentrations (False)

co2_emissions = True

# pre-industrial CO2 mixing ratio (for CO2 forcing)

co2_pi_ppm = co2_start_ppm 

# N2O mixing ratio (for CO2 forcing)

n2o_ppm = 323.
n2o_refr_ppm = 323.

# time scales and coefficients for 4 different carbon pools (Millar et al., 2017)

tau_co2 = [1.e+6, 394.4, 36.54, 4.304]
coef_co2 = [0.2173, 0.2240, 0.2824, 0.2763]

# temperature response function parameters (Collins et al., 2013; Shindell et al., 2010)

tr_c = [0.541, 0.368]                    # K (W/m^2)^-1    - modified by Shindell
tr_d = [8.4, 409.5]                      # yr

# target equilibrium climate sensitivity

ecs = 3.7/3.7                  # K (W/m^2)^-1 (CMIP6, mean, Meehl et al., 2021).
if realm == 'low_ecs':
  ecs = ecs * 0.7 
if realm == 'high_ecs':
  ecs = ecs * 1.3

# area fractions of latitude bands
area_frc_ar = 0.0669873                           # 60N-90N (Arctic) 
area_frc_ml = 0.1982769                           # 28N-60N (mid latitudes)
area_frc_tr = 0.4694716                           # 28N-28S (tropics)
area_frc_sh = 0.2652642                           # 28S-90S (southern hemisphere)

# regional temperature sensitivity (K (W/m^2)^-1, Shindell et al., 2010)
cs_co2_ar_old = np.array([0.31, 0.17, 0.16, 0.06])    # 60N-90N (Arctic) 
cs_co2_ml_old = np.array([0.06, 0.24, 0.17, 0.07])    # 28N-60N (mid latitudes)
cs_co2_tr_old = np.array([0.02, 0.10, 0.24, 0.09])    # 28N-28S (tropics)
cs_co2_sh_old = np.array([0.  , 0.02, 0.05, 0.19])    # 28S-90S (southern hemisphere)

# weightings for the RTP coefficient rescaling
D = 0.10 # percent increase for the area weighted sum of the Arctic RTP
delta = D*area_frc_ar*cs_co2_ar_old
d1 = delta/(area_frc_ar*cs_co2_ar_old)
d2 = -delta/(3*area_frc_ml*cs_co2_ml_old)
d3 = -delta/(3*area_frc_tr*cs_co2_tr_old)
d4 = -delta/(3*area_frc_sh*cs_co2_sh_old)

cs_co2_ar = list(np.array([0.31, 0.17, 0.16, 0.06])*(1+d1))    # 60N-90N (Arctic) 
cs_co2_ml = list(np.array([0.06, 0.24, 0.17, 0.07])*(1+d2))    # 28N-60N (mid latitudes)
cs_co2_tr = list(np.array([0.02, 0.10, 0.24, 0.09])*(1+d3))    # 28N-28S (tropics)
cs_co2_sh = list(np.array([0.  , 0.02, 0.05, 0.19])*(1+d4))    # 28S-90S (southern hemisphere)

# account for cs_co2_sh[0] = 0 (get inf.)
d4 = -delta[0]/(2*area_frc_ml*cs_co2_ml_old[0])
d5 = -delta[0]/(2*area_frc_tr*cs_co2_tr_old[0])

cs_co2_ml[0] = cs_co2_ml_old[0]*(1+d4)    # 28N-60N (mid latitudes)
cs_co2_tr[0] = cs_co2_tr_old[0]*(1+d5)
cs_co2_sh[0] = 0

# region and sector labels

nreg = 7
sregion = [ 'ASIA', 'ROEUR', 'AWEST', 'AEAST', 'ROW', 'ARCTIC', 'NONARC' ]
nsec = 4
ssector = [ 'SURF', 'FLAR', 'FIRE', 'SHIP' ]

# numerical parameters

fillv = -9999.

# check if input files exist

def check_input_files(file_name_lcl):
  slcf_file = file_name_lcl
  if os.path.isfile(slcf_file):
    chck = True
  else:
    print ('File '+slcf_file+' does not exist')
    chck = False
  return chck

spchk1 = False
spchk2 = False
spchk3 = False
spchk4 = False
spchk = False
ntag_co2 = None
ntag_co2_refr = None
ndat_co2 = None
ndat_co2_refr = None
year_co2_mltyr = None
year_co2_refr_mltyr = None
year_start = None
year_end = None
nstpyr = None
emis_co2_names = None
emis_co2_refr_names = None
ntag_fire_co2 = None
nstpyr_fire = None
if fires:
  ntag_fire_co2_refr = None
  ndat_fire_co2 = None
  ndat_fire_co2_refr = None
  year_fire_co2_mltyr = None
  year_fire_co2_refr_mltyr = None
  year_fire_start = None
  year_fire_end = None
  emis_fire_co2_names = None
  emis_fire_co2_refr_names = None

# read emissions, first step

print('Scenario :'+scenario)
ssrcspec = 'CO2'
emi_dat_file = emidir+ssrcspec+'_base_'+scenario+'.csv'
emi_refr_dat_file = emidir+ssrcspec+'_'+reference+'.csv'
if fires:
  emi_fire_dat_file = emidir+ssrcspec+'_fire_base_'+scenario+'.csv'
  emi_fire_refr_dat_file = emidir+ssrcspec+'_fire_'+reference+'.csv'

print(emi_dat_file, emi_refr_dat_file)
spchk1 = check_input_files(emi_dat_file)
spchk2 = check_input_files(emi_refr_dat_file)
spchk = spchk1 and spchk2
if fires:
  spchk3 = check_input_files(emi_fire_dat_file)
  spchk4 = check_input_files(emi_fire_refr_dat_file)
  spchk = spchk and spchk3 and spchk4

if spchk:

# anthropogenic emissions

  co2_file = emi_dat_file
  co2_refr_file = emi_refr_dat_file
  emis_co2_in = np.genfromtxt(co2_file, dtype=None, names=True, delimiter=',')      # open file containing time series and read
  lenx = len(emis_co2_in.dtype.names)
  emis_co2_names = emis_co2_in.dtype.names[1:lenx]                                  # emission tag names
  ntag_co2 = len(emis_co2_names)                                                    # number of emission tags
  ndat_co2 = emis_co2_in['Year'].size                                               # number of data points in input time series
  emis_co2_refr_in = np.genfromtxt(co2_refr_file, dtype=None, names=True, delimiter=',')   # open file containing time series and read
  lenx = len(emis_co2_refr_in.dtype.names)
  emis_co2_refr_names = emis_co2_refr_in.dtype.names[1:lenx]                        # emission tag names
  ntag_co2_refr = len(emis_co2_refr_names)                                          # number of emission tags
  ndat_co2_refr = emis_co2_refr_in['Year'].size                                     # number of data points in input time series
  emis_co2_mltyr = [[0. for i in range(ndat_co2)] for j in range(ntag_co2)]
  emis_co2_refr_mltyr = [[0. for i in range(ndat_co2_refr)] for j in range(ntag_co2_refr)]

# fire emissions

  if fires:
    co2_fire_file = emi_fire_dat_file
    co2_fire_refr_file = emi_fire_refr_dat_file
    emis_fire_co2_in = np.genfromtxt(co2_fire_file, dtype=None, names=True, delimiter=',')  # open file containing time series and read
    lenx = len(emis_fire_co2_in.dtype.names)
    emis_fire_co2_names = emis_fire_co2_in.dtype.names[1:lenx]                      # emission tag names
    ntag_fire_co2 = len(emis_fire_co2_names)                                        # number of emission tags
    ndat_fire_co2 = emis_fire_co2_in['Year'].size                                   # number of data points in input time series
    emis_fire_co2_refr_in = np.genfromtxt(co2_fire_refr_file, dtype=None, names=True, delimiter=',')   # open file containing time series and read
    lenx = len(emis_fire_co2_refr_in.dtype.names)
    emis_fire_co2_refr_names = emis_fire_co2_refr_in.dtype.names[1:lenx]            # emission tag names
    ntag_fire_co2_refr = len(emis_fire_co2_refr_names)                              # number of emission tags
    ndat_fire_co2_refr = emis_fire_co2_refr_in['Year'].size                         # number of data points in input time series
    emis_fire_co2_mltyr = [[0. for i in range(ndat_fire_co2)] for j in range(ntag_fire_co2)]
    emis_fire_co2_refr_mltyr = [[0. for i in range(ndat_fire_co2_refr)] for j in range(ntag_fire_co2_refr)]
  else:
    ntag_fire_co2 = 0

# read emissions, second step

if spchk:
  co2_file = emi_dat_file
  co2_refr_file = emi_refr_dat_file
  emis_co2_in = np.genfromtxt(co2_file, dtype=None, names=True, delimiter=',')
  year_co2_mltyr = emis_co2_in['Year']                                              # input times
  year_start = year_co2_mltyr[0]                                                    # first year of the time series
  year_end = year_co2_mltyr[ndat_co2-1]                                             # last year of the time series
  nstpyr = year_end - year_start + 1                                                # total number of years
  for k in range(ntag_co2):
    emis_co2_mltyr[k] = emis_co2_in[emis_co2_names[k]]                              # copy emission time series
  emis_co2_refr_in = np.genfromtxt(co2_refr_file, dtype=None, names=True, delimiter=',')
  year_co2_refr_mltyr = emis_co2_refr_in['Year']                                    # input times
  if emis_co2_refr_names != emis_co2_names or ntag_co2_refr != ntag_co2 \
    or ndat_co2_refr != ndat_co2 or not np.array_equal(year_co2_refr_mltyr,year_co2_mltyr):
    print ('Inconsistent emission input files')
    sys.exit()
  for k in range(ntag_co2_refr):
    emis_co2_refr_mltyr[k] = emis_co2_refr_in[emis_co2_refr_names[k]]               # copy emission time series

# fire emissions

  if fires:
    co2_fire_file = emi_fire_dat_file
    co2_fire_refr_file = emi_fire_refr_dat_file
    emis_fire_co2_in = np.genfromtxt(co2_fire_file, dtype=None, names=True, delimiter=',')
    year_fire_co2_mltyr = emis_fire_co2_in['Year']                                  # input times
    year_fire_start = year_fire_co2_mltyr[0]                                        # first year of the time series
    year_fire_end = year_fire_co2_mltyr[ndat_fire_co2-1]                            # last year of the time series
    nstpyr_fire = year_fire_end - year_fire_start + 1                               # total number of years
    for k in range(ntag_fire_co2):
      emis_fire_co2_mltyr[k] = emis_fire_co2_in[emis_fire_co2_names[k]]             # copy emission time series
    emis_fire_co2_refr_in = np.genfromtxt(co2_fire_refr_file, dtype=None, names=True, delimiter=',')
    year_fire_co2_refr_mltyr = emis_fire_co2_refr_in['Year']                        # input times
    if emis_fire_co2_refr_names != emis_fire_co2_names or ntag_fire_co2_refr != ntag_fire_co2 \
      or ndat_fire_co2_refr != ndat_fire_co2 or not np.array_equal(year_fire_co2_refr_mltyr,year_fire_co2_mltyr):
      print ('Inconsistent fire emission input files')
      sys.exit()
    for k in range(ntag_fire_co2_refr):
      emis_fire_co2_refr_mltyr[k] = emis_fire_co2_refr_in[emis_fire_co2_refr_names[k]] # copy emission time series
  else:
    nstpyr_fire = nstpyr

ina=-9
nstpyra=ina
if nstpyra == ina:
  nstpyra = nstpyr
if nstpyra != nstpyr or nstpyra != nstpyr_fire:
  print('All emission input files must produce the same number of years')
  sys.exit()

print('Available regions and sectors in input model data file:')
for nr in range(nreg):
  regionp = sregion[nr]
  for ns in range(nsec):
    sectorp = ssector[ns]
    print(regionp+', '+sectorp)

print('Available regions and sectors in emission data file:')
if spchk:
  for k in range(ntag_co2):
    print(ssrcspec+': '+emis_co2_names[k])
  for k in range(ntag_fire_co2):
    print(ssrcspec+': '+emis_fire_co2_names[k])

if spchk:
  isec = [None for i in range(ntag_co2+ntag_fire_co2)]
  ireg = [None for i in range(ntag_co2+ntag_fire_co2)]
for ns in range(nsec):
  sectorp = ssector[ns]
  for nr in range(nreg):
    regionp = sregion[nr]
    emi_tag1 = '{}_{}'.format(sectorp,regionp)              # the header in the emission input file must
    emi_tag2 = '{}_{}'.format(regionp,sectorp)              # have the format "sector_region" or "region_sector" (capital letters)
    if spchk:
      for k in range(ntag_co2):
        if emis_co2_names[k] == emi_tag1 or emis_co2_names[k] == emi_tag2:
          isec[k] = ns                                 # matching sector index in netcdf input data
          ireg[k] = nr                                 # matching region index in netcdf input data
      for k in range(ntag_fire_co2):
        kx = k + ntag_co2
        if emis_fire_co2_names[k] == emi_tag1 or emis_fire_co2_names[k] == emi_tag2:
          isec[kx] = ns                                # matching sector index in netcdf input data
          ireg[kx] = nr                                # matching region index in netcdf input data

print('Processing the following regions and sectors:')
if spchk:
  ifound = 0
  ifoundf = 0
  emi_dat_file = emidir+ssrcspec+'_base_'+scenario+'.csv'
  if fires:
    emi_fire_dat_file = emidir+ssrcspec+'_fire_base_'+scenario+'.csv'
  for k in range(ntag_co2):
    if isec[k] is not None:
      print(ssrcspec+': '+sregion[ireg[k]]+'_'+ssector[isec[k]])          # print emission tags that are being processed
      ifound = ifound + 1
  for k in range(ntag_fire_co2):
    kx = k + ntag_co2
    if isec[kx] is not None:
      print(ssrcspec+': '+sregion[ireg[kx]]+'_'+ssector[isec[kx]])        # print emission tags that are being processed
      ifoundf = ifoundf + 1
  if ifound == 0:
    print('Regions and sectors in file '+emi_dat_file+' does not match information in the input model data file '+model_dat_file)
  if ifoundf == 0 and fires:
    print('Regions and sectors in file '+emi_fire_dat_file+' does not match information in the input model data file '+model_dat_file)

# annual mean emissions from linear interpolation of emissions in time

dtyr = 1
dtyrf = 1.
def interpol_emis(emis_mltyr_lcl,year_mltyr_lcl,year_start_lcl,ndat_lcl,nstpyr_lcl):
  emis_yr = [0. for i in range(nstpyr_lcl)]
  for nt in range(nstpyr_lcl):
    year = year_start_lcl + nt * dtyr
    for nmltyr in range(ndat_lcl-1):
      nmltyrp = nmltyr + 1
      if year >= year_mltyr_lcl[nmltyr] and year < year_mltyr_lcl[nmltyrp]:
         wgt = float(year - year_mltyr_lcl[nmltyr])/float(year_mltyr_lcl[nmltyrp] - year_mltyr_lcl[nmltyr])
         ne = nmltyr
         nep = nmltyrp
    if year == year_mltyr_lcl[ndat_lcl-1]:
       wgt = 0.
       ne = ndat_lcl-1
       nep = ne
    emis_yr[nt] = (1. - wgt) * emis_mltyr_lcl[ne] + wgt * emis_mltyr_lcl[nep]
  return emis_yr

if spchk:
  emis_co2_yr = [[0. for i in range(nstpyr)] for j in range(ntag_co2+ntag_fire_co2)]
  emis_co2_refr_yr = [[0. for i in range(nstpyr)] for j in range(ntag_co2+ntag_fire_co2)]
if spchk:
  for k in range(ntag_co2):
    if isec[k] is not None:
      emis_co2_yr[k] = interpol_emis(emis_co2_mltyr[k],year_co2_mltyr,year_start,ndat_co2,nstpyr)
      emis_co2_refr_yr[k] = interpol_emis(emis_co2_refr_mltyr[k],year_co2_mltyr,year_start,ndat_co2_refr,nstpyr)
  for k in range(ntag_fire_co2):
    kx = k + ntag_co2
    if isec[kx] is not None:
      emis_co2_yr[kx] = interpol_emis(emis_fire_co2_mltyr[k],year_fire_co2_mltyr,year_start,ndat_fire_co2,nstpyr)
      emis_co2_refr_yr[kx] = interpol_emis(emis_fire_co2_refr_mltyr[k],year_fire_co2_mltyr,year_start,ndat_fire_co2_refr,nstpyr)

if spchk:
  emis_co2_names_tmp = ["" for j in range(ntag_co2+ntag_fire_co2)]
if spchk:
  for k in range(ntag_co2):
    if isec[k] is not None:
      emis_co2_names_tmp[k] = emis_co2_names[k] 
  for k in range(ntag_fire_co2):
    kx = k + ntag_co2
    if isec[kx] is not None:
      emis_co2_names_tmp[kx] = emis_fire_co2_names[k] 
emis_co2_names = emis_co2_names_tmp

ntag_co2 = ntag_co2 + ntag_fire_co2

# unit conversion

if spchk:
  for i in range(nstpyr):
    emisum = 0.
    for k in range(ntag_co2):
      emis_co2_yr[k][i] = emis_co2_yr[k][i] * kt2Tg 
      emis_co2_refr_yr[k][i] = emis_co2_refr_yr[k][i] * kt2Tg 
      emisum = emisum + emis_co2_yr[k][i]

if not co2_emissions:
  print ('Option not implemented')
  sys.exit()

# implicit equilibrium climate sensitivity of underlying model
ecs0 = sum(tr_c) 

# regional temperature sensitivity for different latitude bands, assuming globally uniform radiative forcing

cs_co2_ar_tot = sum(cs_co2_ar)
cs_co2_ml_tot = sum(cs_co2_ml)
cs_co2_tr_tot = sum(cs_co2_tr)
cs_co2_sh_tot = sum(cs_co2_sh)
aaaa = cs_co2_ar_tot + cs_co2_ml_tot + cs_co2_tr_tot + cs_co2_sh_tot

# loop over all years to calculate annual mean CO2 concentrations and radiative forcings

dt = 0.01
time = np.zeros((nstpyra), dtype=np.float64)
time = xr.DataArray(time, name='time', dims=('time'), coords={'time':time}, attrs={'long_name':"Time",'standard_name':"time", 'units':"year"})
if spchk:
  for i in range(nstpyra):
    time[i]= year_start + i * dtyrf

def co2_concp0(co2_lcl,emis_co2_lcl):
  coeff_lcl = coef_co2[0]
  emis_lcl = coeff_lcl * emis_co2_lcl
  conc = co2_lcl + emis_lcl * dt
  return conc

def co2_concpx(co2_lcl,emis_co2_lcl,tau_lcl,coeff_lcl):
  emis_lcl = coeff_lcl * emis_co2_lcl
  method = 1
  if method == 1:       # analytic solution
    conc = ( co2_lcl - emis_lcl * tau_lcl ) * np.exp(-dt/tau_lcl) + emis_lcl * tau_lcl
  else:                 # semi-implicit solution
    conc = ( co2_lcl/dt + emis_lcl ) / ( 1./dt + 1./tau_lcl )
  return conc

def co2_radiative_forcing(co2_ppm_lcl,n2o_ppm_lcl):
  rf = ( -2.4e-7 * (co2_ppm_lcl - co2_pi_ppm)**2 +7.2e-4 * np.abs(co2_ppm_lcl - co2_pi_ppm) -2.1e-4 * n2o_ppm_lcl + 5.36 ) * np.log(co2_ppm_lcl/co2_pi_ppm)
  return rf 

def interpol_lin(emis_yr_lcl,wgt_lcl,year_lcl,yearm_lcl):
  if year_lcl <= 0:
    emis = emis_yr_lcl[year_lcl]
  elif year_lcl > year_end - year_start:
    emis = emis_yr_lcl[yearm_lcl]
  else:
    emis = wgt_lcl * emis_yr_lcl[year_lcl] + (1.-wgt_lcl) * emis_yr_lcl[yearm_lcl]
  return emis

nstpyrt = nstpyr
co2_yr = np.zeros((nstpyrt+1))
co2_refr_yr = np.zeros((nstpyrt+1))
rf_yr = [0. for i in range(nstpyrt+1)]
rf_refr_yr = [0. for i in range(nstpyrt+1)]
rf_pref_yr = [0. for i in range(nstpyrt+1)]
emis_tmp= [0. for i in range(nstpyrt+1)]
emis_yr = [0. for i in range(nstpyrt+1)]
emis_refr_yr = [0. for i in range(nstpyrt+1)]
nstp = 1 + int(round(nstpyrt / dt))
nstppyr = nstp / nstpyrt
yrwgt = 1. / float(nstppyr)
conc_pref_ppm = 1.

co2_start = co2_start_ppm * ppm2Tg
co2p0 = co2_start
co2p1 = 0.
co2p2 = 0.
co2p3 = 0.
co2p0_refr = co2_start
co2p1_refr = 0.
co2p2_refr = 0.
co2p3_refr = 0.
co2_ppm = co2p0 / ppm2Tg
co2_refr_ppm = co2p0_refr / ppm2Tg
co2_pert = np.full((ntag_co2), co2_start)
co2_pertp0 = np.full((ntag_co2), co2_start)
co2_pertp1 = np.full((ntag_co2), 0.)
co2_pertp2 = np.full((ntag_co2), 0.)
co2_pertp3 = np.full((ntag_co2), 0.)
co2_pert_ppm = np.full((ntag_co2), co2_start)
co2_pert_ppm = np.full((ntag_co2), co2_ppm)
co2_pert_yr = np.zeros((ntag_co2,nstpyrt+1))
rf_pert_yr = np.zeros((ntag_co2,nstpyrt+1))
rf_pert_pref_yr = np.zeros((ntag_co2,nstpyrt+1))
emis_pert_yr = np.zeros((ntag_co2,nstpyrt+1))

# CO2 abundance and radiative forcing time series

for n in range(nstp):
  timet = n * dt
  timeyr = year_start + timet
  year = int(round(timet))
  yearc = int(np.floor(timet))
  yearm = year - 1
  yearp = year + 1
  wgt = timet - float(yearm) - 0.5
  emis_co2_all = emis_co2_nat
  emis_co2_refr_all = 0.
  if spchk:
    for kt in range(ntag_co2):                                                            # loop over all sectors and regions
      if isec[kt] is not None:
        emis_tmp = emis_co2_yr[kt]
        emis_co2 = interpol_lin(emis_tmp,wgt,year,yearm)                          # linearly time-interpolated emissions for control scenario
        emis_co2_all = emis_co2_all + emis_co2
        emis_tmp = emis_co2_refr_yr[kt]
        emis_co2_refr =  interpol_lin(emis_tmp,wgt,year,yearm)                # linearly time-interpolated emissions for baseine scenario
        emis_co2_refr_all = emis_co2_refr_all + emis_co2_refr
  emis_yr[yearc] = emis_yr[yearc] + emis_co2_all * yrwgt
  emis_refr_yr[yearc] = emis_refr_yr[yearc] + emis_co2_refr_all * yrwgt
  co2p0 = co2_concpx(co2p0,emis_co2_all,tau_co2[0],coef_co2[0])
  co2p1 = co2_concpx(co2p1,emis_co2_all,tau_co2[1],coef_co2[1])
  co2p2 = co2_concpx(co2p2,emis_co2_all,tau_co2[2],coef_co2[2])
  co2p3 = co2_concpx(co2p3,emis_co2_all,tau_co2[3],coef_co2[3])
  co2 = co2p0 + co2p1 + co2p2 + co2p3
  co2p0_refr = co2_concpx(co2p0_refr,emis_co2_refr_all,tau_co2[0],coef_co2[0])
  co2p1_refr = co2_concpx(co2p1_refr,emis_co2_refr_all,tau_co2[1],coef_co2[1])
  co2p2_refr = co2_concpx(co2p2_refr,emis_co2_refr_all,tau_co2[2],coef_co2[2])
  co2p3_refr = co2_concpx(co2p3_refr,emis_co2_refr_all,tau_co2[3],coef_co2[3])
  co2_refr = co2p0_refr + co2p1_refr + co2p2_refr + co2p3_refr
  co2_ppm = co2 / ppm2Tg
  co2_refr_ppm = co2_refr / ppm2Tg
  co2_yr[yearc] = co2_yr[yearc] + co2_ppm * yrwgt
  co2_refr_yr[yearc] = co2_refr_yr[yearc] + co2_refr_ppm * yrwgt
  co2_rf = co2_radiative_forcing(co2_ppm,n2o_ppm)
  co2_rf_refr = co2_radiative_forcing(co2_refr_ppm,n2o_refr_ppm)
  co2_pref_ppm = co2_ppm + conc_pref_ppm
  co2_rf_pref = co2_radiative_forcing(co2_pref_ppm,n2o_refr_ppm)
  rf_yr[yearc] = rf_yr[yearc] + co2_rf * yrwgt
  rf_refr_yr[yearc] = rf_refr_yr[yearc] + co2_rf_refr * yrwgt
  rf_pref_yr[yearc] = rf_pref_yr[yearc] + co2_rf_pref * yrwgt
  if spchk:
    for kt in range(ntag_co2):                                                            # loop over all sectors and regions
      if isec[kt] is not None:
        for i in range(nstpyr):
          emis_tmp[i] = emis_co2_yr[kt][i] - emis_co2_refr_yr[kt][i]
        emis_co2 =  interpol_lin(emis_tmp,wgt,year,yearm)                          # linearly time-interpolated emissions for control scenario
        emis_pert_yr[kt,yearc] = emis_pert_yr[kt,yearc] + emis_co2 * yrwgt
        emis_co2_pert = emis_co2_all - emis_co2
        emis_co2_sel = emis_co2_pert
        co2_pertp0[kt] = co2_concpx(co2_pertp0[kt],emis_co2_sel,tau_co2[0],coef_co2[0])
        co2_pertp1[kt] = co2_concpx(co2_pertp1[kt],emis_co2_sel,tau_co2[1],coef_co2[1])
        co2_pertp2[kt] = co2_concpx(co2_pertp2[kt],emis_co2_sel,tau_co2[2],coef_co2[2])
        co2_pertp3[kt] = co2_concpx(co2_pertp3[kt],emis_co2_sel,tau_co2[3],coef_co2[3])
        co2_pert[kt] = co2_pertp0[kt] + co2_pertp1[kt] + co2_pertp2[kt] + co2_pertp3[kt] 
        co2_pert_ppm[kt] = co2_pert[kt] / ppm2Tg
        co2_pert_yr[kt,yearc] = co2_pert_yr[kt,yearc] + co2_pert_ppm[kt] * yrwgt
        co2_rf = co2_radiative_forcing(co2_pert_ppm[kt],n2o_ppm)
        rf_pert_yr[kt,yearc] = rf_pert_yr[kt,yearc] + co2_rf * yrwgt
        co2_pref_ppm = co2_pert_ppm[kt] + conc_pref_ppm
        co2_rf_pref = co2_radiative_forcing(co2_pref_ppm,n2o_refr_ppm)
        rf_pert_pref_yr[kt,yearc] = rf_pert_pref_yr[kt,yearc] + co2_rf_pref * yrwgt

nhor = len(year_ref)
agtp_hor = np.zeros((nhor))
lats = 4
artp_hor = np.zeros((lats,nhor))
scl = np.full((nstpyrt), 1.)
temp_co2_ar = np.zeros((nstpyrt))
temp_co2_ml = np.zeros((nstpyrt))
temp_co2_tr = np.zeros((nstpyrt))
temp_co2_sh = np.zeros((nstpyrt))
demis_pert_i = np.zeros((ntag_co2))
temp_co2_pert_ar = np.zeros((ntag_co2,nstpyrt))
temp_co2_pert_ml = np.zeros((ntag_co2,nstpyrt))
temp_co2_pert_tr = np.zeros((ntag_co2,nstpyrt))
temp_co2_pert_sh = np.zeros((ntag_co2,nstpyrt))

# ARTP model time-integration

for i in range(nstpyrt):                                                                 # loop over years of pulse emission
  timeyr_i = year_start + i * dtyrf
  demis_i = (emis_yr[i] - emis_refr_yr[i]) * dtyrf                                       # emitted amount of CO2 in one year
  drf = rf_yr[i] - rf_pref_yr[i]
  dconc = -conc_pref_ppm * ppm2Tg
  frceff = drf/dconc                                                                     # forcing efficiency of co2
  drf_ar = drf                                                                           # forcing in the Arctic
  drf_ml = drf                                                                           # forcing at mid latitudes
  drf_tr = drf                                                                           # forcing in the tropics
  drf_sh = drf                                                                           # forcing in the southern hemisphere
  csh_co2_ar = cs_co2_ar[0] * (drf_ar/drf) + cs_co2_ar[1] * (drf_ml/drf) \
             + cs_co2_ar[2] * (drf_tr/drf) + cs_co2_ar[3] * (drf_sh/drf)
  csh_co2_ml = cs_co2_ml[0] * (drf_ar/drf) + cs_co2_ml[1] * (drf_ml/drf) \
             + cs_co2_ml[2] * (drf_tr/drf) + cs_co2_ml[3] * (drf_sh/drf)
  csh_co2_tr = cs_co2_tr[0] * (drf_ar/drf) + cs_co2_tr[1] * (drf_ml/drf) \
             + cs_co2_tr[2] * (drf_tr/drf) + cs_co2_tr[3] * (drf_sh/drf)
  csh_co2_sh = cs_co2_sh[0] * (drf_ar/drf) + cs_co2_sh[1] * (drf_ml/drf) \
             + cs_co2_sh[2] * (drf_tr/drf) + cs_co2_sh[3] * (drf_sh/drf)
  agtp_co2_est = csh_co2_ar * area_frc_ar + csh_co2_ml * area_frc_ml \
               + csh_co2_tr * area_frc_tr + csh_co2_sh * area_frc_sh 
  if scaling:
    scl[i] = ecs0 / agtp_co2_est                                                          # scale regional ARTPs in order to match global ARTP
  if spchk:
    for kt in range(ntag_co2):                                                            # loop over all sectors and regions
      if isec[kt] is not None:
        demis_pert_i[kt] = demis_i                                                        # emitted amount of CO2 in one year
        demis_pert_i[kt] = demis_pert_i[kt] - emis_pert_yr[kt,i] * dtyrf
  for j in range(i+1,nstpyrt):                                                            # loop over years of temperature response to pulse emission
    timeyr_j = year_start + j * dtyrf
    dt_emis = (j - i) * dtyrf
    tr = 0.
    for n in range(2):
      for nj in range(4):
        tr = tr + (coef_co2[nj] * tr_c[n] * tau_co2[nj])/(tau_co2[nj] - tr_d[n]) \
                    * ( np.exp(-dt_emis/tau_co2[nj]) - np.exp(-dt_emis/tr_d[n]) )
    tr = tr * (ecs/ecs0)                                                                 # scale in order to match specified equilibrium climate sensitivity 
    agtp_co2 = frceff * tr                                                               # absolute global temperature potential of CO2
    artp_co2_ar = (agtp_co2/ecs0) * csh_co2_ar * scl[i]                                  # arctic absolute regional temperature potential of CO2
    artp_co2_ml = (agtp_co2/ecs0) * csh_co2_ml * scl[i]                                  # mid-latitude absolute regional temperature potential of CO2
    artp_co2_tr = (agtp_co2/ecs0) * csh_co2_tr * scl[i]                                  # tropics absolute regional temperature potential of CO2
    artp_co2_sh = (agtp_co2/ecs0) * csh_co2_sh * scl[i]                                  # southern hemisphere absolute regional temperature potential of CO2
    temp_co2_ar[j] = temp_co2_ar[j] + artp_co2_ar * demis_i                              # arctic temperature response to pulse emission
    temp_co2_ml[j] = temp_co2_ml[j] + artp_co2_ml * demis_i                              # mid-latitude temperature response to pulse emission
    temp_co2_tr[j] = temp_co2_tr[j] + artp_co2_tr * demis_i                              # tropics temperature response to pulse emission
    temp_co2_sh[j] = temp_co2_sh[j] + artp_co2_sh * demis_i                              # southern hemisphere temperature response to pulse emission
    if spchk:
      for kt in range(ntag_co2):                                                           # loop over all sectors and regions
        if isec[kt] is not None:
          temp_co2_pert_ar[kt,j] = temp_co2_pert_ar[kt,j] + artp_co2_ar * demis_pert_i[kt] # arctic temperature response to pulse emission
          temp_co2_pert_ml[kt,j] = temp_co2_pert_ml[kt,j] + artp_co2_ml * demis_pert_i[kt] # mid-latitude temperature response to pulse emission
          temp_co2_pert_tr[kt,j] = temp_co2_pert_tr[kt,j] + artp_co2_tr * demis_pert_i[kt] # tropics temperature response to pulse emission
          temp_co2_pert_sh[kt,j] = temp_co2_pert_sh[kt,j] + artp_co2_sh * demis_pert_i[kt] # southern hemisphere temperature response to pulse emission

  for ny in range(nhor):                                                                   # diagnose results for different specified time horizons
    if year_ref[ny] == timeyr_i:
      dt_emis = float(thorizon[ny])
      tr = 0.
      for n in range(2):
        for nj in range(4):
          tr = tr + (coef_co2[nj] * tr_c[n] * tau_co2[nj])/(tau_co2[nj] - tr_d[n]) \
                      * ( np.exp(-dt_emis/tau_co2[nj]) - np.exp(-dt_emis/tr_d[n]) )
      tr = tr * (ecs/ecs0)                                                                 # scale in order to match specified equilibrium climate sensitivity 
      agtp_hor[ny] = frceff * tr
      artp_hor[0,ny] = (agtp_hor[ny]/ecs0) * csh_co2_ar * scl[i]
      artp_hor[1,ny] = (agtp_hor[ny]/ecs0) * csh_co2_ml * scl[i]
      artp_hor[2,ny] = (agtp_hor[ny]/ecs0) * csh_co2_tr * scl[i]
      artp_hor[3,ny] = (agtp_hor[ny]/ecs0) * csh_co2_sh * scl[i]

# print ARTPs and AGTPs
      
slats = ['ARCTIC', 'NHMDLAT', 'TROPICS', 'SH']
print('AGTP and ARTP at different time horizons (mK/Tg):')
for ny in range(nhor):
  print(year_ref[ny],'AGTP',thorizon[ny],agtp_hor[ny]*1000.)
  print('ARTP')
  for nl in range(lats):
    print(year_ref[ny],thorizon[ny],slats[nl],artp_hor[nl,ny]*1000.)

# equilibrium temperatures
    
temp_co2_eq_ar = [0. for i in range(nstpyrt)]
temp_co2_eq_ml = [0. for i in range(nstpyrt)]
temp_co2_eq_tr = [0. for i in range(nstpyrt)]
temp_co2_eq_sh = [0. for i in range(nstpyrt)]
temp_co2_eq_glb = [0. for i in range(nstpyrt)]
temp_co2_glb = [0. for i in range(nstpyrt)]
for i in range(nstpyrt):
  drf = rf_yr[i] - rf_refr_yr[i]
  drf_ar = drf                                                                           # forcing in the Arctic
  drf_ml = drf                                                                           # forcing at mid latitudes
  drf_tr = drf                                                                           # forcing in the tropics
  drf_sh = drf                                                                           # forcing in the southern hemisphere
  csh_co2_ar = cs_co2_ar[0] * (drf_ar) + cs_co2_ar[1] * (drf_ml) \
             + cs_co2_ar[2] * (drf_tr) + cs_co2_ar[3] * (drf_sh)
  csh_co2_ml = cs_co2_ml[0] * (drf_ar) + cs_co2_ml[1] * (drf_ml) \
             + cs_co2_ml[2] * (drf_tr) + cs_co2_ml[3] * (drf_sh)
  csh_co2_tr = cs_co2_tr[0] * (drf_ar) + cs_co2_tr[1] * (drf_ml) \
             + cs_co2_tr[2] * (drf_tr) + cs_co2_tr[3] * (drf_sh)
  csh_co2_sh = cs_co2_sh[0] * (drf_ar) + cs_co2_sh[1] * (drf_ml) \
             + cs_co2_sh[2] * (drf_tr) + cs_co2_sh[3] * (drf_sh)
  tr = ecs/ecs0
  temp_co2_eq_ar[i] = csh_co2_ar * scl[i] * tr
  temp_co2_eq_ml[i] = csh_co2_ml * scl[i] * tr
  temp_co2_eq_tr[i] = csh_co2_tr * scl[i] * tr
  temp_co2_eq_sh[i] = csh_co2_sh * scl[i] * tr
  temp_co2_eq_glb[i] = temp_co2_eq_ar[i] * area_frc_ar + temp_co2_eq_ml[i] * area_frc_ml \
                     + temp_co2_eq_tr[i] * area_frc_tr + temp_co2_eq_sh[i] * area_frc_sh 
  temp_co2_glb[i] = temp_co2_ar[i] * area_frc_ar + temp_co2_ml[i] * area_frc_ml \
                  + temp_co2_tr[i] * area_frc_tr + temp_co2_sh[i] * area_frc_sh 
  year = year_start + i * dtyr

# output fields
  
temp_ARCTIC = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
temp_NHML = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
temp_TROPICS = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
temp_SH = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
temp_glb = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
emissions = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
emissions_nat = np.full((nstpyra), fillv, dtype=np.float64)
conc_co2 = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
forc_co2 = np.full((nstpyra,nreg,nsec), fillv, dtype=np.float64)
conc_co2_nat = np.full((nstpyra), fillv, dtype=np.float64)
conc_co2_base = np.full((nstpyra), fillv, dtype=np.float64)
forc_co2_base = np.full((nstpyra), fillv, dtype=np.float64)
emissions_co2_base = np.full((nstpyra), fillv, dtype=np.float64)
artp_ARCTIC = np.full((nhor), fillv, dtype=np.float64)
artp_NHML = np.full((nhor), fillv, dtype=np.float64)
artp_TROPICS = np.full((nhor), fillv, dtype=np.float64)
artp_SH = np.full((nhor), fillv, dtype=np.float64)
agtp = np.full((nhor), fillv, dtype=np.float64)
temp_eq_ARCTIC = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_NHML = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_TROPICS = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_SH = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_glb = np.full((nstpyra), fillv, dtype=np.float64)
for kt in range(ntag_co2):                                                       # loop over all sectors and regions
  if isec[kt] is not None:
    ns = isec[kt]
    nr = ireg[kt]
    for i in range(nstpyra):
      temp_ARCTIC[i,nr,ns]  = temp_co2_ar[i] - temp_co2_pert_ar[kt,i]       # temperature changes due to dir. radiative forcing
      temp_NHML[i,nr,ns]    = temp_co2_ml[i] - temp_co2_pert_ml[kt,i]
      temp_TROPICS[i,nr,ns] = temp_co2_tr[i] - temp_co2_pert_tr[kt,i]
      temp_SH[i,nr,ns]      = temp_co2_sh[i] - temp_co2_pert_sh[kt,i]
      temp_glb[i,nr,ns]     = temp_ARCTIC[i,nr,ns] * area_frc_ar \
                                    + temp_NHML[i,nr,ns] * area_frc_ml \
                                    + temp_TROPICS[i,nr,ns] * area_frc_tr \
                                    + temp_SH[i,nr,ns] * area_frc_sh
      conc_co2[i,nr,ns]        = co2_yr[i] - co2_pert_yr[kt,i]
      forc_co2[i,nr,ns]        = rf_yr[i] - rf_pert_yr[kt,i]
      emissions[i,nr,ns]       = emis_pert_yr[kt,i] * Tg2Gg
      conc_co2_base[i]  = co2_yr[i]
      forc_co2_base[i]  = rf_yr[i]
      emissions_co2_base[i] = emis_yr[i] * Tg2Gg
      emissions_nat[i] = emis_co2_nat
      conc_co2_nat[i] = co2_start_ppm
      artp_ARCTIC[:] = 1000. * artp_hor[0,:]
      artp_NHML[:] = 1000. * artp_hor[1,:]
      artp_TROPICS[:] = 1000. * artp_hor[2,:]
      artp_SH[:] = 1000. * artp_hor[3,:]
      agtp[:] = 1000. * agtp_hor[:]
for i in range(nstpyra):
  temp_eq_ARCTIC[i]  = temp_co2_eq_ar[i]       # equilibrium temperature changes due to dir. radiative forcing
  temp_eq_NHML[i]    = temp_co2_eq_ml[i]
  temp_eq_TROPICS[i] = temp_co2_eq_tr[i]
  temp_eq_SH[i]      = temp_co2_eq_sh[i]
  temp_eq_glb[i]     = temp_co2_eq_glb[i]

# write output files

region  = np.arange((nreg), dtype=np.float64)
sector  = np.arange((nsec), dtype=np.float64)
region[:] = region[:] + 1
sector[:] = sector[:] + 1
region = xr.DataArray(region,name='region', dims=('region'), coords={'region':region}, \
                      attrs={'long_name':"Source regions: 1 - ASIA, 2 - ROEUR, 3 - AWEST, 4 - AEAST, 5 - ROW, 6 - ARCTIC, 7 - NONARC", 'standard_name':"region"})
sector = xr.DataArray(sector,name='sector', dims=('sector'), coords={'sector':sector}, \
                      attrs={'long_name':"Source sectors: 1 - SURF, 2 - FLAR, 3 - FIRE, 4 - SHIP", 'standard_name':"sector"})
thorizon = xr.DataArray(thorizon, name='thorizon', dims=('thorizon'), coords={'thorizon':thorizon}, attrs={'long_name':"Time horizon",'standard_name':"time_horizon", 'units':"year"})
ds1 = xr.Dataset({'temp_ARCTIC': (['time','region','sector'], temp_ARCTIC)}, coords={'time':time,'region':region,'sector':sector})
ds1.temp_ARCTIC.attrs = {('long_name', 'Arctic temperature perturbation associated with regional CO2 emissions'),
                 ('standard_name', 'Arctic_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds2 = xr.Dataset({'temp_NHML': (['time','region','sector'], temp_NHML)}, coords={'time':time,'region':region,'sector':sector})
ds2.temp_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude temperature perturbation associated with regional CO2 emissions'),
                 ('standard_name', 'NH_midlat_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds3 = xr.Dataset({'temp_TROPICS': (['time','region','sector'], temp_TROPICS)}, coords={'time':time,'region':region,'sector':sector})
ds3.temp_TROPICS.attrs = {('long_name', 'Tropical temperature perturbation associated with regional CO2 emissions'),
                 ('standard_name', 'Tropics_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds4 = xr.Dataset({'temp_SH': (['time','region','sector'], temp_SH)}, coords={'time':time,'region':region,'sector':sector})
ds4.temp_SH.attrs = {('long_name', 'Southern Hemisphere temperature perturbation associated with regional CO2 emissions'),
                 ('standard_name', 'SH_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds5 = xr.Dataset({'temp_glb': (['time','region','sector'], temp_glb)}, coords={'time':time,'region':region,'sector':sector})
ds5.temp_glb.attrs = {('long_name', 'Global temperature perturbation associated with regional CO2 emissions'),
                 ('standard_name', 'global_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds6 = xr.Dataset({'emissions': (['time','region','sector'], emissions)}, coords={'time':time,'region':region,'sector':sector})
ds6.emissions.attrs = {('long_name', 'Regional CO2 emissions'),
                 ('standard_name', 'regional_CO2_emissions'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds7 = xr.Dataset({'emissions_nat': (['time'], emissions_nat)}, coords={'time':time})
ds7.emissions_nat.attrs = {('long_name', 'Natural CO2 emissions'),
                 ('standard_name', 'natural_CO2_emissions'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds8 = xr.Dataset({'conc_co2': (['time','region','sector'], conc_co2)}, coords={'time':time,'region':region,'sector':sector})
ds8.conc_co2.attrs = {('long_name', 'Change in concentration of CO2'),
                 ('standard_name', 'conc_co2'),
                 ('units', 'ppm'),
                 ('_FillValue', fillv)}
ds9 = xr.Dataset({'forc_co2': (['time','region','sector'], forc_co2)}, coords={'time':time,'region':region,'sector':sector})
ds9.forc_co2.attrs = {('long_name', 'Change in direct radiative forcing of CO2'),
                 ('standard_name', 'forc_co2'),
                 ('units', 'W/m2'),
                 ('_FillValue', fillv)}
ds10 = xr.Dataset({'artp_ARCTIC': (['thorizon'], artp_ARCTIC)}, coords={'thorizon':thorizon})
ds10.artp_ARCTIC.attrs = {('long_name', 'Arctic absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_ARCTIC'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds11 = xr.Dataset({'artp_NHML': (['thorizon'], artp_NHML)}, coords={'thorizon':thorizon})
ds11.artp_NHML.attrs = {('long_name', 'NH midlatitude absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_NHML'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds12 = xr.Dataset({'artp_TROPICS': (['thorizon'], artp_TROPICS)}, coords={'thorizon':thorizon})
ds12.artp_TROPICS.attrs = {('long_name', 'Tropics absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_TROPICS'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds13 = xr.Dataset({'artp_SH': (['thorizon'], artp_SH)}, coords={'thorizon':thorizon})
ds13.artp_SH.attrs = {('long_name', 'Southern Hemisphere absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_SH'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds14 = xr.Dataset({'agtp': (['thorizon'], agtp)}, coords={'thorizon':thorizon})
ds14.agtp.attrs = {('long_name', 'Absolute global temperarture change potential'),
                 ('standard_name', 'AGTP'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds15 = xr.Dataset({'conc_co2_nat': (['time'], conc_co2_nat)}, coords={'time':time})
ds15.conc_co2_nat.attrs = {('long_name', 'Natural CO2 concentration'),
                 ('standard_name', 'natural_CO2_conc'),
                 ('units', 'ppm'),
                 ('_FillValue', fillv)}
ds16 = xr.Dataset({'conc_co2_base': (['time'], conc_co2_base)}, coords={'time':time})
ds16.conc_co2_base.attrs = {('long_name', 'Concentration of CO2'),
                 ('standard_name', 'conc_co2'),
                 ('units', 'ppm'),
                 ('_FillValue', fillv)}
ds17 = xr.Dataset({'forc_co2_base': (['time'], forc_co2_base)}, coords={'time':time})
ds17.forc_co2_base.attrs = {('long_name', 'Direct radiative forcing of CO2'),
                 ('standard_name', 'forcing_co2'),
                 ('units', 'W/m2'),
                 ('_FillValue', fillv)}
ds18 = xr.Dataset({'emissions_co2_base': (['time'], emissions_co2_base)}, coords={'time':time})
ds18.emissions_co2_base.attrs = {('long_name', 'Emissions of CO2'),
                 ('standard_name', 'emissions_co2'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds19 = xr.Dataset({'temp_eq_ARCTIC': (['time'], temp_eq_ARCTIC)}, coords={'time':time})
ds19.temp_eq_ARCTIC.attrs = {('long_name', 'Arctic equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_Arctic'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds20 = xr.Dataset({'temp_eq_NHML': (['time'], temp_eq_NHML)}, coords={'time':time})
ds20.temp_eq_NHML.attrs = {('long_name', 'Northern Hemisphere ML equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_NHML'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds21 = xr.Dataset({'temp_eq_TROPICS': (['time'], temp_eq_TROPICS)}, coords={'time':time})
ds21.temp_eq_TROPICS.attrs = {('long_name', 'Tropics equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_tropics'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds22 = xr.Dataset({'temp_eq_SH': (['time'], temp_eq_SH)}, coords={'time':time})
ds22.temp_eq_SH.attrs = {('long_name', 'Southern Hemisphere equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_SH'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds23 = xr.Dataset({'temp_eq_glb': (['time'], temp_eq_glb)}, coords={'time':time})
ds23.temp_eq_glb.attrs = {('long_name', 'Global equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_global'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds = xr.merge([ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9, ds10, ds11, ds12, ds13, ds14, ds15, ds16, ds17, ds18, ds19, ds20, ds21, ds22, ds23])
ds.attrs = {('comment', 'contact: knut.vonsalzen@canada.ca'),
            ('data_licence', '1) GRANT OF LICENCE - The Government of Canada (Environment Canada) is the \
owner of all intellectual property rights (including copyright) that may exist in this Data \
product. You (as \"The Licensee\") are hereby granted a non-exclusive, non-assignable, \
non-transferable unrestricted licence to use this data product for any purpose including \
the right to share these data with others and to make value-added and derivative \
products from it. This licence is not a sale of any or all of the owner\'s rights.\
2) NO WARRANTY - This Data product is provided \"as-is\"; it has not been designed or \
prepared to meet the Licensee\'s particular requirements. Environment Canada makes no \
warranty, either express or implied, including but not limited to, warranties of \
merchantability and fitness for a particular purpose. In no event will Environment Canada \
be liable for any indirect, special, consequential or other damages attributed to the \
Licensee\'s use of the Data product.')}
outfile = './temp_co2_'+scenario+'_'+reference+'_emulator_'+realm+'.nc'
ds.to_netcdf(outfile, unlimited_dims=["time"])
