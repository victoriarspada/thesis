# python3 source code for CH4/climate emulator
# author: Knut von Salzen, email: knut.vonsalzen@canada.ca
# Edited by Victoria Spada, 2023
# RTP coefficients were modified from the emulator default settings for a 
# set increase to the Arctic RTP coefficients, and corresponding decreases 
# to the other regional coefficients.

import numpy as np
import os 
import sys
import xarray as xr
#import warnings
#warnings.simplefilter('error')

### read scenario information
# possible choices are:
#   scenario = SSP126 SSP245 SSP370 SSP585 CLE MFR MFR_SDS CFM
#   realm = mean low_ecs high_ecs low_aci high_aci low_ari_bc high_ari_bc low_ari_s high_ari_s low_alb high_alb

data_in = np.genfromtxt('emulator.csv', dtype="<U100,<U100",names=True, delimiter=',')
scenario = str(data_in['scenario'])
realm = str(data_in['realm'])

# check if SSP or AMAP scenario

cmip6_scenario = [ 'SSP119', 'SSP245', 'SSP585', 'SSP126', 'SSP370' ]
cmip6 = False
for nsc in range(len(cmip6_scenario)):
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
tunef = 1.                        # tuning factor for CH4 lifetime

# natural background emissions of CH4 (Tg/year)

emis_ch4_nat = 202.     # Prather et al. (2012)
emis_co_nat = 300.      # Olivie
emis_nox_nat = 15.      # Olivie
emis_voc_nat = 60.      # Olivie

# initial CH4 mixing ratio at the beginning of the time series (ppb)

ch4_start_ppb = 731.41              # Meinshausen (2017) estimate for 1750

# diagnostic parameters

thorizon = [20, 100]                # time horizon for ARTP and AGTP diagnostic calculation (years)
year_ref = [2015, 2015]             # reference year for ARTP and AGTP diagnostic calculation (years)

# factors for unit conversion

ppb2Tg = 2.78                                 # convert ppb to tera gram (Tg)
kt2Tg = 1.e-03                                # convert ktons to tera gram (Tg)
Tg2Gg = 1.e+03                                # convert tera gram (Tg) to giga gram (Gg)
mm_O = 16.                                    # molecular weight of oxygen
mm_N = 14.01                                  # molecular weight of nitrogen
no22n = mm_N / ( mm_N + 2. * mm_O )           # convert NO2 to N

# specified (True) or parameterized CH4 decay time (False)

constant_tau = False 
if constant_tau:
  tau_const = 9.5
  method_lifetime = 0
else:
  method_lifetime = 1                         # method for calculating CH4 lifetime
                                              # 1 - Olivie et al. (AMAP, 2021)
                                              # 2 - Arora et al. (AMAP, 2015)

# specified CH4 emissions (True) or concentrations (False)

ch4_emissions = True

# Averaging of temperature differences (True) for output

avg_temp = False
if avg_temp:
  avg_start_year = 2036
  avg_end_year = 2050

# lifetimes of CH4 due to various processes (years)

if method_lifetime == 1:
  tau_strat = 120.                            # destruction of CH4 by tropospheric OH radicals
  tau_soil = 150.                             # uptake by soils
  tau_tropcl = 200.                           # reaction with tropospheric chlorine
  tau_oh0 = 11.17                             # destruction of CH4 by tropospheric OH radicals at present-day

# Sensitivity coefficients in parameterization of OH lifetime (IPCC, 2001 and Stocker et al., 2013)

if method_lifetime == 1:
  ch40 = 1783.36 * ppb2Tg                     # reference CH4 mixing ratio for calculation of OH lifetime
  s_oh_ch4 = -0.31                            # dimensionless
  s_oh_nox = 0.0042                           # (Tg N/yr)^-1
  s_oh_co = -0.000105                         # (Tg CO/yr)^-1
  s_oh_voc = -0.000315                        # (Tg VOC/yr)^-1
  s_oh_tt = 0.0316                            # K^-1

# temperature trend (K/yr)

if method_lifetime == 1:
  dtemp0 = 0.
  dtemp0_refr = 0.

# pre-industrial CH4 mixing ratio (for CH4 forcing)

ch4_pi_ppb = 750.

# N2O mixing ratio (for CH4 forcing)

n2o_ppb = 323.
n2o_refr_ppb = 323.

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
cs_ch4_ar_old = np.array([0.31, 0.17, 0.16, 0.06])    # 60N-90N (Arctic) 
cs_ch4_ml_old = np.array([0.06, 0.24, 0.17, 0.07])    # 28N-60N (mid latitudes)
cs_ch4_tr_old = np.array([0.02, 0.10, 0.24, 0.09])    # 28N-28S (tropics)
cs_ch4_sh_old = np.array([0.  , 0.02, 0.05, 0.19])    # 28S-90S (southern hemisphere)

# weightings for the RTP coefficient rescaling
D = 0.10 # percent increase for the area weighted sum of the Arctic RTP
delta = D*area_frc_ar*cs_ch4_ar_old
d1 = delta/(area_frc_ar*cs_ch4_ar_old)
d2 = -delta/(3*area_frc_ml*cs_ch4_ml_old)
d3 = -delta/(3*area_frc_tr*cs_ch4_tr_old)
d4 = -delta/(3*area_frc_sh*cs_ch4_sh_old)

cs_ch4_ar = list(np.array([0.31, 0.17, 0.16, 0.06])*(1+d1))    # 60N-90N (Arctic) 
cs_ch4_ml = list(np.array([0.06, 0.24, 0.17, 0.07])*(1+d2))    # 28N-60N (mid latitudes)
cs_ch4_tr = list(np.array([0.02, 0.10, 0.24, 0.09])*(1+d3))    # 28N-28S (tropics)
cs_ch4_sh = list(np.array([0.  , 0.02, 0.05, 0.19])*(1+d4))    # 28S-90S (southern hemisphere)

# account for cs_ch4_sh[0] = 0 (get inf.)
d4 = -delta[0]/(2*area_frc_ml*cs_ch4_ml_old[0])
d5 = -delta[0]/(2*area_frc_tr*cs_ch4_tr_old[0])

cs_ch4_ml[0] = cs_ch4_ml_old[0]*(1+d4)    # 28N-60N (mid latitudes)
cs_ch4_tr[0] = cs_ch4_tr_old[0]*(1+d5)
cs_ch4_sh[0] = 0

# scaling factors for CH4 indirect forcings associated with O3, SO4, and stratospheric water vapour

ch4_o3_scl = 0.
ch4_so4_scl = 0.
ch4_h2o_scl = 0.

nch4 = 0
if constant_tau:
  print ('Using contstant CH4 lifetime')
  nssp = 1
  ssrcspec = [ 'CH4' ]
else:
  if method_lifetime < 0 or method_lifetime > 2:
    print ('Invalid choice, calculation of CH4 lifetime')
    sys.exit()
  print ('Using parameterized CH4 lifetime, method',str(method_lifetime))
  if method_lifetime == 1:
    nssp = 4
    ssrcspec = [ 'CH4', 'CO', 'NOx', 'VOC' ]
    nco  = 1
    nnox = 2
    nvoc = 3

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

spchk1 = [False for n in range(nssp)]
spchk2 = [False for n in range(nssp)]
spchk3 = [False for n in range(nssp)]
spchk4 = [False for n in range(nssp)]
spchk = [False for n in range(nssp)]
ntag_ch4 = [None for n in range(nssp)]
ntag_ch4_refr = [None for n in range(nssp)]
ndat_ch4 = [None for n in range(nssp)]
ndat_ch4_refr = [None for n in range(nssp)]
year_ch4_mltyr = [None for n in range(nssp)]
year_ch4_refr_mltyr = [None for n in range(nssp)]
year_start = [None for n in range(nssp)]
year_end = [None for n in range(nssp)]
nstpyr = [None for n in range(nssp)]
emis_ch4_names = [None for n in range(nssp)]
emis_ch4_refr_names = [None for n in range(nssp)]
ntag_fire_ch4 = [None for n in range(nssp)]
nstpyr_fire = [None for n in range(nssp)]
if fires:
  ntag_fire_ch4_refr = [None for n in range(nssp)]
  ndat_fire_ch4 = [None for n in range(nssp)]
  ndat_fire_ch4_refr = [None for n in range(nssp)]
  year_fire_ch4_mltyr = [None for n in range(nssp)]
  year_fire_ch4_refr_mltyr = [None for n in range(nssp)]
  year_fire_start = [None for n in range(nssp)]
  year_fire_end = [None for n in range(nssp)]
  emis_fire_ch4_names = [None for n in range(nssp)]
  emis_fire_ch4_refr_names = [None for n in range(nssp)]

# read emissions, first step

print('Scenario :'+scenario)
for n in range(nssp):
  emi_dat_file = emidir+ssrcspec[n]+'_base_'+scenario+'.csv'
  emi_refr_dat_file = emidir+ssrcspec[n]+'_'+reference+'.csv'
  if fires:
    if not cmip6:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
      emi_fire_refr_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
    else:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base_'+scenario+'.csv'
      emi_fire_refr_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
  spchk1[n] = check_input_files(emi_dat_file)
  spchk2[n] = check_input_files(emi_refr_dat_file)
  spchk[n] = spchk1[n] and spchk2[n]
  if fires:
    spchk3[n] = check_input_files(emi_fire_dat_file)
    spchk4[n] = check_input_files(emi_fire_refr_dat_file)
    spchk[n] = spchk[n] and spchk3[n] and spchk4[n]

  if spchk[n]:

#   anthropogenic emissions

    ch4_file = emi_dat_file
    ch4_refr_file = emi_refr_dat_file
    emis_ch4_in = np.genfromtxt(ch4_file, dtype=None, names=True, delimiter=',')      # open file containing time series and read
    lenx = len(emis_ch4_in.dtype.names)
    emis_ch4_names[n] = emis_ch4_in.dtype.names[1:lenx]                               # emission tag names
    ntag_ch4[n] = len(emis_ch4_names[n])                                              # number of emission tags
    ndat_ch4[n] = emis_ch4_in['Year'].size                                            # number of data points in input time series
    emis_ch4_refr_in = np.genfromtxt(ch4_refr_file, dtype=None, names=True, delimiter=',')   # open file containing time series and read
    lenx = len(emis_ch4_refr_in.dtype.names)
    emis_ch4_refr_names[n] = emis_ch4_refr_in.dtype.names[1:lenx]                     # emission tag names
    ntag_ch4_refr[n] = len(emis_ch4_refr_names[n])                                    # number of emission tags
    ndat_ch4_refr[n] = emis_ch4_refr_in['Year'].size                                  # number of data points in input time series
    emis_ch4_mltyr = [[[0. for i in range(ndat_ch4[n])] for j in range(ntag_ch4[n])] for ne in range(nssp)]
    emis_ch4_refr_mltyr = [[[0. for i in range(ndat_ch4_refr[n])] for j in range(ntag_ch4_refr[n])] for ne in range(nssp)]

#   fire emissions

    if fires:
      ch4_fire_file = emi_fire_dat_file
      ch4_fire_refr_file = emi_fire_refr_dat_file
      emis_fire_ch4_in = np.genfromtxt(ch4_fire_file, dtype=None, names=True, delimiter=',')  # open file containing time series and read
      lenx = len(emis_fire_ch4_in.dtype.names)
      emis_fire_ch4_names[n] = emis_fire_ch4_in.dtype.names[1:lenx]                     # emission tag names
      ntag_fire_ch4[n] = len(emis_fire_ch4_names[n])                                    # number of emission tags
      ndat_fire_ch4[n] = emis_fire_ch4_in['Year'].size                                  # number of data points in input time series
      emis_fire_ch4_refr_in = np.genfromtxt(ch4_fire_refr_file, dtype=None, names=True, delimiter=',')   # open file containing time series and read
      lenx = len(emis_fire_ch4_refr_in.dtype.names)
      emis_fire_ch4_refr_names[n] = emis_fire_ch4_refr_in.dtype.names[1:lenx]           # emission tag names
      ntag_fire_ch4_refr[n] = len(emis_fire_ch4_refr_names[n])                          # number of emission tags
      ndat_fire_ch4_refr[n] = emis_fire_ch4_refr_in['Year'].size                        # number of data points in input time series
      emis_fire_ch4_mltyr = [[[0. for i in range(ndat_fire_ch4[n])] for j in range(ntag_fire_ch4[n])] for ne in range(nssp)]
      emis_fire_ch4_refr_mltyr = [[[0. for i in range(ndat_fire_ch4_refr[n])] for j in range(ntag_fire_ch4_refr[n])] for ne in range(nssp)]
    else:
      ntag_fire_ch4[n] = 0

# read emissions, second step

for n in range(nssp):
  emi_dat_file = emidir+ssrcspec[n]+'_base_'+scenario+'.csv'
  emi_refr_dat_file = emidir+ssrcspec[n]+'_'+reference+'.csv'
  if fires:
    if not cmip6:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
      emi_fire_refr_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
    else:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base_'+scenario+'.csv'
      emi_fire_refr_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
  if spchk[n]:
    ch4_file = emi_dat_file
    ch4_refr_file = emi_refr_dat_file
    emis_ch4_in = np.genfromtxt(ch4_file, dtype=None, names=True, delimiter=',')
    year_ch4_mltyr[n] = emis_ch4_in['Year']                                           # input times
    year_start[n] = year_ch4_mltyr[n][0]                                              # first year of the time series
    year_end[n] = year_ch4_mltyr[n][ndat_ch4[n]-1]                                    # last year of the time series
    nstpyr[n] = year_end[n] - year_start[n] + 1                                       # total number of years
    for k in range(ntag_ch4[n]):
      emis_ch4_mltyr[n][k] = emis_ch4_in[emis_ch4_names[n][k]]                        # copy emission time series
    emis_ch4_refr_in = np.genfromtxt(ch4_refr_file, dtype=None, names=True, delimiter=',')
    year_ch4_refr_mltyr[n] = emis_ch4_refr_in['Year']                                 # input times
    if emis_ch4_refr_names[n] != emis_ch4_names[n] or ntag_ch4_refr[n] != ntag_ch4[n] \
      or ndat_ch4_refr[n] != ndat_ch4[n] or not np.array_equal(year_ch4_refr_mltyr[n],year_ch4_mltyr[n]):
      print ('Inconsistent emission input files 1')
      sys.exit()
    for k in range(ntag_ch4_refr[n]):
      emis_ch4_refr_mltyr[n][k] = emis_ch4_refr_in[emis_ch4_refr_names[n][k]]         # copy emission time series

#   fire emissions

    if fires:
      ch4_fire_file = emi_fire_dat_file
      ch4_fire_refr_file = emi_fire_refr_dat_file
      emis_fire_ch4_in = np.genfromtxt(ch4_fire_file, dtype=None, names=True, delimiter=',')
      year_fire_ch4_mltyr[n] = emis_fire_ch4_in['Year']                                 # input times
      year_fire_start[n] = year_fire_ch4_mltyr[n][0]                                    # first year of the time series
      year_fire_end[n] = year_fire_ch4_mltyr[n][ndat_fire_ch4[n]-1]                     # last year of the time series
      nstpyr_fire[n] = year_fire_end[n] - year_fire_start[n] + 1                        # total number of years
      for k in range(ntag_fire_ch4[n]):
        emis_fire_ch4_mltyr[n][k] = emis_fire_ch4_in[emis_fire_ch4_names[n][k]]         # copy emission time series
      emis_fire_ch4_refr_in = np.genfromtxt(ch4_fire_refr_file, dtype=None, names=True, delimiter=',')
      year_fire_ch4_refr_mltyr[n] = emis_fire_ch4_refr_in['Year']                       # input times
      if emis_fire_ch4_refr_names[n] != emis_fire_ch4_names[n] or ntag_fire_ch4_refr[n] != ntag_fire_ch4[n] \
        or ndat_fire_ch4_refr[n] != ndat_fire_ch4[n] or not np.array_equal(year_fire_ch4_refr_mltyr[n],year_fire_ch4_mltyr[n]):
        print ('Inconsistent fire emission input files 2')
        sys.exit()
      for k in range(ntag_fire_ch4_refr[n]):
        emis_fire_ch4_refr_mltyr[n][k] = emis_fire_ch4_refr_in[emis_fire_ch4_refr_names[n][k]] # copy emission time series
    else:
      nstpyr_fire[n] = nstpyr[n]

ina=-9
nstpyra=ina
for n in range(nssp):
  if nstpyra == ina:
    nstpyra = nstpyr[n]
  if nstpyra != nstpyr[n] or nstpyra != nstpyr_fire[n]:
    print('All emission input files must produce the same number of years')
    sys.exit()

print('Available regions and sectors in input model data file:')
for nr in range(nreg):
  regionp = sregion[nr]
  for ns in range(nsec):
    sectorp = ssector[ns]
    print(regionp+', '+sectorp)

print('Available regions and sectors in emission data file:')
for n in range(nssp):
  if spchk[n]:
    for k in range(ntag_ch4[n]):
      print(ssrcspec[n]+': '+emis_ch4_names[n][k])
    for k in range(ntag_fire_ch4[n]):
      print(ssrcspec[n]+': '+emis_fire_ch4_names[n][k])

for n in range(nssp):
  if spchk[n]:
    isec = [[None for i in range(ntag_ch4[n]+ntag_fire_ch4[n])] for ne in range(nssp)]
    ireg = [[None for i in range(ntag_ch4[n]+ntag_fire_ch4[n])] for ne in range(nssp)]
for ns in range(nsec):
  sectorp = ssector[ns]
  for nr in range(nreg):
    regionp = sregion[nr]
    emi_tag1 = '{}_{}'.format(sectorp,regionp)              # the header in the emission input file must
    emi_tag2 = '{}_{}'.format(regionp,sectorp)              # have the format "sector_region" or "region_sector" (capital letters)
    for n in range(nssp):
      if spchk[n]:
        for k in range(ntag_ch4[n]):
          if emis_ch4_names[n][k] == emi_tag1 or emis_ch4_names[n][k] == emi_tag2:
            isec[n][k] = ns                                 # matching sector index in netcdf input data
            ireg[n][k] = nr                                 # matching region index in netcdf input data
        for k in range(ntag_fire_ch4[n]):
          kx = k + ntag_ch4[n]
          if emis_fire_ch4_names[n][k] == emi_tag1 or emis_fire_ch4_names[n][k] == emi_tag2:
            isec[n][kx] = ns                                # matching sector index in netcdf input data
            ireg[n][kx] = nr                                # matching region index in netcdf input data

print('Processing the following regions and sectors:')
for n in range(nssp):
  if spchk[n]:
    ifound = 0
    ifoundf = 0
    emi_dat_file = emidir+ssrcspec[n]+'_base_'+scenario+'.csv'
    if fires:
      if not cmip6:
        emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
      else:
        emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base_'+scenario+'.csv'
    for k in range(ntag_ch4[n]):
      if isec[n][k] is not None:
        print(ssrcspec[n]+': '+sregion[ireg[n][k]]+'_'+ssector[isec[n][k]])          # print emission tags that are being processed
        ifound = ifound + 1
    for k in range(ntag_fire_ch4[n]):
      kx = k + ntag_ch4[n]
      if isec[n][kx] is not None:
        print(ssrcspec[n]+': '+sregion[ireg[n][kx]]+'_'+ssector[isec[n][kx]])        # print emission tags that are being processed
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

for n in range(nssp):
  if spchk[n]:
    emis_ch4_yr = [[[0. for i in range(nstpyr[n])] for j in range(ntag_ch4[n]+ntag_fire_ch4[n])] for ne in range(nssp)]
    emis_ch4_refr_yr = [[[0. for i in range(nstpyr[n])] for j in range(ntag_ch4[n]+ntag_fire_ch4[n])] for ne in range(nssp)]
for n in range(nssp):
  if spchk[n]:
    for k in range(ntag_ch4[n]):
      if isec[n][k] is not None:
        emis_ch4_yr[n][k] = interpol_emis(emis_ch4_mltyr[n][k],year_ch4_mltyr[n],year_start[n],ndat_ch4[n],nstpyr[n])
        emis_ch4_refr_yr[n][k] = interpol_emis(emis_ch4_refr_mltyr[n][k],year_ch4_mltyr[n],year_start[n],ndat_ch4_refr[n],nstpyr[n])
    for k in range(ntag_fire_ch4[n]):
      kx = k + ntag_ch4[n]
      if isec[n][kx] is not None:
        emis_ch4_yr[n][kx] = interpol_emis(emis_fire_ch4_mltyr[n][k],year_fire_ch4_mltyr[n],year_start[n],ndat_fire_ch4[n],nstpyr[n])
        emis_ch4_refr_yr[n][kx] = interpol_emis(emis_fire_ch4_refr_mltyr[n][k],year_fire_ch4_mltyr[n],year_start[n],ndat_fire_ch4_refr[n],nstpyr[n])

for n in range(nssp):
  if spchk[n]:
    emis_ch4_names_tmp = [["" for j in range(ntag_ch4[n]+ntag_fire_ch4[n])] for ne in range(nssp)]
for n in range(nssp):
  if spchk[n]:
    for k in range(ntag_ch4[n]):
      if isec[n][k] is not None:
        emis_ch4_names_tmp[n][k] = emis_ch4_names[n][k] 
    for k in range(ntag_fire_ch4[n]):
      kx = k + ntag_ch4[n]
      if isec[n][kx] is not None:
        emis_ch4_names_tmp[n][kx] = emis_fire_ch4_names[n][k] 
emis_ch4_names = emis_ch4_names_tmp

for n in range(nssp):
  ntag_ch4[n] = ntag_ch4[n] + ntag_fire_ch4[n]

# unit conversion

for n in range(nssp):
  if spchk[n]:
    for i in range(nstpyr[n]):
      emisum = 0.
      for k in range(ntag_ch4[n]):
        emis_ch4_yr[n][k][i] = emis_ch4_yr[n][k][i] * kt2Tg 
        emis_ch4_refr_yr[n][k][i] = emis_ch4_refr_yr[n][k][i] * kt2Tg 

if not ch4_emissions:
  print ('Option not implemented')
  sys.exit()

# implicit equilibrium climate sensitivity of underlying model

ecs0 = sum(tr_c)

# regional temperature sensitivity for different latitude bands, assuming globally uniform radiative forcing

cs_ch4_ar_tot = sum(cs_ch4_ar)
cs_ch4_ml_tot = sum(cs_ch4_ml)
cs_ch4_tr_tot = sum(cs_ch4_tr)
cs_ch4_sh_tot = sum(cs_ch4_sh)

dt = 0.01
time = np.zeros((nstpyra), dtype=np.float64)
time = xr.DataArray(time, name='time', dims=('time'), coords={'time':time}, attrs={'long_name':"Time",'standard_name':"time", 'units':"year"})
if spchk[nch4]:
  for i in range(nstpyra):
    time[i]= year_start[nch4] + i * dtyrf

ysmall = 1.e-20
def methane_lifetime_olivie(ch4_lcl,emis_nox_lcl,emis_co_lcl,emis_voc_lcl,dtemp_lcl):
#  lifetime_oh = 1. / ( (1./tau_oh0) * ( (ch4_lcl/ch40)**s_oh_ch4 * np.exp(s_oh_nox * emis_nox_lcl + s_oh_co * emis_co_lcl + s_oh_voc * emis_voc_lcl) + s_oh_tt * dtemp_lcl ) )
  lifetime_oh = (ch4_lcl/ch40)**(-s_oh_ch4) / ( (1./tau_oh0) * np.exp(s_oh_nox * emis_nox_lcl + s_oh_co * emis_co_lcl + s_oh_voc * emis_voc_lcl) + s_oh_tt * dtemp_lcl )
  if lifetime_oh > ysmall:
    lifetime = 1. / ( 1./lifetime_oh + 1./tau_strat + 1./tau_soil + 1./tau_tropcl )
  else:
    lifetime = 0.
  return lifetime

def methane_lifetime_arora(ch4_lcl):
  lifetime = 1./((0.0895*np.exp(0.32*np.log(1754./(ch4_lcl/ppb2Tg))))+1./150.+1./120.+1./200.)
  return lifetime

def methane_conc(ch4_lcl,emis_ch4_lcl,tau_lcl):
  method = 1
  if method == 1:       # analytic solution
    if tau_lcl > ysmall:
      conc = ( ch4_lcl - emis_ch4_lcl * tau_lcl ) * np.exp(-dt/tau_lcl) + emis_ch4_lcl * tau_lcl
    else:
      conc = 0.
  else:                 # semi-implicit solution
    conc = ( ch4_lcl/dt + emis_ch4_lcl ) / ( 1./dt + 1./tau_lcl )
  return conc

def methane_radiative_forcing(ch4_ppb_lcl,n2o_ppb_lcl):
  method = 1
  if method == 1:       # Etmnian et al. (2016)
    ch4_tmp = 0.5 * (ch4_ppb_lcl + ch4_pi_ppb)
    rf = ( 0.043 -1.3e-6 * ch4_tmp -8.2e-6 * n2o_ppb_lcl ) * ( np.sqrt(ch4_ppb_lcl) - np.sqrt(ch4_pi_ppb) )
  else:                 # Myhre et al. (1998)
    fun = 0.47 * np.log(1. + 2.01e-5 * (ch4_ppb_lcl * n2o_ppb_lcl)**0.75 + 5.31e-15 * ch4_ppb_lcl * (ch4_ppb_lcl * n2o_ppb_lcl)**1.52)
    fun0 = 0.47 * np.log(1. + 2.01e-5 * (ch4_pi_ppb * n2o_ppb_lcl)**0.75 + 5.31e-15 * ch4_pi_ppb * (ch4_pi_ppb * n2o_ppb_lcl)**1.52)
    rf = 0.036 * ( np.sqrt(ch4_ppb_lcl) - np.sqrt(ch4_pi_ppb) ) - (fun - fun0)
  return rf 

def interpol_lin(emis_yr_lcl,wgt_lcl,year_lcl,yearm_lcl):
  if year_lcl <= 0:
    emis = emis_yr_lcl[year_lcl]
  elif year_lcl > year_end[nch4] - year_start[nch4]:
    emis = emis_yr_lcl[yearm_lcl]
  else:
    emis = wgt_lcl * emis_yr_lcl[year_lcl] + (1.-wgt_lcl) * emis_yr_lcl[yearm_lcl]
  return emis

# arrays

nstpyrt = nstpyr[nch4]
ch4_yr = np.zeros((nstpyrt+1))
ch4_refr_yr = np.zeros((nstpyrt+1))
tau_yr = np.zeros((nstpyrt+1))
tau_refr_yr = np.zeros((nstpyrt+1))
rf_yr = [0. for i in range(nstpyrt+1)]
rf_refr_yr = [0. for i in range(nstpyrt+1)]
rf_pref_yr = [0. for i in range(nstpyrt+1)]
emis_tmp= [0. for i in range(nstpyrt+1)]
emis_yr = [0. for i in range(nstpyrt+1)]
emis_refr_yr = [0. for i in range(nstpyrt+1)]
nstp = 1 + int(round(nstpyrt / dt))
nstppyr = nstp / nstpyrt
yrwgt = 1. / float(nstppyr)
conc_pref_ppb = 1.

ch4_start = ch4_start_ppb * ppb2Tg
ch4 = ch4_start
ch4_refr = ch4_start
ch4_ppb = ch4 / ppb2Tg
ch4_refr_ppb = ch4_refr / ppb2Tg
ch4_pert = np.full((nssp,ntag_ch4[nch4]), ch4_start)
ch4_pert_ppb = np.full((nssp,ntag_ch4[nch4]), ch4_start)
ch4_pert_ppb = np.full((nssp,ntag_ch4[nch4]), ch4_ppb)
ch4_pert_yr = np.zeros((nssp,ntag_ch4[nch4],nstpyrt+1))
rf_pert_yr = np.zeros((nssp,ntag_ch4[nch4],nstpyrt+1))
rf_pert_pref_yr = np.zeros((nssp,ntag_ch4[nch4],nstpyrt+1))
emis_pert_yr = np.zeros((nssp,ntag_ch4[nch4],nstpyrt+1))
tau_pert_yr = np.zeros((nssp,ntag_ch4[nch4],nstpyrt+1))

# loop over all years to calculate annual mean CH4 concentrations and radiative forcings

print('Calculating concentration and radiative forcing perturbations...')
for n in range(nstp):
  timet = n * dt
  timeyr = year_start[nch4] + timet
  year = int(round(timet))
  yearc = int(np.floor(timet))
  yearm = year - 1
  yearp = year + 1
  wgt = timet - float(yearm) - 0.5
  emis_ch4_all = emis_ch4_nat
  emis_ch4_refr_all = 0.
  emis_co_all = emis_co_nat
  emis_co_refr_all = 0.
  emis_nox_all = emis_nox_nat
  emis_nox_refr_all = 0.
  emis_voc_all = emis_voc_nat
  emis_voc_refr_all = 0.
  if spchk[nch4]:
    for kt in range(ntag_ch4[nch4]):                                          # loop over all sectors and regions
      if isec[nch4][kt] is not None:
        emis_tmp = emis_ch4_yr[nch4][kt]
        emis_ch4 =  interpol_lin(emis_tmp,wgt,year,yearm)                     # linearly time-interpolated emissions for control scenario
        emis_ch4_all = emis_ch4_all + emis_ch4
        emis_tmp = emis_ch4_refr_yr[nch4][kt]
        emis_ch4_refr =  interpol_lin(emis_tmp,wgt,year,yearm)                # linearly time-interpolated emissions for baseine scenario
        emis_ch4_refr_all = emis_ch4_refr_all + emis_ch4_refr
        if method_lifetime == 1:
          emis_tmp = emis_ch4_yr[nco][kt]
          emis_co =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_co_all = emis_co_all + emis_co
          emis_tmp = emis_ch4_refr_yr[nco][kt]
          emis_co_refr =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_co_refr_all = emis_co_refr_all + emis_co_refr

          emis_tmp = emis_ch4_yr[nnox][kt]
          emis_nox =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_nox = emis_nox * no22n
          emis_nox_all = emis_nox_all + emis_nox
          emis_tmp = emis_ch4_refr_yr[nnox][kt]
          emis_nox_refr =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_nox_refr = emis_nox_refr * no22n
          emis_nox_refr_all = emis_nox_refr_all + emis_nox_refr

          emis_tmp = emis_ch4_yr[nvoc][kt]
          emis_voc =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_voc_all = emis_voc_all + emis_voc
          emis_tmp = emis_ch4_refr_yr[nvoc][kt]
          emis_voc_refr =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_voc_refr_all = emis_voc_refr_all + emis_voc_refr

  if method_lifetime == 1:
    if timeyr >= 2015:
      dtemp = dtemp0
      dtemp_refr = dtemp0_refr
    else:
      dtemp = 0.
      dtemp_refr = 0.
  if constant_tau:
    tau = tau_const
  else:
    if method_lifetime == 1:
      tau = methane_lifetime_olivie(ch4,emis_nox_all,emis_co_all,emis_voc_all,dtemp)
    elif method_lifetime == 2:
      tau = methane_lifetime_arora(ch4)
    tau = tau * tunef
  tau_yr[yearc] = tau_yr[yearc] + tau * yrwgt
  emis_yr[yearc] = emis_yr[yearc] + emis_ch4_all * yrwgt
  emis_refr_yr[yearc] = emis_refr_yr[yearc] + emis_ch4_refr_all * yrwgt
  if constant_tau:
    tau_refr = tau_const
  else:
    if method_lifetime == 1:
      tau_refr = methane_lifetime_olivie(ch4_refr,emis_nox_refr_all,emis_co_refr_all,emis_voc_refr_all,dtemp_refr)
    elif method_lifetime == 2:
      tau_refr = methane_lifetime_arora(ch4_refr)
    tau_refr = tau_refr * tunef
  tau_refr_yr[yearc] = tau_refr_yr[yearc] + tau_refr * yrwgt
  ch4 = methane_conc(ch4,emis_ch4_all,tau)
  ch4_refr = methane_conc(ch4_refr,emis_ch4_refr_all,tau_refr)
  ch4_ppb = ch4 / ppb2Tg
  ch4_refr_ppb = ch4_refr / ppb2Tg
  ch4_yr[yearc] = ch4_yr[yearc] + ch4_ppb * yrwgt
  ch4_refr_yr[yearc] = ch4_refr_yr[yearc] + ch4_refr_ppb * yrwgt
  ch4_rf = methane_radiative_forcing(ch4_ppb,n2o_ppb)
  ch4_rf_refr = methane_radiative_forcing(ch4_refr_ppb,n2o_refr_ppb)
  ch4_pref_ppb = ch4_ppb + conc_pref_ppb
  ch4_rf_pref = methane_radiative_forcing(ch4_pref_ppb,n2o_refr_ppb)
  rf_yr[yearc] = rf_yr[yearc] + ch4_rf * yrwgt
  rf_refr_yr[yearc] = rf_refr_yr[yearc] + ch4_rf_refr * yrwgt
  rf_pref_yr[yearc] = rf_pref_yr[yearc] + ch4_rf_pref * yrwgt

  if spchk[nch4]:
    for kt in range(ntag_ch4[nch4]):                                                            # loop over all sectors and regions
      if isec[nch4][kt] is not None:
        for i in range(nstpyr[nch4]):
          emis_tmp[i] = emis_ch4_yr[nch4][kt][i] - emis_ch4_refr_yr[nch4][kt][i]
        emis_ch4 =  interpol_lin(emis_tmp,wgt,year,yearm)                          # linearly time-interpolated emissions for control scenario
        emis_pert_yr[nch4,kt,yearc] = emis_pert_yr[nch4,kt,yearc] + emis_ch4 * yrwgt
        emis_ch4_pert = emis_ch4_all - emis_ch4
        if method_lifetime == 1:
          for i in range(nstpyr[nch4]):
            emis_tmp[i] = emis_ch4_yr[nco][kt][i] - emis_ch4_refr_yr[nco][kt][i]
          emis_co =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_pert_yr[nco,kt,yearc] = emis_pert_yr[nco,kt,yearc] + emis_co * yrwgt
          emis_co_pert = emis_co_all - emis_co

          for i in range(nstpyr[nch4]):
            emis_tmp[i] = emis_ch4_yr[nnox][kt][i] - emis_ch4_refr_yr[nnox][kt][i]
          emis_nox =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_nox = emis_nox * no22n
          emis_pert_yr[nnox,kt,yearc] = emis_pert_yr[nnox,kt,yearc] + emis_nox * yrwgt
          emis_nox_pert = emis_nox_all - emis_nox

          for i in range(nstpyr[nch4]):
            emis_tmp[i] = emis_ch4_yr[nvoc][kt][i] - emis_ch4_refr_yr[nvoc][kt][i]
          emis_voc =  interpol_lin(emis_tmp,wgt,year,yearm)
          emis_pert_yr[nvoc,kt,yearc] = emis_pert_yr[nvoc,kt,yearc] + emis_voc * yrwgt
          emis_voc_pert = emis_voc_all - emis_voc
        for ne in range(nssp):
          if ne == nch4:
            emis_ch4_sel = emis_ch4_pert
          else:
            emis_ch4_sel = emis_ch4_all
          if constant_tau:
            tau_pert = tau_const
          else:
            if method_lifetime == 1:
              if ne == nco:
                emis_co_sel  = emis_co_pert
              else:
                emis_co_sel  = emis_co_all
              if ne == nnox:
                emis_nox_sel = emis_nox_pert
              else:
                emis_nox_sel = emis_nox_all
              if ne == nvoc:
                emis_voc_sel = emis_voc_pert
              else:
                emis_voc_sel = emis_voc_all
              tau_pert = methane_lifetime_olivie(ch4_pert[ne,kt],emis_nox_sel,emis_co_sel,emis_voc_sel,dtemp)
            elif method_lifetime == 2:
              tau_pert = methane_lifetime_arora(ch4_pert[ne,kt])
            tau_pert = tau_pert * tunef
          tau_pert_yr[ne,kt,yearc] = tau_pert_yr[ne,kt,yearc] + tau_pert * yrwgt
          ch4_pert[ne,kt] = methane_conc(ch4_pert[ne,kt],emis_ch4_sel,tau_pert)
          ch4_pert_ppb[ne,kt] = ch4_pert[ne,kt] / ppb2Tg
          ch4_pert_yr[ne,kt,yearc] = ch4_pert_yr[ne,kt,yearc] + ch4_pert_ppb[ne,kt] * yrwgt
          ch4_rf = methane_radiative_forcing(ch4_pert_ppb[ne,kt],n2o_ppb)
          rf_pert_yr[ne,kt,yearc] = rf_pert_yr[ne,kt,yearc] + ch4_rf * yrwgt
          ch4_pref_ppb = ch4_pert_ppb[ne,kt] + conc_pref_ppb
          ch4_rf_pref = methane_radiative_forcing(ch4_pref_ppb,n2o_refr_ppb)
          rf_pert_pref_yr[ne,kt,yearc] = rf_pert_pref_yr[ne,kt,yearc] + ch4_rf_pref * yrwgt

nhor = len(year_ref)
agtp_hor = np.zeros((nhor))
lats = 4
artp_hor = np.zeros((lats,nhor))
scl = np.full((nstpyrt), 1.)
temp_ch4_ar = np.zeros((nstpyrt))
temp_ch4_ml = np.zeros((nstpyrt))
temp_ch4_tr = np.zeros((nstpyrt))
temp_ch4_sh = np.zeros((nstpyrt))
demis_pert_i = np.zeros((nssp,ntag_ch4[nch4]))
temp_ch4_pert_ar = np.zeros((nssp,ntag_ch4[nch4],nstpyrt))
temp_ch4_pert_ml = np.zeros((nssp,ntag_ch4[nch4],nstpyrt))
temp_ch4_pert_tr = np.zeros((nssp,ntag_ch4[nch4],nstpyrt))
temp_ch4_pert_sh = np.zeros((nssp,ntag_ch4[nch4],nstpyrt))

# ARTP model time-integration to generate temperature perturbations

print('Calculating temperature perturbations...')
for i in range(nstpyrt):                                                                 # loop over years of pulse emission
  timeyr_i = year_start[nch4] + i * dtyrf
  demis_i = (emis_yr[i] - emis_refr_yr[i]) * dtyrf                                       # emitted amount of CH4 in one year
  drf = rf_yr[i] - rf_pref_yr[i]
  dconc = -conc_pref_ppb * ppb2Tg
  frceff = (drf/dconc) * (1. + ch4_o3_scl + ch4_so4_scl + ch4_h2o_scl)                   # forcing efficiency of methane, accounting for indirect effects
  drf_ar = drf                                                                           # forcing in the Arctic
  drf_ml = drf                                                                           # forcing at mid latitudes
  drf_tr = drf                                                                           # forcing in the tropics
  drf_sh = drf                                                                           # forcing in the southern hemisphere
  csh_ch4_ar = cs_ch4_ar[0] * (drf_ar/drf) + cs_ch4_ar[1] * (drf_ml/drf) \
             + cs_ch4_ar[2] * (drf_tr/drf) + cs_ch4_ar[3] * (drf_sh/drf)
  csh_ch4_ml = cs_ch4_ml[0] * (drf_ar/drf) + cs_ch4_ml[1] * (drf_ml/drf) \
             + cs_ch4_ml[2] * (drf_tr/drf) + cs_ch4_ml[3] * (drf_sh/drf)
  csh_ch4_tr = cs_ch4_tr[0] * (drf_ar/drf) + cs_ch4_tr[1] * (drf_ml/drf) \
             + cs_ch4_tr[2] * (drf_tr/drf) + cs_ch4_tr[3] * (drf_sh/drf)
  csh_ch4_sh = cs_ch4_sh[0] * (drf_ar/drf) + cs_ch4_sh[1] * (drf_ml/drf) \
             + cs_ch4_sh[2] * (drf_tr/drf) + cs_ch4_sh[3] * (drf_sh/drf)
  agtp_ch4_est = csh_ch4_ar * area_frc_ar + csh_ch4_ml * area_frc_ml \
               + csh_ch4_tr * area_frc_tr + csh_ch4_sh * area_frc_sh 
  if scaling:
    scl[i] = ecs0 / agtp_ch4_est                                                           # scale regional ARTPs in order to match global ARTP

  if spchk[nch4]:
    for kt in range(ntag_ch4[nch4]):                                                            # loop over all sectors and regions
      if isec[nch4][kt] is not None:
        for ne in range(nssp):
          demis_pert_i[ne,kt] = demis_i                                                      # emitted amount of CH4 in one year
          if ne == nch4:
            demis_pert_i[ne,kt] = demis_pert_i[ne,kt] - emis_pert_yr[nch4,kt,i] * dtyrf

  for j in range(i+1,nstpyrt):                                                           # loop over years of temperature response to pulse emission
    timeyr_j = year_start[nch4] + j * dtyrf
    dt_emis = (j - i) * dtyrf
    tr = 0.
    for n in range(2):
      tr = tr + (tau_yr[i] * tr_c[n])/(tau_yr[i] - tr_d[n]) * ( np.exp(-dt_emis/tau_yr[i]) - np.exp(-dt_emis/tr_d[n]) )
    tr = tr * (ecs/ecs0)                                                                 # scale in order to match specified equilibrium climate sensitivity 
    agtp_ch4 = frceff * tr                                                               # absolute global temperature potential of CH4
    artp_ch4_ar = (agtp_ch4/ecs0) * csh_ch4_ar * scl[i]                                  # arctic absolute regional temperature potential of CH4
    artp_ch4_ml = (agtp_ch4/ecs0) * csh_ch4_ml * scl[i]                                  # mid-latitude absolute regional temperature potential of CH4
    artp_ch4_tr = (agtp_ch4/ecs0) * csh_ch4_tr * scl[i]                                  # tropics absolute regional temperature potential of CH4
    artp_ch4_sh = (agtp_ch4/ecs0) * csh_ch4_sh * scl[i]                                  # southern hemisphere absolute regional temperature potential of CH4
    temp_ch4_ar[j] = temp_ch4_ar[j] + artp_ch4_ar * demis_i                              # arctic temperature response to pulse emission
    temp_ch4_ml[j] = temp_ch4_ml[j] + artp_ch4_ml * demis_i                              # mid-latitude temperature response to pulse emission
    temp_ch4_tr[j] = temp_ch4_tr[j] + artp_ch4_tr * demis_i                              # tropics temperature response to pulse emission
    temp_ch4_sh[j] = temp_ch4_sh[j] + artp_ch4_sh * demis_i                              # southern hemisphere temperature response to pulse emission

    if spchk[nch4]:
      for kt in range(ntag_ch4[nch4]):                                            # loop over all sectors and regions
        if isec[nch4][kt] is not None:
          for ne in range(nssp):
            tr = 0.
            for n in range(2):
              tr = tr + (tau_pert_yr[ne,kt,i] * tr_c[n])/(tau_pert_yr[ne,kt,i] - tr_d[n]) \
                      * ( np.exp(-dt_emis/tau_pert_yr[ne,kt,i]) - np.exp(-dt_emis/tr_d[n]) )
            tr = tr * (ecs/ecs0)                                                  # scale in order to match specified equilibrium climate sensitivity 
            agtp_ch4 = frceff * tr                                                # absolute global temperature potential of CH4
            artp_ch4_ar = (agtp_ch4/ecs0) * csh_ch4_ar * scl[i]                   # arctic absolute regional temperature potential of CH4
            artp_ch4_ml = (agtp_ch4/ecs0) * csh_ch4_ml * scl[i]                   # mid-latitude absolute regional temperature potential of CH4
            artp_ch4_tr = (agtp_ch4/ecs0) * csh_ch4_tr * scl[i]                   # tropics absolute regional temperature potential of CH4
            artp_ch4_sh = (agtp_ch4/ecs0) * csh_ch4_sh * scl[i]                   # southern hemisphere absolute regional temperature potential of CH4

            temp_ch4_pert_ar[ne,kt,j] = temp_ch4_pert_ar[ne,kt,j] + artp_ch4_ar * demis_pert_i[ne,kt] # arctic temperature response to pulse emission
            temp_ch4_pert_ml[ne,kt,j] = temp_ch4_pert_ml[ne,kt,j] + artp_ch4_ml * demis_pert_i[ne,kt] # mid-latitude temperature response to pulse emission
            temp_ch4_pert_tr[ne,kt,j] = temp_ch4_pert_tr[ne,kt,j] + artp_ch4_tr * demis_pert_i[ne,kt] # tropics temperature response to pulse emission
            temp_ch4_pert_sh[ne,kt,j] = temp_ch4_pert_sh[ne,kt,j] + artp_ch4_sh * demis_pert_i[ne,kt] # southern hemisphere temperature response to pulse emission

  for ny in range(nhor):                                                                 # diagnose results for different specified time horizons
    if year_ref[ny] == timeyr_i:
      dt_emis = float(thorizon[ny])
      tr = 0.
      for n in range(2):
        tr = tr + (tau_yr[i] * tr_c[n])/(tau_yr[i] - tr_d[n]) * ( np.exp(-dt_emis/tau_yr[i]) - np.exp(-dt_emis/tr_d[n]) )
      tr = tr * (ecs/ecs0)
      agtp_hor[ny] = frceff * tr
      artp_hor[0,ny] = (agtp_hor[ny]/ecs0) * csh_ch4_ar * scl[i]
      artp_hor[1,ny] = (agtp_hor[ny]/ecs0) * csh_ch4_ml * scl[i]
      artp_hor[2,ny] = (agtp_hor[ny]/ecs0) * csh_ch4_tr * scl[i]
      artp_hor[3,ny] = (agtp_hor[ny]/ecs0) * csh_ch4_sh * scl[i]

# print AGTP and ARTPs
      
slats = ['ARCTIC', 'NHMDLAT', 'TROPICS', 'SH']
print('AGTP and ARTP at different time horizons (mK/Tg)')
for ny in range(nhor):
  print(year_ref[ny],'AGTP',thorizon[ny],agtp_hor[ny]*1000.)
  print('ARTP')
  for nl in range(lats):
    print(year_ref[ny],thorizon[ny],slats[nl],artp_hor[nl,ny]*1000.)

# equilibrium temperatures
    
temp_ch4_eq_ar = [0. for i in range(nstpyrt)]
temp_ch4_eq_ml = [0. for i in range(nstpyrt)]
temp_ch4_eq_tr = [0. for i in range(nstpyrt)]
temp_ch4_eq_sh = [0. for i in range(nstpyrt)]
temp_ch4_eq_glb = [0. for i in range(nstpyrt)]
temp_ch4_glb = [0. for i in range(nstpyrt)]
if avg_temp:
  temp_avg_ar = 0.
  temp_avg_glb = 0.
  twgt = 1./float(avg_end_year - avg_start_year + 1)
for i in range(nstpyrt):
  drf = rf_yr[i] - rf_refr_yr[i]
  drf_ar = drf                                                                           # forcing in the Arctic
  drf_ml = drf                                                                           # forcing at mid latitudes
  drf_tr = drf                                                                           # forcing in the tropics
  drf_sh = drf                                                                           # forcing in the southern hemisphere
  csh_ch4_ar = cs_ch4_ar[0] * (drf_ar) + cs_ch4_ar[1] * (drf_ml) \
             + cs_ch4_ar[2] * (drf_tr) + cs_ch4_ar[3] * (drf_sh)
  csh_ch4_ml = cs_ch4_ml[0] * (drf_ar) + cs_ch4_ml[1] * (drf_ml) \
             + cs_ch4_ml[2] * (drf_tr) + cs_ch4_ml[3] * (drf_sh)
  csh_ch4_tr = cs_ch4_tr[0] * (drf_ar) + cs_ch4_tr[1] * (drf_ml) \
             + cs_ch4_tr[2] * (drf_tr) + cs_ch4_tr[3] * (drf_sh)
  csh_ch4_sh = cs_ch4_sh[0] * (drf_ar) + cs_ch4_sh[1] * (drf_ml) \
             + cs_ch4_sh[2] * (drf_tr) + cs_ch4_sh[3] * (drf_sh)
  tr = ecs/ecs0 * (1. + ch4_o3_scl + ch4_so4_scl + ch4_h2o_scl)
  temp_ch4_eq_ar[i] = csh_ch4_ar * scl[i] * tr
  temp_ch4_eq_ml[i] = csh_ch4_ml * scl[i] * tr
  temp_ch4_eq_tr[i] = csh_ch4_tr * scl[i] * tr
  temp_ch4_eq_sh[i] = csh_ch4_sh * scl[i] * tr
  temp_ch4_eq_glb[i] = temp_ch4_eq_ar[i] * area_frc_ar + temp_ch4_eq_ml[i] * area_frc_ml \
                     + temp_ch4_eq_tr[i] * area_frc_tr + temp_ch4_eq_sh[i] * area_frc_sh 
  temp_ch4_glb[i] = temp_ch4_ar[i] * area_frc_ar + temp_ch4_ml[i] * area_frc_ml \
                  + temp_ch4_tr[i] * area_frc_tr + temp_ch4_sh[i] * area_frc_sh 
  year = year_start[nch4] + i * dtyr
  if avg_temp:
    if year >= avg_start_year and year <= avg_end_year:
      temp_avg_ar = temp_avg_ar + temp_ch4_ar[i] * twgt
      temp_avg_glb = temp_avg_glb + temp_ch4_glb[i] * twgt

nfsp = 1
srctofrc = [[None for n in range(nfsp)] for ne in range(nssp)]
for ne in range(nssp):                                                                       # loop over all emitted species
  if spchk[ne]:
    srctofrc[ne][nch4] = nch4 

temp_ARCTIC = np.full((nstpyra,nreg,nsec,nfsp,nssp), fillv, dtype=np.float64)
temp_NHML = np.full((nstpyra,nreg,nsec,nfsp,nssp), fillv, dtype=np.float64)
temp_TROPICS = np.full((nstpyra,nreg,nsec,nfsp,nssp), fillv, dtype=np.float64)
temp_SH = np.full((nstpyra,nreg,nsec,nfsp,nssp), fillv, dtype=np.float64)
temp_glb = np.full((nstpyra,nreg,nsec,nfsp,nssp), fillv, dtype=np.float64)
emissions = np.full((nstpyra,nreg,nsec,nssp), fillv, dtype=np.float64)
emissions_nat = np.full((nstpyra,nssp), fillv, dtype=np.float64)
lifetime = np.full((nstpyra,nreg,nsec,nssp), fillv, dtype=np.float64)
conc_ch4 = np.full((nstpyra,nreg,nsec,nssp), fillv, dtype=np.float64)
forc_ch4 = np.full((nstpyra,nreg,nsec,nssp), fillv, dtype=np.float64)
conc_ch4_nat = np.full((nstpyra), fillv, dtype=np.float64)
lifetime_ch4_base = np.full((nstpyra), fillv, dtype=np.float64)
conc_ch4_base = np.full((nstpyra), fillv, dtype=np.float64)
forc_ch4_base = np.full((nstpyra), fillv, dtype=np.float64)
emissions_ch4_base = np.full((nstpyra), fillv, dtype=np.float64)
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

# final temperatures from temperature perturbations

print('Calculating final temperatures...')
for ne in range(nssp):                                                                       # loop over all emitted species
  if spchk[ne]:
    try:
      nfspt = len(srctofrc[ne][:])
    except:
      nfspt = 0
    if nfspt > 0:
      for nft in range(nfspt):                                                               # loop over all forcing species
        nf = srctofrc[ne][nft]
        for kt in range(ntag_ch4[ne]):                                                       # loop over all sectors and regions
          if isec[ne][kt] is not None:
            ns = isec[ne][kt]
            nr = ireg[ne][kt]
            for i in range(nstpyra):
              temp_ARCTIC[i,nr,ns,nf,ne]  = temp_ch4_ar[i] - temp_ch4_pert_ar[ne,kt,i]       # temperature changes due to dir. radiative forcing
              temp_NHML[i,nr,ns,nf,ne]    = temp_ch4_ml[i] - temp_ch4_pert_ml[ne,kt,i]
              temp_TROPICS[i,nr,ns,nf,ne] = temp_ch4_tr[i] - temp_ch4_pert_tr[ne,kt,i]
              temp_SH[i,nr,ns,nf,ne]      = temp_ch4_sh[i] - temp_ch4_pert_sh[ne,kt,i]
              temp_glb[i,nr,ns,nf,ne]     = temp_ARCTIC[i,nr,ns,nf,ne] * area_frc_ar \
                                          + temp_NHML[i,nr,ns,nf,ne] * area_frc_ml \
                                          + temp_TROPICS[i,nr,ns,nf,ne] * area_frc_tr \
                                          + temp_SH[i,nr,ns,nf,ne] * area_frc_sh
              lifetime[i,nr,ns,ne]        = tau_yr[i] - tau_pert_yr[ne,kt,i]
              conc_ch4[i,nr,ns,ne]        = ch4_yr[i] - ch4_pert_yr[ne,kt,i]
              forc_ch4[i,nr,ns,ne]        = rf_yr[i] - rf_pert_yr[ne,kt,i]
              emissions[i,nr,ns,ne]       = emis_pert_yr[ne,kt,i] * Tg2Gg
              lifetime_ch4_base[i]  = tau_yr[i]
              conc_ch4_base[i]  = ch4_yr[i]
              forc_ch4_base[i]  = rf_yr[i]
              emissions_ch4_base[i] = emis_yr[i] * Tg2Gg
              if ne == nch4:
                emissions_nat[i,ne] = emis_ch4_nat
                conc_ch4_nat[i] = ch4_start_ppb
              if not constant_tau:
                if ne == nnox:
                  emissions[i,nr,ns,ne] = emissions[i,nr,ns,ne]/no22n
                  emissions_nat[i,ne] = emis_nox_nat
                if ne == nvoc:
                  emissions_nat[i,ne] = emis_voc_nat
                if ne == nco:
                  emissions_nat[i,ne] = emis_co_nat
              artp_ARCTIC[:] = 1000. * artp_hor[0,:]
              artp_NHML[:] = 1000. * artp_hor[1,:]
              artp_TROPICS[:] = 1000. * artp_hor[2,:]
              artp_SH[:] = 1000. * artp_hor[3,:]
              agtp[:] = 1000. * agtp_hor[:]

for i in range(nstpyra):
  temp_eq_ARCTIC[i]  = temp_ch4_eq_ar[i]       # equilibrium temperature changes due to dir. radiative forcing
  temp_eq_NHML[i]    = temp_ch4_eq_ml[i]
  temp_eq_TROPICS[i] = temp_ch4_eq_tr[i]
  temp_eq_SH[i]      = temp_ch4_eq_sh[i]
  temp_eq_glb[i]     = temp_ch4_eq_glb[i]

# write output files

region  = np.arange((nreg), dtype=np.float64)
sector  = np.arange((nsec), dtype=np.float64)
srcspec = np.arange((nssp), dtype=np.float64)
frcspec = np.arange((nfsp), dtype=np.float64)
region[:] = region[:] + 1
sector[:] = sector[:] + 1
srcspec[:] = srcspec[:] + 1
frcspec[:] = frcspec[:] + 1
region = xr.DataArray(region,name='region', dims=('region'), coords={'region':region}, \
                      attrs={'long_name':"Source regions: 1 - ASIA, 2 - ROEUR, 3 - AWEST, 4 - AEAST, 5 - ROW, 6 - ARCTIC, 7 - NONARC", 'standard_name':"region"})
sector = xr.DataArray(sector,name='sector', dims=('sector'), coords={'sector':sector}, \
                      attrs={'long_name':"Source sectors: 1 - SURF, 2 - FLAR, 3 - FIRE, 4 - SHIP", 'standard_name':"sector"})
if nssp == 4:
  srcspec = xr.DataArray(srcspec,name='srcspec', dims=('srcspec'), coords={'srcspec':srcspec}, \
                        attrs={'long_name':"Source species: 1 - CH4, 2 - CO, 3 - NOX, 4 - VOC", 'standard_name':"srcspec"})
else:
  srcspec = xr.DataArray(srcspec,name='srcspec', dims=('srcspec'), coords={'srcspec':srcspec}, \
                        attrs={'long_name':"Source species: 1 - CH4", 'standard_name':"srcspec"})
frcspec = xr.DataArray(frcspec,name='frcspec', dims=('frcspec'), coords={'frcspec':frcspec}, \
                      attrs={'long_name':"Forcing species: 1 - CH4", 'standard_name':"frcspec"})
thorizon = xr.DataArray(thorizon, name='thorizon', dims=('thorizon'), coords={'thorizon':thorizon}, attrs={'long_name':"Time horizon",'standard_name':"time_horizon", 'units':"year"})
ds1 = xr.Dataset({'temp_ARCTIC': (['time','region','sector','frcspec','srcspec'], temp_ARCTIC)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec})
ds1.temp_ARCTIC.attrs = {('long_name', 'Arctic temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds2 = xr.Dataset({'temp_NHML': (['time','region','sector','frcspec','srcspec'], temp_NHML)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec})
ds2.temp_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds3 = xr.Dataset({'temp_TROPICS': (['time','region','sector','frcspec','srcspec'], temp_TROPICS)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec})
ds3.temp_TROPICS.attrs = {('long_name', 'Tropical temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds4 = xr.Dataset({'temp_SH': (['time','region','sector','frcspec','srcspec'], temp_SH)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec})
ds4.temp_SH.attrs = {('long_name', 'Southern Hemisphere temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds5 = xr.Dataset({'temp_glb': (['time','region','sector','frcspec','srcspec'], temp_glb)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec})
ds5.temp_glb.attrs = {('long_name', 'Global temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'global_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds6 = xr.Dataset({'emissions': (['time','region','sector','srcspec'], emissions)}, coords={'time':time,'region':region,'sector':sector,'srcspec':srcspec})
ds6.emissions.attrs = {('long_name', 'Regional SLCF emissions'),
                 ('standard_name', 'regional_SLCF_emissions'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds7 = xr.Dataset({'emissions_nat': (['time','srcspec'], emissions_nat)}, coords={'time':time,'srcspec':srcspec})
ds7.emissions_nat.attrs = {('long_name', 'Natural SLCF emissions'),
                 ('standard_name', 'natural_SLCF_emissions'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds8 = xr.Dataset({'lifetime': (['time','region','sector','srcspec'], lifetime)}, coords={'time':time,'region':region,'sector':sector,'srcspec':srcspec})
ds8.lifetime.attrs = {('long_name', 'Change in lifetime of CH4 oxidation'),
                 ('standard_name', 'lifetime'),
                 ('units', 'year'),
                 ('_FillValue', fillv)}
ds9 = xr.Dataset({'conc_ch4': (['time','region','sector','srcspec'], conc_ch4)}, coords={'time':time,'region':region,'sector':sector,'srcspec':srcspec})
ds9.conc_ch4.attrs = {('long_name', 'Change in concentration of CH4'),
                 ('standard_name', 'conc_ch4'),
                 ('units', 'ppb'),
                 ('_FillValue', fillv)}
ds10 = xr.Dataset({'forc_ch4': (['time','region','sector','srcspec'], forc_ch4)}, coords={'time':time,'region':region,'sector':sector,'srcspec':srcspec})
ds10.forc_ch4.attrs = {('long_name', 'Change in direct radiative forcing of CH4'),
                 ('standard_name', 'forc_ch4'),
                 ('units', 'W/m2'),
                 ('_FillValue', fillv)}
ds11 = xr.Dataset({'artp_ARCTIC': (['thorizon'], artp_ARCTIC)}, coords={'thorizon':thorizon})
ds11.artp_ARCTIC.attrs = {('long_name', 'Arctic absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_ARCTIC'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds12 = xr.Dataset({'artp_NHML': (['thorizon'], artp_NHML)}, coords={'thorizon':thorizon})
ds12.artp_NHML.attrs = {('long_name', 'NH midlatitude absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_NHML'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds13 = xr.Dataset({'artp_TROPICS': (['thorizon'], artp_TROPICS)}, coords={'thorizon':thorizon})
ds13.artp_TROPICS.attrs = {('long_name', 'Tropics absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_TROPICS'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds14 = xr.Dataset({'artp_SH': (['thorizon'], artp_SH)}, coords={'thorizon':thorizon})
ds14.artp_SH.attrs = {('long_name', 'Southern Hemisphere absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_SH'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds15 = xr.Dataset({'agtp': (['thorizon'], agtp)}, coords={'thorizon':thorizon})
ds15.agtp.attrs = {('long_name', 'Absolute global temperarture change potential'),
                 ('standard_name', 'AGTP'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds16 = xr.Dataset({'conc_ch4_nat': (['time'], conc_ch4_nat)}, coords={'time':time})
ds16.conc_ch4_nat.attrs = {('long_name', 'Natural CH4 concentration'),
                 ('standard_name', 'natural_CH4_conc'),
                 ('units', 'ppb'),
                 ('_FillValue', fillv)}
ds17 = xr.Dataset({'lifetime_ch4_base': (['time'], lifetime_ch4_base)}, coords={'time':time})
ds17.lifetime_ch4_base.attrs = {('long_name', 'Lifetime of CH4'),
                 ('standard_name', 'lifetime_ch4'),
                 ('units', 'year'),
                 ('_FillValue', fillv)}
ds18 = xr.Dataset({'conc_ch4_base': (['time'], conc_ch4_base)}, coords={'time':time})
ds18.conc_ch4_base.attrs = {('long_name', 'Concentration of CH4'),
                 ('standard_name', 'conc_ch4'),
                 ('units', 'ppb'),
                 ('_FillValue', fillv)}
ds19 = xr.Dataset({'forc_ch4_base': (['time'], forc_ch4_base)}, coords={'time':time})
ds19.forc_ch4_base.attrs = {('long_name', 'Direct radiative forcing of CH4'),
                 ('standard_name', 'forcing_ch4'),
                 ('units', 'W/m2'),
                 ('_FillValue', fillv)}
ds20 = xr.Dataset({'emissions_ch4_base': (['time'], emissions_ch4_base)}, coords={'time':time})
ds20.emissions_ch4_base.attrs = {('long_name', 'Emissions of CH4'),
                 ('standard_name', 'emissions_ch4'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds21 = xr.Dataset({'temp_eq_ARCTIC': (['time'], temp_eq_ARCTIC)}, coords={'time':time})
ds21.temp_eq_ARCTIC.attrs = {('long_name', 'Arctic equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_Arctic'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds22 = xr.Dataset({'temp_eq_NHML': (['time'], temp_eq_NHML)}, coords={'time':time})
ds22.temp_eq_NHML.attrs = {('long_name', 'Northern Hemisphere ML equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_NHML'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds23 = xr.Dataset({'temp_eq_TROPICS': (['time'], temp_eq_TROPICS)}, coords={'time':time})
ds23.temp_eq_TROPICS.attrs = {('long_name', 'Tropics equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_tropics'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds24 = xr.Dataset({'temp_eq_SH': (['time'], temp_eq_SH)}, coords={'time':time})
ds24.temp_eq_SH.attrs = {('long_name', 'Southern Hemisphere equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_SH'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds25 = xr.Dataset({'temp_eq_glb': (['time'], temp_eq_glb)}, coords={'time':time})
ds25.temp_eq_glb.attrs = {('long_name', 'Global equilibrium temperature change'),
                 ('standard_name', 'Eq_temp_global'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds = xr.merge([ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9, ds10, ds11, ds12, ds13, ds14, ds15, ds16, ds17, ds18, ds19, ds20, ds21, ds22, ds23, ds24, ds25])
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
outfile = './temp_ch4_'+scenario+'_'+reference+'_emulator_'+realm+'.nc'
ds.to_netcdf(outfile, unlimited_dims=["time"])
