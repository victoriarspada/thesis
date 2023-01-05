# python3 source code for aerosol/climate emulator
# author: Knut von Salzen, email: knut.vonsalzen@canada.ca

import numpy as np
import os 
import sys
import xarray as xr
import re

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

### emulator type
    
case = 5            # 1 - CanAM, 2 - MRI-ESM, 3 - UKESM, 4 - CESM, 5 - multi-model emulator, 6 - AMAP (2015)
model = ['CanAM', 'MRI-ESM', 'UKESM', 'CESM', 'emulator', 'AMAP']
iamap = 6

### input files
# emulator input 1: file for radiative forcing sensitivity from climate models

datdir = './netcdf/'
if realm == 'low_ecs' or realm == 'high_ecs':
  realmx = 'mean'
else:
  realmx = realm
if case == 1:
  model_dat_file = datdir+'CanAM5PAM_data_'+realmx+'.nc'
elif case == 2:
  model_dat_file = datdir+'MRI-ESM_v3_data_'+realmx+'.nc'
elif case == 3:
  model_dat_file = datdir+'UKESM1_data_'+realmx+'.nc'
elif case == 4:
  model_dat_file = datdir+'CESM_data_'+realmx+'.nc'
elif case == 5:
  model_dat_file = datdir+'emulator_data_'+realmx+'.nc'
elif case == 6:
  model_dat_file = datdir+'amap_2015_refmt_models_data.nc'

# emulator input 2: files for emissions, name of the directory that contains the files

datdir = './emissions/'
if not cmip6:
  emidir = datdir+'AMAP/'
else:
  emidir = datdir+'CMIP6/'
    
### emulator configuration
# radiative forcing options

scaling_rad = [False, True, False]   # scale RTPs to match equilibrium climate sensitivity
scaling_cld = [False, False, False]  # scale RTPs to match equilibrium climate sensitivity

frcchk = False                       # this selects direct radiative forcings. Otherwise effective radiative forcing is used
achk   = True                        # this selects albedo forcings
cchk   = True                        # this selects cloud forcings
fires  = True                        # simulate fire emission impacts

reference = '1750'  # reference year

# diagnostic parameters

thorizon = np.zeros((2), dtype=np.float64)
thorizon = [20., 100.]               # time horizon for ARTP and AGTP diagnostic calculation (years)
year_ref = [2015, 2015]              # reference year for ARTP and AGTP diagnostic calculation (years)

# temperature response function parameters (Collins et al., 2013; Shindell et al., 2010)

tr_c = [0.541, 0.368]                # in K (W/m^2)^-1 - modified by Shindell
tr_d = [8.4, 409.5]                  # in yr

# target equilibrium climate sensitivity

ecs = 3.7/3.7                        # in K (W/m^2)^-1 (CMIP6, Meehl et al., 2021).
if realm == 'low_ecs':
  ecs = ecs * 0.7 
if realm == 'high_ecs':
  ecs = ecs * 1.3

# regional and species information
  
lats0 = 4
slats0 = [' ' for j in range(lats0)]
slats0 = ['ARCTIC', 'NHMDLAT', 'TROPICS', 'SH']
nfsp0 = 6
sfrcspec0 =  [' ' for i in range(nfsp0)]
sfrcspec0 = ['BC', 'BCSN', 'S', 'OC', 'O3', 'CLOUD']

# numerical parameters

fillv = -9999.

# regional temperature sensitivity (K (W/m^2)^-1, Shindell et al., 2010)

cs = [[[0. for k in range(lats0)] for j in range(lats0) ] for i in range(nfsp0)]

# regional temperature sensitivities for black carbon (AMAP, 2015)
# note that the sensitivity in the Arctic is omitted here (zero coefficient)
# - it will be provided for each region and sector using the model data input file

cs[0][0] = [0.  , 0.15, 0.31, 0.06]                             # 60N-90N (Arctic) 
cs[0][1] = [0.08, 0.14, 0.24, 0.07]                             # 28N-60N (mid latitudes)
cs[0][2] = [0.02, 0.07, 0.17, 0.09]                             # 28N-28S (tropics)
cs[0][3] = [0.  , 0.02, 0.06, 0.19]                             # 28S-90S (southern hemisphere)

# regional temperature sensitivities for black carbon in snow (AMAP, 2015).

if case == iamap:
  cs[1][0] = [1.46, 0.45, 0.93, 0.18]                           # 60N-90N (Arctic) 
  cs[1][1] = [0.24, 0.42, 0.62, 0.21]                           # 28N-60N (mid latitudes)
  cs[1][2] = [0.06, 0.21, 0.51, 0.27]                           # 28N-28S (tropics)
  cs[1][3] = [0.  , 0.06, 0.18, 0.57]                           # 28S-90S (southern hemisphere)
else:
  cs[1][0] = [0.49, 0.15, 0.31, 0.06]                           # 60N-90N (Arctic) 
  cs[1][1] = [0.08, 0.14, 0.24, 0.07]                           # 28N-60N (mid latitudes)
  cs[1][2] = [0.02, 0.07, 0.17, 0.09]                           # 28N-28S (tropics)
  cs[1][3] = [0.  , 0.02, 0.06, 0.19]                           # 28S-90S (southern hemisphere)

# regional temperature sensitivities for sulfate (AMAP, 2015)

cs[2][0] = [0.31, 0.17, 0.16, 0.06]                             # 60N-90N (Arctic) 
cs[2][1] = [0.06, 0.24, 0.17, 0.07]                             # 28N-60N (mid latitudes)
cs[2][2] = [0.02, 0.10, 0.24, 0.09]                             # 28N-28S (tropics)
cs[2][3] = [0.  , 0.02, 0.05, 0.19]                             # 28S-90S (southern hemisphere)

# regional temperature sensitivities for organic aerosol (AMAP, 2015)

cs[3][0] = [0.31, 0.17, 0.16, 0.06]                             # 60N-90N (Arctic) 
cs[3][1] = [0.06, 0.24, 0.17, 0.07]                             # 28N-60N (mid latitudes)
cs[3][2] = [0.02, 0.10, 0.24, 0.09]                             # 28N-28S (tropics)
cs[3][3] = [0.  , 0.02, 0.05, 0.19]                             # 28S-90S (southern hemisphere)

# regional temperature sensitivities for ozone (AMAP, 2015)

cs[4][0] = [0.07, 0.05, 0.13, 0.06]                             # 60N-90N (Arctic) 
cs[4][1] = [0.06, 0.20, 0.15, 0.07]                             # 28N-60N (mid latitudes)
cs[4][2] = [0.02, 0.09, 0.26, 0.09]                             # 28N-28S (tropics)
cs[4][3] = [-0.03,-0.06,0.13, 0.19]                             # 28S-90S (southern hemisphere)

# regional temperature sensitivities for clouds (AMAP, 2015)

cs[5][0] = [0.31, 0.17, 0.16, 0.06]                             # 60N-90N (Arctic) 
cs[5][1] = [0.06, 0.24, 0.17, 0.07]                             # 28N-60N (mid latitudes)
cs[5][2] = [0.02, 0.10, 0.24, 0.09]                             # 28N-28S (tropics)
cs[5][3] = [0.  , 0.02, 0.05, 0.19]                             # 28S-90S (southern hemisphere)

# area fractions of latitude bands

area_frc = [0. for j in range(lats0)]
area_frc[0] = 0.0669873                                         # 60N-90N (Arctic) 
area_frc[1] = 0.1982769                                         # 28N-60N (mid latitudes)
area_frc[2] = 0.4694716                                         # 28N-28S (tropics)
area_frc[3] = 0.2652642                                         # 28S-90S (southern hemisphere)

# functions

def get_coord_size(varn,ncfile):                                # get size of the coordinate array
  coord = ncfile.variables[varn][:]
  ncoord = len(coord)
  return ncoord

def get_coord_attr(varn,ncoord,ncfile):                         # get names associated with coordinate indices
  coord = ncfile.variables[varn][:]
  coord_attr = str(coord.attrs)
  coord_attspl = re.split('\d - ',coord_attr)                   # split long_name using integers as separator
  scoord = [" " for i in range(ncoord)]
  for n in range(ncoord):
    npt = n + 1
#    if npt < ncoord:                                             # remove redundant characters at the end of the string
    coord_attspl_len = len(coord_attspl[npt]) - 2
#    else:
#      coord_attspl_len = len(coord_attspl[npt]) - 4
    form = "".join(["{:.","{}".format(coord_attspl_len),"}"])   # determine format of the string
    scoord[n] = form.format(coord_attspl[npt])                   # read the string using the format
    scoord_attspl = re.split('\'',scoord[n])                    # split long_name using integers as separator
    coord_attspl_len = len(scoord_attspl) 
    if coord_attspl_len > 1:
      scoord[n] = scoord_attspl[0]
  return scoord

### read radiative forcing sensitivity data, based on climate model simulations
# read parameter and coordinate information
    
print('Scenario :'+scenario)
print('Radiative forcing sensitivity file: '+model_dat_file)
ncfile = xr.open_dataset(model_dat_file)
lats     = get_coord_size('lats',ncfile)
slats    = get_coord_attr('lats',lats,ncfile)
nfsp     = get_coord_size('frcspec',ncfile)
sfrcspec = get_coord_attr('frcspec',nfsp,ncfile)
nreg     = get_coord_size('region',ncfile)
sregion  = get_coord_attr('region',nreg,ncfile)
nsec     = get_coord_size('sector',ncfile)
ssector  = get_coord_attr('sector',nsec,ncfile)
nssp     = get_coord_size('srcspec',ncfile)
ssrcspec = get_coord_attr('srcspec',nssp,ncfile)

drfrc    = ncfile.variables['drfrc'][:]
srfrc    = ncfile.variables['srfrc'][:]
scfrc    = ncfile.variables['scfrc'][:]
albfrc   = ncfile.variables['albfrc'][:]
dteqdr   = ncfile.variables['dteqdr'][:]
dteqsr   = ncfile.variables['dteqsr'][:]
dteqsc   = ncfile.variables['dteqsc'][:]
dteqalb  = ncfile.variables['dteqalb'][:]
emis     = ncfile.variables['emis'][:]
arbcsens = ncfile.variables['arbcsens'][:]
region   = ncfile.variables['region'][:]
sector   = ncfile.variables['sector'][:]
srcspec  = ncfile.variables['srcspec'][:]
frcspec  = ncfile.variables['frcspec'][:]

# analyze and establish connection between source and forcing species by checking whether
# direct or effective radiative forcings are available in the model data input file

srctofrc = [None for ne in range(nssp)]
for ne in range(nssp):
  ic = 0
  for nf in range(nfsp):
    if not drfrc[:,nf,:,:,ne].isnull().all() or not srfrc[:,nf,:,:,ne].isnull().all() :
      ic = ic + 1
  if ic > 0:
    srctofrc[ne] = [None for i in range(ic)]
  ic = 0
  for nf in range(nfsp):
    if not drfrc[:,nf,:,:,ne].isnull().all() or not srfrc[:,nf,:,:,ne].isnull().all():
      srctofrc[ne][ic] = nf 
      ic = ic + 1
for ne in range(nssp):
  try:
    nfspt = len(srctofrc[ne][:])
  except:
    nfspt = 0
  if nfspt > 0:
    print('Emitted species '+ssrcspec[ne])
    for i in range(nfspt):
      nf = srctofrc[ne][i]
      print('-> Forcing species '+sfrcspec[nf])
albchk = [False for ne in range(nssp)]
for ne in range(nssp):
  if not albfrc[:,:,:,ne].isnull().all() and achk :
    albchk[ne] = True
cldchk = [False for ne in range(nssp)]
for ne in range(nssp):
  if not scfrc[:,:,:,ne].isnull().all() and cchk :
    cldchk[ne] = True

# find and save specified RTP coefficients for given forcing species from model data input file

csidx_fsp = [None for nf in range(nfsp)]
for nf0 in range(nfsp0):
  for nf in range (nfsp):
    if sfrcspec[nf] == sfrcspec0[nf0]:
       csidx_fsp[nf] = nf0
csidx_lat = [None for nl in range(lats)]
for nl0 in range(lats0):
  for nl in range (lats):
    if slats[nl] == slats0[nl0]:
       csidx_lat[nl] = nl0
csidx_alb = sfrcspec0.index('BCSN')
csidx_cld = sfrcspec0.index('CLOUD')

# check if input files exist

def check_input_files(file_name_lcl):
  slcf_file = file_name_lcl
  if os.path.isfile(slcf_file):
    chck = True
  else:
    print ('File '+slcf_file+' does not exist')
    chck = False
  return chck

# configure arrays

spchk1 = [False for n in range(nssp)]
spchk2 = [False for n in range(nssp)]
spchk3 = [False for n in range(nssp)]
spchk4 = [False for n in range(nssp)]
spchk = [False for n in range(nssp)]
ntag_aer = [None for n in range(nssp)]
ntag_aer_base = [None for n in range(nssp)]
ndat_aer = [None for n in range(nssp)]
ndat_aer_base = [None for n in range(nssp)]
year_aer_mltyr = [None for n in range(nssp)]
year_aer_base_mltyr = [None for n in range(nssp)]
year_start = [None for n in range(nssp)]
year_end = [None for n in range(nssp)]
nstpyr = [None for n in range(nssp)]
emis_aer_names = [None for n in range(nssp)]
emis_aer_base_names = [None for n in range(nssp)]
ntag_fire_aer = [None for n in range(nssp)]
nstpyr_fire = [None for n in range(nssp)]
if fires:
  ntag_fire_aer_base = [None for n in range(nssp)]
  ndat_fire_aer = [None for n in range(nssp)]
  ndat_fire_aer_base = [None for n in range(nssp)]
  year_fire_aer_mltyr = [None for n in range(nssp)]
  year_fire_aer_base_mltyr = [None for n in range(nssp)]
  year_fire_start = [None for n in range(nssp)]
  year_fire_end = [None for n in range(nssp)]
  emis_fire_aer_names = [None for n in range(nssp)]
  emis_fire_aer_base_names = [None for n in range(nssp)]

# read emissions, first step

for n in range(nssp):
  emi_dat_file = emidir+ssrcspec[n]+'_base_'+scenario+'.csv'
  emi_base_dat_file = emidir+ssrcspec[n]+'_'+reference+'.csv'
  if fires:
    if not cmip6:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
      emi_fire_base_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
    else:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base_'+scenario+'.csv'
      emi_fire_base_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
  spchk1[n] = check_input_files(emi_dat_file)
  spchk2[n] = check_input_files(emi_base_dat_file)
  spchk[n] = spchk1[n] and spchk2[n]
  if fires:
    spchk3[n] = check_input_files(emi_fire_dat_file)
    spchk4[n] = check_input_files(emi_fire_base_dat_file)
    spchk[n] = spchk[n] and spchk3[n] and spchk4[n]

  if spchk[n]:

#   anthropogenic emissions

    aer_file = emi_dat_file
    aer_base_file = emi_base_dat_file
    emis_aer_in = np.genfromtxt(aer_file, dtype=None, names=True, delimiter=',')      # open file containing time series and read
    lenx = len(emis_aer_in.dtype.names)
    emis_aer_names[n] = emis_aer_in.dtype.names[1:lenx]                               # emission tag names
    ntag_aer[n] = len(emis_aer_names[n])                                              # number of emission tags
    ndat_aer[n] = emis_aer_in['Year'].size                                            # number of data points in input time series
    emis_aer_base_in = np.genfromtxt(aer_base_file, dtype=None, names=True, delimiter=',')   # open file containing time series and read
    lenx = len(emis_aer_base_in.dtype.names)
    emis_aer_base_names[n] = emis_aer_base_in.dtype.names[1:lenx]                     # emission tag names
    ntag_aer_base[n] = len(emis_aer_base_names[n])                                    # number of emission tags
    ndat_aer_base[n] = emis_aer_base_in['Year'].size                                  # number of data points in input time series
    emis_aer_mltyr = [[[0. for i in range(ndat_aer[n])] for j in range(ntag_aer[n])] for ne in range(nssp)]
    emis_aer_base_mltyr = [[[0. for i in range(ndat_aer_base[n])] for j in range(ntag_aer_base[n])] for ne in range(nssp)]

#   fire emissions

    if fires:
      aer_fire_file = emi_fire_dat_file
      aer_fire_base_file = emi_fire_base_dat_file
      emis_fire_aer_in = np.genfromtxt(aer_fire_file, dtype=None, names=True, delimiter=',')  # open file containing time series and read
      lenx = len(emis_fire_aer_in.dtype.names)
      emis_fire_aer_names[n] = emis_fire_aer_in.dtype.names[1:lenx]                     # emission tag names
      ntag_fire_aer[n] = len(emis_fire_aer_names[n])                                    # number of emission tags
      ndat_fire_aer[n] = emis_fire_aer_in['Year'].size                                  # number of data points in input time series
      emis_fire_aer_base_in = np.genfromtxt(aer_fire_base_file, dtype=None, names=True, delimiter=',')   # open file containing time series and read
      lenx = len(emis_fire_aer_base_in.dtype.names)
      emis_fire_aer_base_names[n] = emis_fire_aer_base_in.dtype.names[1:lenx]           # emission tag names
      ntag_fire_aer_base[n] = len(emis_fire_aer_base_names[n])                          # number of emission tags
      ndat_fire_aer_base[n] = emis_fire_aer_base_in['Year'].size                        # number of data points in input time series
      emis_fire_aer_mltyr = [[[0. for i in range(ndat_fire_aer[n])] for j in range(ntag_fire_aer[n])] for ne in range(nssp)]
      emis_fire_aer_base_mltyr = [[[0. for i in range(ndat_fire_aer_base[n])] for j in range(ntag_fire_aer_base[n])] for ne in range(nssp)]
    else:
      ntag_fire_aer[n] = 0

# read emissions, second step

for n in range(nssp):
  emi_dat_file = emidir+ssrcspec[n]+'_base_'+scenario+'.csv'
  emi_base_dat_file = emidir+ssrcspec[n]+'_'+reference+'.csv'
  if fires:
    if not cmip6:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
      emi_fire_base_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
    else:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base_'+scenario+'.csv'
      emi_fire_base_dat_file = emidir+ssrcspec[n]+'_fire_'+reference+'.csv'
  if spchk[n]:
    aer_file = emi_dat_file
    aer_base_file = emi_base_dat_file
    emis_aer_in = np.genfromtxt(aer_file, dtype=None, names=True, delimiter=',')
    year_aer_mltyr[n] = emis_aer_in['Year']                                           # input times
    year_start[n] = year_aer_mltyr[n][0]                                              # first year of the time series
    year_end[n] = year_aer_mltyr[n][ndat_aer[n]-1]                                    # last year of the time series
    nstpyr[n] = year_end[n] - year_start[n] + 1                                       # total number of years
    for k in range(ntag_aer[n]):
      emis_aer_mltyr[n][k] = emis_aer_in[emis_aer_names[n][k]]                        # copy emission time series
    emis_aer_base_in = np.genfromtxt(aer_base_file, dtype=None, names=True, delimiter=',')
    year_aer_base_mltyr[n] = emis_aer_base_in['Year']                                 # input times
    if emis_aer_base_names[n] != emis_aer_names[n] or ntag_aer_base[n] != ntag_aer[n] \
      or ndat_aer_base[n] != ndat_aer[n] or not np.array_equal(year_aer_base_mltyr[n],year_aer_mltyr[n]):
      print ('Inconsistent emission input files')
      sys.exit()
    for k in range(ntag_aer_base[n]):
      emis_aer_base_mltyr[n][k] = emis_aer_base_in[emis_aer_base_names[n][k]]         # copy emission time series

#   fire emissions

    if fires:
      aer_fire_file = emi_fire_dat_file
      aer_fire_base_file = emi_fire_base_dat_file
      emis_fire_aer_in = np.genfromtxt(aer_fire_file, dtype=None, names=True, delimiter=',')
      year_fire_aer_mltyr[n] = emis_fire_aer_in['Year']                                 # input times
      year_fire_start[n] = year_fire_aer_mltyr[n][0]                                    # first year of the time series
      year_fire_end[n] = year_fire_aer_mltyr[n][ndat_fire_aer[n]-1]                     # last year of the time series
      nstpyr_fire[n] = year_fire_end[n] - year_fire_start[n] + 1                        # total number of years
      for k in range(ntag_fire_aer[n]):
        emis_fire_aer_mltyr[n][k] = emis_fire_aer_in[emis_fire_aer_names[n][k]]         # copy emission time series
      emis_fire_aer_base_in = np.genfromtxt(aer_fire_base_file, dtype=None, names=True, delimiter=',')
      year_fire_aer_base_mltyr[n] = emis_fire_aer_base_in['Year']                       # input times
      if emis_fire_aer_base_names[n] != emis_fire_aer_names[n] or ntag_fire_aer_base[n] != ntag_fire_aer[n] \
        or ndat_fire_aer_base[n] != ndat_fire_aer[n] or not np.array_equal(year_fire_aer_base_mltyr[n],year_fire_aer_mltyr[n]):
        print ('Inconsistent fire emission input files')
        sys.exit()
      for k in range(ntag_fire_aer_base[n]):
        emis_fire_aer_base_mltyr[n][k] = emis_fire_aer_base_in[emis_fire_aer_base_names[n][k]] # copy emission time series
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
    for k in range(ntag_aer[n]):
      print(ssrcspec[n]+': '+emis_aer_names[n][k])
    for k in range(ntag_fire_aer[n]):
      print(ssrcspec[n]+': '+emis_fire_aer_names[n][k])

for n in range(nssp):
  if spchk[n]:
    isec = [[None for i in range(ntag_aer[n]+ntag_fire_aer[n])] for ne in range(nssp)]
    ireg = [[None for i in range(ntag_aer[n]+ntag_fire_aer[n])] for ne in range(nssp)]
for ns in range(nsec):
  sectorp = ssector[ns]
  for nr in range(nreg):
    regionp = sregion[nr]
    emi_tag1 = '{}_{}'.format(sectorp,regionp)              # the header in the emission input file must
    emi_tag2 = '{}_{}'.format(regionp,sectorp)              # have the format "sector_region" or "region_sector" (capital letters)
    for n in range(nssp):
      if spchk[n]:
        for k in range(ntag_aer[n]):
          if emis_aer_names[n][k] == emi_tag1 or emis_aer_names[n][k] == emi_tag2:
            isec[n][k] = ns                                 # matching sector index in netcdf input data
            ireg[n][k] = nr                                 # matching region index in netcdf input data
        for k in range(ntag_fire_aer[n]):
          kx = k + ntag_aer[n]
          if emis_fire_aer_names[n][k] == emi_tag1 or emis_fire_aer_names[n][k] == emi_tag2:
            isec[n][kx] = ns                                # matching sector index in netcdf input data
            ireg[n][kx] = nr                                # matching region index in netcdf input data

print('Processing the following regions and sectors:')
for n in range(nssp):
  if spchk[n]:
    ifound = 0
    ifoundf = 0
    emi_dat_file = emidir+ssrcspec[n]+'_base_'+scenario+'.csv'
    if fires:
      emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
    for k in range(ntag_aer[n]):
      if isec[n][k] is not None:
        print(ssrcspec[n]+': '+sregion[ireg[n][k]]+'_'+ssector[isec[n][k]])          # print emission tags that are being processed
        ifound = ifound + 1
    for k in range(ntag_fire_aer[n]):
      kx = k + ntag_aer[n]
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
    emis_aer_yr = [[[0. for i in range(nstpyr[n])] for j in range(ntag_aer[n]+ntag_fire_aer[n])] for ne in range(nssp)]
    emis_aer_base_yr = [[[0. for i in range(nstpyr[n])] for j in range(ntag_aer[n]+ntag_fire_aer[n])] for ne in range(nssp)]
for n in range(nssp):
  if spchk[n]:
    for k in range(ntag_aer[n]):
      if isec[n][k] is not None:
        emis_aer_yr[n][k] = interpol_emis(emis_aer_mltyr[n][k],year_aer_mltyr[n],year_start[n],ndat_aer[n],nstpyr[n])
        emis_aer_base_yr[n][k] = interpol_emis(emis_aer_base_mltyr[n][k],year_aer_mltyr[n],year_start[n],ndat_aer_base[n],nstpyr[n])
    for k in range(ntag_fire_aer[n]):
      kx = k + ntag_aer[n]
      if isec[n][kx] is not None:
        emis_aer_yr[n][kx] = interpol_emis(emis_fire_aer_mltyr[n][k],year_fire_aer_mltyr[n],year_start[n],ndat_fire_aer[n],nstpyr[n])
        emis_aer_base_yr[n][kx] = interpol_emis(emis_fire_aer_base_mltyr[n][k],year_fire_aer_mltyr[n],year_start[n],ndat_fire_aer_base[n],nstpyr[n])

for n in range(nssp):
  if spchk[n]:
    emis_aer_names_tmp = [["" for j in range(ntag_aer[n]+ntag_fire_aer[n])] for ne in range(nssp)]
for n in range(nssp):
  if spchk[n]:
    for k in range(ntag_aer[n]):
      if isec[n][k] is not None:
        emis_aer_names_tmp[n][k] = emis_aer_names[n][k] 
    for k in range(ntag_fire_aer[n]):
      kx = k + ntag_aer[n]
      if isec[n][kx] is not None:
        emis_aer_names_tmp[n][kx] = emis_fire_aer_names[n][k] 
emis_aer_names = emis_aer_names_tmp

for n in range(nssp):
  ntag_aer[n] = ntag_aer[n] + ntag_fire_aer[n]

# implicit equilibrium climate sensitivity of underlying RTP model

ecs0 = sum(tr_c)

fact = [0. for n in range(2)]
for n in range(2):
  fact[n] = dtyrf/tr_d[n]

# arrays
  
nhor = len(year_ref)
for n in range(nssp):
  if spchk[n]:
    agtp_hor = np.zeros((nssp,ntag_aer[n],nhor))
    artp_hor = np.zeros((nssp,ntag_aer[n],lats,nhor))
    temp_rs = np.zeros((4,nssp,nfsp,lats,ntag_aer[n],nstpyr[n]))
    temp = np.zeros((nssp,lats,nstpyr[n]))
    temp_glb = np.zeros((nssp,nstpyr[n]))
    try:
      nfspt = len(srctofrc[n][:])
    except:
      nfspt = 0
    if nfspt > 0:
      scl_dir = np.zeros((nssp,ntag_aer[n],nfspt))
      scl_rad = np.zeros((nssp,ntag_aer[n],nfspt))
      scl_cld = np.zeros((nssp,ntag_aer[n],nfspt))
frceff_dir = np.zeros((lats))
frceff_rad = np.zeros((lats))
frceff_alb = np.zeros((lats))
frceff_cld = np.zeros((lats))
artp_dir = np.zeros((lats))
artp_rad = np.zeros((lats))
artp_alb = np.zeros((lats))
artp_cld = np.zeros((lats))
artph_dir = np.zeros((lats))
artph_rad = np.zeros((lats))
artph_alb = np.zeros((lats))
artph_cld = np.zeros((lats))
csh_dir = np.zeros((lats,lats))
csh_rad = np.zeros((lats,lats))
csh_alb = np.zeros((lats,lats))
csh_cld = np.zeros((lats,lats))

# time coordinate 

time = np.zeros((nstpyra), dtype=np.float64)
time = xr.DataArray(time, name='time', dims=('time'), coords={'time':time}, attrs={'long_name':"Time",'standard_name':"time", 'units':"year"})

# ARTP model time-integration

print('Calculating temperatures...')
for ne in range(nssp):
  if spchk[ne]:
    for i in range(nstpyra):
      time[i]= year_start[ne] + i * dtyrf
for ne in range(nssp):                                                                       # loop over all emitted species
  if spchk[ne]:
    print('Source species: '+ssrcspec[ne])
    temp_chk1 = [0. for j in range(nstpyr[ne])]
    try:
      nfspt = len(srctofrc[ne][:])
    except:
      nfspt = 0
    if nfspt > 0:
      for nft in range(nfspt):                                                               # loop over all forcing species
        nf = srctofrc[ne][nft]
        for kt in range(ntag_aer[ne]):                                                       # loop over all sectors and regions
          if isec[ne][kt] is not None:
            ns = isec[ne][kt]
            nr = ireg[ne][kt]
            if not emis[nr,ns,ne].isnull():
              frceffg_dir = 0.
              frceffg_rad = 0.
              frceffg_alb = 0.
              frceffg_cld = 0.
              rchk1 = True
              rchk2 = True
              rchk3 = True
              rchk4 = True
              for nl in range(lats):
                if not drfrc[nl,nf,nr,ns,ne].isnull():
                  frceff_dir[nl] = drfrc[nl,nf,nr,ns,ne] / (1.E-06*emis[nr,ns,ne]*dtyrf)       # forcing efficiency from direct radiative effect for each latitude band
                else:
                  frceff_dir[nl] = 0.
                  rchk1 = False
                if not srfrc[nl,nf,nr,ns,ne].isnull():
                  frceff_rad[nl] = srfrc[nl,nf,nr,ns,ne] / (1.E-06*emis[nr,ns,ne]*dtyrf)       # forcing efficiency from SLCF/radiation interactions for each latitude band
                else:
                  frceff_rad[nl] = 0.
                  rchk2 = False
                frceffg_dir = frceffg_dir + frceff_dir[nl] * area_frc[csidx_lat[nl]]           # corresponding global mean radiative forcing efficiency, direct forcing
                frceffg_rad = frceffg_rad + frceff_rad[nl] * area_frc[csidx_lat[nl]]           # corresponding global mean radiative forcing efficiency, effective forcing
                if albchk[ne] and not albfrc[nl,nr,ns,ne].isnull():
                  frceff_alb[nl] = albfrc[nl,nr,ns,ne] / (1.E-06*emis[nr,ns,ne]*dtyrf)         # forcing efficiency from impacts on surface albedos for each latitude band
                  frceffg_alb = frceffg_alb + frceff_alb[nl] * area_frc[csidx_lat[nl]]         # corresponding global mean radiative forcing efficiency
                else:
                  frceffg_alb = 0.
                  rchk3 = False
                if cldchk[ne] and not scfrc[nl,nr,ns,ne].isnull():
                  frceff_cld[nl] = scfrc[nl,nr,ns,ne] / (1.E-06*emis[nr,ns,ne]*dtyrf)          # forcing efficiency from SLCF/cloud interactions for each latitude band
                  frceffg_cld = frceffg_cld + frceff_cld[nl] * area_frc[csidx_lat[nl]]         # corresponding global mean radiative forcing efficiency
                else:
                  frceffg_cld = 0.
                  rchk4 = False
              if rchk1 or rchk2 or rchk3 or rchk4:
                agtph_dir = 0.
                agtph_rad = 0.
                for nl in range(lats):
                  for k in range(lats):
                    if k == 0 and nl == 0 and sfrcspec[nf] == 'BC':
                      cst = arbcsens[nr,ns]                                                      # Arctic BC RTP coefficient from model data input file
                      cst = cst.values
                    else:
                      cst = cs[csidx_fsp[nf]][csidx_lat[nl]][csidx_lat[k]]                       # specified RTP coefficient
                    if frceffg_dir != 0.:
                      csh_dir[nl,k] = cst * frceff_dir[k] / frceffg_dir                         # RTP coefficient times forcing efficiency fraction
                    else:
                      csh_dir[nl,k] = 0.
                    if frceffg_rad != 0.:
                      csh_rad[nl,k] = cst * frceff_rad[k] / frceffg_rad
                    else:
                      csh_rad[nl,k] = 0.
                    if sfrcspec[nf] != 'BC' or k > 0:
                      agtph_dir = agtph_dir + csh_dir[nl,k] * area_frc[csidx_lat[nl]]
                      agtph_rad = agtph_rad + csh_rad[nl,k] * area_frc[csidx_lat[nl]]
                if agtph_dir != 0. and scaling_rad[ne]:
                  scl_dir[ne,kt,nft] = ecs0 / agtph_dir                                          # scaling factor for RTP coefficient to ensure
                else:                                                                            # that the sum of ARTPs matches the AGTP
                                                                                                 # if the Arctic forcing is nil. Inverse of efficacy
                  scl_dir[ne,kt,nft] = 1.
                if agtph_rad != 0. and scaling_rad[ne]:
                  scl_rad[ne,kt,nft] = ecs0 / agtph_rad                               
                else:
                  scl_rad[ne,kt,nft] = 1.
                if albchk[ne]:                                                                   # repeat for albedo forcing
                  for nl in range(lats):
                    for k in range(lats):
                      cst = cs[csidx_alb][csidx_lat[nl]][csidx_lat[k]] 
                      if frceffg_alb != 0:
                        csh_alb[nl,k] = cst * frceff_alb[k] / frceffg_alb
                      else:
                        csh_alb[nl,k] = 0.
                if cldchk[ne]:                                                                   # repeat for cloud forcing
                  agtph_cld = 0.
                  for nl in range(lats):
                    for k in range(lats):
                      cst = cs[csidx_cld][csidx_lat[nl]][csidx_lat[k]]
                      if frceffg_cld != 0:
                        csh_cld[nl,k] = cst * frceff_cld[k] / frceffg_cld
                      else:
                        csh_cld[nl,k] = 0.
                      agtph_cld = agtph_cld + csh_cld[nl,k] * area_frc[csidx_lat[nl]]
                  if agtph_cld != 0. and scaling_cld[ne]:
                    scl_cld[ne,kt,nft] = ecs0 / agtph_cld
                  else:
                    scl_cld[ne,kt,nft] = 1.
                for nl in range(lats):
                  artpt_dir = 0.
                  artpt_rad = 0.
                  for k in range(lats):
                    if sfrcspec[nf] != 'BC' or k > 0:
                      artpt_dir = artpt_dir + csh_dir[nl,k] * scl_dir[ne,kt,nft]
                      artpt_rad = artpt_rad + csh_rad[nl,k] * scl_rad[ne,kt,nft]
                    else:
                      artpt_dir = artpt_dir + csh_dir[nl,k]
                      artpt_rad = artpt_rad + csh_rad[nl,k]
                  artph_dir[nl] = artpt_dir
                  artph_rad[nl] = artpt_rad
                  if albchk[ne]:                                                                 # repeat for albedo forcing
                    artpt_alb = 0.
                    for k in range(lats):
                      artpt_alb = artpt_alb + csh_alb[nl,k]
                    artph_alb[nl] = artpt_alb
                  else:
                    artph_alb[nl] = 0.
                  if cldchk[ne]:                                                                 # repeat for cloud forcing
                    artpt_cld = 0.
                    for k in range(lats):
                      artpt_cld = artpt_cld + csh_cld[nl,k] * scl_cld[ne,kt,nft]
                    artph_cld[nl] = artpt_cld
                  else:
                    artph_cld[nl] = 0.
                for i in range(nstpyr[ne]):                                                      # loop over years of pulse emission
                  timeyr_i = year_start[ne] + i * dtyrf
                  demis_i = (emis_aer_yr[ne][kt][i] - emis_aer_base_yr[ne][kt][i]) * dtyrf       # emitted amount of SLCF in one year
                  for j in range(i+1,nstpyr[ne]):                                                # loop over years of temperature response to pulse emission
                    timeyr_j = year_start[ne] + j * dtyrf
                    dt_emis = (j - i) * dtyrf
                    tr = 0.
                    for n in range(2):
                      tr = tr + tr_c[n] * np.exp(-dt_emis/tr_d[n]) * fact[n]
                    tr = tr * (ecs/ecs0)                                                         # scale in order to match specified equilibrium climate sensitivity 
                    agtpt_dir = frceffg_dir * tr                                                 # absolute global temperature potential 
                    agtpt_rad = frceffg_rad * tr
                    agtp_dir = agtpt_dir
                    agtp_rad = agtpt_rad
                    for nl in range(lats):
                      artp_dir[nl] = artph_dir[nl] * (agtp_dir/ecs0)                             # absolute regional temperature potential for each latitude band
                      artp_rad[nl] = artph_rad[nl] * (agtp_rad/ecs0)
                      temp_rs[0,ne,nf,nl,kt,j] = temp_rs[0,ne,nf,nl,kt,j] \
                                                     + artp_dir[nl] * demis_i                    # temperature change in each band for each year, region, and sector
                      temp_rs[1,ne,nf,nl,kt,j] = temp_rs[1,ne,nf,nl,kt,j] \
                                                     + artp_rad[nl] * demis_i                    # temperature change in each band for each year, region, and sector
                      if frcchk:
                        temp[ne,nl,j] = temp[ne,nl,j] + artp_dir[nl] * demis_i                   # temperature change in each band for each year
                      else:
                        temp[ne,nl,j] = temp[ne,nl,j] + artp_rad[nl] * demis_i
                    if frcchk:
                      temp_chk1[j] = temp_chk1[j] + agtp_dir * demis_i
                    else:
                      temp_chk1[j] = temp_chk1[j] + agtp_rad * demis_i
                    if albchk[ne]:                                                               # repeat for albedo forcing
                      agtpt_alb = frceffg_alb * tr
                      agtp_alb = agtpt_alb
                      for nl in range(lats):
                        artp_alb[nl] = artph_alb[nl] * (agtp_alb/ecs0)
                        temp_rs[2,ne,nf,nl,kt,j] = temp_rs[2,ne,nf,nl,kt,j] \
                                                       + artp_alb[nl] * demis_i
                        temp[ne,nl,j] = temp[ne,nl,j] + artp_alb[nl] * demis_i
                      temp_chk1[j] = temp_chk1[j] + agtp_alb * demis_i
                    else:
                      for nl in range(lats):
                        artp_alb[nl] = 0.
                    if cldchk[ne]:                                                               # repeat for cloud forcing
                      agtpt_cld = frceffg_cld * tr
                      agtp_cld = agtpt_cld
                      for nl in range(lats):
                        artp_cld[nl] = artph_cld[nl] * (agtp_cld/ecs0)
                        temp_rs[3,ne,nf,nl,kt,j] = temp_rs[3,ne,nf,nl,kt,j] \
                                                     + artp_cld[nl] * demis_i
                        temp[ne,nl,j] = temp[ne,nl,j] + artp_cld[nl] * demis_i
                      temp_chk1[j] = temp_chk1[j] + agtp_cld * demis_i
                    else:
                      for nl in range(lats):
                        artp_cld[nl] = 0.
                    agtpf = 0.
                    for nl in range(lats):
                      if frcchk:
                        agtpf = agtpf + (artp_dir[nl] + artp_alb[nl] + artp_cld[nl]) \
                                     * area_frc[csidx_lat[nl]]                                 # absolute global temperarture potential diagnosed from ARTPs
                      else:
                        agtpf = agtpf + (artp_rad[nl] + artp_alb[nl] + artp_cld[nl]) \
                                     * area_frc[csidx_lat[nl]]                                 # absolute global temperarture potential diagnosed from ARTPs
                  for ny in range(nhor):                                                       # diagnose results for different specified time horizons
                    if year_ref[ny] == timeyr_i:
                      dt_emis = float(thorizon[ny])
                      tr = 0.
                      for n in range(2):
                        tr = tr + tr_c[n] * np.exp(-dt_emis/tr_d[n]) * fact[n]
                      tr = tr * (ecs/ecs0)
                      if frcchk:
#                        agtpt = (frceffg_dir + frceffg_alb + frceffg_cld) * tr
                        agtpt = frceffg_dir * tr
                      else:
#                        agtpt = (frceffg_rad + frceffg_alb + frceffg_cld) * tr
                        agtpt = frceffg_rad * tr
                      agtp = agtpt
                      agtp_hor[ne,kt,ny] = 0.
                      for nl in range(lats):
                        if frcchk:
                          artp_hor[ne,kt,nl,ny] = (artph_dir[nl] + artph_alb[nl] + artph_cld[nl])* (agtp/ecs0)
                        else:
                          artp_hor[ne,kt,nl,ny] = (artph_rad[nl] + artph_alb[nl] + artph_cld[nl])* (agtp/ecs0)
                        agtp_hor[ne,kt,ny] = agtp_hor[ne,kt,ny] + artp_hor[ne,kt,nl,ny] * area_frc[csidx_lat[nl]]

# global temperature, from regional temperatures
          
for ne in range(nssp):
  if spchk[ne]:
    for j in range(nstpyr[ne]):
      temp_glb[ne,j] = 0.
      for nl in range(lats):
        temp_glb[ne,j] = temp_glb[ne,j] + temp[ne,nl,j] * area_frc[csidx_lat[nl]]

# arrays
        
for n in range(nssp):
  if spchk[n]:
    emis_diff = np.zeros((nssp,ntag_aer[n],nstpyr[n]))
    frc_diff = np.zeros((4,nssp,nfsp,lats,ntag_aer[n],nstpyr[n]))
    temp_rs_eq = np.zeros((4,nssp,nfsp,lats,ntag_aer[n],nstpyr[n]))
    temp_eq = np.zeros((nssp,lats,nstpyr[n]))
    temp_eq_glb = np.zeros((nssp,nstpyr[n]))
frc_dir = np.zeros((lats))
frc_rad = np.zeros((lats))
frc_alb = np.zeros((lats))
frc_cld = np.zeros((lats))

# equilibrium temperatures, RTPs

print('Calculating equilibrium temperatures...')
for ne in range(nssp):
  if spchk[ne]:
    print('Source species: '+ssrcspec[ne])
    try:
      nfspt = len(srctofrc[ne][:])
    except:
      nfspt = 0
    if nfspt > 0:
      for nft in range(nfspt):                                                                    # loop over all forcing species
        nf = srctofrc[ne][nft]
        for i in range(nstpyr[ne]):
          temp_chk2 = 0.
          for kt in range(ntag_aer[ne]):                                                          # loop over all sectors and regions
            if isec[ne][kt] is not None:
              ns = isec[ne][kt]
              nr = ireg[ne][kt]
              if not emis[nr,ns,ne].isnull():
                frcg_dir = 0.
                frcg_rad = 0.
                frcg_alb = 0.
                frcg_cld = 0.
                rchk1 = True
                rchk2 = True
                rchk3 = True
                rchk4 = True
                for nl in range(lats):
                  if not drfrc[nl,nf,nr,ns,ne].isnull():
                    frct_dir = drfrc[nl,nf,nr,ns,ne]
                  else:
                    frct_dir = 0.
                    rchk1 = False
                  if not srfrc[nl,nf,nr,ns,ne].isnull():
                    frct_rad = srfrc[nl,nf,nr,ns,ne]
                  else:
                    frct_rad = 0.
                    rchk2 = False
                  frc_dir[nl] = frct_dir                                                           # forcing for each latitude band
                  frc_rad[nl] = frct_rad                                                           # forcing for each latitude band
                  if albchk[ne] and not albfrc[nl,nr,ns,ne].isnull():                              # repeat for albedo forcing
                    frct = albfrc[nl,nr,ns,ne]
                    frc_alb[nl] = frct.values
                  else:
                    frc_alb[nl] = 0.
                    rchk3 = False
                  if cldchk[ne] and not scfrc[nl,nr,ns,ne].isnull():                               # repeat for cloud forcing
                    frct = scfrc[nl,nr,ns,ne]
                    frc_cld[nl] = frct.values
                  else:
                    frc_cld[nl] = 0.
                    rchk4 = False
                  frcg_dir = frcg_dir + frc_dir[nl] * area_frc[csidx_lat[nl]]                      # global mean radiative forcing
                  frcg_rad = frcg_rad + frc_rad[nl] * area_frc[csidx_lat[nl]]
                  frcg_alb = frcg_alb + frc_alb[nl] * area_frc[csidx_lat[nl]]
                  frcg_cld = frcg_cld + frc_cld[nl] * area_frc[csidx_lat[nl]]
                if rchk1 or rchk2 or rchk3 or rchk4:
                  demis_i = emis_aer_yr[ne][kt][i] - emis_aer_base_yr[ne][kt][i]                     # SLCF emission flux change
                  emis_diff[ne,kt,i] = emis_aer_yr[ne][kt][i]
                  emifrc = demis_i / (1.E-06*emis[nr,ns,ne])
                  for nl in range(lats):
                    frc_diff[0,ne,nf,nl,kt,i] = frc_dir[nl] * emifrc
                    frc_diff[1,ne,nf,nl,kt,i] = frc_rad[nl] * emifrc
                    frc_diff[2,ne,nf,nl,kt,i] = frc_alb[nl] * emifrc
                    frc_diff[3,ne,nf,nl,kt,i] = frc_cld[nl] * emifrc
                    for k in range(lats):
                      if k == 0 and nl == 0 and sfrcspec[nf] == 'BC':
                        cst = arbcsens[nr,ns]                                                         # Arctic BC RTP coefficient from model data input file
                        cst = cst.values
                      else:
                        cst = cs[csidx_fsp[nf]][csidx_lat[nl]][csidx_lat[k]]                          # specified RTP coefficient
                      if sfrcspec[nf] != 'BC' or k > 0:
                        cst_dir = cst * scl_dir[ne,kt,nft]
                        cst_rad = cst * scl_rad[ne,kt,nft]
                      else:
                        cst_dir = cst
                        cst_rad = cst
                      temp_rs_eq[0,ne,nf,nl,kt,i] = temp_rs_eq[0,ne,nf,nl,kt,i] + (ecs/ecs0) * frc_dir[k] * cst_dir * emifrc.values
                      temp_rs_eq[1,ne,nf,nl,kt,i] = temp_rs_eq[1,ne,nf,nl,kt,i] + (ecs/ecs0) * frc_rad[k] * cst_rad * emifrc.values
                      if frcchk:
                        temp_eq[ne,nl,i] = temp_eq[ne,nl,i] + (ecs/ecs0) * frc_dir[k] * cst_dir * emifrc.values
                      else:
                        temp_eq[ne,nl,i] = temp_eq[ne,nl,i] + (ecs/ecs0) * frc_rad[k] * cst_rad * emifrc.values
                      if albchk[ne]:                                                                   # repeat for albedo forcing
                        cst_alb = cs[csidx_alb][csidx_lat[nl]][csidx_lat[k]] 
                        temp_rs_eq[2,ne,nf,nl,kt,i] = temp_rs_eq[2,ne,nf,nl,kt,i] + (ecs/ecs0) * frc_alb[k] * cst_alb * emifrc.values
                        temp_eq[ne,nl,i] = temp_eq[ne,nl,i] + (ecs/ecs0) * frc_alb[k] * cst_alb * emifrc.values
                      if cldchk[ne]:                                                                   # repeat for cloud forcing
                        cst_cld = cs[csidx_cld][csidx_lat[nl]][csidx_lat[k]] *scl_cld[ne,kt,nft]
                        temp_rs_eq[3,ne,nf,nl,kt,i] = temp_rs_eq[3,ne,nf,nl,kt,i] + (ecs/ecs0) * frc_cld[k] * cst_cld * emifrc.values
                        temp_eq[ne,nl,i] = temp_eq[ne,nl,i] + (ecs/ecs0) * frc_cld[k] * cst_cld * emifrc.values
                  if frcchk:
                    temp_chk2 = temp_chk2 + (frcg_dir + frcg_alb + frcg_cld) * ecs * emifrc.values
                  else:
                    temp_chk2 = temp_chk2 + (frcg_rad + frcg_alb + frcg_cld) * ecs * emifrc.values
                  temp_eq_glb[ne,i] = 0.
          for nl in range(lats):
            temp_eq_glb[ne,i] = temp_eq_glb[ne,i] + temp_eq[ne,nl,i] * area_frc[csidx_lat[nl]]

# output arrays
            
temp_ARCTIC = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_NHML = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_TROPICS = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_SH = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_ARCTIC = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_NHML = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_TROPICS = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_SH = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_raw_ARCTIC = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_raw_NHML = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_raw_TROPICS = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
temp_eq_raw_SH = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_ARCTIC = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_NHML = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_TROPICS = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_SH = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_raw_ARCTIC = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_raw_NHML = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_raw_TROPICS = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
forcing_raw_SH = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
emissions = np.full((nstpyra,nreg,nsec,nssp), fillv, dtype=np.float64)
artp_ARCTIC = np.full((nhor,nreg,nsec,nssp), fillv, dtype=np.float64)
artp_NHML = np.full((nhor,nreg,nsec,nssp), fillv, dtype=np.float64)
artp_TROPICS = np.full((nhor,nreg,nsec,nssp), fillv, dtype=np.float64)
artp_SH = np.full((nhor,nreg,nsec,nssp), fillv, dtype=np.float64)

# copy emulator results to output arrays

for ne in range(nssp):                                                                       # loop over all emitted species
  if spchk[ne]:
    try:
      nfspt = len(srctofrc[ne][:])
    except:
      nfspt = 0
    if nfspt > 0:
      for nft in range(nfspt):                                                               # loop over all forcing species
        nf = srctofrc[ne][nft]
        for kt in range(ntag_aer[ne]):                                                       # loop over all sectors and regions
          if isec[ne][kt] is not None:
            ns = isec[ne][kt]
            nr = ireg[ne][kt]
            nx = 0                                                                           # temperature changes due to dir. radiative forcing
            if not drfrc[:,nf,nr,ns,ne].isnull().all():
              temp_ARCTIC[:,nr,ns,nf,ne,nx]     = temp_rs[nx,ne,nf,0,kt,:]
              temp_NHML[:,nr,ns,nf,ne,nx]       = temp_rs[nx,ne,nf,1,kt,:]
              temp_TROPICS[:,nr,ns,nf,ne,nx]    = temp_rs[nx,ne,nf,2,kt,:]
              temp_SH[:,nr,ns,nf,ne,nx]         = temp_rs[nx,ne,nf,3,kt,:]
              temp_eq_ARCTIC[:,nr,ns,nf,ne,nx]  = temp_rs_eq[nx,ne,nf,0,kt,:]
              temp_eq_NHML[:,nr,ns,nf,ne,nx]    = temp_rs_eq[nx,ne,nf,1,kt,:]
              temp_eq_TROPICS[:,nr,ns,nf,ne,nx] = temp_rs_eq[nx,ne,nf,2,kt,:]
              temp_eq_SH[:,nr,ns,nf,ne,nx]      = temp_rs_eq[nx,ne,nf,3,kt,:]
              forcing_ARCTIC[:,nr,ns,nf,ne,nx]  = frc_diff[nx,ne,nf,0,kt,:]
              forcing_NHML[:,nr,ns,nf,ne,nx]    = frc_diff[nx,ne,nf,1,kt,:]
              forcing_TROPICS[:,nr,ns,nf,ne,nx] = frc_diff[nx,ne,nf,2,kt,:]
              forcing_SH[:,nr,ns,nf,ne,nx]      = frc_diff[nx,ne,nf,3,kt,:]
            if not dteqdr[:,nf,nr,ns,ne].isnull().all():
              temp_eq_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = dteqdr[0,nf,nr,ns,ne]
              temp_eq_raw_NHML[:,nr,ns,nf,ne,nx]    = dteqdr[1,nf,nr,ns,ne]
              temp_eq_raw_TROPICS[:,nr,ns,nf,ne,nx] = dteqdr[2,nf,nr,ns,ne]
              temp_eq_raw_SH[:,nr,ns,nf,ne,nx]      = dteqdr[3,nf,nr,ns,ne]
            if not drfrc[:,nf,nr,ns,ne].isnull().all():
              forcing_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = drfrc[0,nf,nr,ns,ne]
              forcing_raw_NHML[:,nr,ns,nf,ne,nx]    = drfrc[1,nf,nr,ns,ne]
              forcing_raw_TROPICS[:,nr,ns,nf,ne,nx] = drfrc[2,nf,nr,ns,ne]
              forcing_raw_SH[:,nr,ns,nf,ne,nx]      = drfrc[3,nf,nr,ns,ne]
            nx = 1                                                                           # temperature changes due to SLCF/radiation interaction
            if not srfrc[:,nf,nr,ns,ne].isnull().all():
              temp_ARCTIC[:,nr,ns,nf,ne,nx]     = temp_rs[nx,ne,nf,0,kt,:]
              temp_NHML[:,nr,ns,nf,ne,nx]       = temp_rs[nx,ne,nf,1,kt,:]
              temp_TROPICS[:,nr,ns,nf,ne,nx]    = temp_rs[nx,ne,nf,2,kt,:]
              temp_SH[:,nr,ns,nf,ne,nx]         = temp_rs[nx,ne,nf,3,kt,:]
              temp_eq_ARCTIC[:,nr,ns,nf,ne,nx]  = temp_rs_eq[nx,ne,nf,0,kt,:]
              temp_eq_NHML[:,nr,ns,nf,ne,nx]    = temp_rs_eq[nx,ne,nf,1,kt,:]
              temp_eq_TROPICS[:,nr,ns,nf,ne,nx] = temp_rs_eq[nx,ne,nf,2,kt,:]
              temp_eq_SH[:,nr,ns,nf,ne,nx]      = temp_rs_eq[nx,ne,nf,3,kt,:]
              forcing_ARCTIC[:,nr,ns,nf,ne,nx]  = frc_diff[nx,ne,nf,0,kt,:]
              forcing_NHML[:,nr,ns,nf,ne,nx]    = frc_diff[nx,ne,nf,1,kt,:]
              forcing_TROPICS[:,nr,ns,nf,ne,nx] = frc_diff[nx,ne,nf,2,kt,:]
              forcing_SH[:,nr,ns,nf,ne,nx]      = frc_diff[nx,ne,nf,3,kt,:]
            if not dteqsr[:,nf,nr,ns,ne].isnull().all():
              temp_eq_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = dteqsr[0,nf,nr,ns,ne]
              temp_eq_raw_NHML[:,nr,ns,nf,ne,nx]    = dteqsr[1,nf,nr,ns,ne]
              temp_eq_raw_TROPICS[:,nr,ns,nf,ne,nx] = dteqsr[2,nf,nr,ns,ne]
              temp_eq_raw_SH[:,nr,ns,nf,ne,nx]      = dteqsr[3,nf,nr,ns,ne]
            if not srfrc[:,nf,nr,ns,ne].isnull().all():
              forcing_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = srfrc[0,nf,nr,ns,ne]
              forcing_raw_NHML[:,nr,ns,nf,ne,nx]    = srfrc[1,nf,nr,ns,ne]
              forcing_raw_TROPICS[:,nr,ns,nf,ne,nx] = srfrc[2,nf,nr,ns,ne]
              forcing_raw_SH[:,nr,ns,nf,ne,nx]      = srfrc[3,nf,nr,ns,ne]
            nx = 2                                                                           # temperature changes due to SLCF/albedo interaction
            if not albfrc[:,nr,ns,ne].isnull().all():
              temp_ARCTIC[:,nr,ns,nf,ne,nx]     = temp_rs[nx,ne,nf,0,kt,:]
              temp_NHML[:,nr,ns,nf,ne,nx]       = temp_rs[nx,ne,nf,1,kt,:]
              temp_TROPICS[:,nr,ns,nf,ne,nx]    = temp_rs[nx,ne,nf,2,kt,:]
              temp_SH[:,nr,ns,nf,ne,nx]         = temp_rs[nx,ne,nf,3,kt,:]
              temp_eq_ARCTIC[:,nr,ns,nf,ne,nx]  = temp_rs_eq[nx,ne,nf,0,kt,:]
              temp_eq_NHML[:,nr,ns,nf,ne,nx]    = temp_rs_eq[nx,ne,nf,1,kt,:]
              temp_eq_TROPICS[:,nr,ns,nf,ne,nx] = temp_rs_eq[nx,ne,nf,2,kt,:]
              temp_eq_SH[:,nr,ns,nf,ne,nx]      = temp_rs_eq[nx,ne,nf,3,kt,:]
              forcing_ARCTIC[:,nr,ns,nf,ne,nx]  = frc_diff[nx,ne,nf,0,kt,:]
              forcing_NHML[:,nr,ns,nf,ne,nx]    = frc_diff[nx,ne,nf,1,kt,:]
              forcing_TROPICS[:,nr,ns,nf,ne,nx] = frc_diff[nx,ne,nf,2,kt,:]
              forcing_SH[:,nr,ns,nf,ne,nx]      = frc_diff[nx,ne,nf,3,kt,:]
            if not dteqalb[:,nr,ns,ne].isnull().all():
              temp_eq_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = dteqalb[0,nr,ns,ne]
              temp_eq_raw_NHML[:,nr,ns,nf,ne,nx]    = dteqalb[1,nr,ns,ne]
              temp_eq_raw_TROPICS[:,nr,ns,nf,ne,nx] = dteqalb[2,nr,ns,ne]
              temp_eq_raw_SH[:,nr,ns,nf,ne,nx]      = dteqalb[3,nr,ns,ne]
            if not albfrc[:,nr,ns,ne].isnull().all() and nf == ne:
              forcing_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = albfrc[0,nr,ns,ne]
              forcing_raw_NHML[:,nr,ns,nf,ne,nx]    = albfrc[1,nr,ns,ne]
              forcing_raw_TROPICS[:,nr,ns,nf,ne,nx] = albfrc[2,nr,ns,ne]
              forcing_raw_SH[:,nr,ns,nf,ne,nx]      = albfrc[3,nr,ns,ne]
            nx = 3                                                                           # temperature changes due to SLCF/cloud interaction
            if not scfrc[:,nr,ns,ne].isnull().all():
              temp_ARCTIC[:,nr,ns,nf,ne,nx]     = temp_rs[nx,ne,nf,0,kt,:]
              temp_NHML[:,nr,ns,nf,ne,nx]       = temp_rs[nx,ne,nf,1,kt,:]
              temp_TROPICS[:,nr,ns,nf,ne,nx]    = temp_rs[nx,ne,nf,2,kt,:]
              temp_SH[:,nr,ns,nf,ne,nx]         = temp_rs[nx,ne,nf,3,kt,:]
              temp_eq_ARCTIC[:,nr,ns,nf,ne,nx]  = temp_rs_eq[nx,ne,nf,0,kt,:]
              temp_eq_NHML[:,nr,ns,nf,ne,nx]    = temp_rs_eq[nx,ne,nf,1,kt,:]
              temp_eq_TROPICS[:,nr,ns,nf,ne,nx] = temp_rs_eq[nx,ne,nf,2,kt,:]
              temp_eq_SH[:,nr,ns,nf,ne,nx]      = temp_rs_eq[nx,ne,nf,3,kt,:]
              forcing_ARCTIC[:,nr,ns,nf,ne,nx]  = frc_diff[nx,ne,nf,0,kt,:]
              forcing_NHML[:,nr,ns,nf,ne,nx]    = frc_diff[nx,ne,nf,1,kt,:]
              forcing_TROPICS[:,nr,ns,nf,ne,nx] = frc_diff[nx,ne,nf,2,kt,:]
              forcing_SH[:,nr,ns,nf,ne,nx]      = frc_diff[nx,ne,nf,3,kt,:]
            if not dteqsc[:,nr,ns,ne].isnull().all():
              temp_eq_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = dteqsc[0,nr,ns,ne]
              temp_eq_raw_NHML[:,nr,ns,nf,ne,nx]    = dteqsc[1,nr,ns,ne]
              temp_eq_raw_TROPICS[:,nr,ns,nf,ne,nx] = dteqsc[2,nr,ns,ne]
              temp_eq_raw_SH[:,nr,ns,nf,ne,nx]      = dteqsc[3,nr,ns,ne]
            if not scfrc[:,nr,ns,ne].isnull().all() and nf == ne:
              forcing_raw_ARCTIC[:,nr,ns,nf,ne,nx]  = scfrc[0,nr,ns,ne]
              forcing_raw_NHML[:,nr,ns,nf,ne,nx]    = scfrc[1,nr,ns,ne]
              forcing_raw_TROPICS[:,nr,ns,nf,ne,nx] = scfrc[2,nr,ns,ne]
              forcing_raw_SH[:,nr,ns,nf,ne,nx]      = scfrc[3,nr,ns,ne]
            emissions[:,nr,ns,ne] = emis_diff[ne,kt,:]
            artp_ARCTIC[:,nr,ns,ne] = 1000. * artp_hor[ne,kt,0,:]
            artp_NHML[:,nr,ns,ne] = 1000. * artp_hor[ne,kt,1,:]
            artp_TROPICS[:,nr,ns,ne] = 1000. * artp_hor[ne,kt,2,:]
            artp_SH[:,nr,ns,ne] = 1000. * artp_hor[ne,kt,3,:]

# global temperature

temp_glb = np.full((nstpyra,nreg,nsec,nfsp,nssp,4), fillv, dtype=np.float64)
for nt in range(nstpyra):
  for nr in range(nreg):
    for ns in range(nsec):
      for nf in range(nfsp):
        for ne in range(nssp):
          for npr in range(4):
            if not temp_ARCTIC[nt,nr,ns,nf,ne,npr] == fillv:
              temp_glb[nt,nr,ns,nf,ne,npr] = area_frc[0]*temp_ARCTIC[nt,nr,ns,nf,ne,npr] \
                                           + area_frc[1]*temp_NHML[nt,nr,ns,nf,ne,npr] \
                                           + area_frc[2]*temp_TROPICS[nt,nr,ns,nf,ne,npr] \
                                           + area_frc[3]*temp_SH[nt,nr,ns,nf,ne,npr]
# write output files

region  = np.arange((nreg), dtype=np.float64)
sector  = np.arange((nsec), dtype=np.float64)
srcspec = np.arange((nssp), dtype=np.float64)
frcspec = np.arange((nfsp), dtype=np.float64)
process = np.arange((4), dtype=np.float64)
region[:] = region[:] + 1
sector[:] = sector[:] + 1
srcspec[:] = srcspec[:] + 1
frcspec[:] = frcspec[:] + 1
process[:] = process[:] + 1
region = xr.DataArray(region,name='region', dims=('region'), coords={'region':region}, \
                      attrs={'long_name':"Source regions: 1 - ASIA, 2 - ROEUR, 3 - AWEST, 4 - AEAST, 5 - ROW, 6 - ARCTIC, 7 - NONARC", 'standard_name':"region"})
sector = xr.DataArray(sector,name='sector', dims=('sector'), coords={'sector':sector}, \
                      attrs={'long_name':"Source sectors: 1 - SURF, 2 - FLAR, 3 - FIRE, 4 - SHIP", 'standard_name':"sector"})
srcspec = xr.DataArray(srcspec,name='srcspec', dims=('srcspec'), coords={'srcspec':srcspec}, \
                      attrs={'long_name':"Source species: 1 - BC, 2 - S, 3 - OC", 'standard_name':"srcspec"})
frcspec = xr.DataArray(frcspec,name='frcspec', dims=('frcspec'), coords={'frcspec':frcspec}, \
                      attrs={'long_name':"Forcing species: 1 - BC, 2 - S, 3 - OC", 'standard_name':"frcspec"})
process = xr.DataArray(process,name='process', dims=('process'), coords={'process':process}, \
                      attrs={'long_name':"Forcing process: 1 - direct radiative forcing, 2 - SLCF/radiation, 3 - SLCF/albedo, 4 - SLCF/cloud", 'standard_name':"process"})

thorizon = xr.DataArray(thorizon, name='thorizon', dims=('thorizon'), coords={'thorizon':thorizon}, attrs={'long_name':"Time horizon",'standard_name':"time_horizon", 'units':"year"})
ds0 = xr.Dataset({'temp_glb': (['time','region','sector','frcspec','srcspec','process'], temp_glb)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds0.temp_glb.attrs = {('long_name', 'Global temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Global_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds1 = xr.Dataset({'temp_ARCTIC': (['time','region','sector','frcspec','srcspec','process'], temp_ARCTIC)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds1.temp_ARCTIC.attrs = {('long_name', 'Arctic temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds2 = xr.Dataset({'temp_NHML': (['time','region','sector','frcspec','srcspec','process'], temp_NHML)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds2.temp_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds3 = xr.Dataset({'temp_TROPICS': (['time','region','sector','frcspec','srcspec','process'], temp_TROPICS)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds3.temp_TROPICS.attrs = {('long_name', 'Tropical temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds4 = xr.Dataset({'temp_SH': (['time','region','sector','frcspec','srcspec','process'], temp_SH)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds4.temp_SH.attrs = {('long_name', 'Southern Hemisphere temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds5 = xr.Dataset({'temp_eq_ARCTIC': (['time','region','sector','frcspec','srcspec','process'], temp_eq_ARCTIC)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds5.temp_eq_ARCTIC.attrs = {('long_name', 'Arctic equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds6 = xr.Dataset({'temp_eq_NHML': (['time','region','sector','frcspec','srcspec','process'], temp_eq_NHML)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds6.temp_eq_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds7 = xr.Dataset({'temp_eq_TROPICS': (['time','region','sector','frcspec','srcspec','process'], temp_eq_TROPICS)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds7.temp_eq_TROPICS.attrs = {('long_name', 'Tropical equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds8 = xr.Dataset({'temp_eq_SH': (['time','region','sector','frcspec','srcspec','process'], temp_eq_SH)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds8.temp_eq_SH.attrs = {('long_name', 'Southern Hemisphere equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds9 = xr.Dataset({'forcing_ARCTIC': (['time','region','sector','frcspec','srcspec','process'], forcing_ARCTIC)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds9.forcing_ARCTIC.attrs = {('long_name', 'Arctic radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds10 = xr.Dataset({'forcing_NHML': (['time','region','sector','frcspec','srcspec','process'], forcing_NHML)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds10.forcing_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds11 = xr.Dataset({'forcing_TROPICS': (['time','region','sector','frcspec','srcspec','process'], forcing_TROPICS)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds11.forcing_TROPICS.attrs = {('long_name', 'Tropical radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds12 = xr.Dataset({'forcing_SH': (['time','region','sector','frcspec','srcspec','process'], forcing_SH)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds12.forcing_SH.attrs = {('long_name', 'Southern Hemisphere radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds13 = xr.Dataset({'emissions': (['time','region','sector','srcspec'], emissions)}, coords={'time':time,'region':region,'sector':sector,'srcspec':srcspec})
ds13.emissions.attrs = {('long_name', 'Regional SLCF emissions'),
                 ('standard_name', 'regional_SLCF_emissions'),
                 ('units', 'Gg/year'),
                 ('_FillValue', fillv)}
ds14 = xr.Dataset({'temp_eq_raw_ARCTIC': (['time','region','sector','frcspec','srcspec','process'], temp_eq_raw_ARCTIC)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds14.temp_eq_raw_ARCTIC.attrs = {('long_name', 'Arctic equilibrium temperature perturbation associated with regional reference SLCF emissions'),
                 ('standard_name', 'Arctic_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds15 = xr.Dataset({'temp_eq_raw_NHML': (['time','region','sector','frcspec','srcspec','process'], temp_eq_raw_NHML)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds15.temp_eq_raw_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude equilibrium temperature perturbation associated with regional reference SLCF emissions'),
                 ('standard_name', 'NH_midlat_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds16 = xr.Dataset({'temp_eq_raw_TROPICS': (['time','region','sector','frcspec','srcspec','process'], temp_eq_raw_TROPICS)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds16.temp_eq_raw_TROPICS.attrs = {('long_name', 'Tropical equilibrium temperature perturbation associated with regional reference SLCF emissions'),
                 ('standard_name', 'Tropics_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds17 = xr.Dataset({'temp_eq_raw_SH': (['time','region','sector','frcspec','srcspec','process'], temp_eq_raw_SH)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds17.temp_eq_raw_SH.attrs = {('long_name', 'Southern Hemisphere equilibrium temperature perturbation associated with regional reference SLCF emissions'),
                 ('standard_name', 'SH_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds18 = xr.Dataset({'artp_ARCTIC': (['thorizon','region','sector','srcspec'], artp_ARCTIC)}, coords={'thorizon':thorizon,'region':region,'sector':sector,'srcspec':srcspec})
ds18.artp_ARCTIC.attrs = {('long_name', 'Arctic absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_ARCTIC'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds19 = xr.Dataset({'artp_NHML': (['thorizon','region','sector','srcspec'], artp_NHML)}, coords={'thorizon':thorizon,'region':region,'sector':sector,'srcspec':srcspec})
ds19.artp_NHML.attrs = {('long_name', 'NH midlatitude absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_NHML'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds20 = xr.Dataset({'artp_TROPICS': (['thorizon','region','sector','srcspec'], artp_TROPICS)}, coords={'thorizon':thorizon,'region':region,'sector':sector,'srcspec':srcspec})
ds20.artp_TROPICS.attrs = {('long_name', 'Tropics absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_TROPICS'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds21 = xr.Dataset({'artp_SH': (['thorizon','region','sector','srcspec'], artp_SH)}, coords={'thorizon':thorizon,'region':region,'sector':sector,'srcspec':srcspec})
ds21.artp_SH.attrs = {('long_name', 'Southern Hemisphere absolute regional temperarture change potential'),
                 ('standard_name', 'ARTP_SH'),
                 ('units', 'mK/Tg'),
                 ('_FillValue', fillv)}
ds22 = xr.Dataset({'forcing_raw_ARCTIC': (['time','region','sector','frcspec','srcspec','process'], forcing_raw_ARCTIC)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds22.forcing_raw_ARCTIC.attrs = {('long_name', 'Arctic equilibrium radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_radiative_forcing_raw'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds23 = xr.Dataset({'forcing_raw_NHML': (['time','region','sector','frcspec','srcspec','process'], forcing_raw_NHML)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds23.forcing_raw_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude equilibrium radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_radiative_forcing_raw'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds24 = xr.Dataset({'forcing_raw_TROPICS': (['time','region','sector','frcspec','srcspec','process'], forcing_raw_TROPICS)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds24.forcing_raw_TROPICS.attrs = {('long_name', 'Tropical radiative equilibrium forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_radiative_forcing_raw'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds25 = xr.Dataset({'forcing_raw_SH': (['time','region','sector','frcspec','srcspec','process'], forcing_raw_SH)}, coords={'time':time,'region':region,'sector':sector,'frcspec':frcspec,'srcspec':srcspec,'process':process})
ds25.forcing_raw_SH.attrs = {('long_name', 'Southern Hemisphere equilibrium radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_radiative_forcing_raw'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds = xr.merge([ds0,ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9, ds10, ds11, ds12, ds13, ds14, ds15, ds16, ds17, ds18, ds19, ds20, ds21, ds22, ds23, ds24, ds25])
ds.attrs = {('comment', 'Model: '+model[case-1]+'; contact: knut.vonsalzen@canada.ca'),
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
outfile = './temp_aerosol_'+scenario+'_'+reference+'_'+model[case-1]+'_'+realm+'.nc'
ds.to_netcdf(outfile, unlimited_dims=["time"])
