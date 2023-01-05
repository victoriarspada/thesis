# python3 source code for pm25 emulator
# author: Knut von Salzen, email: knut.vonsalzen@canada.ca

import numpy as np
import os 
import sys
import xarray as xr
import re

### read scenario information
# possible choices are:
#   scenario = SSP126 SSP245 SSP370 SSP585 CLE MFR MFR_SDS CFM

# data_in = np.genfromtxt('emulator.csv', dtype="<U100,<U100",names=True, delimiter=',')
# scenario = str(data_in['scenario'])
scenario = 'SSP126'

# check if SSP or AMAP scenario

cmip6_scenario = [ 'SSP245', 'SSP585', 'SSP126', 'SSP370' ]
cmip6 = False
for nsc in range(len(cmip6_scenario)):
  if scenario == cmip6_scenario[nsc]:
    cmip6 = True

# emulator configuration

case = 5                           # 1 - CanAM, 2 - MRI-ESM, 3 - UKESM, 4 - CESM, 5 - emulator
cmip6 = False

# emulator output directory

outdir = './'

### input files
# emulator input 1: file for concentration sensitivity from climate models

datdir = './netcdf/'
if case == 1:
  model_dat_file = datdir+'CanAM5PAM_aq_data.nc'
elif case == 2: 
  model_dat_file = datdir+'MRI-ESM_aq_data.nc'
elif case == 3: 
  model_dat_file = datdir+'UKESM1_aq_data.nc'
elif case == 4: 
  model_dat_file = datdir+'CESM_aq_data.nc'
elif case == 5: 
  model_dat_file = datdir+'emulator_red_aq_data.nc'

# emulator input 2: files for emissions, name of the directory that contains the files
  
datdir = './emissions/'
if not cmip6:
  emidir = datdir+'AMAP/'
else:
  emidir = datdir+'CMIP6/'

# diagnostic parameters

reference = '2015'
fires = False
extra = True
sel = True
if sel:
  selyear = [1990, 2015, 2030, 2050]
  nsel = len(selyear)

# other parameters and i/o routines
  
model = ['CanAM', 'MRI-ESM', 'UKESM', 'CESM', 'emulator_red']
scl_so4 = 132.1369/96.06 # scaling factor SO4 -> (NH4)SO4
fillv = -9999.

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

# read coordinate information from netcdf model data file

ncfile = xr.open_dataset(model_dat_file)
nreg     = get_coord_size('region',ncfile)
sregion  = get_coord_attr('region',nreg,ncfile)
nsec     = get_coord_size('sector',ncfile)
ssector  = get_coord_attr('sector',nsec,ncfile)
nssp     = get_coord_size('srcspec',ncfile)
ssrcspec = get_coord_attr('srcspec',nssp,ncfile)
nlat     = get_coord_size('lat',ncfile)
nlon     = get_coord_size('lon',ncfile)
pm25     = ncfile.variables['pm25'][:]
emis     = ncfile.variables['emis'][:]
lat      = ncfile.variables['lat'][:]
lon      = ncfile.variables['lon'][:]
lat.attrs={'long_name':"Latitude",'standard_name':"latitude", 'units':"degrees_north", 'axis':"Y"}
lon.attrs={'long_name':"Longitude",'standard_name':"longitude", 'units':"degrees_east", 'axis':"X"}

# analyze and establish connection between source and forcing species by checking whether
# direct or effective radiative forcings are available in the model data input file

valid = [None for ne in range(nssp)]
for ne in range(nssp):
  ic = 0
  if not pm25[:,:,ne,:,:].isnull().all() and not emis[:,:,ne].isnull().all() :
    ic = ic + 1
  if ic > 0:
    valid[ne] = True
for ne in range(nssp):
  if valid[ne]:
    print('Emitted species '+ssrcspec[ne])

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
#
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

#   anthropogenic emissions

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
  region = sregion[nr]
  for ns in range(nsec):
    sector = ssector[ns]
    print(region+', '+sector)

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
  sector = ssector[ns]
  for nr in range(nreg):
    region = sregion[nr]
    emi_tag1 = '{}_{}'.format(sector,region)              # the header in the emission input file must
    emi_tag2 = '{}_{}'.format(region,sector)              #   have the format "sector_region" or "region_sector" (capital letters)
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
      if not cmip6:
        emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base.csv'
      else:
        emi_fire_dat_file = emidir+ssrcspec[n]+'_fire_base_'+scenario+'.csv'
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

# output arrays
  
time = np.zeros((nstpyra), dtype=np.float64)
time = xr.DataArray(time, name='time', dims=('time'), coords={'time':time}, attrs={'long_name':"Time",'standard_name':"time", 'units':"year"})
region = np.arange((nreg), dtype=np.float64)
region[:] = region[:] + 1
region = xr.DataArray(region,name='region', dims=('region'), coords={'region':region}, \
                      attrs={'long_name':"Response regions: 1 - ASIA, 2 - ROEUR, 3 - AWEST, 4 - AEAST, 5 - ROW, 6 - ARCTIC, 7 - NONARC", 'standard_name':"region"})
sector = np.arange((nsec), dtype=np.float64)
sector[:] = sector[:] + 1
sector = xr.DataArray(sector,name='sector', dims=('sector'), coords={'sector':sector}, \
                      attrs={'long_name':"Response sectors: 1 - SURF, 2 - FLAR, 3 - FIRE, 4 - SHIP", 'standard_name':"sector"})
srcspec = np.arange((nssp), dtype=np.float64)
srcspec[:] = srcspec[:] + 1
srcspec = xr.DataArray(srcspec,name='srcspec', dims=('srcspec'), coords={'srcspec':srcspec}, \
                      attrs={'long_name':"Source species: 1 - BC, 2 - S, 3 - OC", 'standard_name':"srcspec"})

# PM2.5 scenarios from multiplication of region/sector contributions to PM2.5 by emissions

pm25ker = np.zeros((nlat,nlon), dtype=np.float32)
pm25tmp = np.zeros((nlat,nlon), dtype=np.float32)
pm25tmpt = np.zeros((nlat,nlon), dtype=np.float32)
if not sel:
  pm25sim = np.zeros((nstpyra,nlat,nlon), dtype=np.float32)
  pm25sims = np.zeros((nstpyra,nlat,nlon,nssp), dtype=np.float32)
  pm25simt = np.zeros((nstpyra,nlat,nlon), dtype=np.float32)
  pm25simst = np.zeros((nstpyra,nlat,nlon,nssp), dtype=np.float32)
  pm25tag = np.zeros((nstpyra,nlat,nlon,nreg,nsec,nssp), dtype=np.float32)
  pm25tagx = np.zeros((nstpyra,nlat,nlon,nreg,nsec,nssp), dtype=np.float32)
else:
  pm25sim = np.zeros((nsel,nlat,nlon), dtype=np.float32)
  pm25sims = np.zeros((nsel,nlat,nlon,nssp), dtype=np.float32)
  pm25simt = np.zeros((nsel,nlat,nlon), dtype=np.float32)
  pm25simst = np.zeros((nsel,nlat,nlon,nssp), dtype=np.float32)
  pm25tag = np.zeros((nsel,nlat,nlon,nreg,nsec,nssp), dtype=np.float32)
  pm25tagx = np.zeros((nsel,nlat,nlon,nreg,nsec,nssp), dtype=np.float32)
pm25tag[:,:,:,:,:,:] = fillv
pm25tagx[:,:,:,:,:,:] = fillv
if not sel:
  for ne in range(nssp):                                                                     # loop over all emitted species
    if spchk[ne] and valid[ne]:
      for i in range(nstpyra):
        time[i]= year_start[ne] + i * dtyrf
else:
  time = np.float32(selyear)
for ne in range(nssp):                                                                       # loop over all emitted species
  if spchk[ne] and valid[ne]:
    for kt in range(ntag_aer[ne]):                                                           # loop over all sectors and regions
      if isec[ne][kt] is not None:
        ns = isec[ne][kt]
        nr = ireg[ne][kt]
        if not pm25[nr,ns,ne,:,:].isnull().all() :
          pm25ker = 1.e+9 * pm25[nr,ns,ne,:,:] / (1.E-06*emis[nr,ns,ne]*dtyrf)               # kernel (efficiency factor) for PM2.5
          if ssrcspec[ne] == 'S':
            pm25ker = pm25ker * scl_so4
          ix = 0
          times = time[ix]
          for i in range(nstpyra):                                                           # loop over years
            timet= year_start[ne] + i * dtyrf
            if not sel or (int(times) == int(timet)):
              demis_i = (emis_aer_yr[ne][kt][i] - emis_aer_base_yr[ne][kt][i]) * dtyrf       # emitted amount of SCLF in one year
              demis_it = emis_aer_yr[ne][kt][i] * dtyrf                                      # emitted amount of SCLF in one year
              pm25tmp = pm25ker * demis_i
              pm25tmpt = pm25ker * demis_it
              pm25sim[ix,:,:] = np.add(pm25sim[ix,:,:],pm25tmp)                              # simulated PM2.5 scenario, accounting for all tags
              pm25sims[ix,:,:,ne] = np.add(pm25sims[ix,:,:,ne],pm25tmp)                      # simulated PM2.5 scenario, accounting for all tags
              pm25simt[ix,:,:] = np.add(pm25simt[ix,:,:],pm25tmpt)                           # simulated PM2.5 scenario, accounting for all tags
              pm25simst[ix,:,:,ne] = np.add(pm25simst[ix,:,:,ne],pm25tmpt)                   # simulated PM2.5 scenario, accounting for all tags
              if extra:
                pm25tagx[ix,:,:,nr,ns,ne] = pm25tmp[:,:]                                     # PM2.5 contribution to total PM2.5 for each tag 
                demis_i = emis_aer_yr[ne][kt][i] * dtyrf                                     # emitted amount of SCLF in one year
                pm25tmp = pm25ker * demis_i
                pm25tag[ix,:,:,nr,ns,ne] = pm25tmp[:,:]                                      # absolute PM2.5 contribution to total PM2.5 for each tag
              ix = ix + 1
              if sel:
                if ix < nsel:
                  times = time[ix]

# output
                  
ds1 = xr.Dataset({'pm25': (['time','lat','lon'], pm25sim)}, coords={'time':time,'lat':lat,'lon':lon})
ds1.pm25.attrs = {('long_name', 'PM2.5 concentration anomaly'),
                  ('standard_name', 'pm25_concentration_anomaly'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds2 = xr.Dataset({'pm25bc': (['time','lat','lon'], pm25sims[:,:,:,0])}, coords={'time':time,'lat':lat,'lon':lon})
ds2.pm25bc.attrs = {('long_name', 'PM2.5 concentration anomaly due to BC'),
                  ('standard_name', 'pm25_concentration_anomaly_BC'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds3 = xr.Dataset({'pm25s': (['time','lat','lon'], pm25sims[:,:,:,1])}, coords={'time':time,'lat':lat,'lon':lon})
ds3.pm25s.attrs = {('long_name', 'PM2.5 concentration anomaly due to S'),
                  ('standard_name', 'pm25_concentration_anomaly_S'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds4 = xr.Dataset({'pm25oc': (['time','lat','lon'], pm25sims[:,:,:,2])}, coords={'time':time,'lat':lat,'lon':lon})
ds4.pm25oc.attrs = {('long_name', 'PM2.5 concentration anomaly due to OC'),
                  ('standard_name', 'pm25_concentration_anomaly_OC'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds = xr.merge([ds1,ds2,ds3,ds4])
ds.attrs = {('comment', 'Model: CanAM5/PAM; contact: knut.vonsalzen@canada.ca'),
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
outfile = outdir+'pm25_anomaly_'+scenario+'_ref_'+reference+'_'+model[case-1]+'.nc'
ds.to_netcdf(outfile, unlimited_dims=["time"])

ds1 = xr.Dataset({'pm25': (['time','lat','lon'], pm25simt)}, coords={'time':time,'lat':lat,'lon':lon})
ds1.pm25.attrs = {('long_name', 'PM2.5 concentration anomaly'),
                  ('standard_name', 'pm25_concentration_anomaly'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds2 = xr.Dataset({'pm25bc': (['time','lat','lon'], pm25simst[:,:,:,0])}, coords={'time':time,'lat':lat,'lon':lon})
ds2.pm25bc.attrs = {('long_name', 'PM2.5 concentration anomaly due to BC'),
                  ('standard_name', 'pm25_concentration_anomaly_BC'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds3 = xr.Dataset({'pm25s': (['time','lat','lon'], pm25simst[:,:,:,1])}, coords={'time':time,'lat':lat,'lon':lon})
ds3.pm25s.attrs = {('long_name', 'PM2.5 concentration anomaly due to S'),
                  ('standard_name', 'pm25_concentration_anomaly_S'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds4 = xr.Dataset({'pm25oc': (['time','lat','lon'], pm25simst[:,:,:,2])}, coords={'time':time,'lat':lat,'lon':lon})
ds4.pm25oc.attrs = {('long_name', 'PM2.5 concentration anomaly due to OC'),
                  ('standard_name', 'pm25_concentration_anomaly_OC'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
ds = xr.merge([ds1,ds2,ds3,ds4])
ds.attrs = {('comment', 'Model: CanAM5/PAM; contact: knut.vonsalzen@canada.ca'),
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
outfile = outdir+'pm25_total_'+scenario+'_ref_'+reference+'_'+model[case-1]+'.nc'
ds.to_netcdf(outfile, unlimited_dims=["time"])

if extra:
  ds2 = xr.Dataset({'pm25': (['time','lat','lon','region','sector','srcspec'], pm25tag)}, coords={'time':time,'lat':lat,'lon':lon,'region':region,'sector':sector,'srcspec':srcspec})
  ds2.pm25.attrs = {('long_name', 'PM2.5 concentration associated with emissions from individual regions and sectors'),
                  ('standard_name', 'pm25_concentration_tagged'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
  ds2.attrs = ds.attrs
  outfile = outdir+'pm25_tagged_abs_'+scenario+'_ref_'+reference+'_'+model[case-1]+'.nc'
  ds2.to_netcdf(outfile, unlimited_dims=["time"])

  ds3 = xr.Dataset({'pm25': (['time','lat','lon','region','sector','srcspec'], pm25tagx)}, coords={'time':time,'lat':lat,'lon':lon,'region':region,'sector':sector,'srcspec':srcspec})
  ds3.pm25.attrs = {('long_name', 'PM2.5 concentration associated with emissions from individual regions and sectors relative to base scenario'),
                  ('standard_name', 'pm25_concentration_tagged'),
                  ('units', 'kg m-3'),
                  ('_FillValue', fillv)}
  ds3.attrs = ds.attrs
  outfile = outdir+'pm25_tagged_rel_'+scenario+'_ref_'+reference+'_'+model[case-1]+'.nc'
  ds3.to_netcdf(outfile, unlimited_dims=["time"])
