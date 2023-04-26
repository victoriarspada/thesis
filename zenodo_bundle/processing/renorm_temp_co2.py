"""
Author: Knut von Salzen
Edited by: Victoria Spada

Uses the AMAP emulator CO2 output as input to generate normalized temperatures, 
relative to 1995-2014, which is the reference time period of the ESM simulations, 
allowing direct comparisons between the AMAP emulator and CMIP6 datasets.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date # http://unidata.github.io/netcdf4-python/
import re
import sys

year_starty = 1995
year_endy = 2014

# routines for reading netcdf data

fillv = -9999.
def get_coord_size(varn,ncfile):                                # get size of the coordinate array
  coord = ncfile.variables[varn][:]
  ncoord = len(coord)
  return ncoord

def get_coord_attr(varn,ncoord,ncfile):                         # get names associated with coordinate indices
  coord_attr = ncfile.variables[varn].long_name
  coord_attspl = re.split('\d - ',coord_attr)                   # split long_name using integers as separator
  scoord = [" " for i in range(ncoord)]
  for n in range(ncoord):
    npt = n + 1
    if npt < ncoord:                                             # remove redundant characters at the end of the string
      coord_attspl_len = len(coord_attspl[npt]) - 2
    else:
      coord_attspl_len = len(coord_attspl[npt])
    form = "".join(["{:.","{}".format(coord_attspl_len),"}"])   # determine format of the string
    scoord[n] = form.format(coord_attspl[npt])                   # read the string using the format
    scoord_attspl = re.split('\'',scoord[n])                    # split long_name using integers as separator
    coord_attspl_len = len(scoord_attspl) 
    if coord_attspl_len > 1:
      scoord[n] = scoord_attspl[0]
  return scoord

# simulated temperature change due to CO2
datadir = 'C:/Users/victo/Downloads/EngSci Year 4 Sem 1/ESC499 - Thesis/zenodo_bundle/AMAP Projections/'

modtemp = [ \
'temp_co2_SSP126_zero_emulator_high_ecs', \
'temp_co2_SSP126_zero_emulator_low_ecs', \
'temp_co2_SSP126_zero_emulator_mean', \
'temp_co2_SSP245_zero_emulator_high_ecs', \
'temp_co2_SSP245_zero_emulator_low_ecs', \
'temp_co2_SSP245_zero_emulator_mean', \
'temp_co2_SSP370_zero_emulator_high_ecs', \
'temp_co2_SSP370_zero_emulator_low_ecs', \
'temp_co2_SSP370_zero_emulator_mean', \
'temp_co2_SSP585_zero_emulator_high_ecs', \
'temp_co2_SSP585_zero_emulator_low_ecs', \
'temp_co2_SSP585_zero_emulator_mean' ]
# modtemp = [ \
# 'temp_co2_CLE_zero_emulator_high_ecs', \
# 'temp_co2_CLE_zero_emulator_low_ecs', \
# 'temp_co2_CLE_zero_emulator_mean', \
# 'temp_co2_MFR_zero_emulator_high_ecs', \
# 'temp_co2_MFR_zero_emulator_low_ecs', \
# 'temp_co2_MFR_zero_emulator_mean', \
# 'temp_co2_MFR_SDS_zero_emulator_high_ecs', \
# 'temp_co2_MFR_SDS_zero_emulator_low_ecs', \
# 'temp_co2_MFR_SDS_zero_emulator_mean', \
# 'temp_co2_CFM_zero_emulator_high_ecs', \
# 'temp_co2_CFM_zero_emulator_low_ecs', \
# 'temp_co2_CFM_zero_emulator_mean' ]
nszen = len(modtemp)

for nz in range(nszen):
  # filein = './'+modtemp[nz]+'.nc'
  filein = datadir+modtemp[nz]+'.nc'
  ncfile = Dataset(filein,"r",format="NETCDF4")
  str_year = ncfile.variables["time"][:]
  ntime = len(str_year)
  nreg     = get_coord_size('region',ncfile)
  sregion  = get_coord_attr("region",nreg,ncfile)
  sector = ncfile.variables["sector"][:]
  nsec = len(sector)
  if nz==0:
    temp_mod_in_co2_ARCTIC = np.zeros((nszen,ntime,nreg,nsec))
    temp_mod_in_co2_glb = np.zeros((nszen,ntime,nreg,nsec))
  temp_mod_in_co2_ARCTIC[nz,:,:] = ncfile.variables["temp_ARCTIC"][:]
  temp_mod_in_co2_glb[nz,:,:] = ncfile.variables["temp_glb"][:]
  fillvin = ncfile.variables["temp_glb"]._FillValue
  mask = np.array(temp_mod_in_co2_ARCTIC[nz,:,:])
  temp_mod_in_co2_ARCTIC[nz,:,:] = np.where(mask==fillvin, 0., mask).copy()
  mask = np.array(temp_mod_in_co2_glb[nz,:,:])
  temp_mod_in_co2_glb[nz,:,:] = np.where(mask==fillvin, 0., mask).copy()

temp_mod_avg_ARCTIC = np.zeros((nszen,nreg,nsec))
temp_mod_avg_glb = np.zeros((nszen,nreg,nsec))
ntt = 0
for nt in range(ntime):
  if (int(str_year[nt]) >= year_starty) and (int(str_year[nt]) <= year_endy):
    temp_mod_avg_ARCTIC[:,:,:] = temp_mod_avg_ARCTIC[:,:,:] + temp_mod_in_co2_ARCTIC[:,nt,:,:]
    temp_mod_avg_glb[:,:,:] = temp_mod_avg_glb[:,:,:] + temp_mod_in_co2_glb[:,nt,:,:]
    ntt = ntt + 1
temp_mod_avg_ARCTIC = temp_mod_avg_ARCTIC/ntt
temp_mod_avg_glb = temp_mod_avg_glb/ntt

temp_ARCTIC = np.zeros((nszen,ntime,nreg,nsec))
temp_glb = np.zeros((nszen,ntime,nreg,nsec))

for nt in range(ntime):
  temp_ARCTIC[:,nt,:,:] = temp_mod_in_co2_ARCTIC[:,nt,:,:] - temp_mod_avg_ARCTIC[:,:,:]
  temp_glb[:,nt,:,:] = temp_mod_in_co2_glb[:,nt,:,:] - temp_mod_avg_glb[:,:,:]

# output

region  = np.arange((nreg), dtype=np.float64)
sector  = np.arange((nsec), dtype=np.float64)
region[:] = region[:] + 1
sector[:] = sector[:] + 1
region = xr.DataArray(region,name='region', dims=('region'), coords={'region':region}, \
                      attrs={'long_name':"Source regions: 1 - ASIA, 2 - ROEUR, 3 - AWEST, 4 - AEAST, 5 - ROW, 6 - ARCTIC, 7 - NONARC", 'standard_name':"region"})
sector = xr.DataArray(sector,name='sector', dims=('sector'), coords={'sector':sector}, \
                      attrs={'long_name':"Source sectors: 1 - SURF, 2 - FLAR, 3 - FIRE, 4 - SHIP", 'standard_name':"sector"})
for nz in range(nszen):
  ds1 = xr.Dataset({'temp_ARCTIC': (['time','region','sector'], temp_ARCTIC[nz,:,:,:])}, coords={'time':str_year,'region':region,'sector':sector})
  ds1.temp_ARCTIC.attrs = {('long_name', \
                            'Arctic temperature perturbation associated with regional CO2 emissions, relative to '+str(year_starty)+'-'+str(year_endy)), \
                 ('standard_name', 'Arctic_temperature'), \
                 ('units', 'K'), \
                 ('_FillValue', fillv)}
  ds2 = xr.Dataset({'temp_glb': (['time','region','sector'], temp_glb[nz,:,:,:])}, coords={'time':str_year,'region':region,'sector':sector})
  ds2.temp_glb.attrs = {('long_name', \
                            'Global temperature perturbation associated with regional CO2 emissions, relative to '+str(year_starty)+'-'+str(year_endy)), \
                 ('standard_name', 'Arctic_temperature'), \
                 ('units', 'K'), \
                 ('_FillValue', fillv)}
  ds = xr.merge([ds1, ds2])
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
  outfile = './'+modtemp[nz]+'_'+str(year_starty)+'_'+str(year_endy)+'.nc'
  ds.to_netcdf(outfile, unlimited_dims=["time"])
