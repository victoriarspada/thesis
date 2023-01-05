# python3 postprocessor for emulator-simulated temperature changes due to O3
# author: Knut von Salzen, email: knut.vonsalzen@canada.ca
#
# call this code after emulator_O3.py and ensure that "indir"
# is pointing to the directory that contains the output generated by
# emulator_O3.py.

import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date # http://unidata.github.io/netcdf4-python/
import re
import sys

# postprocessor configuration parameters and input file names

modtemp = [ 'temp_o3_CLE_emulator_mean' ]  #  input file names (without file extension)
year_starty = 1995
year_endy = 2014

# input files with preliminary emulator-simulated temperatures, name of the directory that contains the files
    
indir = '../'

# routines for reading netcdf data

fillv = -9999.

# simulated temperature change due to CO2

nszen = len(modtemp)
for nz in range(nszen):
  filein = indir+modtemp[nz]+'.nc'
  ncfile = Dataset(filein,"r",format="NETCDF4")
  str_year = ncfile.variables["time"][:]
  ntime = len(str_year)
  if nz==0:
    temp_mod_in_o3_ARCTIC = np.zeros((nszen,ntime))
    temp_mod_in_o3_glb = np.zeros((nszen,ntime))
  temp_mod_in_o3_ARCTIC[nz,:] = ncfile.variables["temp_ARCTIC"][:]
  temp_mod_in_o3_glb[nz,:] = ncfile.variables["temp_global"][:]
  fillvin = ncfile.variables["temp_global"]._FillValue
  mask = np.array(temp_mod_in_o3_ARCTIC[nz,:])
  temp_mod_in_o3_ARCTIC[nz,:] = np.where(mask==fillvin, 0., mask).copy()
  mask = np.array(temp_mod_in_o3_glb[nz,:])
  temp_mod_in_o3_glb[nz,:] = np.where(mask==fillvin, 0., mask).copy()

temp_mod_avg_ARCTIC = np.zeros((nszen))
temp_mod_avg_glb = np.zeros((nszen))
ntt = 0
for nt in range(ntime):
  if (int(str_year[nt]) >= year_starty) and (int(str_year[nt]) <= year_endy):
    temp_mod_avg_ARCTIC[:] = temp_mod_avg_ARCTIC[:] + temp_mod_in_o3_ARCTIC[:,nt]
    temp_mod_avg_glb[:] = temp_mod_avg_glb[:] + temp_mod_in_o3_glb[:,nt]
    ntt = ntt + 1
temp_mod_avg_ARCTIC = temp_mod_avg_ARCTIC/ntt
temp_mod_avg_glb = temp_mod_avg_glb/ntt

temp_ARCTIC = np.zeros((nszen,ntime))
temp_glb = np.zeros((nszen,ntime))

for nt in range(ntime):
  temp_ARCTIC[:,nt] = temp_mod_in_o3_ARCTIC[:,nt] - temp_mod_avg_ARCTIC[:]
  temp_glb[:,nt] = temp_mod_in_o3_glb[:,nt] - temp_mod_avg_glb[:]

# output

for nz in range(nszen):
  ds1 = xr.Dataset({'temp_ARCTIC': (['time'], temp_ARCTIC[nz,:])}, coords={'time':str_year})
  ds1.temp_ARCTIC.attrs = {('long_name', \
                            'Arctic temperature perturbation associated with regional O3 precursor emissions, relative to '+str(year_starty)+'-'+str(year_endy)), \
                 ('standard_name', 'Arctic_temperature'), \
                 ('units', 'K'), \
                 ('_FillValue', fillv)}
  ds2 = xr.Dataset({'temp_glb': (['time'], temp_glb[nz,:])}, coords={'time':str_year})
  ds2.temp_glb.attrs = {('long_name', \
                            'Global temperature perturbation associated with regional O3 precursor emissions, relative to '+str(year_starty)+'-'+str(year_endy)), \
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
