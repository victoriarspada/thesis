# python3 source code for O3/climate emulator
# author: Knut von Salzen, email: knut.vonsalzen@canada.ca

import numpy as np
import os 
import sys
import xarray as xr
import re

### read scenario information
# possible choices are:
#   scenario = CLE MFR MFR_SDS CFM

# data_in = np.genfromtxt('emulator.csv', dtype="<U100,<U100",names=True, delimiter=',')
# scenario = str(data_in['scenario'])
scenario = 'CFM'

amap_scenario = [ 'CLE', 'MFR', 'MFR_SDS', 'CFM' ]
amap = False
for nsc in range(len(amap_scenario)):
  if scenario == amap_scenario[nsc]:
    amap = True
if amap == False:
  print('Scenario '+scenario+' is not available')
  sys.exit()

# emulator input: file for radiative forcings from climate models

datdir = './netcdf/'
  
# temperature response function parameters (Collins et al., 2013; Shindell et al., 2010)

tr_c = [0.541, 0.368]                    # K (W/m^2)^-1    - modified by Shindell
tr_d = [8.4, 409.5]                      # yr

# target equilibrium climate sensitivity

ecs = 3.9/3.7                  # K (W/m^2)^-1 (CMIP6, mean results, Zelinka et al., 2020).

# numerical parameters

fillv = -9999.

# regional temperature sensitivity (K (W/m^2)^-1, Shindell et al., 2010)

lats0 = 4
slats0 = [' ' for j in range(lats0)]
slats0 = ['ARCTIC', 'NHMDLAT', 'TROPICS', 'SH']
cs = [[0. for k in range(lats0)] for j in range(lats0) ]

# Regional temperature sensitivities for ozone (AMAP, 2015)

cs[0] = [0.07, 0.05, 0.13, 0.06]                             # 60N-90N (Arctic) 
cs[1] = [0.06, 0.20, 0.15, 0.07]                             # 28N-60N (mid latitudes)
cs[2] = [0.02, 0.09, 0.26, 0.09]                             # 28N-28S (tropics)
cs[3] = [-0.03,-0.06,0.13, 0.19]                             # 28S-90S (southern hemisphere)

# area fractions of latitude bands

area_frc = [0. for j in range(lats0)]
area_frc[0] = 0.0669873                                         # 60N-90N (Arctic) 
area_frc[1] = 0.1982769                                         # 28N-60N (mid latitudes)
area_frc[2] = 0.4694716                                         # 28N-28S (tropics)
area_frc[3] = 0.2652642                                         # 28S-90S (southern hemisphere)

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
    if npt < ncoord:                                             # remove redundant characters at the end of the string
      coord_attspl_len = len(coord_attspl[npt]) - 2
    else:
      coord_attspl_len = len(coord_attspl[npt]) - 4
    form = "".join(["{:.","{}".format(coord_attspl_len),"}"])   # determine format of the string
    scoord[n] = form.format(coord_attspl[npt])                   # read the string using the format
    scoord_attspl = re.split('\'',scoord[n])                    # split long_name using integers as separator
    coord_attspl_len = len(scoord_attspl) 
    if coord_attspl_len > 1:
      scoord[n] = scoord_attspl[0]
  return scoord

model_dat_file = datdir+'ozone_rf_'+scenario+'_1990_2050.nc'

ncfile = xr.open_dataset(model_dat_file)
lats     = get_coord_size('lats',ncfile)
slats    = get_coord_attr('lats',lats,ncfile)
ntime    = get_coord_size('time',ncfile)
drfrci   = ncfile.variables['drfrc'][:]
dteqdri  = ncfile.variables['dteqdr'][:]
timei    = ncfile.variables['time'][:]

year_start = int(timei[0].values)
year_end = int(timei[ntime-1].values)
nstpyra = year_end - year_start + 1
dtyr = 1
if lats != lats0:
  print("inconsistent number of latitude bands")
  sys.exit()

# all forcings relative to start of time series

drfrcit = np.zeros((ntime,lats))
for n in range(ntime):
  drfrcit[n,:] = drfrci[n,:] - drfrci[0,:]
  
# annual mean forcing time series from linear interpolation in time

ina = -9
drfrc = np.zeros((nstpyra,lats))
dteqdr = np.zeros((nstpyra,lats))
for n in range(nstpyra):
  iyear = year_start + n * dtyr
  nxr = ina
  for nx in range(ntime-1):
    if iyear >= timei[nx] and iyear <= timei[nx+1]:
      nxr = nx
  if nxr == ina:
    print("incorrect time range")
    sys.exit()
  wgt = (iyear - int(timei[nxr].values))/(int(timei[nxr+1].values) - int(timei[nxr].values))
  drfrc[n,:] = (1. - wgt) * drfrcit[nxr,:] + wgt * drfrcit[nxr+1,:]  
  dteqdr[n,:] = (1. - wgt) * dteqdri[nxr,:] + wgt * dteqdri[nxr+1,:]  
  
# global mean radiative forcing

drfrc_glb = np.zeros((nstpyra))
for n in range(lats):
  drfrc_glb[:] = drfrc_glb[:] + drfrc[:,n] * area_frc[n]

# implicit equilibrium climate sensitivity of underlying model

ecs0 = sum(tr_c)

dtyrf = 1.
fact = [0. for n in range(2)]
for n in range(2):
  fact[n] = dtyrf/tr_d[n]
temp_rs = np.zeros((lats,nstpyra))
temp_rs_eq = np.zeros((lats,nstpyra))
temp_eq = np.zeros((lats,nstpyra))
temp = np.zeros((lats,nstpyra))
rtpfrc = np.zeros((lats))
time = np.zeros((nstpyra), dtype=np.float64)
time = xr.DataArray(time, name='time', dims=('time'), coords={'time':time}, attrs={'long_name':"Time",'standard_name':"time", 'units':"year"})
for i in range(nstpyra):
  time[i]= year_start + i * dtyrf
for i in range(nstpyra):                                                     # loop over years of pulse emission
  timeyr_i = year_start + i * dtyrf
  for nl in range(lats):
    rtpfrc[nl] = 0.
    for k in range(lats):
      rtpfrc[nl]  = rtpfrc[nl] + cs[nl][k] * drfrc[i,k] * (ecs/ecs0)
    temp_rs_eq[nl,i] = rtpfrc[nl]                                            # equilibrium temperature change, adjusted for climate sensitivity
    temp_eq[nl,i] = dteqdr[i,nl]                                             # raw equilibrium temperature change, from input
  for j in range(i+1,nstpyra):                                               # loop over years of temperature response to pulse emission
    timeyr_j = year_start + j * dtyrf
    dt_emis = (j - i) * dtyrf
    tr = 0.
    for n in range(2):
      tr = tr + tr_c[n] * np.exp(-dt_emis/tr_d[n]) * fact[n]
    tr = tr/ecs0                                                             # scale in order to match specified equilibrium climate sensitivity 
    for nl in range(lats):
      temp_rs[nl,j] = temp_rs[nl,j] + rtpfrc[nl] * tr                        # temperature change in each band for each year, region, and sector

# output fields
      
temp_ARCTIC = np.full((nstpyra), fillv, dtype=np.float64)
temp_NHML = np.full((nstpyra), fillv, dtype=np.float64)
temp_TROPICS = np.full((nstpyra), fillv, dtype=np.float64)
temp_SH = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_ARCTIC = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_NHML = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_TROPICS = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_SH = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_raw_ARCTIC = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_raw_NHML = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_raw_TROPICS = np.full((nstpyra), fillv, dtype=np.float64)
temp_eq_raw_SH = np.full((nstpyra), fillv, dtype=np.float64)
forcing_ARCTIC = np.full((nstpyra), fillv, dtype=np.float64)
forcing_NHML = np.full((nstpyra), fillv, dtype=np.float64)
forcing_TROPICS = np.full((nstpyra), fillv, dtype=np.float64)
forcing_SH = np.full((nstpyra), fillv, dtype=np.float64)

temp_ARCTIC[:]     = temp_rs[0,:]
temp_NHML[:]       = temp_rs[1,:]
temp_TROPICS[:]    = temp_rs[2,:]
temp_SH[:]         = temp_rs[3,:]
temp_eq_ARCTIC[:]  = temp_rs_eq[0,:]
temp_eq_NHML[:]    = temp_rs_eq[1,:]
temp_eq_TROPICS[:] = temp_rs_eq[2,:]
temp_eq_SH[:]      = temp_rs_eq[3,:]
temp_eq_raw_ARCTIC[:]  = temp_eq[0,:]
temp_eq_raw_NHML[:]    = temp_eq[1,:]
temp_eq_raw_TROPICS[:] = temp_eq[2,:]
temp_eq_raw_SH[:]      = temp_eq[3,:]
forcing_ARCTIC[:]  = drfrc[:,0]
forcing_NHML[:]    = drfrc[:,1]
forcing_TROPICS[:] = drfrc[:,2]
forcing_SH[:]      = drfrc[:,3]
temp_global = np.zeros((nstpyra))
temp_tot_ARCTIC = np.zeros((nstpyra))
temp_tot_NHML = np.zeros((nstpyra))
temp_tot_TROPICS = np.zeros((nstpyra))
temp_tot_SH = np.zeros((nstpyra))
for n in range(lats):
  temp_global[:] = temp_global[:] + temp_rs[n,:] * area_frc[n]
temp_tot_ARCTIC[:] = temp_rs[0,:]
temp_tot_NHML[:] = temp_rs[1,:]
temp_tot_TROPICS[:] = temp_rs[2,:]
temp_tot_SH[:] = temp_rs[3,:]

# write output files

ds0 = xr.Dataset({'temp_global': (['time'], temp_global)}, coords={'time':time})
ds0.temp_global.attrs = {('long_name', 'Global temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Global_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds1 = xr.Dataset({'temp_ARCTIC': (['time'], temp_ARCTIC)}, coords={'time':time})
ds1.temp_ARCTIC.attrs = {('long_name', 'Arctic temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds2 = xr.Dataset({'temp_NHML': (['time'], temp_NHML)}, coords={'time':time})
ds2.temp_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds3 = xr.Dataset({'temp_TROPICS': (['time'], temp_TROPICS)}, coords={'time':time})
ds3.temp_TROPICS.attrs = {('long_name', 'Tropical temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds4 = xr.Dataset({'temp_SH': (['time'], temp_SH)}, coords={'time':time})
ds4.temp_SH.attrs = {('long_name', 'Southern Hemisphere temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds5 = xr.Dataset({'temp_eq_ARCTIC': (['time'], temp_eq_ARCTIC)}, coords={'time':time})
ds5.temp_eq_ARCTIC.attrs = {('long_name', 'Arctic equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds6 = xr.Dataset({'temp_eq_NHML': (['time'], temp_eq_NHML)}, coords={'time':time})
ds6.temp_eq_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds7 = xr.Dataset({'temp_eq_TROPICS': (['time'], temp_eq_TROPICS)}, coords={'time':time})
ds7.temp_eq_TROPICS.attrs = {('long_name', 'Tropical equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds8 = xr.Dataset({'temp_eq_SH': (['time'], temp_eq_SH)}, coords={'time':time})
ds8.temp_eq_SH.attrs = {('long_name', 'Southern Hemisphere equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_equilibrium_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds9 = xr.Dataset({'temp_eq_raw_ARCTIC': (['time'], temp_eq_raw_ARCTIC)}, coords={'time':time})
ds9.temp_eq_raw_ARCTIC.attrs = {('long_name', 'Raw Arctic equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds10 = xr.Dataset({'temp_eq_raw_NHML': (['time'], temp_eq_raw_NHML)}, coords={'time':time})
ds10.temp_eq_raw_NHML.attrs = {('long_name', 'Raw Northern Hemisphere midlatitude equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds11 = xr.Dataset({'temp_eq_raw_TROPICS': (['time'], temp_eq_raw_TROPICS)}, coords={'time':time})
ds11.temp_eq_raw_TROPICS.attrs = {('long_name', 'Raw Tropical equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds12 = xr.Dataset({'temp_eq_raw_SH': (['time'], temp_eq_raw_SH)}, coords={'time':time})
ds12.temp_eq_raw_SH.attrs = {('long_name', 'Raw Southern Hemisphere equilibrium temperature perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_equilibrium_temperature_raw'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds13 = xr.Dataset({'forcing_ARCTIC': (['time'], forcing_ARCTIC)}, coords={'time':time})
ds13.forcing_ARCTIC.attrs = {('long_name', 'Arctic radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Arctic_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds14 = xr.Dataset({'forcing_NHML': (['time'], forcing_NHML)}, coords={'time':time})
ds14.forcing_NHML.attrs = {('long_name', 'Northern Hemisphere midlatitude radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'NH_midlat_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds15 = xr.Dataset({'forcing_TROPICS': (['time'], forcing_TROPICS)}, coords={'time':time})
ds15.forcing_TROPICS.attrs = {('long_name', 'Tropical radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'Tropics_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds16 = xr.Dataset({'forcing_SH': (['time'], forcing_SH)}, coords={'time':time})
ds16.forcing_SH.attrs = {('long_name', 'Southern Hemisphere radiative forcing perturbation associated with regional SLCF emissions'),
                 ('standard_name', 'SH_radiative_forcing'),
                 ('units', 'Wm-2'),
                 ('_FillValue', fillv)}
ds17 = xr.Dataset({'temp_tot_ARCTIC': (['time'], temp_tot_ARCTIC)}, coords={'time':time})
ds17.temp_tot_ARCTIC.attrs = {('long_name', 'Net Arctic temperature perturbation from all regional SLCF emissions'),
                 ('standard_name', 'Net_Arctic_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds18 = xr.Dataset({'temp_tot_NHML': (['time'], temp_tot_NHML)}, coords={'time':time})
ds18.temp_tot_NHML.attrs = {('long_name', 'Net Northern Hemisphere midlatitude temperature perturbation from all regional SLCF emissions'),
                 ('standard_name', 'Net_NH_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds19 = xr.Dataset({'temp_tot_TROPICS': (['time'], temp_tot_TROPICS)}, coords={'time':time})
ds19.temp_tot_TROPICS.attrs = {('long_name', 'Net Tropical temperature perturbation from all regional SLCF emissions'),
                 ('standard_name', 'Net_Tropics_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds20 = xr.Dataset({'temp_tot_SH': (['time'], temp_tot_SH)}, coords={'time':time})
ds20.temp_tot_SH.attrs = {('long_name', 'Net Southern Hemisphere temperature perturbation from all regional SLCF emissions'),
                 ('standard_name', 'Net_SH_temperature'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
ds = xr.merge([ds0,ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9, ds10, ds11, ds12, ds13, ds14, ds15, ds16, ds17, ds18, ds19, ds20])
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
outfile = './temp_o3_'+scenario+'_emulator_mean.nc'
ds.to_netcdf(outfile, unlimited_dims=["time"])
