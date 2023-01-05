import numpy as np
import xarray as xr
from netCDF4 import Dataset # http://unidata.github.io/netcdf4-python/
#from netCDF4 import Dataset,netcdftime,num2date # http://unidata.github.io/netcdf4-python/
import re

netcdf = 1

#data_in = np.genfromtxt('aerosol.in', dtype="<U100",names=True, delimiter=',')
#scenario = str(data_in['scenario'])
#realm = str(data_in['realm'])

scenarioi = [  'SSP370', 'SSP585', 'SSP119', 'SSP126', 'SSP245']
scenarioxi = scenarioi

# scenarioi = [ 'CLE', 'MFR', 'MFR_SDS', 'CFM' ]
# scenarioxi = [ 'CLE', 'MFR', 'MFR_SDS', 'CFM' ]
# scenarioi = scenarioxi

scen = 4
scenario = scenarioi[scen]
scenariox = scenarioxi[scen]
case = 0   # 0 - Arctic, 4 - global                                     

realm_gas = [ 'mean', 'low_ecs', 'high_ecs' ]
nrgas = len(realm_gas)
realm_aer = [ 'mean', 'low_ecs', 'high_ecs', 'low_aci_s', 'high_aci_s', 'low_aci_bc', 'high_aci_bc', \
              'low_aci_oc', 'high_aci_oc', 'low_ari_s', 'high_ari_s', 'low_ari_bc', 'high_ari_bc', \
              'low_ari_oc', 'high_ari_oc', 'low_alb_s', 'high_alb_s', 'low_alb_bc', 'high_alb_bc', \
              'low_alb_oc', 'high_alb_oc' ]
realm_aer = realm_gas
realm_aerx = realm_aer
#realm_aerx = [ 'mean', 'low', 'high', 'low', 'high', 'low', 'high', \
#              'low', 'high', 'low', 'high', 'low', 'high', \
#              'low', 'high', 'low', 'high', 'low', 'high', \
#              'low', 'high' ]
combined = [ True, True, False, False, False, False, \
             False, False, False, False, False, False, \
             False, False, False, False, False, False, \
             False, False ]
nraer = len(realm_aer) # Number of scenarios for aerosols
ncase = int((nraer - 1)/2) - 1 # 1/2 the length of the aerosols list
nper = 3**nraer  #ncase
indx = np.full((ncase,3),0)
for ncs in range(ncase):
  indx[ncs,0] = 3 + ncs*2
  indx[ncs,1] = 0
  indx[ncs,2] = indx[ncs,0] + 1

if case == 0:
  temp_lab = 'temp_ARCTIC'
elif case == 4:
  temp_lab = 'temp_glb'

precur = True

#rmarship = False  # remove Arctic BC shipping emissions
rmarship = True  # remove Arctic BC shipping emissions
arship = [0, 5, 3]  # Arctic BC shipping emissions: species, region, sector

year_start = 1990
year_end = 2050
year_start_avg = 1995
year_end_avg = 2014
refyr = 2015

fillv = -9999.

# routines for reading netcdf data

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
    coord_attspl_len = len(coord_attspl[npt]) - 2
    form = "".join(["{:.","{}".format(coord_attspl_len),"}"])   # determine format of the string
    scoord[n] = form.format(coord_attspl[npt])                   # read the string using the format
    scoord_attspl = re.split('\'',scoord[n])                    # split long_name using integers as separator
    coord_attspl_len = len(scoord_attspl) 
    if coord_attspl_len > 1:
      scoord[n] = scoord_attspl[0]
  return scoord

# simulated temperature change due to CO2

datadir = 'C:/Users/victo/Downloads/EngSci Year 4 Sem 1/ESC499 - Thesis/zenodo_bundle/AMAP Projections/'

for nre in range(nrgas):
  print(scenariox)
  # modconc = datadir+'temp_co2_'+scenariox+'_zero_'+realm_gas[nre]+'.nc'
  modconc = datadir+'temp_co2_'+scenariox+'_zero_emulator_'+realm_gas[nre]+'.nc'
  
  # modconc = datadir+'temp_co2_'+scenariox+'_zero_emulator_'+realm_gas[nre]+'_1995_2014.nc'
  ncfile = Dataset(modconc,"r",format="NETCDF4")
  str_year = ncfile.variables["time"][:]
  ntime = len(str_year)
  nreg     = get_coord_size('region',ncfile)
  sregion  = get_coord_attr("region",nreg,ncfile)
  sector = ncfile.variables["sector"][:]
  nsec = len(sector)
  if nre==0:
    temp_mod_in_co2 = np.zeros((nrgas,ntime,nreg,nsec))
  temp_mod_in_co2[nre,:,:,:] = ncfile.variables[temp_lab][:]

# simulated temperature change due to CH4

temp_mod_in_ch4 = np.zeros((nrgas,ntime,nreg,nsec))
for nre in range(nrgas):
  #modconc = './temp_ch4_'+scenariox+'_zero_'+realm_gas[nre]+'.nc'
  modconc = datadir+'temp_ch4_'+scenariox+'_zero_emulator_'+realm_gas[nre]+'.nc'
  
  # modconc = datadir+'temp_ch4_'+scenariox+'_zero_emulator_'+realm_gas[nre]+'_1995_2014.nc'
#  ncfile = xr.open_dataset(modconc)
  ncfile = Dataset(modconc,"r",format="NETCDF4")
  nfsp     = get_coord_size('frcspec',ncfile)
  nssp     = get_coord_size('srcspec',ncfile)
  temp_mod_in_ch4_tmp = ncfile.variables[temp_lab][:]
  temp_mod_in_ch4[nre,:,:,:] = temp_mod_in_ch4_tmp[:,:,:,0,0]
  if precur:
    for nps in range(nssp):
      temp_mod_in_ch4[nre,:,:,:] = temp_mod_in_ch4[nre,:,:,:] + temp_mod_in_ch4_tmp[:,:,:,0,nps]

ntx = 0
for nt in range(ntime):
  if (int(str_year[nt]) >= year_start) and (int(str_year[nt]) <= year_end):
    ntx = ntx + 1
ntimem = ntx

# simulated temperature change due aerosol and ozone
    
nspc = 3
pltall = True
#pltall = False
if pltall:
  proc_plt = [ 2, 3, 4 ]
else:
  proc_plt = 2  # 1 - direct (CO2, CH4), 2 - ARI, 3 - albedo, 4 - cloud
#
temp_mod_in_aer = np.zeros((nraer,nspc,ntime,nreg,nsec))
temp_mod_in_aer_prc = np.zeros((nraer,nspc,ntime,nreg,nsec,4))
time_mod = np.full((ntimem),0)
if pltall:
  process_sum = proc_plt
else:
  process_sum = [proc_plt]
# datadir = './'

for nre in range(nraer):
  modconc = datadir+'temp_aerosol_'+scenario+'_1750_emulator_'+realm_aer[nre]+'.nc'
  
  # modconc = datadir+'temp_aerosol_'+scenario+'_1750_emulator_'+realm_aer[nre]+'_1995_2014.nc'
  ncfile = Dataset(modconc,"r",format="NETCDF4")
  nfsp     = get_coord_size('frcspec',ncfile)
  nssp     = get_coord_size('srcspec',ncfile)
  nproc    = get_coord_size('process',ncfile)
  process  = ncfile.variables['process'][:]
  sproc    = get_coord_attr('process',nproc,ncfile)
  species  = ncfile.variables['srcspec'][:]
  ssrcspec = get_coord_attr('srcspec',nssp,ncfile)
  timei    = ncfile.variables['time'][:]
  ntimei   = len(timei)
  temp_mod_in_aer_tmp = ncfile.variables[temp_lab][:]
  fillvin = ncfile.variables[temp_lab]._FillValue
  mask = np.array(temp_mod_in_aer_tmp)
  temp_mod_in_aer_tmp = np.where(mask==fillvin, 0., mask).copy()
  if rmarship:
    temp_mod_in_aer_tmp[:,arship[1],arship[2],arship[0],arship[0],:] = 0.
  for npt in range(len(process_sum)):
    for npr in range(nproc):
      for nfs in range(nfsp):
        for nss in range(nssp):
          if nss == nfs: 
            temp_mod_in_aer_prc[nre,nss,:,:,:,npr] = temp_mod_in_aer_tmp[:,:,:,nfs,nss,npr]
      if process[npr] == process_sum[npt]:
        for nfs in range(nfsp):
          for nss in range(nssp):
            if nss == nfs: 
              temp_mod_in_aer[nre,nss,:,:,:] = temp_mod_in_aer[nre,nss,:,:,:] + temp_mod_in_aer_tmp[:,:,:,nfs,nss,npr]
  ntx = 0
  for nt in range(ntime):
    if (int(timei[nt]) >= year_start) and (int(timei[nt]) <= year_end):
      time_mod[ntx] = int(timei[nt])
      ntx = ntx + 1
        
# total simulated temperature change

nspc = 5
temp_mod_in = np.zeros((nraer,nspc,ntime,nreg,nsec))
temp_mod_in[:,0,:,:,:] = temp_mod_in_aer[:,0,:,:,:]
temp_mod_in[:,2,:,:,:] = temp_mod_in_aer[:,1,:,:,:]
temp_mod_in[:,4,:,:,:] = temp_mod_in_aer[:,2,:,:,:]
for nre in range(nraer):
  temp_mod_in[nre,3,:,:,:] = temp_mod_in_ch4[0,:,:,:]
  temp_mod_in[nre,1,:,:,:] = temp_mod_in_co2[0,:,:,:]
for nre in range(nrgas):
  temp_mod_in[nre,3,:,:,:] = temp_mod_in_ch4[nre,:,:,:]
  temp_mod_in[nre,1,:,:,:] = temp_mod_in_co2[nre,:,:,:]
temp_mod_in_prc = np.zeros((nraer,nspc,ntime,nreg,nsec,nproc))
temp_mod_in_prc[:,0,:,:,:,:] = temp_mod_in_aer_prc[:,0,:,:,:,:]
temp_mod_in_prc[:,2,:,:,:,:] = temp_mod_in_aer_prc[:,1,:,:,:,:]
temp_mod_in_prc[:,4,:,:,:,:] = temp_mod_in_aer_prc[:,2,:,:,:,:]
for nre in range(nraer):
  temp_mod_in_prc[nre,3,:,:,:,1] = temp_mod_in_ch4[0,:,:,:]
  temp_mod_in_prc[nre,1,:,:,:,1] = temp_mod_in_co2[0,:,:,:]
for nre in range(nrgas):
  temp_mod_in_prc[nre,3,:,:,:,1] = temp_mod_in_ch4[nre,:,:,:]
  temp_mod_in_prc[nre,1,:,:,:,1] = temp_mod_in_co2[nre,:,:,:]

# sum over all regions and sectors

temp_mod_tmp = np.full((nraer,nspc,ntime),0.)
temp_mod_tmp_prc = np.full((nraer,nspc,ntime,nproc),0.)
temp_mod_tmp_co2 = np.full((nrgas,ntime),0.)
for nr in range(nreg):
  for ns in range(nsec):
    temp_mod_tmp[:,:,:] = temp_mod_tmp[:,:,:] + temp_mod_in[:,:,:,nr,ns]
    temp_mod_tmp_prc[:,:,:,:] = temp_mod_tmp_prc[:,:,:,:] + temp_mod_in_prc[:,:,:,nr,ns,:]
    temp_mod_tmp_co2[:,:] = temp_mod_tmp_co2[:,:] + temp_mod_in[1:nrgas-1,1,:,nr,ns]

# save temperature time series and adjust to match reference time period
    
temp_mod = np.full((nraer,ntimem,nspc), 0.)
temp_mod_prc = np.full((nraer,ntimem,nspc,nproc), 0.)
temp_mod_time = np.full((ntimem), fillv)
temp_mod_ref = np.full((nraer,nspc), 0.)
temp_mod_prc_ref = np.full((nraer,nspc,nproc), 0.)
temp_mod_co2 = np.full((nrgas,ntimem), 0.)
temp_mod_co2_ref = np.full((nrgas), 0.)
ntx = 0
nta = 0
for nt in range(ntime):
  if (int(str_year[nt]) >= year_start) and (int(str_year[nt]) <= year_end):
    temp_mod[:,ntx,:] = temp_mod_tmp[:,:,nt]    
    temp_mod_prc[:,ntx,:,:] = temp_mod_tmp_prc[:,:,nt,:]    
    temp_mod_co2[:,ntx] = temp_mod_tmp_co2[:,nt]    
    temp_mod_time[ntx] = float(str_year[nt])
    if int(str_year[nt]) >= year_start_avg and int(str_year[nt]) <= year_end_avg :
      temp_mod_ref[:,:] = temp_mod_ref[:,:] + temp_mod_tmp[:,:,nt]
      temp_mod_prc_ref[:,:,:] = temp_mod_prc_ref[:,:,:] + temp_mod_tmp_prc[:,:,nt,:]
      temp_mod_co2_ref[:] = temp_mod_co2_ref[:] + temp_mod_tmp_co2[:,nt]
      nta = nta + 1
    ntx = ntx + 1

temp_mod_ref[:,:] = temp_mod_ref[:,:]/float(nta)
temp_mod_prc_ref[:,:,:] = temp_mod_prc_ref[:,:,:]/float(nta)
temp_mod_co2_ref[:] = temp_mod_co2_ref[:]/float(nta)

for nt in range(ntimem):
  temp_mod[:,nt,:] = temp_mod[:,nt,:] - temp_mod_ref[:,:]
  temp_mod_prc[:,nt,:,:] = temp_mod_prc[:,nt,:,:] - temp_mod_prc_ref[:,:,:]
  temp_mod_co2[:,nt] = temp_mod_co2[:,nt] - temp_mod_co2_ref[:]

temp_mod_all = np.full((nraer,ntimem), 0.)
temp_mod_so4 = np.full((nraer,ntimem), 0.)
temp_mod_slcf = np.full((nraer,ntimem), 0.)
temp_mod_co2 = np.full((nraer,ntimem), 0.)
temp_mod_all[:,:] = temp_mod[:,:,0]
for nss in range(nspc-1):
  temp_mod_all[:,:] = temp_mod_all[:,:] + temp_mod[:,:,nss+1]
temp_mod_slcf[:,:] = temp_mod[:,:,0] + temp_mod[:,:,2] + temp_mod[:,:,3] + temp_mod[:,:,4]
temp_mod_so4[:,:] = temp_mod[:,:,2]
temp_mod_co2[:,:] = temp_mod[:,:,1]

# temperature uncertainty from simple error propagation. dt = sqrt(dt_ecs**2 + sum(dt_frc**2)).
# for separate SLCF and CO2 temperature impacts, the uncertainties are partly
# caused by dependent/correlated (via ECS uncertainy) and independent/uncorrelated (via forcing uncertainty) processes,
# yielding an uncertain estimate that falls between these two cases.
  
temp_uncert_pos = np.full((ntimem), 0.)
temp_uncert_neg = np.full((ntimem), 0.)
temp_uncert_pos_slcf = np.full((ntimem), 0.)
temp_uncert_neg_slcf = np.full((ntimem), 0.)
temp_uncert_pos_so4 = np.full((ntimem), 0.)
temp_uncert_neg_so4 = np.full((ntimem), 0.)
temp_uncert_pos_co2 = np.full((ntimem), 0.)
temp_uncert_neg_co2 = np.full((ntimem), 0.)
temp_mean = np.full((ntimem), 0.)
temp_mean_slcf = np.full((ntimem), 0.)
temp_mean_so4 = np.full((ntimem), 0.)
temp_mean_co2 = np.full((ntimem), 0.)
for nre in range(nraer-1):
#for nre in range(2):
  if not combined[nre]:
    for nss in range(nspc):
      for nt in range(ntimem):
        diff = temp_mod[nre+1,nt,nss] - temp_mod[0,nt,nss]
        if diff >= 0 :
          temp_uncert_pos[nt] = temp_uncert_pos[nt] + diff*diff
          if nss == 1:
            temp_uncert_pos_co2[nt] = temp_uncert_pos_co2[nt] + diff*diff
          if nss == 0 or nss == 2 or nss == 4 or nss ==3 :
            temp_uncert_pos_slcf[nt] = temp_uncert_pos_slcf[nt] + diff*diff
            if nt == ntimem-1:
              print(str(nss)+' '+str(diff*diff))
          if nss == 2 :
            temp_uncert_pos_so4[nt] = temp_uncert_pos_so4[nt] + diff*diff
        else:
          temp_uncert_neg[nt] = temp_uncert_neg[nt] + diff*diff
          if nss == 1:
            temp_uncert_neg_co2[nt] = temp_uncert_neg_co2[nt] + diff*diff
          if nss == 0 or nss == 2 or nss == 4 or nss ==3 :
            temp_uncert_neg_slcf[nt] = temp_uncert_neg_slcf[nt] + diff*diff
          if nss == 2 :
            temp_uncert_neg_so4[nt] = temp_uncert_neg_so4[nt] + diff*diff
  else:
    for nt in range(ntimem):
      diff = temp_mod_all[nre+1,nt] - temp_mod_all[0,nt]
      temp_mean[nt] = temp_mod_all[0,nt] 
      if diff >= 0 :
        temp_uncert_pos[nt] = temp_uncert_pos[nt] + diff*diff 
      else:
        temp_uncert_neg[nt] = temp_uncert_neg[nt] + diff*diff
        
      diff = temp_mod_slcf[nre+1,nt] - temp_mod_slcf[0,nt]
      temp_mean_slcf[nt] = temp_mod_slcf[0,nt] 
      if diff >= 0 :
        temp_uncert_pos_slcf[nt] = temp_uncert_pos_slcf[nt] + diff*diff 
      else:
        temp_uncert_neg_slcf[nt] = temp_uncert_neg_slcf[nt] + diff*diff 

      diff = temp_mod_so4[nre+1,nt] - temp_mod_so4[0,nt]
      temp_mean_so4[nt] = temp_mod_so4[0,nt] 
      if diff >= 0 :
        temp_uncert_pos_so4[nt] = temp_uncert_pos_so4[nt] + diff*diff 
      else:
        temp_uncert_neg_so4[nt] = temp_uncert_neg_so4[nt] + diff*diff 
        
      diff = temp_mod_co2[nre+1,nt] - temp_mod_co2[0,nt]
      temp_mean_co2[nt] = temp_mod_co2[0,nt] 
      if diff >= 0 :
        temp_uncert_pos_co2[nt] = temp_uncert_pos_co2[nt] + diff*diff 
      else:
        temp_uncert_neg_co2[nt] = temp_uncert_neg_co2[nt] + diff*diff 

for nt in range(ntimem):
  temp_uncert_pos[nt] = np.sqrt(temp_uncert_pos[nt])
  temp_uncert_neg[nt] = np.sqrt(temp_uncert_neg[nt])
  temp_uncert_pos_slcf[nt] = np.sqrt(temp_uncert_pos_slcf[nt])
  temp_uncert_neg_slcf[nt] = np.sqrt(temp_uncert_neg_slcf[nt])
  temp_uncert_pos_so4[nt] = np.sqrt(temp_uncert_pos_so4[nt])
  temp_uncert_neg_so4[nt] = np.sqrt(temp_uncert_neg_so4[nt])
  temp_uncert_pos_co2[nt] = np.sqrt(temp_uncert_pos_co2[nt])
  temp_uncert_neg_co2[nt] = np.sqrt(temp_uncert_neg_co2[nt])
  if nt == ntimem-1:
    print(str(temp_mod_time[nt])+": all "+str(temp_mean[nt])+"+"+str(temp_uncert_pos[nt])+"-"+str(temp_uncert_neg[nt]))
    print(str(temp_mod_time[nt])+": slcf "+str(temp_mean_slcf[nt])+"+"+str(temp_uncert_pos_slcf[nt])+"-"+str(temp_uncert_neg_slcf[nt]))
    print(str(temp_mod_time[nt])+": so4 "+str(temp_mean_so4[nt])+"+"+str(temp_uncert_pos_so4[nt])+"-"+str(temp_uncert_neg_so4[nt]))
    print(str(temp_mod_time[nt])+": co2 "+str(temp_mean_co2[nt])+"+"+str(temp_uncert_pos_co2[nt])+"-"+str(temp_uncert_neg_co2[nt]))

# check all possible combinations of forcings and climate sensitivity
# and compute temperatures for these cases
    
temp_per = np.full((nper,ntimem), 0.)
info = np.full((nper), '                                                            ')
npt = 0
# for ip1 in range(3):
#   print(indx)
#   in1 = indx[0,ip1]
#   npr = 3  # ACI
#   for nss in range(nspc):
#     for ntx in range(int(nper/3)):
#       temp_per[npt+ntx,:] = temp_per[npt+ntx,:] + temp_mod_prc[in1,:,nss,npr]
#   for ip2 in range(3):
#     in2 = indx[1,ip2]
#     npr = 1  # ARI_BC
#     for nss in range(nspc):
#       if nss == 0:
#         for ntx in range(int(nper/(3**2))):
#           temp_per[npt+ntx,:] = temp_per[npt+ntx,:] + temp_mod_prc[in2,:,nss,npr]
#     for ip3 in range(3):
#       in3 = indx[2,ip3]
#       npr = 1  # ARI_S
#       for nss in range(nspc):
#         if nss == 2:
#           for ntx in range(int(nper/(3**3))):
#             temp_per[npt+ntx,:] = temp_per[npt+ntx,:] + temp_mod_prc[in3,:,nss,npr]
#       for ip4 in range(3):
#         in4 = indx[3,ip4]
#         npr = 2 # ALB
#         for nss in range(nspc):
#           if nss == 0:
#             temp_per[npt,:] = temp_per[npt,:] + temp_mod_prc[in4,:,nss,npr]
#           else:
#             temp_per[npt,:] = temp_per[npt,:] + temp_mod_prc[0,:,nss,npr]
#         npr = 1 # ARI
#         for nss in range(nspc):
#           if not (nss == 0 or nss == 2):
#             temp_per[npt,:] = temp_per[npt,:] + temp_mod_prc[0,:,nss,npr]

#         info[npt] = ','+realm_aerx[in1]+','+realm_aerx[in2]+','+realm_aerx[in3]+','+realm_aerx[in4]
# ##        print(str(temp_per[npt,ntimem-1])+','+realm_aerx[in1]+','+realm_aerx[in2]+','+realm_aerx[in3]+','+realm_aerx[in4])
#         npt = npt + 1

#indxs = np.argsort(temp_per[:,ntimem-1])
#for ip in range(nper):
#  print(str(temp_per[indxs[ip],ntimem-1])+','+str(temp_per[indxs[ip],ntimem-1]-temp_mod_co2[0,ntimem-1])+info[indxs[ip]])

# output
  
time = xr.DataArray(time_mod, name='time', dims=('time'), coords={'time':time_mod}, attrs={'long_name':"Time",'standard_name':"time", 'units':"year"})
if case == 0:
  ds1 = xr.Dataset({'dtemp_ARCTIC_pos': (['time'], temp_uncert_pos)}, coords={'time':time_mod})
  ds1.dtemp_ARCTIC_pos.attrs = {('long_name', 'Arctic temperature process uncertainty'),
                 ('standard_name', 'Arctic_temperature_uncertainty'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds2 = xr.Dataset({'dtemp_ARCTIC_SLCF_pos': (['time'], temp_uncert_pos_slcf)}, coords={'time':time_mod})
  ds2.dtemp_ARCTIC_SLCF_pos.attrs = {('long_name', 'Arctic temperature process uncertainty due to SLCFs'),
                 ('standard_name', 'Arctic_temperature_uncertainty_SLCF'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds3 = xr.Dataset({'dtemp_ARCTIC_CO2_pos': (['time'], temp_uncert_pos_co2)}, coords={'time':time_mod})
  ds3.dtemp_ARCTIC_CO2_pos.attrs = {('long_name', 'Arctic temperature process uncertainty due to CO2'),
                 ('standard_name', 'Arctic_temperature_uncertainty_CO2'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds4 = xr.Dataset({'dtemp_ARCTIC_neg': (['time'], temp_uncert_neg)}, coords={'time':time_mod})
  ds4.dtemp_ARCTIC_neg.attrs = {('long_name', 'Arctic temperature process uncertainty'),
                 ('standard_name', 'Arctic_temperature_uncertainty'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds5 = xr.Dataset({'dtemp_ARCTIC_SLCF_neg': (['time'], temp_uncert_neg_slcf)}, coords={'time':time_mod})
  ds5.dtemp_ARCTIC_SLCF_neg.attrs = {('long_name', 'Arctic temperature process uncertainty due to SLCFs'),
                 ('standard_name', 'Arctic_temperature_uncertainty_SLCF'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds6 = xr.Dataset({'dtemp_ARCTIC_CO2_neg': (['time'], temp_uncert_neg_co2)}, coords={'time':time_mod})
  ds6.dtemp_ARCTIC_CO2_neg.attrs = {('long_name', 'Arctic temperature process uncertainty due to CO2'),
                 ('standard_name', 'Arctic_temperature_uncertainty_CO2'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds7 = xr.Dataset({'dtemp_ARCTIC_SO4_neg': (['time'], temp_uncert_neg_so4)}, coords={'time':time_mod})
  ds7.dtemp_ARCTIC_SO4_neg.attrs = {('long_name', 'Arctic temperature process uncertainty due to SO4'),
                 ('standard_name', 'Arctic_temperature_uncertainty_SO4'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds8 = xr.Dataset({'dtemp_ARCTIC_SO4_pos': (['time'], temp_uncert_pos_so4)}, coords={'time':time_mod})
  ds8.dtemp_ARCTIC_SO4_pos.attrs = {('long_name', 'Arctic temperature process uncertainty due to SO4'),
                 ('standard_name', 'Arctic_temperature_uncertainty_SO4'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  #outfile = './dtemp_ARCTIC_'+scenario+'.nc'
  outfile = datadir+'Ensemble/dtemp_ARCTIC_'+scenario+'.nc'
  print(outfile)
if case == 4:
  ds1 = xr.Dataset({'dtemp_GLOBAL_pos': (['time'], temp_uncert_pos)}, coords={'time':time_mod})
  ds1.dtemp_GLOBAL_pos.attrs = {('long_name', 'Global temperature process uncertainty'),
                 ('standard_name', 'Global_temperature_uncertainty'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds2 = xr.Dataset({'dtemp_GLOBAL_SLCF_pos': (['time'], temp_uncert_pos_slcf)}, coords={'time':time_mod})
  ds2.dtemp_GLOBAL_SLCF_pos.attrs = {('long_name', 'Global temperature process uncertainty due to SLCFs'),
                 ('standard_name', 'Global_temperature_uncertainty_SLCF'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds3 = xr.Dataset({'dtemp_GLOBAL_CO2_pos': (['time'], temp_uncert_pos_co2)}, coords={'time':time_mod})
  ds3.dtemp_GLOBAL_CO2_pos.attrs = {('long_name', 'Global temperature process uncertainty due to CO2'),
                 ('standard_name', 'Global_temperature_uncertainty_CO2'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds4 = xr.Dataset({'dtemp_GLOBAL_neg': (['time'], temp_uncert_neg)}, coords={'time':time_mod})
  ds4.dtemp_GLOBAL_neg.attrs = {('long_name', 'Global temperature process uncertainty'),
                 ('standard_name', 'Global_temperature_uncertainty'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds5 = xr.Dataset({'dtemp_GLOBAL_SLCF_neg': (['time'], temp_uncert_neg_slcf)}, coords={'time':time_mod})
  ds5.dtemp_GLOBAL_SLCF_neg.attrs = {('long_name', 'Global temperature process uncertainty due to SLCFs'),
                 ('standard_name', 'Global_temperature_uncertainty_SLCF'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds6 = xr.Dataset({'dtemp_GLOBAL_CO2_neg': (['time'], temp_uncert_neg_co2)}, coords={'time':time_mod})
  ds6.dtemp_GLOBAL_CO2_neg.attrs = {('long_name', 'Global temperature process uncertainty due to CO2'),
                 ('standard_name', 'Global_temperature_uncertainty_CO2'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds7 = xr.Dataset({'dtemp_GLOBAL_SO4_neg': (['time'], temp_uncert_neg_so4)}, coords={'time':time_mod})
  ds7.dtemp_GLOBAL_SO4_neg.attrs = {('long_name', 'Global temperature process uncertainty due to SO4'),
                 ('standard_name', 'Arctic_temperature_uncertainty_SO4'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  ds8 = xr.Dataset({'dtemp_GLOBAL_SO4_pos': (['time'], temp_uncert_pos_so4)}, coords={'time':time_mod})
  ds8.dtemp_GLOBAL_SO4_pos.attrs = {('long_name', 'Global temperature process uncertainty due to SO4'),
                 ('standard_name', 'Arctic_temperature_uncertainty_SO4'),
                 ('units', 'K'),
                 ('_FillValue', fillv)}
  #outfile = './dtemp_GLOBAL_'+scenario+'.nc'
  outfile = datadir+'Ensemble/dtemp_GLOBAL_'+scenario+'.nc'
  print(outfile)
ds = xr.merge([ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8])
ds.attrs = {('comment', 'Scenario: '+scenario+'; contact: knut.vonsalzen@canada.ca'),
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
if netcdf == 1:
  ds.to_netcdf(outfile, unlimited_dims=["time"])
