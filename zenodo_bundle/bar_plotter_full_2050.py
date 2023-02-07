"""
Author: Knut von Salzen
Edits by: Victoria Spada
"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date # http://unidata.github.io/netcdf4-python/
import re
import sys


# Routines for reading netcdf data

def get_coord_size(varn,ncfile):  # Get size of the coordinate array
  coord = ncfile.variables[varn][:]
  ncoord = len(coord)
  return ncoord

def get_coord_attr(varn,ncoord,ncfile):  # get names associated with coordinate indices
  coord_attr = ncfile.variables[varn].long_name
  coord_attspl = re.split('\d - ',coord_attr)                  # Split long_name using integers as separator
  scoord = [" " for i in range(ncoord)]
  for n in range(ncoord):
    npt = n + 1
    if npt < ncoord:                                           # Remove redundant characters at the end of the string
      coord_attspl_len = len(coord_attspl[npt]) - 2
    else:
      coord_attspl_len = len(coord_attspl[npt])
    form = "".join(["{:.","{}".format(coord_attspl_len),"}"])  # Determine format of the string
    scoord[n] = form.format(coord_attspl[npt])                 # Read the string using the format
    scoord_attspl = re.split('\'',scoord[n])                   # Split long_name using integers as separator
    coord_attspl_len = len(scoord_attspl) 
    if coord_attspl_len > 1:
      scoord[n] = scoord_attspl[0]
  return scoord

clean = False

ylabel = r"Relative to 1995 - 2014 ($^{\circ}\mathrm{C}$)"
ylabelx = r"Relative to 1880 - 1920 ($^{\circ}\mathrm{C}$)"
title = ['Global (AMAP)', 'Arctic (AMAP)', 'Global (SSP)', 'Arctic (SSP)']
scenario = [ 'SSP119', 'SSP126', 'SSP245', 'SSP370', 'SSP585', 'CLE', 'MFR', 'CFM', 'MFR_SDS', 'extra']
scenariox = [  'SSP1-1.9', 'SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5', 'CLE', 'MFR', 'CFM', 'MFR_SDS', 'extra']
cmip6 = [ True, True, True, True, True, False, False, False, False, False ]

refyrs = '1995_2014'
year_start = 1990 #2026 #1990
year_end = 2050 #2035 #2050
ptimei = [ year_end ]

ici = [ 0, 2, 1, 3 ]
cmip6_plt = [ False, False, True, True ]
cmip6_type = 1 # 0 - mean, 1 - median

fillv = -9999.
nszen = len(scenario) # Number of Scenarioes
nplt = nszen
nplt_cmip6 = 0
nplt_amap = 0
for nz in range(nszen):
  if cmip6[nz]:
    nplt_cmip6 = nplt_cmip6 + 1
  else:
    nplt_amap = nplt_amap + 1
    
#nplt_cmip6 += 1

casei = [ 4, 0 ] # 0 - Arctic, 4 - global
ic = 0
models = [ 'emulator' ]
cases = ['ARCTIC', 'GLOBAL']

#rmarship = False  # remove Arctic BC shipping emissions
rmarship = True  # remove Arctic BC shipping emissions
arship = [0, 5, 3]  # Arctic BC shipping emissions: species, region, sector

indxbc = 0  # BC emissions index
indxs = 1  # Sulfur emissions index
indxoc = 2  # OC emissions index

prange = True
pextra = True
# prange = False
# pextra = False

ylarge = 10.
nsspy = 5#6 # Exclude ozone
region_values_pos = np.zeros((4,nplt))
region_splits_pos = np.zeros((4,nplt,nsspy))
region_values_neg = np.zeros((4,nplt))
region_splits_neg = np.zeros((4,nplt,nsspy))
region_values_net = np.zeros((4,nplt))
region_values_min1 = np.zeros((4,nplt))
region_values_max1 = np.zeros((4,nplt))
region_values_min2 = np.zeros((4,nplt))
region_values_max2 = np.zeros((4,nplt))
region_values_netx = np.zeros((4,nplt))
region_values_minx = np.zeros((4,nplt))
region_values_maxx = np.zeros((4,nplt))
region_values_nety = np.zeros((4,nplt))
region_values_miny = np.zeros((4,nplt))
region_values_maxy = np.zeros((4,nplt))


ymin = np.zeros((4))
ymax = np.zeros((4))
temp_off = np.zeros((4))
yticks = np.full((4,12),ylarge)
yticksx = np.full((4,12),ylarge)

partial = False
legend = True
titleinside = False
tightl = True
tightlx = False

esm_data_dir = 'FW__Decision,__Science_Advances_/'
amap_data_dir = 'AMAP Projections/Renormed/'
amap_data_dir_oz = 'AMAP Projections Ozone/'
amap_data_dir_ensemble = 'AMAP Projections/Renormed Ensemble/'

for cs in range(2): # For Arctic and global cases
  case = casei[cs]
  yrsesm = np.arange(year_start, year_start, 1)
  print(cs, case)
  if case==0: # If Arctic scenario
    # Open CMIP6 files
    print('## ARCTIC CASE ##')
    esmdat_median_arc = esm_data_dir+'temp_CMIP6_Arctic_median.csv'
    esmdat_min_arc = esm_data_dir+'temp_CMIP6_Arctic_low.csv'
    esmdat_max_arc = esm_data_dir+'temp_CMIP6_Arctic_high.csv'
    
    esmdat_median_arc = esm_data_dir+'temp_CMIP6_Arctic_median_best.csv'
    esmdat_min_arc = esm_data_dir+'temp_CMIP6_Arctic_low_best.csv'
    esmdat_max_arc = esm_data_dir+'temp_CMIP6_Arctic_high_best.csv'    
    
    data_esm_min_arc = np.genfromtxt(esmdat_min_arc, delimiter=",")[:,1]
    data_esm_median_arc = np.genfromtxt(esmdat_median_arc, delimiter=",")[:,1]
    data_esm_max_arc = np.genfromtxt(esmdat_max_arc, delimiter=",")[:,1]
    
    print('data_esm_min_arc',data_esm_min_arc)
    print('data_esm_median_arc',data_esm_median_arc)
    print('data_esm_max_arc',data_esm_max_arc)
    
  else: # If Global scenario 
    # Open CMIP6 files
    print('## GLOBAL CASE ##')
    esmdat_median_glbl = esm_data_dir+'temp_CMIP6_global_median.csv'
    esmdat_min_glbl = esm_data_dir+'temp_CMIP6_global_low.csv'
    esmdat_max_glbl = esm_data_dir+'temp_CMIP6_global_high.csv'

    esmdat_median_glbl = esm_data_dir+'temp_CMIP6_global_median_best.csv'
    esmdat_min_glbl = esm_data_dir+'temp_CMIP6_global_low_best.csv'
    esmdat_max_glbl = esm_data_dir+'temp_CMIP6_global_high_best.csv'
    
    data_esm_min_glbl = np.genfromtxt(esmdat_min_glbl, delimiter=",")[:,1]
    data_esm_median_glbl = np.genfromtxt(esmdat_median_glbl, delimiter=",")[:,1]
    data_esm_max_glbl = np.genfromtxt(esmdat_max_glbl, delimiter=",")[:,1]

    print('data_esm_min_glbl',data_esm_min_glbl)
    print('data_esm_median_glbl',data_esm_median_glbl)
    print('data_esm_max_glbl',data_esm_max_glbl)
    
  # data_esm_min = np.genfromtxt(esmdat_min)
  # data_esm_median = np.genfromtxt(esmdat_median)
  # data_esm_max = np.genfromtxt(esmdat_max)
  
  # print(data_esm_min)
  # print(data_esm_median)
  # print(data_esm_max)
   
   # Three columns in each AMAP file
   # data_esm_median = np.genfromtxt(esmdat_median, dtype='f,f,f', names=True, delimiter=',')
   # data_esm_min = np.genfromtxt(esmdat_min, dtype='f,f,f', names=True, delimiter=',')
   # data_esm_max = np.genfromtxt(esmdat_max, dtype='f,f,f', names=True, delimiter=',')
   # yrsesm = data_esm_median['Year']
   # temp_esm_median=np.zeros((nszen,len(yrsesm)))
   # temp_esm_min=np.zeros((nszen,len(yrsesm)))
   # temp_esm_max=np.zeros((nszen,len(yrsesm)))

  if case!=0: # GLOBAL
    esmdat_median = esm_data_dir+'temp_CMIP6_global_median.csv'
    esmdat_min = esm_data_dir+'temp_CMIP6_global_low.csv'
    esmdat_max = esm_data_dir+'temp_CMIP6_global_high.csv'
  else: # ARCTIC
    esmdat_median = esm_data_dir+'temp_CMIP6_Arctic_median.csv'
    esmdat_min = esm_data_dir+'temp_CMIP6_Arctic_low.csv'
    esmdat_max = esm_data_dir+'temp_CMIP6_Arctic_high.csv'
  # # Six rows in each CMIP csv file
  # data_esm_cmip6_median = np.genfromtxt(esmdat_median, dtype='f,f,f,f,f,f', names=True, delimiter=',')
  # data_esm_cmip6_min = np.genfromtxt(esmdat_min, dtype='f,f,f,f,f,f', names=True, delimiter=',')
  # data_esm_cmip6_max = np.genfromtxt(esmdat_max, dtype='f,f,f,f,f,f', names=True, delimiter=',')

  for nz in range(nszen-1): # For each scenario
    # if not cmip6[nz]: # and scenario[nz] != 'MFR_SDS': # If not an SSP pathway & not MFR_SDS
    #   #if scenario[nz] != 'CFM' : 
          
    #   # Each column of the temp_esm_median is the median temp 
    #   temp_esm_median[nz,:] = data_esm_median[scenario[nz]]
    #   temp_esm_min[nz,:] = data_esm_min[scenario[nz]]
    #   temp_esm_max[nz,:] = data_esm_max[scenario[nz]]

    # if cmip6[nz]: # If a CMIP6 scenario (SSP Pathways)
    #    temp_esm_median[nz,:] = data_esm_cmip6_median[scenario[nz]]
    #    temp_esm_min[nz,:] = data_esm_cmip6_min[scenario[nz]]
    #    temp_esm_max[nz,:] = data_esm_cmip6_max[scenario[nz]]

# Simulated temperature change due to CO2 (take MEAN ECS output)
    modconc = amap_data_dir+'temp_co2_'+scenario[nz]+'_zero_emulator_mean_'+refyrs+'.nc'
    ncfile = Dataset(modconc,"r",format="NETCDF4")
    str_year = ncfile.variables["time"][:]
    ntime = len(str_year)
    nreg     = get_coord_size('region',ncfile)
    sregion  = get_coord_attr("region",nreg,ncfile)
    sector = ncfile.variables["sector"][:]
    nsec = len(sector)
    if nz==0:
      temp_mod_in_co2 = np.zeros((nszen,ntime,nreg,nsec))
    if case == 0:
      temp_mod_in_co2[nz,:,:] = ncfile.variables["temp_ARCTIC"][:]
    elif case == 4:
      temp_mod_in_co2[nz,:,:] = ncfile.variables["temp_glb"][:]
    fillvin = ncfile.variables["temp_glb"]._FillValue
    mask = np.array(temp_mod_in_co2[nz,:,:])
    temp_mod_in_co2[nz,:,:] = np.where(mask==fillvin, 0., mask).copy()

# Simulated temperature change due to CH4
  nmod = len(models)
  ch4_flag = True
  if ch4_flag:
    temp_mod_in_ch4 = np.zeros((nszen,ntime,nreg,nsec))
    for nz in range(nszen-1):
      modconc = amap_data_dir+'temp_ch4_'+scenario[nz]+'_zero_emulator_mean_'+refyrs+'.nc'
      ncfile = Dataset(modconc,"r",format="NETCDF4")
      nfsp     = get_coord_size('frcspec',ncfile)
      nssp     = get_coord_size('srcspec',ncfile)
      if case == 0:
        temp_mod_in_ch4_tmp = ncfile.variables["temp_ARCTIC"][:]
      elif case == 4:
        temp_mod_in_ch4_tmp = ncfile.variables["temp_glb"][:]
      fillvin = ncfile.variables["temp_glb"]._FillValue
      for nfs in range(nfsp):
        for nss in range(nssp):
          mask = np.array(temp_mod_in_ch4_tmp[:,:,:,nfs,nss])
          temp_mod_in_ch4_tmp[:,:,:,nfs,nss] = np.where(mask==fillvin, 0., mask).copy()
          temp_mod_in_ch4[nz,:,:,:] = temp_mod_in_ch4[nz,:,:,:] + temp_mod_in_ch4_tmp[:,:,:,nfs,nss] 
        
  ntx = 0
  for nt in range(ntime):
    if (int(str_year[nt]) >= year_start) and (int(str_year[nt]) <= year_end):
      ntx = ntx + 1
  ntimem = ntx

# Simulated temperature change due to ozone
  # nfirst = True
  # for nz in range(nszen):
  #   temp_mod_oz_plt = np.zeros((nszen,ntimem))
  #   if not cmip6[nz]:
  #     # modconc = amap_data_dir_oz+'temp_ozone_'+scenario[nz]+'_1990_2050_1995_2014.nc'
  #     modconc = amap_data_dir_oz+'temp_o3_'+scenario[nz]+'_emulator_mean.nc'
  #     ncfile = Dataset(modconc,"r",format="NETCDF4")
  #     timeoz = ncfile.variables['time'][:]
  #     ntimeoz = len(timeoz)
  #     if nfirst:
  #       temp_mod_in_oz = np.zeros((nszen,ntimeoz))
  #       nfirst = False
  #     if case == 0:
  #       temp_mod_in = ncfile.variables['temp_ARCTIC'][:]
  #     elif case == 4:
  #       temp_mod_in = ncfile.variables['temp_glb'][:]
  #     fillvin = ncfile.variables["temp_ARCTIC"]._FillValue
  #     mask = np.array(temp_mod_in)
  #     temp_mod_in_oz[nz,:] = np.where(mask==fillvin, 0., mask).copy()

  #     temp_mod_oz = np.full((nszen,ntimem), 0.)
  #     ntx = 0
  #     for nt in range(ntimeoz):
  #       if (int(timeoz[nt]) >= year_start) and (int(timeoz[nt]) <= year_end):
  #         for nz in range(nszen):
  #           temp_mod_oz[nz,ntx] = temp_mod_in_oz[nz,nt]
  #         ntx = ntx + 1
  #     temp_mod_oz_plt = temp_mod_oz

# Simulated temperature change due to aerosol

  temp_mod_in_aer = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_aerl = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_aerh = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_aer_s = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_aer_bc = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_aer_oc = np.zeros((nszen,nmod,ntime,nreg,nsec))
  time_mod = np.full((nmod,ntime),False)
  process_sum = [2, 3, 4]
  bcdiffall = True
  if bcdiffall:
    process_sum_bc = [2, 3, 4]
  else:
    process_sum_bc = [2]
  species_sum = [1, 2, 3]
  # species_sum = [3]
  nsspx = len(species_sum)

  for nm in range(nmod):
    if models[nm] == 'CanAM' or models[nm] == 'MRI-ESM' or models[nm] == 'emulator':
      nsec1 = nsec
      nreg1 = nreg
      for nz in range(nszen-1):
        modconc = amap_data_dir+'temp_aerosol_'+scenario[nz]+'_1750_'+models[nm]+'_mean_'+refyrs+'.nc'
        # modconc = amap_data_dir+'temp_aerosol_'+scenario[nz]+'_1750_'+models[nm]+'_mean_tst_'+refyrs+'.nc'
        ncfile = Dataset(modconc,"r",format="NETCDF4")
        nfsp     = get_coord_size('frcspec',ncfile)
        nssp     = get_coord_size('srcspec',ncfile)
        nproc    = get_coord_size('process',ncfile)
        process  = ncfile.variables['process'][:]
        species  = ncfile.variables['srcspec'][:]
        nreg     = get_coord_size('region',ncfile)
        nsec     = get_coord_size('sector',ncfile)
        timei    = ncfile.variables['time'][:]
        ntimei   = len(timei)
        if case == 0:
          temp_mod_in_aer_tmp = ncfile.variables["temp_ARCTIC"][:]
        elif case == 4:
          temp_mod_in_aer_tmp = ncfile.variables["temp_glb"][:]
        fillvin = ncfile.variables["temp_glb"]._FillValue
        mask = np.array(temp_mod_in_aer_tmp)
        temp_mod_in_aer_tmp = np.where(mask==fillvin, 0., mask).copy()
        if rmarship:
          temp_mod_in_aer_tmp[:,arship[1],arship[2],arship[0],arship[0],:] = 0.
        for npt in range(len(process_sum)):
          for npr in range(nproc):
            if process[npr] == process_sum[npt]:
              temp_mod_in_aer_s[nz,nm,:,:,:] = temp_mod_in_aer_s[nz,nm,:,:,:] + temp_mod_in_aer_tmp[:,:,:,indxs,indxs,npr]
              temp_mod_in_aer_bc[nz,nm,:,:,:] = temp_mod_in_aer_bc[nz,nm,:,:,:] + temp_mod_in_aer_tmp[:,:,:,indxbc,indxbc,npr]
              temp_mod_in_aer_oc[nz,nm,:,:,:] = temp_mod_in_aer_oc[nz,nm,:,:,:] + temp_mod_in_aer_tmp[:,:,:,indxoc,indxoc,npr]
              for nfs in range(nfsp):
                for nss in range(nssp):
                  for nsx in range(nsspx):
                    if species_sum[nsx] == species[nss] and nss == nfs: 
                      temp_mod_in_aer[nz,nm,:,:,:] = temp_mod_in_aer[nz,nm,:,:,:] + temp_mod_in_aer_tmp[:,:,:,nfs,nss,npr]
      for nt in range(ntime):
        for ntx in range(ntimei):
          if timei[ntx] == str_year[nt]:
            time_mod[nm,nt] = True

# Total simulated temperature change

  temp_mod_in = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_s = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_bc = np.zeros((nszen,nmod,ntime,nreg,nsec))
  temp_mod_in_oc = np.zeros((nszen,nmod,ntime,nreg,nsec))
  for nz in range(nszen):
    for nm in range(nmod):
      # if models[nm] == 'CanAM' or models[nm] == 'MRI-ESM' or models[nm] == 'emulator' and not cmip6[nz]:
      #   temp_mod_in[nz,nm,0,0,0] = temp_mod_in_oz[nz,nm]
      temp_mod_in[nz,nm,:,:,:] = temp_mod_in[nz,nm,:,:,:] + temp_mod_in_aer[nz,nm,:,:,:]
      temp_mod_in_s[nz,nm,:,:,:] = temp_mod_in_aer_s[nz,nm,:,:,:]
      temp_mod_in_bc[nz,nm,:,:,:] = temp_mod_in_aer_bc[nz,nm,:,:,:]
      temp_mod_in_oc[nz,nm,:,:,:] = temp_mod_in_aer_oc[nz,nm,:,:,:]
      if models[nm] == 'CanAM' or models[nm] == 'MRI-ESM' or models[nm] == 'emulator':
        temp_mod_in[nz,nm,:,:,:] = temp_mod_in[nz,nm,:,:,:] + temp_mod_in_co2[nz,:,:,:] + temp_mod_in_ch4[nz,:,:,:]

  temp_mod_tmp = np.full((nszen,nmod,ntime),0.)
  temp_mod_s_tmp = np.full((nszen,nmod,ntime),0.)
  temp_mod_bc_tmp = np.full((nszen,nmod,ntime),0.)
  temp_mod_oc_tmp = np.full((nszen,nmod,ntime),0.)
  temp_mod_co2_tmp = np.full((nszen,nmod,ntime),0.)
  temp_mod_ch4_tmp = np.full((nszen,nmod,ntime),0.)
  for nz in range(nszen):
    for nm in range(nmod):
      for nr in range(nreg):
        for ns in range(nsec):
          temp_mod_tmp[nz,nm,:] = temp_mod_tmp[nz,nm,:] + temp_mod_in[nz,nm,:,nr,ns]
          temp_mod_s_tmp[nz,nm,:] = temp_mod_s_tmp[nz,nm,:] + temp_mod_in_s[nz,nm,:,nr,ns]
          temp_mod_bc_tmp[nz,nm,:] = temp_mod_bc_tmp[nz,nm,:] + temp_mod_in_bc[nz,nm,:,nr,ns]
          temp_mod_oc_tmp[nz,nm,:] = temp_mod_oc_tmp[nz,nm,:] + temp_mod_in_oc[nz,nm,:,nr,ns]

  for nz in range(nszen):
    for nr in range(nreg):
      for ns in range(nsec):
        temp_mod_co2_tmp[nz,0,:] = temp_mod_co2_tmp[nz,0,:] + temp_mod_in_co2[nz,:,nr,ns]
        temp_mod_ch4_tmp[nz,0,:] = temp_mod_ch4_tmp[nz,0,:] + temp_mod_in_ch4[nz,:,nr,ns]


# Select time period for plotting and adjust model results to match observations 
# in reference year

  temp_mod = np.full((nszen,nmod,ntimem), fillv)
  temp_mod_co2 = np.full((nszen,nmod,ntimem), fillv)
  temp_mod_ch4 = np.full((nszen,nmod,ntimem), fillv)
  temp_mod_s = np.full((nszen,nmod,ntimem), fillv)
  temp_mod_bc = np.full((nszen,nmod,ntimem), fillv)
  temp_mod_oc = np.full((nszen,nmod,ntimem), fillv)
  temp_mod_time = np.full((ntimem), fillv)
  ntx = 0
  for nt in range(ntime):
    if (int(str_year[nt]) >= year_start) and (int(str_year[nt]) <= year_end):
      for nz in range(nszen):
        temp_mod_co2[nz,:,ntx] = temp_mod_co2_tmp[nz,:,nt]
        temp_mod_ch4[nz,:,ntx] = temp_mod_ch4_tmp[nz,:,nt]
        temp_mod_s[nz,:,ntx] = temp_mod_s_tmp[nz,:,nt]
        temp_mod_bc[nz,:,ntx] = temp_mod_bc_tmp[nz,:,nt]
        temp_mod_oc[nz,:,ntx] = temp_mod_oc_tmp[nz,:,nt]
        temp_mod[nz,:,ntx] = temp_mod_tmp[nz,:,nt]
      temp_mod_time[ntx] = float(str_year[nt])
      ntx = ntx + 1

# Simulated temperature change and its range for multi-model ensemble

  mean = False
  temp_mod_plt = np.zeros((nszen,ntimem))
  temp_mod_pltl = np.zeros((nszen,ntimem))
  temp_mod_plth = np.zeros((nszen,ntimem))
  # temp_mod_all_plt = np.zeros((nszen,ntimem))
  temp_mod_spec_plt = np.zeros((nszen,ntimem))
  temp_mod_co2_plt = np.zeros((nszen,ntimem))
  temp_mod_ch4_plt = np.zeros((nszen,ntimem))
  temp_mod_s_plt = np.zeros((nszen,ntimem))
  temp_mod_bc_plt = np.zeros((nszen,ntimem))
  temp_mod_oc_plt = np.zeros((nszen,ntimem))
  temp_mod_plt = temp_mod[:,0,:]
  temp_mod_co2_plt = temp_mod_co2[:,0,:]
  temp_mod_ch4_plt = temp_mod_ch4[:,0,:]
  temp_mod_s_plt = temp_mod_s[:,0,:]
  temp_mod_bc_plt = temp_mod_bc[:,0,:]
  temp_mod_oc_plt = temp_mod_oc[:,0,:]

# Uncertainty ranges
  for nz in range(nszen-1):
    if case == 0:
      modconc = amap_data_dir_ensemble+'dtemp_ARCTIC_'+scenario[nz]+'.nc'
    else:
      modconc = amap_data_dir_ensemble+'dtemp_GLOBAL_'+scenario[nz]+'.nc'
    ncfile = Dataset(modconc,"r",format="NETCDF4")
    timex    = ncfile.variables['time'][:]
    ntimex   = len(timex)
    if nz == 0:
      dtemp_pos = np.zeros((nszen,ntimex))
      dtemp_neg = np.zeros((nszen,ntimex))
      dtemp_SLCF_pos = np.zeros((nszen,ntimex))
      dtemp_SLCF_neg = np.zeros((nszen,ntimex))
    if case == 0:
      dtemp_pos[nz,:] = ncfile.variables['dtemp_ARCTIC_pos'][:]
      dtemp_neg[nz,:] = ncfile.variables['dtemp_ARCTIC_neg'][:]
      dtemp_SLCF_pos[nz,:] = ncfile.variables['dtemp_ARCTIC_SLCF_pos'][:]
      dtemp_SLCF_neg[nz,:] = ncfile.variables['dtemp_ARCTIC_SLCF_neg'][:]
    else:
      dtemp_pos[nz,:] = ncfile.variables['dtemp_GLOBAL_pos'][:]
      dtemp_neg[nz,:] = ncfile.variables['dtemp_GLOBAL_neg'][:]
      dtemp_SLCF_pos[nz,:] = ncfile.variables['dtemp_GLOBAL_SLCF_pos'][:]
      dtemp_SLCF_neg[nz,:] = ncfile.variables['dtemp_GLOBAL_SLCF_neg'][:]

# Total   
  for nt in range(ntimem):
    temp_mod_spec_plt[:,nt] = temp_mod_co2_plt[:,nt] + temp_mod_ch4_plt[:,nt] + temp_mod_s_plt[:,nt] + temp_mod_bc_plt[:,nt] + temp_mod_oc_plt[:,nt] #+ temp_mod_oz_plt[:,nt]

# PRINT RESULTS
  # for nz in range(nszen): # For each of the AMAP and CMIP scenarios
  #   print('\n', cases[cs], scenario[nz])
  #   for nt in range(ntimem):
  #     if temp_mod_time[nt] == year_end: #2030 or temp_mod_time[nt] == 2050:
  #       print(str(temp_mod_time[nt])+' '+str(temp_mod_plt[nz,nt]))
  # for nz in range(nszen):
  #   for nt in range(ntimem):
  #     if temp_mod_time[nt] == year_end: #2030 or temp_mod_time[nt] == 2050:
  #       print(str(temp_mod_time[nt])+' '+str(temp_mod_plt[nz,nt])+' '+str(temp_mod_co2_plt[nz,nt])+' '+str(temp_mod_ch4_plt[nz,nt])+' '+str(temp_mod_s_plt[nz,nt])+' '+str(temp_mod_bc_plt[nz,nt])+' '+str(temp_mod_oc_plt[nz,nt])+' ')#+str(temp_mod_all_plt[nz,nt]))

  pt = 0
  ptime = ptimei[pt] # ptime is start year
  for ptx in range(2):
    ict = ici[ic]
    if case==4:
      temp_off0 = 0.81867805  # global temperature in 2015, relative to 1951-1980
      temp_off1 = 0.57535783  # global temperature in 1995 to 2014, relative to 1951-1980
      temp_off2 = 0.24855012  # global temperature anomaly 1951-1980, relative to 1880-1920
    else:
      temp_off0 = 2.10117541  # arctic temperature in 2015, relative to 1951-1980
      temp_off1 = 1.34416833  # arctic temperature in 1995 to 2014, relative to 1951-1980
      temp_off2 = 0.56775893  # arctic temperature anomaly 1951-1980, relative to 1880-1920
#    temp0 = temp_off0 - temp_off1  # temperature relative to 1995-2014, for reference year 2015
#    temp_off[ict] = temp0 + temp_off1 + temp_off2   # temperature relative to 1880-1920, for reference year 1995-2014
    temp_off[ict] = temp_off1 + temp_off2   # temperature relative to 1880-1920, for reference year 1995-2014

#    yticks[ict,:] = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]
#    yticksx[ict,:] = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]
    if case == 4: 
      ymin[ict] = -.25
      ymax[ict] = 2.5
      yticks[ict,:8] = [-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
      yticksx[ict,:7] = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]
    else:
      ymin[ict] = -1.0
      if not clean:
        ymax[ict] = 5.0
        yticks[ict,:8] = [-1, 0, 1, 2, 3, 4, 5, 6]
        yticksx[ict,:8] = [1, 2, 3, 4, 5, 6, 7, 8]
      else:
        ymax[ict] = 3.8
        yticks[ict,:8] = [-4, -3, -2, -1, 0, 1, 2, 3]
        yticksx[ict,:8] = [-3, -2, -1, 1, 2, 3, 4, 5]

    if case == 0:
      if cmip6_plt[ict]:
        title[ict] = "Arctic (SSP)"
      else:
        title[ict] = "Arctic (AMAP)"
    else:
      if cmip6_plt[ict]:
        title[ict] = "Global (SSP)"
      else:
        title[ict] = "Global (AMAP)"
    print('+++ '+title[ict])
    ntplt = -9
    for nt in range(ntimem):
      if temp_mod_time[nt] == ptime: # If we have reached the start year 
        ntplt = nt
    ntpltx = -9
    for nt in range(ntimex):
      if timex[nt] == ptime:
        ntpltx = nt
    ntesm = -9
    for nt in range(len(yrsesm)):
      if int(yrsesm[nt]) == ptime:
        ntesm = nt

# Select emulator temperatures
    pltval = np.zeros((nszen,nsspy))
    pltval_net = np.zeros((nszen))
    pltval_min1 = np.zeros((nszen))
    pltval_max1 = np.zeros((nszen))
    pltval_min2 = np.zeros((nszen))
    pltval_max2 = np.zeros((nszen))
    pltval_netx = np.zeros((nszen))
    pltval_minx = np.zeros((nszen))
    pltval_maxx = np.zeros((nszen))
    pltval_nety = np.zeros((nszen))
    pltval_miny = np.zeros((nszen))
    pltval_maxy = np.zeros((nszen))
    for ip in range(nplt):
      pltval[ip,0] = temp_mod_co2_plt[ip,ntplt]
      pltval[ip,1] = temp_mod_bc_plt[ip,ntplt]
      pltval[ip,2] = temp_mod_s_plt[ip,ntplt]
      pltval[ip,3] = temp_mod_oc_plt[ip,ntplt]
      pltval[ip,4] = temp_mod_ch4_plt[ip,ntplt]
      # pltval[ip,5] = temp_mod_oz_plt[ip,ntplt]
    if prange:
      pltval_net[:] = temp_mod_plt[:,ntplt]
      pltval_min1[:] = pltval_net[:] - dtemp_SLCF_neg[:,ntpltx]
      pltval_max1[:] = pltval_net[:] + dtemp_SLCF_pos[:,ntpltx]
      pltval_min2[:] = pltval_net[:] - dtemp_neg[:,ntpltx]
      pltval_max2[:] = pltval_net[:] + dtemp_pos[:,ntpltx]
    # if pextra:
    #   pltval_netx[:] = temp_esm_median[:,ntesm]
    #   pltval_minx[:] = temp_esm_min[:,ntesm]
    #   pltval_maxx[:] = temp_esm_max[:,ntesm]
    # for ip in range(nplt):
    #   print(str(ip)+' '+scenario[ip]+' net T = '+str(pltval_net[ip])+', range = -'+str(dtemp_SLCF_neg[ip,ntpltx])+' to +'+str(dtemp_SLCF_pos[ip,ntpltx]))


# Convert to plotting input
    for ip in range(nplt):
      for isp in range(nsspy):
        if pltval[ip,isp] >= 0.:
          region_splits_pos[ict,ip,isp] = pltval[ip,isp]
        else:
          region_splits_neg[ict,ip,isp] = pltval[ip,isp]

    for isp in range(nsspy):
      region_values_pos[ict,:] = np.add(region_values_pos[ict,:], region_splits_pos[ict,:,isp])
      region_values_neg[ict,:] = np.add(region_values_neg[ict,:], region_splits_neg[ict,:,isp])
    if prange:
      region_values_net[ict,:] = pltval_net[:]
      region_values_min1[ict,:] = pltval_min1[:]
      region_values_max1[ict,:] = pltval_max1[:]
      region_values_min2[ict,:] = pltval_min2[:]
      region_values_max2[ict,:] = pltval_max2[:]
    else:
      region_values_net[ict,:] = np.add(region_values_net[ict,:], np.add(region_values_pos[ict,:], region_values_neg[ict,:]))
    if pextra: # Assign rows for each variable
      region_values_netx[ict,:] = pltval_netx[:]
      region_values_nety[ict,:] = pltval_nety[:]
      region_values_minx[ict,:] = pltval_minx[:]
      region_values_maxx[ict,:] = pltval_maxx[:]
      region_values_miny[ict,:] = pltval_miny[:]
      region_values_maxy[ict,:] = pltval_maxy[:]
    ic = ic + 1

xlabelsx = np.full((nplt),'        ')
x = np.zeros((nplt))
dx = 1.
for i in range(nplt):
  x[i] = i * dx
  xlabelsx[i] = scenariox[i]
xlabelsx_amap = np.full((nplt_amap),'        ')
xlabelsx_cmip6 = np.full((nplt_cmip6),'        ')
x_amap = np.zeros((nplt_amap))
x_cmip6 = np.zeros((nplt_cmip6))
ica = 0
icc = 0
for i in range(nplt):
  if not cmip6[i]:
    x_amap[ica] = ica * dx
    xlabelsx_amap[ica] = scenariox[i]
    ica = ica + 1
  else:
    x_cmip6[icc] = icc * dx
    xlabelsx_cmip6[icc] = scenariox[i]
    icc = icc + 1
    
# Define labels
species = [ "CO$_{2}$", "BC", "SO$_{4}$", "OC", "CH$_{4}$"] #, "O$_{3}$" ]

spltlabel = np.full((nsspy),'        ')
for isp in range(nsspy):
  spltlabel[isp] = species[isp]

# Define colours
icper = [ 10, 0, 2, 1, 4]#, 19 ]
cmap1 = matplotlib.cm.get_cmap('Accent')
dcol1 = 1./8.
icol1 = dcol1/2
cmap2 = matplotlib.cm.get_cmap('Paired')
dcol2 = 1./12.
icol2 = dcol2/2
spltcolor = np.zeros((nsspy,4))
for ic in range(nsspy):
  icx = icper[ic]
  if icx < 10:
    spltcolor[ic,:] = cmap1(icol1 + icx * dcol1)
  elif icx == -1:
#    spltcolor[ic,:] = [0, 0, 0, 0.2]    # grey
#    spltcolor[ic,:] = [1, 1, 1, 1]    # white
#    spltcolor[ic,:] = [0, 0, 0, 1]    # black
#    spltcolor[ic,:] = [1, 0, 0, 1]    # red
    spltcolor[ic,:] = cmap1(icol1 + 6 * dcol1)
  else:
    spltcolor[ic,:] = cmap2(icol2 + (icx-10) * dcol2)


###############################################################################
def plot_bars(region_values_pos,region_values_neg,region_splits_pos,region_splits_neg,spltcolor,spltlabel,xlabelsx_amap,x_amap,xlabelsx_cmip6,x_cmip6,ylabel,title,ymin,ymax,yticks,yticksx,nplt,nsspy,titleinside,region_values_min1,region_values_min2,region_values_net,region_values_max1,region_values_max2,region_values_netx,region_values_minx,region_values_maxx,region_values_nety,region_values_miny,region_values_maxy,temp_off,pextra,prange,legend,cmip6,cmip6_plt,nplt_amap):
  global2050, arctic2050 = [], []
  g5ae, a5ae = [], []
    
  nfram = 4
  splthgt_pos       = np.zeros((nfram,nplt_amap,nsspy+1))
  splthgt_neg       = np.zeros((nfram,nplt_amap,nsspy+1))

  splthgt_pos[:,:,0] = 0.
  splthgt_neg[:,:,0] = 0.
  nzx = 0
  ina = -9
  iss = np.full((nfram),ina)
  ise = np.full((nfram),ina)
  ftsize_ax = 10
  for nfr in range(nfram):  # For each of the 4 subplots
    print(nfram, '### '+title[nfr]+' ###')
    for nz in range(nplt):
      if (cmip6_plt[nfr] and cmip6[nz]) or (not cmip6_plt[nfr] and not cmip6[nz]):
        if iss[nfr] == ina:
          iss[nfr] = nz
        ise[nfr] = nz
    ib = 0
    for i2 in range(nsspy): # Cycle through each species (CO2, BC, S0)
      # if nfr >= 2:
      #     iss[0], iss[1] = 4,4
      # print(str(iss[nfr])+' '+str(ise[nfr]))
      # for i2x in range(nplt):
      #   if i2x >= iss[nfr] and i2x <= ise[nfr]:
      #     print('pos '+str(i2x)+' '+str(region_splits_pos[nfr,i2x,i2]))
      # for i2x in range(nplt):
      #   if i2x >= iss[nfr] and i2x <= ise[nfr]:
      #     print('neg '+str(i2x)+' '+str(region_splits_neg[nfr,i2x,i2]))

      splthgt_pos[nfr,:,ib+1] = np.add(splthgt_pos[nfr,:,ib], region_splits_pos[nfr,iss[nfr]:ise[nfr]+1,i2]).tolist()
      splthgt_neg[nfr,:,ib+1] = np.add(splthgt_neg[nfr,:,ib], region_splits_neg[nfr,iss[nfr]:ise[nfr]+1,i2]).tolist()
      ib = ib + 1

#  x = np.arange(nplt)  # the label locations
  width2 = 0.7  # the width of the bars

  ylegend1  = 0.1
  xlegend   = 0.
  ylegincr1 = 0.1
  xlegbox   = 0.1
  ylegbox   = 0.05 
  yoff      = -0.003
  xoff      = 0.03

  fig = plt.figure(figsize=(10,7), tight_layout=True)
  gs = gridspec.GridSpec(2, 3)

  axa = fig.add_subplot(gs[0, 0])
  axb = fig.add_subplot(gs[0, 1])
  axc = fig.add_subplot(gs[1, 0])
  axd = fig.add_subplot(gs[1, 1])
  axl = fig.add_subplot(gs[0, 2])
  
  axa2 = axa.twinx()  # instantiate a second axes that shares the same x-axis
  axb2 = axb.twinx()  # instantiate a second axes that shares the same x-axis
  axc2 = axc.twinx()  # instantiate a second axes that shares the same x-axis
  axd2 = axd.twinx()  # instantiate a second axes that shares the same x-axis

  def abst2relt(temp,temp_off):
    return temp - temp_off
  def relt2abst(temp,temp_off):
    return temp + temp_off
  def convert_axa_to_pi(axa):
    y1, y2 = axa.get_ylim()
    axa2.set_ylim(relt2abst(y1,temp_off[0]), relt2abst(y2,temp_off[0]))
    axa2.figure.canvas.draw()
  def convert_axb_to_pi(axb):
    y1, y2 = axb.get_ylim()
    axb2.set_ylim(relt2abst(y1,temp_off[1]), relt2abst(y2,temp_off[1]))
    axb2.figure.canvas.draw()
  def convert_axc_to_pi(axc):
    y1, y2 = axc.get_ylim()
    axc2.set_ylim(relt2abst(y1,temp_off[2]), relt2abst(y2,temp_off[2]))
    axc2.figure.canvas.draw()
  def convert_axd_to_pi(axd):
    y1, y2 = axd.get_ylim()
    axd2.set_ylim(relt2abst(y1,temp_off[3]), relt2abst(y2,temp_off[3]))
    axd2.figure.canvas.draw()

  # automatically update ylim of axa2 when ylim of ax1 changes.
  axa.callbacks.connect("ylim_changed", convert_axa_to_pi)
  axb.callbacks.connect("ylim_changed", convert_axb_to_pi)
  axc.callbacks.connect("ylim_changed", convert_axc_to_pi)
  axd.callbacks.connect("ylim_changed", convert_axd_to_pi)
  axa.set_ylabel(ylabel, fontsize=ftsize_ax)
  axc.set_ylabel(ylabel, fontsize=ftsize_ax)

  xmin = 1.e+20
  xmax = -1.e+20
  for i in range(len(x_amap)):
    xmin = min(xmin,x_amap[i])
    xmax = max(xmax,x_amap[i])
  if not partial:
    if not titleinside:
      axa.set_title(title[0], fontsize=10)
      axb.set_title(title[1], fontsize=10)
      axc.set_title(title[2], fontsize=10)
      axd.set_title(title[3], fontsize=10)
    else:
      xpos = xmin + 0.5 * (xmax - xmin)
      dyza = 0.05*(ymax[0]-ymin[0])
      dyzb = 0.05*(ymax[1]-ymin[1])
      dyzc = 0.05*(ymax[2]-ymin[2])
      dyzd = 0.05*(ymax[3]-ymin[3])
      axa.text(xpos,ymax[0]-dyza,title[0],verticalalignment='center',\
               horizontalalignment='center',color='black',fontsize=10)
      axb.text(xpos,ymax[1]-dyzb,title[1],verticalalignment='center',\
               horizontalalignment='center',color='black',fontsize=10)
      axc.text(xpos,ymax[2]-dyzc,title[2],verticalalignment='center',\
               horizontalalignment='center',color='black',fontsize=10)
      axd.text(xpos,ymax[3]-dyzd,title[3],verticalalignment='center',\
               horizontalalignment='center',color='black',fontsize=10)
  ftsize=7
  for i1 in range(2):
    for i2 in range(3):
      axa2.set_yticks(yticksx[0,:])
      axb2.set_yticks(yticksx[1,:])
      axc2.set_yticks(yticksx[2,:])
      axd2.set_yticks(yticksx[3,:])
      axb2.set_ylabel(ylabelx, color='grey', fontsize=ftsize_ax)  # we already handled the x-label with ax1
      axd2.set_ylabel(ylabelx, color='grey', fontsize=ftsize_ax)  # we already handled the x-label with ax1
      axa.tick_params(axis='y', labelsize=ftsize )
      axb.tick_params(axis='y', labelsize=ftsize )
      axc.tick_params(axis='y', labelsize=ftsize )
      axd.tick_params(axis='y', labelsize=ftsize )
      axa2.tick_params(axis='y', labelcolor='grey', labelsize=ftsize)
      axb2.tick_params(axis='y', labelcolor='grey', labelsize=ftsize)
      axc2.tick_params(axis='y', labelcolor='grey', labelsize=ftsize)
      axd2.tick_params(axis='y', labelcolor='grey', labelsize=ftsize)
      axa2.tick_params(colors='grey', which='both')
      axb2.tick_params(colors='grey', which='both')
      axc2.tick_params(colors='grey', which='both')
      axd2.tick_params(colors='grey', which='both')

  xtol = 0.8

  hlinex = [ xmin-xtol, xmax+xtol ]
  hliney = [ abst2relt(1.5,temp_off[0]), abst2relt(1.5,temp_off[0]) ] 
  rect = axa.plot(hlinex, hliney, linewidth=1, color='grey', \
                  linestyle='--', zorder=0)
  hliney = [ abst2relt(1.5,temp_off[2]), abst2relt(1.5,temp_off[2]) ] 
  rect = axc.plot(hlinex, hliney, linewidth=1, color='grey', \
                  linestyle='--', zorder=0)
  hlinex = [ xmin-xtol, xmax+xtol ]
  hliney = [ abst2relt(2.,temp_off[0]), abst2relt(2.,temp_off[0]) ] 
  rect = axa.plot(hlinex, hliney, linewidth=1, color='grey', \
                  linestyle='--', zorder=0)
  hliney = [ abst2relt(2.,temp_off[2]), abst2relt(2.,temp_off[2]) ] 
  rect = axc.plot(hlinex, hliney, linewidth=1, color='grey', \
                  linestyle='--', zorder=0)

  ib = 0
  for i2 in range(nsspy):
    # Bars corresponding to species contributions
    # print(region_splits_pos[0,iss[0]:ise[0]+1,i2].size)
    # print(splthgt_pos[0,:,ib].size)
    rectsa = axa.bar(x_amap, region_splits_pos[0,iss[0]:ise[0]+1,i2] , width2, bottom=splthgt_pos[0,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsa = axa.bar(x_amap, region_splits_neg[0,iss[0]:ise[0]+1,i2] , width2, bottom=splthgt_neg[0,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsb = axb.bar(x_amap, region_splits_pos[1,iss[1]:ise[1]+1,i2] , width2, bottom=splthgt_pos[1,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsb = axb.bar(x_amap, region_splits_neg[1,iss[1]:ise[1]+1,i2] , width2, bottom=splthgt_neg[1,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsc = axc.bar(x_amap, region_splits_pos[2,iss[2]:ise[2]+1,i2] , width2, bottom=splthgt_pos[2,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsc = axc.bar(x_amap, region_splits_neg[2,iss[2]:ise[2]+1,i2] , width2, bottom=splthgt_neg[2,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsd = axd.bar(x_amap, region_splits_pos[3,iss[3]:ise[3]+1,i2] , width2, bottom=splthgt_pos[3,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    rectsd = axd.bar(x_amap, region_splits_neg[3,iss[3]:ise[3]+1,i2] , width2, bottom=splthgt_neg[3,:,ib], \
                     color=spltcolor[i2,:], label=spltlabel[i2], edgecolor='black', linewidth=0)
    ib = ib + 1
  # Create a rectangle patch and text for the legend
  for i2 in range(nsspy):
    if legend:
      if ymin[0] * ymax[0] < 0 and ymin[0] + ymax[0] < 0:
        i2x = nsspy - i2 - 1
      else:
        i2x = i2
      ylegend = ylegend1 + ylegincr1 * i2
      # axl.scatter(0, 0.6, c='k', edgecolors='black', s=40, marker='D',linewidths=0)
      # axl.scatter(0, 0.65, c='k', edgecolors='black', s=40, marker='.',linewidths=0)
      rect = patches.Rectangle((xlegend,ylegend),xlegbox,ylegbox,facecolor=spltcolor[i2x,:], edgecolor='black', linewidth=0)
      axl.add_patch(rect)
      xpos = xlegend + xlegbox + xoff
      ypos = ylegend + ylegbox/2 + yoff
      axl.text(xpos,ypos,spltlabel[i2x],verticalalignment='center',horizontalalignment='left',color='black',fontsize=ftsize_ax)

  # Add some text for labels, title and custom x-axis tick labels, etc.

  if prange:
    if pextra:
      dxb = 0.2
    else:
      dxb = 0.
#    dbx = 0.1
    dbx = 0.12
  if clean:
    dxb = 0.
  if prange:
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if cmip6[nr]:
        if cmip6_plt[0]:
          rects3 = axa.scatter(x[nra]-dxb, region_values_net[0,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nra = nra + 1
        if cmip6_plt[1]:
          rects3 = axb.scatter(x[nrb]-dxb, region_values_net[1,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nrb = nrb + 1
        if cmip6_plt[2]:
          print('not cmip6 plot [2]', region_values_net[2,nr])
          global2050 += [region_values_net[2,nr]]
          rects3 = axc.scatter(x[nrc]-dxb, region_values_net[2,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nrc = nrc + 1
        if cmip6_plt[3]:
          print('not cmip6 plot [3]', region_values_net[3,nr])
          arctic2050 += [region_values_net[3,nr]]
          rects3 = axd.scatter(x[nrd]-dxb, region_values_net[3,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nrd = nrd + 1
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if not cmip6[nr]:
        if not cmip6_plt[0]:
          rects3 = axa.scatter(x[nra]-dxb, region_values_net[0,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nra = nra + 1
        if not cmip6_plt[1]:
          rects3 = axb.scatter(x[nrb]-dxb, region_values_net[1,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nrb = nrb + 1
        if not cmip6_plt[2]:
          print('not cmip6 plot [2]', region_values_net[2,nr])
          rects3 = axc.scatter(x[nrc]-dxb, region_values_net[2,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nrc = nrc + 1
        if not cmip6_plt[3]:
          print('not cmip6 plot [3]', region_values_net[3,nr])
          rects3 = axd.scatter(x[nrd]-dxb, region_values_net[3,nr], c='black', \
                                 edgecolors='black', s=60, zorder=9999, marker='.',linewidths=0)
          nrd = nrd + 1
  if pextra and not clean:
    dxb1 = 0.5 * dxb
    dxb2 = 1.5 * dxb
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if cmip6[nr]: # Check for CMIP6 Scenario
        # if cmip6_plt[0]:
        #   print('cmip6 plot [0]')
        #   rects3 = axa.scatter(x[nra]+dxb, data_esm_median[0], c='black', \
        #                          edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
        #   nra = nra + 1
        # if cmip6_plt[1]:
        #   print('cmip6 plot [1]')
        #   rects3 = axb.scatter(x[nrb]+dxb, data_esm_median[1], c='black', \
        #                          edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
        #   nrb = nrb + 1
        if cmip6_plt[2]:
          print('cmip6 plot [2]', data_esm_min_glbl[nr-5], data_esm_median_glbl[nr-5], data_esm_max_glbl[nr-5])
          rects3 = axc.scatter(x[nrc]+dxb, data_esm_median_glbl[nr-5], c='k', \
                                 edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
          # add error bars
          axc.plot((x[nrc]+dxb)*np.ones(2), 
                   np.array([data_esm_min_glbl[nr-5], data_esm_max_glbl[nr-5]]),
                   color='k',
                   linewidth=1)
          nrc = nrc + 1
        if cmip6_plt[3]:
          print('cmip6 plot [3]', data_esm_min_arc[nr-5], data_esm_max_arc[nr-5])  
          # print('cmip6 plot [3]', data_esm_min_arc[nr-5], data_esm_median_arc[nr-5], data_esm_max_arc[nr-5])
          rects3 = axd.scatter(x[nrd]+dxb, data_esm_median_arc[nr-5], c='k', \
                                  edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
          axd.plot((x[nrd]+dxb)*np.ones(2), 
              np.array([data_esm_min_arc[nr-5], data_esm_max_arc[nr-5]]),
              color='k',
              linewidth=1)
          # add error bars
          nrd = nrd + 1
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    # for nr in range(len(x)):
    #   if not cmip6[nr]:
    #     if not cmip6_plt[0] and xlabelsx_amap[nra] != 'CFM' and xlabelsx_amap[nra] != 'MFR_SDS':
    #       rects3 = axa.scatter(x[nra]+dxb1, region_values_netx[0,nr], c='black', \
    #                              edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
    #     if not cmip6_plt[1] and xlabelsx_amap[nrb] != 'CFM' and xlabelsx_amap[nrb] != 'MFR_SDS':
    #       rects3 = axb.scatter(x[nrb]+dxb1, region_values_netx[1,nr], c='black', \
    #                              edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
    #     if not cmip6_plt[2] and xlabelsx_amap[nrc] != 'CFM' and xlabelsx_amap[nrc] != 'MFR_SDS':
    #       rects3 = axc.scatter(x[nrc]+dxb1, region_values_netx[2,nr], c='black', \
    #                              edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
    #     if not cmip6_plt[3] and xlabelsx_amap[nrd] != 'CFM' and xlabelsx_amap[nrd] != 'MFR_SDS':
    #       rects3 = axd.scatter(x[nrd]+dxb1, region_values_netx[3,nr], c='black', \
    #                              edgecolors='black', s=20, zorder=9999, marker='D',linewidths=0)
              
    #     if not cmip6_plt[0] and xlabelsx_amap[nra] != 'MFR_SDS':
    #       rects3 = axa.scatter(x[nra]+dxb2, region_values_nety[0,nr], c='white', \
    #                              edgecolors='black', s=30, marker='^', zorder=9999,linewidths=0.8)
    #     if not cmip6_plt[1] and xlabelsx_amap[nrb] != 'MFR_SDS':
    #       rects3 = axb.scatter(x[nrb]+dxb2, region_values_nety[1,nr], c='white', \
    #                              edgecolors='black', s=30, marker='^', zorder=9999,linewidths=0.8)
    #     if not cmip6_plt[2] and xlabelsx_amap[nrc] != 'MFR_SDS':
    #       rects3 = axc.scatter(x[nrc]+dxb2, region_values_nety[2,nr], c='white', \
    #                              edgecolors='black', s=30, marker='^', zorder=9999,linewidths=0.8)
    #     if not cmip6_plt[3] and xlabelsx_amap[nrd] != 'MFR_SDS':
    #       rects3 = axd.scatter(x[nrd]+dxb2, region_values_nety[3,nr], c='white', \
    #                              edgecolors='black', s=30, marker='^', zorder=9999,linewidths=0.8)
    #     if not cmip6_plt[0]:
    #       nra = nra + 1
    #     if not cmip6_plt[1]:
    #       nrb = nrb + 1
    #     if not cmip6_plt[2]:
    #       nrc = nrc + 1
    #     if not cmip6_plt[3]:
    #       nrd = nrd + 1

  if prange:
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if cmip6[nr]:
        if cmip6_plt[0]:
          hlinex = [x[nra]-dxb, x[nra]-dxb]
          hliney = [max(region_values_max1[0,nr],ymin[0]),min(region_values_max2[0,nr],ymax[0])]
          rects3 = axa.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[0,nr],ymin[0]),min(region_values_min2[0,nr],ymax[0])]
          rects3 = axa.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nra]-dxb-dbx,max(region_values_min1[0,nr],ymin[0])),\
                               2.*dbx,region_values_max1[0,nr]-region_values_min1[0,nr],\
                               fill=False, edgecolor='black', linewidth=1)
        
          axa.add_patch(rect)
          nra = nra + 1
        if cmip6_plt[1]:
          hlinex = [x[nrb]-dxb, x[nrb]-dxb]
          hliney = [max(region_values_max1[1,nr],ymin[1]),min(region_values_max2[1,nr],ymax[1])]
          rects3 = axb.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[1,nr],ymin[1]),min(region_values_min2[1,nr],ymax[1])]
          rects3 = axb.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nrb]-dxb-dbx,max(region_values_min1[1,nr],ymin[1])),\
                               2.*dbx,region_values_max1[1,nr]-region_values_min1[1,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          axb.add_patch(rect)
          nrb = nrb + 1
        if cmip6_plt[2]:
          hlinex = [x[nrc]-dxb, x[nrc]-dxb]
          hliney = [max(region_values_max1[2,nr],ymin[2]),min(region_values_max2[2,nr],ymax[2])]
          rects3 = axc.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[2,nr],ymin[2]),min(region_values_min2[2,nr],ymax[2])]
          rects3 = axc.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nrc]-dxb-dbx,max(region_values_min1[2,nr],ymin[2])),\
                               2.*dbx,region_values_max1[2,nr]-region_values_min1[2,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          print('cmip6 plot [2] uncertainty',rect.get_height()/2)
          g5ae += [rect.get_height()/2]
          axc.add_patch(rect)
          nrc = nrc + 1
        if cmip6_plt[3]:
          hlinex = [x[nrd]-dxb, x[nrd]-dxb]
          hliney = [max(region_values_max1[3,nr],ymin[3]),min(region_values_max2[3,nr],ymax[3])]
          rects3 = axd.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[3,nr],ymin[3]),min(region_values_min2[3,nr],ymax[3])]
          rects3 = axd.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nrd]-dxb-dbx,max(region_values_min1[3,nr],ymin[3])),\
                               2.*dbx,region_values_max1[3,nr]-region_values_min1[3,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          print('cmip6 plot [3] uncertainty',rect.get_height()/2)
          a5ae += [rect.get_height()/2]
          axd.add_patch(rect)
          nrd = nrd + 1
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if not cmip6[nr]:
        if not cmip6_plt[0]:
          hlinex = [x[nra]-dxb, x[nra]-dxb]
          hliney = [max(region_values_max1[0,nr],ymin[0]),min(region_values_max2[0,nr],ymax[0])]
          rects3 = axa.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[0,nr],ymin[0]),min(region_values_min2[0,nr],ymax[0])]
          rects3 = axa.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nra]-dxb-dbx,max(region_values_min1[0,nr],ymin[0])),\
                               2.*dbx,region_values_max1[0,nr]-region_values_min1[0,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          axa.add_patch(rect)
          nra = nra + 1
        if not cmip6_plt[1]:
          hlinex = [x[nrb]-dxb, x[nrb]-dxb]
          hliney = [max(region_values_max1[1,nr],ymin[1]),min(region_values_max2[1,nr],ymax[1])]
          rects3 = axb.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[1,nr],ymin[1]),min(region_values_min2[1,nr],ymax[1])]
          rects3 = axb.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nrb]-dxb-dbx,max(region_values_min1[1,nr],ymin[1])),\
                               2.*dbx,region_values_max1[1,nr]-region_values_min1[1,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          axb.add_patch(rect)
          nrb = nrb + 1
        if not cmip6_plt[2]:
          hlinex = [x[nrc]-dxb, x[nrc]-dxb]
          hliney = [max(region_values_max1[2,nr],ymin[2]),min(region_values_max2[2,nr],ymax[2])]
          rects3 = axc.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[2,nr],ymin[2]),min(region_values_min2[2,nr],ymax[2])]
          rects3 = axc.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nrc]-dxb-dbx,max(region_values_min1[2,nr],ymin[2])),\
                               2.*dbx,region_values_max1[2,nr]-region_values_min1[2,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          axc.add_patch(rect)
          nrc = nrc + 1
        if not cmip6_plt[3]:
          hlinex = [x[nrd]-dxb, x[nrd]-dxb]
          hliney = [max(region_values_max1[3,nr],ymin[3]),min(region_values_max2[3,nr],ymax[3])]
          rects3 = axd.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          hliney = [max(region_values_min1[3,nr],ymin[3]),min(region_values_min2[3,nr],ymax[3])]
          rects3 = axd.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          rect = patches.Rectangle((x[nrd]-dxb-dbx,max(region_values_min1[3,nr],ymin[3])),\
                               2.*dbx,region_values_max1[3,nr]-region_values_min1[3,nr],\
                               fill=False, edgecolor='black', linewidth=1)
          axd.add_patch(rect)
          nrd = nrd + 1
  if pextra and not clean:
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if cmip6[nr]:
        if cmip6_plt[0]:
          hlinex = [x[nra]+dxb, x[nra]+dxb]
          hliney = [max(region_values_minx[0,nr],ymin[0]),min(region_values_maxx[0,nr],ymax[0])]
          rects3 = axa.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nra = nra + 1
        if cmip6_plt[1]:
          hlinex = [x[nrb]+dxb, x[nrb]+dxb]
          hliney = [max(region_values_minx[1,nr],ymin[1]),min(region_values_maxx[1,nr],ymax[1])]
          rects3 = axb.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nrb = nrb + 1
        if cmip6_plt[2]:
          hlinex = [x[nrc]+dxb, x[nrc]+dxb]
          hliney = [max(region_values_minx[2,nr],ymin[2]),min(region_values_maxx[2,nr],ymax[2])]
          rects3 = axc.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nrc = nrc + 1
        if cmip6_plt[3]:
          hlinex = [x[nrd]+dxb, x[nrd]+dxb]
          hliney = [max(region_values_minx[3,nr],ymin[3]),min(region_values_maxx[3,nr],ymax[3])]
          rects3 = axd.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nrd = nrd + 1
    nra = 0
    nrb = 0
    nrc = 0
    nrd = 0
    for nr in range(len(x)):
      if not cmip6[nr]:
        if not cmip6_plt[0]:
          hlinex = [x[nra]+dxb1, x[nra]+dxb1]
          hliney = [max(region_values_minx[0,nr],ymin[0]),min(region_values_maxx[0,nr],ymax[0])]
          rects3 = axa.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nra = nra + 1
        if not cmip6_plt[1]:
          hlinex = [x[nrb]+dxb1, x[nrb]+dxb1]
          hliney = [max(region_values_minx[1,nr],ymin[1]),min(region_values_maxx[1,nr],ymax[1])]
          rects3 = axb.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nrb = nrb + 1
        if not cmip6_plt[2]:
          hlinex = [x[nrc]+dxb1, x[nrc]+dxb1]
          hliney = [max(region_values_minx[2,nr],ymin[2]),min(region_values_maxx[2,nr],ymax[2])]
          rects3 = axc.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nrc = nrc + 1
        if not cmip6_plt[3]:
          hlinex = [x[nrd]+dxb1, x[nrd]+dxb1]
          hliney = [max(region_values_minx[3,nr],ymin[3]),min(region_values_maxx[3,nr],ymax[3])]
          rects3 = axd.plot(hlinex, hliney, linewidth=1, color='black', linestyle='-')
          nrd = nrd + 1

  if cmip6_plt[0]:
    axa.set_xticks(x_cmip6)
  else:
    axa.set_xticks(x_amap)
  if cmip6_plt[1]:
    axb.set_xticks(x_cmip6)
  else:
    axb.set_xticks(x_amap)
  if cmip6_plt[2]:
    axc.set_xticks(x_cmip6)
  else:
    axc.set_xticks(x_amap)
  if cmip6_plt[3]:
    axd.set_xticks(x_cmip6)
  else:
    axd.set_xticks(x_amap)
  xtolx = xtol
  axa.set_xlim([xmin-xtolx, xmax-1+xtolx])
  axb.set_xlim([xmin-xtolx, xmax-1+xtolx])
  axc.set_xlim([xmin-xtolx, xmax+xtolx])
  axd.set_xlim([xmin-xtolx, xmax+xtolx])
  ftsize = 9
  ftsize_ax = 10
  if not tightlx:
    if cmip6_plt[0]:
      axa.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=45)
    else:
      axa.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=45)
    if cmip6_plt[1]:
      axb.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=45)
    else:
      axb.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=45)
    if cmip6_plt[2]:
      axc.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=45)
    else:
      axc.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=45)
    if cmip6_plt[3]:
      axd.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=45)
    else:
      axd.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=45)
  else:
    if cmip6_plt[0]:
      axa.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=90)
    else:
      axa.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=90)
    if cmip6_plt[1]:
      axb.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=90)
    else:
      axb.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=90)
    if cmip6_plt[2]:
      axc.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=90)
    else:
      axc.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=90)
    if cmip6_plt[3]:
      axd.set_xticklabels(xlabelsx_cmip6, fontsize=ftsize, rotation=90)
    else:
      axd.set_xticklabels(xlabelsx_amap, fontsize=ftsize, rotation=90)

  axa.set_yticks(yticks[0,:])
  axb.set_yticks(yticks[1,:])
  axc.set_yticks(yticks[2,:])
  axd.set_yticks(yticks[3,:])

  for n, labela in enumerate(axa.yaxis.get_ticklabels()):
    if abs(yticks[0,n]-int(yticks[0,n])) >= 0.001:
      labela.set_visible(False)
  new_labels = np.full((len(yticks[0,:])),0)
  for n in range(len(yticks[0,:])):
    new_labels[n] = int(yticks[0,n])
  axa.set_yticklabels(new_labels)

  for n, labelb in enumerate(axc.yaxis.get_ticklabels()):
    if abs(yticks[2,n]-int(yticks[2,n])) >= 0.001:
      labelb.set_visible(False)
  new_labels = np.full((len(yticks[2,:])),0)
  for n in range(len(yticks[2,:])):
    new_labels[n] = int(yticks[2,n])
  axc.set_yticklabels(new_labels)
  
  axa.set_ylim([ymin[0], ymax[0]])
  axb.set_ylim([ymin[1], ymax[1]])
  axc.set_ylim([ymin[2], ymax[2]])
  axd.set_ylim([ymin[3], ymax[3]])
  axa.spines['top'].set_visible(False)
  axb.spines['top'].set_visible(False)
  axa.spines['right'].set_visible(False)
  axb.spines['right'].set_visible(False)
  axc.spines['right'].set_visible(False)
  axd.spines['right'].set_visible(False)
  
  if ymin[0] * ymax[0] < 0:
    axa.axhline(0, linewidth=1, color='black')
  if ymin[1] * ymax[1] < 0:
    axb.axhline(0, linewidth=1, color='black')
  if ymin[2] * ymax[2] < 0:
    axc.axhline(0, linewidth=1, color='black')
  if ymin[3] * ymax[3] < 0:
    axd.axhline(0, linewidth=1, color='black')

  # Empty plot for the subplot that contains the legend
  if legend:
    axl.spines['top'].set_visible(False)
    axl.spines['right'].set_visible(False)
    axl.spines['bottom'].set_visible(False)
    axl.spines['left'].set_visible(False)
    axl.set_xticks([])
    axl.set_yticks([])
    axl.set_xticklabels([])
    axl.set_yticklabels([])
  if legend:
    if tightl:
      fig.tight_layout(rect=[0, 0, 0.67, 1])
    else:
      fig.tight_layout()
  else:
    if tightl:
      fig.tight_layout(rect=[0, 0, 0.53, 1])
    else:
      fig.tight_layout(rect=[0, 0, 0.76, 1])
  
  if not partial:
    axa.text(xmin-2,ymax[0],'(a)',verticalalignment='center',horizontalalignment='center',color='black',fontweight='bold', fontsize=10)
    axb.text(xmin-2,ymax[1],'(b)',verticalalignment='center',horizontalalignment='center',color='black',fontweight='bold', fontsize=10)
    axc.text(xmin-2,ymax[2],'(c)',verticalalignment='center',horizontalalignment='center',color='black',fontweight='bold', fontsize=10)
    axd.text(xmin-2,ymax[3],'(d)',verticalalignment='center',horizontalalignment='center',color='black',fontweight='bold', fontsize=10)
  if partial:
    if aoff:
      axa.set_visible(False)
      axa.spines['top'].set_visible(False)
      axa.spines['right'].set_visible(False)
      axa.spines['bottom'].set_visible(False)
      axa.spines['left'].set_visible(False)
      axa.set_xticks([])
      axa.set_yticks([])
      axa.set_xticklabels([])
      axa.set_yticklabels([])
      axa2.set_visible(False)
      axa2.spines['top'].set_visible(False)
      axa2.spines['right'].set_visible(False)
      axa2.spines['bottom'].set_visible(False)
      axa2.spines['left'].set_visible(False)
      axa2.set_xticks([])
      axa2.set_yticks([])
      axa2.set_xticklabels([])
      axa2.set_yticklabels([])
    else:
      axa.set_ylabel(ylabel, fontsize=ftsize_ax)
      axa2.set_ylabel(ylabelx, color='grey', fontsize=ftsize_ax)  # we already handled the x-label with ax1
    if boff:
      axb.set_visible(False)
      axb.spines['top'].set_visible(False)
      axb.spines['right'].set_visible(False)
      axb.spines['bottom'].set_visible(False)
      axb.spines['left'].set_visible(False)
      axb.set_xticks([])
      axb.set_yticks([])
      axb.set_xticklabels([])
      axb.set_yticklabels([])
      axb2.set_visible(False)
      axb2.spines['top'].set_visible(False)
      axb2.spines['right'].set_visible(False)
      axb2.spines['bottom'].set_visible(False)
      axb2.spines['left'].set_visible(False)
      axb2.set_xticks([])
      axb2.set_yticks([])
      axb2.set_xticklabels([])
      axb2.set_yticklabels([])
    else:
      axb.set_ylabel(ylabel, fontsize=ftsize_ax)
      axb2.set_ylabel(ylabelx, color='grey', fontsize=ftsize_ax)  # we already handled the x-label with ax1
    if coff:
      axc.set_visible(False)
      axc.spines['top'].set_visible(False)
      axc.spines['right'].set_visible(False)
      axc.spines['bottom'].set_visible(False)
      axc.spines['left'].set_visible(False)
      axc.set_xticks([])
      axc.set_yticks([])
      axc.set_xticklabels([])
      axc.set_yticklabels([])
      axc2.set_visible(False)
      axc2.spines['top'].set_visible(False)
      axc2.spines['right'].set_visible(False)
      axc2.spines['bottom'].set_visible(False)
      axc2.spines['left'].set_visible(False)
      axc2.set_xticks([])
      axc2.set_yticks([])
      axc2.set_xticklabels([])
      axc2.set_yticklabels([])
    else:
      axc.set_ylabel(ylabel, fontsize=ftsize_ax)
      axc2.set_ylabel(ylabelx, color='grey', fontsize=ftsize_ax)  # we already handled the x-label with ax1
    if doff:
      axd.set_visible(False)
      axd.spines['top'].set_visible(False)
      axd.spines['right'].set_visible(False)
      axd.spines['bottom'].set_visible(False)
      axd.spines['left'].set_visible(False)
      axd.set_xticks([])
      axd.set_yticks([])
      axd.set_xticklabels([])
      axd.set_yticklabels([])
      axd2.set_visible(False)
      axd2.spines['top'].set_visible(False)
      axd2.spines['right'].set_visible(False)
      axd2.spines['bottom'].set_visible(False)
      axd2.spines['left'].set_visible(False)
      axd2.set_xticks([])
      axd2.set_yticks([])
      axd2.set_xticklabels([])
      axd2.set_yticklabels([])
    else:
      axd.set_ylabel(ylabel, fontsize=ftsize_ax)
      axd2.set_ylabel(ylabelx, color='grey', fontsize=ftsize_ax)  # we already handled the x-label with ax1
    if loff:
      axl.set_visible(False)
      axl.spines['top'].set_visible(False)
      axl.spines['right'].set_visible(False)
      axl.spines['bottom'].set_visible(False)
      axl.spines['left'].set_visible(False)
      axl.set_xticks([])
      axl.set_yticks([])
      axl.set_xticklabels([])
      axl.set_yticklabels([])
  
  #plt.show()
  plt.savefig('foo.png', bbox_inches='tight', dpi=600)
  plt.savefig('foo.svg', bbox_inches='tight', dpi=600)
  plt.savefig('foo.pdf', bbox_inches='tight', dpi=600)
  
  # Save results to a text file
  results = np.array([arctic2050, a5ae, global2050, g5ae])
  np.savetxt('2050_results_amap.txt', results, fmt='%s')




###############################################################################

# Do the plot

plot_bars(region_values_pos,region_values_neg,region_splits_pos,region_splits_neg,spltcolor,spltlabel,xlabelsx_amap,x_amap,xlabelsx_cmip6,x_cmip6,ylabel,title,ymin,ymax,yticks,yticksx,nplt,nsspy,titleinside,region_values_min1,region_values_min2,region_values_net,region_values_max1,region_values_max2,region_values_netx,region_values_minx,region_values_maxx,region_values_nety,region_values_miny,region_values_maxy,temp_off,pextra,prange,legend,cmip6,cmip6_plt,nplt_amap)
