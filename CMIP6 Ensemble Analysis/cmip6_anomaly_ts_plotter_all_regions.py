# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

Produce 5 subplots: the global, Arctic, NHML, Tropical, and SH
surface air temperature anomalies projected by the CMIP6 ensemble, for 4 SSPs. .
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from netCDF4 import Dataset # http://unidata.github.io/netcdf4-python//
import matplotlib.gridspec as gridspec

regions = ['GLOBAL', 'ARCTIC', 'NHML', 'TROPICS', 'SH']
n_regions = len(regions) # includes global
ssps = ['SSP126', 'SSP245', 'SSP370', 'SSP585']
ssp_labels = ['SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
colors0 = ['rebeccapurple', 'midnightblue', 'darkgoldenrod', 'firebrick', 'seagreen']
colors0 = ['rebeccapurple', 'seagreen', 'darkgoldenrod', 'crimson', 'seagreen']
project_start, project_end = 2015, 2050
ref_start, ref_end = 1995, 2014

region_prefixes = ['gma', 'arc', 'nhm', 'tro', 'she']     
n_regions = len(region_prefixes) # includes global

fig1 = plt.figure(figsize=(8,11)) #, tight_layout=True)
gs = gridspec.GridSpec(10, 1)
gs.update(wspace=0.25, hspace=1.2)
axs1 = fig1.add_subplot(gs[0:2, 0], )
axs2 = fig1.add_subplot(gs[2:4, 0])
axs3 = fig1.add_subplot(gs[4:6, 0])
axs4 = fig1.add_subplot(gs[6:8, 0])
axs5 = fig1.add_subplot(gs[8:10, 0])
ax1 = [axs1, axs2, axs3, axs4, axs5]


fig1.suptitle('Regional and Global SAT Anomalies Projected by CMIP6 Ensemble,\n Relative to 1995-2014',
              fontsize='medium', y=0.92, x=0.5)
ax1_xlabel = 'Year'
ax1_ylabel = ['Global Anomaly (°C)', 'Arctic Anomaly (°C)', 'NHML Anomaly (°C)', 
              'Tropical Anomaly (°C)', 'SH Anomaly (°C)']
minor_locs = [0.25, 0.25, 0.25, 0.25, 0.25]
major_locs = [1.0, 1.0, 1.0, 1.0, 1.0]

for i in range(0, 5, 1):
    ax1[i].set_xlabel(ax1_xlabel, fontsize='small')
    ax1[i].set_ylabel(ax1_ylabel[i], fontsize='small')
    ax1[i].axhline(y=0, linewidth=1, color='grey')
    ax1[i].axvline(x=2005, linestyle='--', linewidth=1, color='grey')
    ax1[i].grid(color='gainsboro')
    ax1[i].xaxis.set_minor_locator(MultipleLocator(1))
    ax1[i].xaxis.set_major_locator(MultipleLocator(5))
    ax1[i].yaxis.set_major_locator(MultipleLocator(major_locs[i]))
    ax1[i].yaxis.set_minor_locator(MultipleLocator(minor_locs[i]))
    ax1[i].set_xlim([ref_start+0.1,2050])

ax1[0].set_ylim([-1, 2.5])
ax1[1].set_ylim([-2, 5])
ax1[2].set_ylim([-1.25, 2.5])
ax1[3].set_ylim([-1, 2])
ax1[4].set_ylim([-1, 1.75])

for i in range(0, len(ssp_labels), 1):
    ax1[0].text(x=2050.75, y=1.65-i*0.47, s=ssp_labels[i], color=colors0[i], 
              fontsize='large', fontweight='bold')
    # ax1[0].text(x=2050.5, y=4.9-i*0.70 , s=ssp_labels[i], color=colors0[i], 
    #           fontsize='large', fontweight='bold')
ax1_years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator

n_ssps = len(ssp_labels)-1

# Arrays for computing region time series
project_start, project_end = 2015, 2050
ref_start, ref_end = 1995, 2014

data_dir = 'cmip6_data/'
ssp_dirs = ['DATA_ssp126/', 'DATA_ssp245/',
            'DATA_ssp370/', 'DATA_ssp585/']
n_years = 451
n_cmip6 = 1
n_ssps = len(ssp_labels)
# Arrays for computing region median time series
years = np.arange(1850, 2301, 1)
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

# Average amplifications for each of the regions (for each of the SSPs)
amplifications = np.zeros((n_regions-1, n_ssps)) 
amplifications_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_ssps)) # std dev of amps in ensemble at year

# Plot the time series for each CMIP6
for s in range(0, 4, 1):
    print(ssp_dirs[s])
    curr_data_dir = data_dir+ssp_dirs[s]
    curr_files = os.listdir(curr_data_dir)
    curr_files.remove('CanESM5') # folder with CanESM files copied
    curr_files.remove('Extras') # folders with duplicates
    n_files = len(curr_files)
    
    tas_count = np.zeros(n_regions)
    
    # Arrays for computing region mean/median time series
    gm_tas_time_series = np.zeros((1, n_years)) # TAS values from following years
    arc_tas_time_series = np.zeros((1, n_years))
    nh_tas_time_series = np.zeros((1, n_years))
    trop_tas_time_series = np.zeros((1, n_years))
    sh_tas_time_series = np.zeros((1, n_years))
    
    for f in range(0, n_files, 1): # For each file of the current SSP
        file = curr_files[f]
        prefix = file[0:3]
        ds = Dataset(curr_data_dir+file,"r",format="NETCDF4")  # Open dataset
        ds_md = ds.__dict__  # Extract metadata
        
        tas = ds['tas'][:][:,0,0] - 273.16
        time = ds['time'][:]
        curr_n = len(time)
        if prefix in region_prefixes:
            region = region_prefixes.index(prefix)
            tas_count[region] +=1
            
            # Add the projection to the regional array
            new_row = np.zeros(n_years)
            new_row[0:curr_n] = tas 
            
            if prefix == 'gma':
                gm_tas_time_series = np.vstack([gm_tas_time_series, new_row])
            if prefix == 'arc':
                arc_tas_time_series = np.vstack([arc_tas_time_series, new_row])
            if prefix == 'nhm':
                nh_tas_time_series = np.vstack([nh_tas_time_series, new_row])
            if prefix == 'tro':
                trop_tas_time_series = np.vstack([trop_tas_time_series, new_row])
            if prefix == 'she':
                sh_tas_time_series = np.vstack([sh_tas_time_series, new_row])
    gm_tas_time_series = np.delete(gm_tas_time_series, 0, 0)
    arc_tas_time_series = np.delete(arc_tas_time_series, 0, 0)
    nh_tas_time_series = np.delete(nh_tas_time_series, 0, 0)
    trop_tas_time_series = np.delete(trop_tas_time_series, 0, 0)
    sh_tas_time_series = np.delete(sh_tas_time_series, 0, 0)
    
    # Take the MEANS, each row is a region
    tas_time_series_mean = np.zeros((n_regions, n_years)) # Define mean temp. matrix
    tas_mean_anomaly = np.zeros((n_regions, n_years)) # Define mean temp. anomaly matrix
    sigma = np.zeros((n_regions, n_years))
    errorbars = np.zeros((n_regions, n_years))

    mean_refs = np.zeros(n_regions) # Define array for regional mean temperatures of ref. period
    for r in range(0, n_regions, 1):
        tas_time_series_mean[0, :] = np.mean(gm_tas_time_series, axis = 0)
        tas_time_series_mean[1, :] = np.mean(arc_tas_time_series, axis = 0)
        tas_time_series_mean[2, :] = np.mean(nh_tas_time_series, axis = 0)
        tas_time_series_mean[3, :] = np.mean(trop_tas_time_series, axis = 0)
        tas_time_series_mean[4, :] = np.mean(sh_tas_time_series, axis = 0)
        
        sigma[0, :] = np.std(gm_tas_time_series, axis = 0)
        sigma[1, :] = np.std(arc_tas_time_series, axis = 0)
        sigma[2, :] = np.std(nh_tas_time_series, axis = 0)
        sigma[3, :] = np.std(trop_tas_time_series, axis = 0)
        sigma[4, :] = np.std(sh_tas_time_series, axis = 0)
        
        mean_refs[r] = np.mean(tas_time_series_mean[r, ref_start_i:ref_end_i]) 
        
        ref = np.mean(tas_time_series_mean[r, ref_start_i:ref_end_i])
        global_ref = np.mean(tas_time_series_mean[0, ref_start_i:ref_end_i])
        print('global ref', global_ref)
        ref_sigma = np.sqrt(np.sum(np.square(sigma[r, ref_start_i:ref_end_i])))/tas_count[r]
        # Compute mean time series SAT anomaly for each region
        tas_mean_anomaly[r,:] = tas_time_series_mean[r,:] - ref
        # compute standard deviation of the anomaly for each region
        errorbars[r,:] = np.sqrt(np.square(sigma[r,:])+ref_sigma**2)
        if r == 0:
            gm_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
        if r == 1:
            arc_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
        if r == 2:
            nh_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
        if r == 3:
            trop_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
        if r == 4:
            sh_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
                
            ax1[0].errorbar(x=years[ref_start_i:project_end_i+1], 
                            y=gm_ensemble_anomaly[ref_start_i:project_end_i+1], 
                            yerr=errorbars[0, ref_start_i:project_end_i+1],
                            color=colors0[s], ecolor=colors0[s], marker='o', 
                            markersize=2, capsize=0, elinewidth=0, linestyle='-'
                            )
            ax1[1].errorbar(x=years[ref_start_i:project_end_i+1], 
                            y=arc_ensemble_anomaly[ref_start_i:project_end_i+1],
                            yerr=errorbars[1, ref_start_i:project_end_i+1],
                            color=colors0[s], ecolor=colors0[s], marker='o', 
                            markersize=2, capsize=0, elinewidth=0, linestyle='-')
            ax1[2].errorbar(x=years[ref_start_i:project_end_i+1], 
                            y=nh_ensemble_anomaly[ref_start_i:project_end_i+1],
                            yerr=errorbars[2, ref_start_i:project_end_i+1],
                            color=colors0[s], ecolor=colors0[s], marker='o', 
                            markersize=2, capsize=0, elinewidth=0, linestyle='-')
            ax1[3].errorbar(x=years[ref_start_i:project_end_i+1], 
                            y=trop_ensemble_anomaly[ref_start_i:project_end_i+1], 
                            yerr=errorbars[3, ref_start_i:project_end_i+1],
                            color=colors0[s], ecolor=colors0[s], marker='o', 
                            markersize=2, capsize=0, elinewidth=0, linestyle='-')
            ax1[4].errorbar(x=years[ref_start_i:project_end_i+1], 
                            y=sh_ensemble_anomaly[ref_start_i:project_end_i+1], 
                            yerr=errorbars[4, ref_start_i:project_end_i+1],
                            color=colors0[s], ecolor=colors0[s], marker='o', 
                            markersize=2, capsize=0, elinewidth=0, linestyle='-')

            ax1[0].fill_between(years[ref_start_i:project_end_i+1], 
                                gm_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                gm_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.12)
            ax1[1].fill_between(years[ref_start_i:project_end_i+1], 
                                arc_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                arc_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.12)
            ax1[2].fill_between(years[ref_start_i:project_end_i+1], 
                                nh_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                nh_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.12)
            ax1[3].fill_between(years[ref_start_i:project_end_i+1], 
                                trop_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                trop_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.12)
            ax1[4].fill_between(years[ref_start_i:project_end_i+1], 
                                sh_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                sh_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.12)
    
fig1.savefig('cmip6_ts_default_all_glb.png')