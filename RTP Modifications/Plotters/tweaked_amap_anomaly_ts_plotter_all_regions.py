# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

Plotter for two figures. 

The first plot has 5 subplots: the global, Arctic, NHML, Tropical, and SH
surface air temperature anomalies projected by the AMAP emulator for various 
RTP adjustments, over time in SSP5-8.5. CMIP6 results are also plotted for each 
region, for comparison.

In the second plot, there is only the global subplot from the first figure.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from netCDF4 import Dataset # http://unidata.github.io/netcdf4-python//
import matplotlib.gridspec as gridspec

regions = ['GLOBAL', 'ARCTIC', 'NHML', 'TROPICS', 'SH']
n_regions = len(regions) # includes global
ssps = ['SSP119', 'SSP126', 'SSP245', 'SSP370', 'SSP585']
ssp = ssps[4]
ssp_labels = ['\u03B4=0.00', '\u03B4=0.10', '\u03B4=0.25', '\u03B4=0.50', 'CMIP6']
colors0 = ['mediumaquamarine', 'maroon', 'slateblue', 'mediumvioletred', 'darkgoldenrod']
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

fig2, ax2 = plt.subplots() 

fig1.suptitle('Regional SAT Anomalies Projected by Tweaked AMAP Ensembles and CMIP6 Ensemble,\n Relative to 1995-2014 (SSP5-8.5)',
              fontsize='medium', y=0.92, x=0.5)
fig2.suptitle('Global SAT Anomalies Projected by Tweaked AMAP Ensembles and CMIP6 Ensemble,\n Relative to 1995-2014 (SSP5-8.5)',
              fontsize='medium', x=0.5)
ax1_xlabel = 'Year'
ax1_ylabel = ['Global Anomaly (°C)', 'Arctic Anomaly (°C)', 'NHML Anomaly (°C)', 
              'Tropical Anomaly (°C)', 'SH Anomaly (°C)']
minor_locs = [0.25, 0.5, 0.25, 0.25, 0.25]
major_locs = [1.0, 2.0, 1.0, 1.0, 1.0]

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
    ax1[i].set_xlim([1995+0.1,project_end+0.1])
    
ax2.set_xlabel(ax1_xlabel, fontsize='small')
ax2.set_ylabel(ax1_ylabel[0], fontsize='small')
ax2.axhline(y=0, linewidth=1, color='grey')
ax2.axvline(x=2005, linestyle='--', linewidth=1, color='grey')
ax2.grid(color='gainsboro')
ax2.xaxis.set_minor_locator(MultipleLocator(1))
ax2.xaxis.set_major_locator(MultipleLocator(5))
ax2.yaxis.set_major_locator(MultipleLocator(major_locs[0]))
ax2.yaxis.set_minor_locator(MultipleLocator(minor_locs[0]))
ax2.set_xlim([1995+0.1,project_end+0.1])    
    
for i in range(0, len(ssp_labels), 1):
    ax1[0].text(x=2050.75, y=2.8-i*0.75, s=ssp_labels[i], color=colors0[i], 
              fontsize='large', fontweight='bold')
    ax2.text(x=2050.75, y=3-i*0.5, s=ssp_labels[i], color=colors0[i], 
              fontsize='large', fontweight='bold')
ax1_years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator


data_dir = 'AMAP Projections/Renormed Mod1 Ensemble/'
data_dirs = ['AMAP Projections/Renormed Ensemble All/',
             'AMAP Projections/Renormed Mod3 Ensemble/',
             'AMAP Projections/Renormed Mod2 Ensemble/',
             'AMAP Projections/Renormed Mod1 Ensemble/']

fname_prefix = 'dtemp_'
n_ssps = len(ssp_labels)-1

# Arrays for computing region time series
n_years = 61
project_start, project_end = 2015, 2050
ref_start, ref_end = 1995, 2014
years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

# Average amplifications for each of the regions (for each of the SSPs)
amplifications = np.zeros((n_regions-1, n_ssps)) 
amplifications_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_ssps)) # std dev of anomalies in ensemble at year

n_ensemble = 24 # number of ensemble members

# Plot the time series for each SSP
for s in range(0, n_ssps, 1): # For each file (current delta setting) 
    file = 'dtemp_' + ssp + '_ensemble.nc'
    print(file)
    ds = Dataset(data_dirs[s]+file,"r",format="NETCDF4")  # Open dataset
    ds_md = ds.__dict__  # Extract 
    
    time = ds['time'][:]
    # arrays of time series anomaly ensembles
    gm_tas_time_series_ensemble = ds['temp_glb'][:]
    arc_tas_time_series_ensemble = ds['temp_ARCTIC'][:]
    nhml_tas_time_series_ensemble = ds['temp_NHML'][:]
    tropics_tas_time_series_ensemble = ds['temp_TROPICS'][:]
    sh_tas_time_series_ensemble = ds['temp_SH'][:]
    
    # anomaly averages
    gm_tas_time_series = np.mean(gm_tas_time_series_ensemble, axis = 0)
    arc_tas_time_series = np.mean(arc_tas_time_series_ensemble, axis = 0)
    nhml_tas_time_series = np.mean(nhml_tas_time_series_ensemble, axis = 0)
    tropics_tas_time_series = np.mean(tropics_tas_time_series_ensemble, axis = 0)
    sh_tas_time_series = np.mean(sh_tas_time_series_ensemble, axis = 0)
    # anomaly standard deviations
    gm_tas_time_series_sigma = np.std(gm_tas_time_series_ensemble, axis = 0)
    arc_tas_time_series_sigma = np.std(arc_tas_time_series_ensemble, axis = 0)
    nhml_tas_time_series_sigma = np.std(nhml_tas_time_series_ensemble, axis = 0)
    tropics_tas_time_series_sigma = np.std(tropics_tas_time_series_ensemble, axis = 0)
    sh_tas_time_series_sigma = np.std(sh_tas_time_series_ensemble, axis = 0)
    
    # Plotting for Figure 1
    time_series = [arc_tas_time_series, nhml_tas_time_series, tropics_tas_time_series, 
                   sh_tas_time_series]

    # Ensemble of amplifications
    curr_arc_amp_ensemble = arc_tas_time_series_ensemble/gm_tas_time_series_ensemble
    curr_nhml_amp_ensemble = nhml_tas_time_series_ensemble/gm_tas_time_series_ensemble
    curr_tropics_amp_ensemble = tropics_tas_time_series_ensemble/gm_tas_time_series_ensemble
    curr_sh_amp_ensemble = sh_tas_time_series_ensemble/gm_tas_time_series_ensemble
    # Take the average
    curr_arc_amp = np.mean(curr_arc_amp_ensemble, axis = 0)
    curr_nhml_amp = np.mean(curr_nhml_amp_ensemble, axis = 0)
    curr_tropics_amp = np.mean(curr_tropics_amp_ensemble, axis = 0)
    curr_sh_amp = np.mean(curr_sh_amp_ensemble, axis = 0)
    # take the mean for the year of interest
    amplifications[0,s] = curr_arc_amp[project_end_i]
    amplifications[1,s] = curr_nhml_amp[project_end_i]
    amplifications[2,s] = curr_tropics_amp[project_end_i]
    amplifications[3,s] = curr_sh_amp[project_end_i]
    
    # standard deviations of ensemble of amplifications
    curr_arc_amp_std = np.std(curr_arc_amp_ensemble, axis = 0)
    curr_nhml_amp_std = np.std(curr_nhml_amp_ensemble, axis = 0)
    curr_tropics_amp_std = np.std(curr_tropics_amp_ensemble, axis = 0)
    curr_sh_amp_std = np.std(curr_sh_amp_ensemble, axis = 0)
    # take the std for the year of interest
    amplifications_std[0,s] = curr_arc_amp_std[project_end_i]
    amplifications_std[1,s] = curr_nhml_amp_std[project_end_i]
    amplifications_std[2,s] = curr_tropics_amp_std[project_end_i]
    amplifications_std[3,s] = curr_sh_amp_std[project_end_i]
    
    # # Individual uncertainties from dividing the two averages for a mean amp
    curr_arc_amp_sigma = np.sqrt( np.square(arc_tas_time_series_sigma/arc_tas_time_series) +\
                                  np.square(gm_tas_time_series_sigma/gm_tas_time_series) )
    curr_nhml_amp_sigma = np.sqrt( np.square(nhml_tas_time_series_sigma/nhml_tas_time_series_sigma) +\
                                  np.square(gm_tas_time_series_sigma/gm_tas_time_series) )
    curr_tropics_amp_sigma = np.sqrt( np.square(tropics_tas_time_series_sigma/tropics_tas_time_series) +\
                                  np.square(gm_tas_time_series_sigma/gm_tas_time_series) )
    curr_sh_amp_sigma = np.sqrt( np.square(sh_tas_time_series_sigma/sh_tas_time_series) +\
                                  np.square(gm_tas_time_series_sigma/gm_tas_time_series_sigma) )
    
    # Plot time series for fig1 (AMAP)
    ax1[0].errorbar(x=ax1_years, y=gm_tas_time_series, 
                    yerr=gm_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    # markersize=0, capsize=2, elinewidth=1, linestyle='-')
                    markersize=2, capsize=2, elinewidth=1, linestyle='-')
    ax2.errorbar(x=ax1_years, y=gm_tas_time_series, 
                    yerr=gm_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=0, elinewidth=0, linestyle='-')
    ax2.fill_between(ax1_years,
                     gm_tas_time_series-gm_tas_time_series_sigma, 
                     gm_tas_time_series+gm_tas_time_series_sigma,
                     color=colors0[s], alpha=0.10)
    
    # ax1[0].scatter(x=ax1_years, y=curr_arc_amp, 
    #             c=colors0[s], marker='o',s=5) 
    # ax1[0].fill_between(ax1_years,
    #                     curr_arc_amp-curr_arc_amp_std, 
    #                     curr_arc_amp+curr_arc_amp_std,
    #                     color=colors0[s], alpha=0.10)
    
    ax1[1].errorbar(x=ax1_years, y=arc_tas_time_series, 
                    yerr=arc_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=2, elinewidth=1, linestyle='-')
    # ax1[1].fill_between(ax1_years,
    #                     curr_nhml_amp-curr_nhml_amp_std, 
    #                     curr_nhml_amp+curr_nhml_amp_std,
    #                     color=colors0[s], alpha=0.15)
    
    ax1[2].errorbar(x=ax1_years, y=nhml_tas_time_series, 
                    yerr=nhml_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=2, elinewidth=1, linestyle='-')
    # ax1[2].fill_between(ax1_years,
    #                     curr_tropics_amp-curr_tropics_amp_std, 
    #                     curr_tropics_amp+curr_tropics_amp_std,
    #                     color=colors0[s], alpha=0.10)
    
    ax1[3].errorbar(x=ax1_years, y=tropics_tas_time_series, 
                    yerr=tropics_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=0, elinewidth=0, linestyle='-')
    # ax1[3].fill_between(ax1_years,
    #                     curr_sh_amp-curr_sh_amp_std, 
    #                     curr_sh_amp+curr_sh_amp_std,
    #                     color=colors0[s], alpha=0.10)
    
    ax1[4].errorbar(x=ax1_years, y=sh_tas_time_series, 
                    yerr=sh_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=0, elinewidth=0, linestyle='-')
    

data_dir = 'cmip6_data/'
ssp_dirs = ['DATA_ssp585/']
colors0 = ['darkgoldenrod']
n_years = 451
n_cmip6 = 1
# Arrays for computing region median time series
years = np.arange(1850, 2301, 1)
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

gm_anomalies = np.zeros((n_cmip6, n_years)) # mean anomalies for each ssp
arc_anomalies = np.zeros((n_cmip6+1, n_years))
nh_anomalies = np.zeros((n_cmip6+1, n_years))
trop_anomalies = np.zeros((n_cmip6+1, n_years))
sh_anomalies = np.zeros((n_cmip6+1, n_years))

# Average amplifications for each of the regions (for each of the SSPs)
amplifications = np.zeros((n_regions-1, n_cmip6)) 
amplifications_sigma = np.zeros((n_regions-1, n_cmip6)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_cmip6)) # std dev of amps in ensemble at year

# Plot the time series for each CMIP6
for s in range(0, 1, 1):
    print(ssp_dirs[s])
    curr_data_dir = data_dir+ssp_dirs[s]
    curr_files = os.listdir(curr_data_dir)
    curr_files.remove('CanESM5') # folder with CanESM files copied
    curr_files.remove('Extras') # folders with duplicates
    n_files = len(curr_files)
    
    base_tas = np.zeros((n_regions, n_years)) # base SAT values from ref. period
    base_tas_count = np.zeros(n_regions)
    
    tas_time_series = np.zeros((n_regions, n_years)) # TAS values from following years
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
            tas_time_series[region, 0:curr_n] += tas 
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
    
    # Take the means, each row is a region
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
        
        mean_refs[r] = np.mean(tas_time_series[r, ref_start_i:ref_end_i]) 
        ref = np.mean(tas_time_series[r, ref_start_i:ref_end_i]/tas_count[r]) 
        ref_sigma = np.sqrt(np.sum(np.square(sigma[r, ref_start_i:ref_end_i])))/tas_count[r]
        # Compute mean time series SAT anomaly for each region
        tas_mean_anomaly[r,:] = tas_time_series_mean[r,:] - ref
        # compute standard deviation of the anomaly for each region
        errorbars[r,:] = np.sqrt(np.square(sigma[r,:])+ref_sigma**2)
        if r == 0:
            gm_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(gm_tas_time_series[:, y] - ref)
        if r == 1:
            arc_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(arc_tas_time_series[:, y] - ref)
        if r == 2:
            nh_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(nh_tas_time_series[:, y] - ref)
        if r == 3:
            trop_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(trop_tas_time_series[:, y] - ref)
        if r == 4:
            sh_ensemble_anomaly =  tas_time_series_mean[r, :] - ref
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(sh_tas_time_series[:, y] - ref)
                
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
            
            ax2.errorbar(x=years[ref_start_i:project_end_i+1], 
                            y=gm_ensemble_anomaly[ref_start_i:project_end_i+1], 
                            yerr=errorbars[0, ref_start_i:project_end_i+1],
                            color=colors0[s], ecolor=colors0[s], marker='o', 
                            markersize=2, capsize=0, elinewidth=0, linestyle='-'
                            )
            
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
                                color=colors0[s], alpha=0.18)
            ax2.fill_between(years[ref_start_i:project_end_i+1], 
                                gm_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                gm_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.18)
            
            ax1[1].fill_between(years[ref_start_i:project_end_i+1], 
                                arc_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                arc_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.18)
            ax1[2].fill_between(years[ref_start_i:project_end_i+1], 
                                nh_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                nh_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.18)
            ax1[3].fill_between(years[ref_start_i:project_end_i+1], 
                                trop_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                trop_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.18)
            ax1[4].fill_between(years[ref_start_i:project_end_i+1], 
                                sh_ensemble_anomaly[ref_start_i:project_end_i+1]-errorbars[0, ref_start_i:project_end_i+1],
                                sh_ensemble_anomaly[ref_start_i:project_end_i+1]+errorbars[0, ref_start_i:project_end_i+1],
                                color=colors0[s], alpha=0.18)
    
fig1.savefig('amap_tweaked_ts_all.png')