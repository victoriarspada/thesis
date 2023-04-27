# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

This code produces a figure with 5 subplots: the global, Arctic, NHML, Tropical, 
and SH surface air temperature anomalies projected by the AMAP emulator for 4 
SSPs. CMIP6 results are also plotted for each region, for comparison.
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


fig1.suptitle('Regional and Global SAT Anomalies Projected by AMAP Ensemble,\n Relative to 1995-2014',
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
    ax1[i].set_xlim([1995+0.1,project_end+0.1])
    
for i in range(0, len(ssp_labels), 1):
    ax1[0].text(x=2050.75, y=2.7-i*0.8, s=ssp_labels[i], color=colors0[i], 
              fontsize='large', fontweight='bold')
ax1_years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator

datadir = 'C:/Users/victo/Downloads/EngSci Year 4 Sem 1/ESC499 - Thesis/zenodo_bundle/AMAP Projections/'
renormdir = 'Renormed All/'
outdir = 'Renormed Ensemble All/'

fname_prefix = 'dtemp_'
n_ssps = len(ssp_labels)

# Arrays for computing region time series
n_years = 61
project_start, project_end = 2015, 2050
ref_start, ref_end = 1995, 2014
years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

# Average amplifications for each of the regions (for each of the SSPs)
amplifications = np.zeros((n_regions-1, n_years)) 
amplifications_sigma = np.zeros((n_regions-1, n_years)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_years)) # std dev of amps in ensemble at year

n_ensemble = 24 # number of ensemble members

# Plot the time series for each SSP
for s in range(0, n_ssps, 1): # For each file (current delta setting) 
    file = 'dtemp_' + ssps[s] + '_ensemble.nc'
    print(file)
    ds = Dataset(datadir+outdir+file,"r",format="NETCDF4")  # Open dataset
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
                    markersize=2, capsize=1, elinewidth=1, linestyle='-')
    
    ax1[1].errorbar(x=ax1_years, y=arc_tas_time_series, 
                    yerr=arc_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=1, elinewidth=1, linestyle='-')
    
    ax1[2].errorbar(x=ax1_years, y=nhml_tas_time_series, 
                    yerr=nhml_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=1, elinewidth=1, linestyle='-')
    
    ax1[3].errorbar(x=ax1_years, y=tropics_tas_time_series,
                    yerr=tropics_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=1, elinewidth=1, linestyle='-')
    
    ax1[4].errorbar(x=ax1_years, y=sh_tas_time_series,
                    yerr=sh_tas_time_series_sigma, 
                    color=colors0[s], ecolor=colors0[s], marker='o', 
                    markersize=2, capsize=1, elinewidth=1, linestyle='-')
    
fig1.savefig('amap_ts_default_all.png')