# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

Plotter for modified AMAP emulator results.
Produces a figure with 4 subplots: the Arctic, NHML, Tropical, and SH
regional amplifications projected by the AMAP emulator for 4 SSPs. 
CMIP6 results are also plotted for each region, for comparison.

The regional amplfication for a given year is modeled as the ratio of the 
regional temperature anomaly to the global anomaly relative to 1995-2014.
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

fig1 = plt.figure(figsize=(8,10)) #, tight_layout=True)
gs = gridspec.GridSpec(8, 1)
gs.update(wspace=0.25, hspace=0.9)
axs1 = fig1.add_subplot(gs[0:2, 0], )
axs2 = fig1.add_subplot(gs[2:4, 0])
axs3 = fig1.add_subplot(gs[4:6, 0])
axs4 = fig1.add_subplot(gs[6:8, 0])
ax1 = [axs1, axs2, axs3, axs4]


fig1.suptitle('Regional Amplifications Projected by AMAP Ensemble,\n Relative to 1995-2014',
              fontsize='large', y=0.93, x=0.5)
ax1_xlabel = 'Year'
ax1_ylabel = ['Arctic Amplification', 'NHML Amplification', 'Tropical Amplification',
              'SH Amplification']
minor_locs = [0.05, 0.01, 0.01, 0.01]
major_locs = [0.25, 0.05, 0.05, 0.05]

for i in range(0, 4, 1):
    ax1[i].set_xlabel(ax1_xlabel, fontsize='large')
    ax1[i].set_ylabel(ax1_ylabel[i], fontsize='large')
    ax1[i].xaxis.set_minor_locator(MultipleLocator(1))
    ax1[i].xaxis.set_major_locator(MultipleLocator(5))
    ax1[i].yaxis.set_major_locator(MultipleLocator(major_locs[i]))
    ax1[i].yaxis.set_minor_locator(MultipleLocator(minor_locs[i]))
    ax1[i].set_xlim([ref_end+0.1,project_end+0.1])
    
ax1[0].set_ylim([1.4, 2.2])
ax1[1].set_ylim([1.28, 1.40])
ax1[2].set_ylim([0.93, 1.06])
ax1[3].set_ylim([0.47, 0.60])
for i in range(0, len(ssp_labels), 1):
    ax1[0].text(x=2050.75, y=2.1 -i*0.15, s=ssp_labels[i], color=colors0[i], 
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
amplifications = np.zeros((n_regions-1, n_ssps)) 
amplifications_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_ssps)) # std dev of anomalies in ensemble at year

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
    ax1[0].errorbar(x=ax1_years, y=curr_arc_amp, yerr=curr_arc_amp_std, 
                   color=colors0[s], ecolor=colors0[s], marker='o', 
                   markersize=2, capsize=1, elinewidth=1, linestyle='-')
    # ax1[0].scatter(x=ax1_years, y=curr_arc_amp, 
    #             c=colors0[s], marker='o',s=5) 
    # ax1[0].fill_between(ax1_years,
    #                     curr_arc_amp-curr_arc_amp_std, 
    #                     curr_arc_amp+curr_arc_amp_std,
    #                     color=colors0[s], alpha=0.10)
    
    ax1[1].errorbar(x=ax1_years, y=curr_nhml_amp, yerr=curr_nhml_amp_std, 
                  color=colors0[s], ecolor=colors0[s], marker='o', 
                  markersize=2, capsize=1, elinewidth=1, linestyle='-')
    # ax1[1].fill_between(ax1_years,
    #                     curr_nhml_amp-curr_nhml_amp_std, 
    #                     curr_nhml_amp+curr_nhml_amp_std,
    #                     color=colors0[s], alpha=0.15)
    
    ax1[2].errorbar(x=ax1_years, y=curr_tropics_amp, yerr=curr_tropics_amp_std, 
                  color=colors0[s], ecolor=colors0[s], marker='o', 
                  markersize=2, capsize=1, elinewidth=1, linestyle='-')
    # ax1[2].fill_between(ax1_years,
    #                     curr_tropics_amp-curr_tropics_amp_std, 
    #                     curr_tropics_amp+curr_tropics_amp_std,
    #                     color=colors0[s], alpha=0.10)
    
    ax1[3].errorbar(x=ax1_years, y=curr_sh_amp, yerr=curr_sh_amp_std, 
                  color=colors0[s], ecolor=colors0[s], marker='o', 
                  markersize=2, capsize=1, elinewidth=1, linestyle='-')
    # ax1[3].fill_between(ax1_years,
    #                     curr_sh_amp-curr_sh_amp_std, 
    #                     curr_sh_amp+curr_sh_amp_std,
    #                     color=colors0[s], alpha=0.10)
    
fig1.savefig('amap_ts_default_all.png')