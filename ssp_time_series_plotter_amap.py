# -*- coding: utf-8 -*-
"""
@author: Victoria Spada
SSP Comparisons (Time series)
"""

import os
import scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import xarray as xr
from netCDF4 import Dataset # http://unidata.github.io/netcdf4-python//
import re
import csv
from regression import *
from scipy.stats import pearsonr

data_dir = 'AMAP Projections/Renormed Ensemble All/'
regions = ['GLOBAL', 'ARCTIC', 'NHML', 'TROPICS', 'SH']
ssps = ['SSP119', 'SSP126', 'SSP245', 'SSP370', 'SSP585']
ssp_labels = ['SSP1-1.9', 'SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
colors0 = ['green', 'blue', 'darkorange', 'red', 'purple']
fname_prefix = 'dtemp_'
n_ssps = len(ssps)

# Arrays for computing region median time series
n_years = 61
project_start, project_end = 2015, 2050
ref_start, ref_end = 1995, 2014
years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

gm_anomalies = np.zeros((n_ssps+1, n_years)) #average anomalies for each of the SSPs
arc_anomalies = np.zeros((n_ssps+1, n_years))

gm_tas_time_series = np.zeros((1, n_years)) # TAS values from following years
arc_tas_time_series = np.zeros((1, n_years))

gm_tas_time_series_err = np.zeros((1, n_years)) # TAS uncertainty values from following years
arc_tas_time_series_err = np.zeros((1, n_years))

n_ensemble = 63 # number of ensemble members

# Plot regional SAT anomalies vs the global SAT anomalies (Figure 1)
fig1, ax1 = plt.subplots(1, 4, figsize=(11.5,3.5))
ax1_ylim = [[5.0, 4.0, 3.0, 2.0], [2.25, 1.25, 0.80, 0.75]]
ax1_xlabels = ['Arctic SAT Anomaly (°C)', 'NH Mid-latitude SAT Anomaly (°C)', 
             'Tropical SAT Anomaly (°C)', 'SH extratropical SAT Anomaly (°C)']
ax1_ylabel = 'Global SAT Anomaly (°C)'
suptitle1 = 'Regional SAT Anomalies vs Global Anomaly of AMAP Ensemble ('+str(project_start)+' - '+str(project_end)+')'
fig1.suptitle(suptitle1)
fig1.subplots_adjust(wspace=0.3)
    
# Set axis titles and ticks for each region (Figure 1)
for n in range(0,len(ax1_xlabels),1):
    ax1[n].set_xlabel(ax1_ylabel)
    ax1[n].set_ylabel(ax1_xlabels[n], fontsize='medium')
    if project_end==2050:
        ax1[n].set_xlim([0,3.0])
        ax1[n].set_ylim([0, ax1_ylim[0][n]])
        ax1[n].xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1[n].xaxis.set_major_locator(MultipleLocator(0.5))
        ax1[n].yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1[n].yaxis.set_major_locator(MultipleLocator(0.5))
    else:
        ax1[n].set_xlim([0,1.0])
        ax1[n].set_ylim([0, ax1_ylim[1][n]])
        ax1[n].xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1[n].xaxis.set_major_locator(MultipleLocator(0.5))
        if n==4:
            ax1[n].yaxis.set_minor_locator(MultipleLocator(0.025))
            ax1[n].yaxis.set_major_locator(MultipleLocator(0.1))
        if n==3:
            ax1[n].yaxis.set_minor_locator(MultipleLocator(0.025))
            ax1[n].yaxis.set_major_locator(MultipleLocator(0.05))
        if n==2:
            ax1[n].yaxis.set_minor_locator(MultipleLocator(0.05))
            ax1[n].yaxis.set_major_locator(MultipleLocator(0.1))
        else:
            ax1[n].yaxis.set_minor_locator(MultipleLocator(0.1))
            ax1[n].yaxis.set_major_locator(MultipleLocator(0.5))
    ax1[n].tick_params(direction="in", which="both", top=True, 
                      bottom=True, left=True, right=True)
    ax1[n].axhline(0, color='grey', linewidth=0.8)    
    
    
# Plot regional amplifications for SSP (Figure 2)
fig2, ax2 = plt.subplots(1, 3, figsize=(13,3.5))
fig3, ax3 = plt.subplots(1, 3, figsize=(13,3.5))
fig4, ax4 = plt.subplots(1, 3, figsize=(13,3.5))
fig5, ax5 = plt.subplots(1, 3, figsize=(13,3.5))

ax2_ylim = [[5.0, 7.0, 1.9], [3.0, 4.0, 2.1]]
ax3_ylim = [[5.0, 6.0, 1.5], [3.0, 3.0, 1.5]]
ax4_ylim = [[5.0, 4.0, 1.1], [3.0, 2.0, 1.1]]
ax5_ylim = [[5.0, 3.0, 0.7], [3.0, 1.5, 0.7]]
 
ax2_nylim = [-0.5, -6, -8]
ax3_nylim = [-0.5, -6, -8]
ax4_nylim = [-0.5 -6, -8]
ax5_nylim = [-0.5, -6, -8]
ax2_xticklabels = ['SSP1-1.9', 'SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
ssp_ticks = [0.5, 1.5, 2.5, 3.5, 4.5]
ax2_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Arctic SAT Anomaly (°C)', 
             'Future Arctic Amplification']
ax3_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Regional SAT Anomaly (°C)', 
             'Future Regional Amplification']
ax4_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Regional SAT Anomaly (°C)', 
             'Future Regional Amplification']
ax5_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Regional SAT Anomaly (°C)', 
             'Future Regional Amplification']

region_names = ['Arctic', 'NH Midlatitudes', 'Tropical', 'SH Extratropical']
suptitle2 = ' SAT Anomalies vs Global Anomaly of CMIP6 Ensemble ('+str(project_start)+' - '+str(project_end)+')'

fig2.suptitle(region_names[0]+suptitle2)
fig3.suptitle(region_names[1]+suptitle2)
fig4.suptitle(region_names[2]+suptitle2)
fig5.suptitle(region_names[3]+suptitle2)
fig2.subplots_adjust(wspace=0.3)
fig3.subplots_adjust(wspace=0.3)
fig4.subplots_adjust(wspace=0.3)
fig5.subplots_adjust(wspace=0.3)

# Set axis titles and ticks for each subplot (Arctic Amplification) (Figure 2)
for n in range(0,len(ax2_ylabels),1):
    if project_end==2050:
        if n!=2:
            ax2[n].set_ylim([-0.5, ax2_ylim[0][n]])
            ax3[n].set_ylim([-0.5, ax3_ylim[0][n]])
            ax4[n].set_ylim([-0.5, ax4_ylim[0][n]])
            ax5[n].set_ylim([-0.5, ax5_ylim[0][n]])
         
            ax2[n].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax2[n].yaxis.set_major_locator(MultipleLocator(1.0))
        # if n==2:
        #     ax2[n].set_ylim([1.3, ax2_ylim[0][n]])
        #     ax3[n].set_ylim([1.2, ax3_ylim[0][n]])
        #     ax4[n].set_ylim([1.0, ax4_ylim[0][n]])
        #     ax5[n].set_ylim([0.4, ax5_ylim[0][n]])
        #     ax2[n].yaxis.set_minor_locator(MultipleLocator(0.01))
        #     ax2[n].yaxis.set_major_locator(MultipleLocator(0.1))
        #     ax3[n].yaxis.set_minor_locator(MultipleLocator(0.01))
        #     ax3[n].yaxis.set_major_locator(MultipleLocator(0.1))
        #     ax4[n].yaxis.set_minor_locator(MultipleLocator(0.01))
        #     ax4[n].yaxis.set_major_locator(MultipleLocator(0.1))
        #     ax2[n].axhline(0, color='grey', linewidth=0.8) 
        #     ax3[n].axhline(0, color='grey', linewidth=0.8) 
        #     ax4[n].axhline(0, color='grey', linewidth=0.8) 
        #     ax5[n].axhline(0, color='grey', linewidth=0.8) 
    else:
        if n!=2:
            ax2[n].set_ylim([-0.5, ax2_ylim[1][n]])
            ax3[n].set_ylim([-0.5, ax3_ylim[1][n]])
            ax4[n].set_ylim([-0.5, ax4_ylim[1][n]])
            ax5[n].set_ylim([-0.5, ax5_ylim[1][n]])
            ax2[n].yaxis.set_minor_locator(MultipleLocator(0.25))
            ax2[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
            ax2[n].axhline(0, color='grey', linewidth=0.8) 
            ax3[n].axhline(0, color='grey', linewidth=0.8) 
            ax4[n].axhline(0, color='grey', linewidth=0.8) 
            ax5[n].axhline(0, color='grey', linewidth=0.8)
        # if n==2:
        #     ax2[n].set_ylim([1.2, ax2_ylim[1][n]])
        #     ax3[n].set_ylim([1.2, ax3_ylim[1][n]])
        #     ax4[n].set_ylim([0.9, ax4_ylim[1][n]])
        #     ax5[n].set_ylim([0.4, ax5_ylim[1][n]])
        #     ax2[n].yaxis.set_minor_locator(MultipleLocator(0.01))
        #     ax2[n].yaxis.set_major_locator(MultipleLocator(0.1))
        #     ax3[n].yaxis.set_minor_locator(MultipleLocator(0.01))
        #     ax3[n].yaxis.set_major_locator(MultipleLocator(0.1))
        #     ax4[n].yaxis.set_minor_locator(MultipleLocator(0.01))
        #     ax4[n].yaxis.set_major_locator(MultipleLocator(0.1))
        #     ax2[n].axhline(0, color='grey', linewidth=0.8) 
        #     ax3[n].axhline(0, color='grey', linewidth=0.8) 
        #     ax4[n].axhline(0, color='grey', linewidth=0.8) 
        #     ax5[n].axhline(0, color='grey', linewidth=0.8) 

    ax2[n].set_xlim([0, 5.0])
    ax2[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                      rotation=20, minor=False)
    ax3[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                  rotation=20, minor=False)
    ax4[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                  rotation=20, minor=False)
    ax5[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                  rotation=20, minor=False)
    
    ax2[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax3[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax3[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax5[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    
    ax2[n].set_ylabel(ax2_ylabels[n])
    ax3[n].set_ylabel(ax3_ylabels[n])
    ax4[n].set_ylabel(ax4_ylabels[n])
    ax5[n].set_ylabel(ax5_ylabels[n])

# Plot the time series for each SSP
for s in range(0, n_ssps, 1): # For each file (current SSP) 
    file = 'dtemp_' + ssps[s] + '_ensemble.nc'
    print(file)
    ds = Dataset(data_dir+file,"r",format="NETCDF4")  # Open dataset
    ds_md = ds.__dict__  # Extract 
    # print(ds.__dict__) # Print metadata
    
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
    xmaxes = [0.19, 0.19, 0.19, 0.19]
    ymaxes = [4.6, 3.7, 2.75, 1.85]
    delta = [0.3, 0.23, 0.18, 0.15]
    for r in range(1, len(regions), 1):
        curr_time_series = time_series[r-1]
        x,y,d = xmaxes[r-1], ymaxes[r-1], delta[r-1]
        ax1[3].text(x=3.15, y=1.8-s*0.15, s=ssp_labels[s], color=colors0[s], 
                    fontsize='medium', fontweight='bold')
        
        ax1[r-1].plot(gm_tas_time_series[project_start_i:project_end_i],
                        curr_time_series[project_start_i:project_end_i],
                        color=colors0[s],
                        linewidth=1)
        lsr = least_squares_fit(gm_tas_time_series[project_start_i:project_end_i],
                                curr_time_series[project_start_i:project_end_i])
        corr, _ = pearsonr(gm_tas_time_series[project_start_i:project_end_i],
                           curr_time_series[project_start_i:project_end_i])
            
        ax1[r-1].text(s='m =  %.2f' % lsr[0] + ' ± %.2f' % lsr[2],
                      x=x, y=y-s*d,
                      color=colors0[s],
                      fontsize=8)
        
        # plot the global anomalies for the amplification plots
        ax2[0].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                       y=gm_tas_time_series_ensemble[:,project_end_i],
                       color=colors0[0], s=10)
        ax2[0].scatter(x=ssp_ticks[s],
                       y=gm_tas_time_series[project_end_i],
                       color='k', s=10)
        
        ax3[0].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                       y=gm_tas_time_series_ensemble[:,project_end_i],
                       color=colors0[0], s=10)
        ax3[0].scatter(x=ssp_ticks[s],
                       y=gm_tas_time_series[project_end_i],
                       color='k', s=10)
        
        ax4[0].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                       y=gm_tas_time_series_ensemble[:,project_end_i],
                       color=colors0[0], s=10)
        ax4[0].scatter(x=ssp_ticks[s],
                       y=gm_tas_time_series[project_end_i],
                       color='k', s=10)
        
        ax5[0].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                       y=gm_tas_time_series_ensemble[:,project_end_i],
                       color=colors0[0], s=10)
        ax5[0].scatter(x=ssp_ticks[s],
                       y=gm_tas_time_series[project_end_i],
                       color='k', s=10)
        
    
    curr_arc_amp = arc_tas_time_series/gm_tas_time_series
    curr_nhml_amp = nhml_tas_time_series/gm_tas_time_series
    curr_tropics_amp = tropics_tas_time_series/gm_tas_time_series
    curr_sh_amp = sh_tas_time_series/gm_tas_time_series
    
    curr_arc_amp_sigma = np.sqrt( np.square(arc_tas_time_series_sigma) +\
                                  np.square(gm_tas_time_series_sigma) )
    curr_nhml_amp_sigma = np.sqrt( np.square(nhml_tas_time_series_sigma) +\
                                  np.square(gm_tas_time_series_sigma) )
    curr_tropics_amp_sigma = np.sqrt( np.square(tropics_tas_time_series_sigma) +\
                                  np.square(gm_tas_time_series_sigma) )
    curr_sh_amp_sigma = np.sqrt( np.square(sh_tas_time_series_sigma) +\
                                  np.square(gm_tas_time_series_sigma) )
        
        
    # Finished cycling through SSP files.
    # Now plot the amplifications
    
    ax2[1].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=arc_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[1], s=10)
    ax2[1].scatter(x=ssp_ticks[s],
                    y=arc_tas_time_series[project_end_i],
                    color='k', s=10)
    ax2[2].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=arc_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[s], s=10)
    ax2[2].scatter(x=ssp_ticks[s],
                    y=arc_tas_time_series[project_end_i]/gm_tas_time_series[project_end_i],
                    color='k', s=10)

    ax3[1].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=nhml_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[2], s=10)
    ax3[1].scatter(x=ssp_ticks[s],
                    y=nhml_tas_time_series[project_end_i],
                    color='k', s=10)
    ax3[2].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=nhml_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[s], s=10)
    ax3[2].scatter(x=ssp_ticks[s],
                    y=nhml_tas_time_series[project_end_i]/gm_tas_time_series[project_end_i],
                    color='k', s=10)    

    ax4[1].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=tropics_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[3], s=10)
    ax4[1].scatter(x=ssp_ticks[s],
                    y=tropics_tas_time_series[project_end_i],
                    color='k', s=10)
    ax4[2].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=tropics_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[s], s=10)
    ax4[2].scatter(x=ssp_ticks[s],
                    y=tropics_tas_time_series[project_end_i]/gm_tas_time_series[project_end_i],
                    color='k', s=10) 
    
    ax5[1].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=sh_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[4], s=10)
    ax5[1].scatter(x=ssp_ticks[s],
                    y=sh_tas_time_series[project_end_i],
                    color='k', s=10)
    ax5[2].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=sh_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[s], s=10)
    ax5[2].scatter(x=ssp_ticks[s],
                    y=sh_tas_time_series[project_end_i]/gm_tas_time_series[project_end_i],
                    color='k', s=10) 
    
fig1.savefig('regional_sat_anomalies_vs_global_amap_ensemble_'+str(project_end)+'.png')
fig2.savefig('arctic_amplification_amap_ensemble_'+str(project_end)+'.png')
fig3.savefig('nhml_amplification_amap_ensemble_'+str(project_end)+'.png')
fig4.savefig('tropics_amplification_amap_ensemble_'+str(project_end)+'.png')
fig5.savefig('sh_amplification_amap_ensemble_'+str(project_end)+'.png')