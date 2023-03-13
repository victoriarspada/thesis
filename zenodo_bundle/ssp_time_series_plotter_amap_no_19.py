# -*- coding: utf-8 -*-
"""
@author: Victoria Spada
SSP Comparisons (Time series)

Same as ssp_time_series_plotter_amap but without SSP1-1.9
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

data_dir = 'AMAP Projections/Renormed Mod1 Ensemble/'
data_dir = 'AMAP Projections/Renormed Ensemble All/'
regions = ['GLOBAL', 'ARCTIC', 'NHML', 'TROPICS', 'SH']
n_regions = len(regions) # includes global

ssps = ['SSP126', 'SSP245', 'SSP370', 'SSP585']
ssp_labels = ['SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
colors0 = ['green', 'blue', 'darkorange', 'red', 'purple']
fname_prefix = 'dtemp_'
n_ssps = len(ssps)

# Arrays for computing region median time series
n_years = 61
project_start, project_end = 2015, 2030
ref_start, ref_end = 1995, 2014
years = np.arange(1990, 2051, 1) # years emulated by the AMAP emulator
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

# Average amplifications for each of the regions (for each of the SSPs)
amplifications = np.zeros((n_regions-1, n_ssps)) 
amplifications_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_ssps)) # std dev of anomalies in ensemble at year

slopes = np.zeros((n_regions-1, n_ssps)) 
slopes_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the slopes

n_ensemble = 24 # number of ensemble members

# # Plot the time series of the global and regional temperatures for each SSP
# fig0, ax0 = plt.subplots(1, 5, figsize=(17,3))
# ax0_ylim = [[2.1, 4.1, 2.1, 2.1, 2.1], [1.1, 2.3, 1.3, 0.8, 0.8]]
# ssp_labels = ['SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
# colors0 = ['green', 'blue', 'darkorange', 'red', 'purple']
# ax0_titles = ['Global mean', 'Arctic mean', 'NH Mid-latitude mean', 
#              'Tropical mean', 'SH extratropical mean']
# ax0_ylabel = 'SAT Anomaly (°C)'
# suptitle0 = 'SAT Anomalies of CMIP6 Ensemble'
# fig0.suptitle(suptitle0, x=0.32, y=1.25)
# fig0.subplots_adjust(left=0.0,
#                     right=0.7,
#                     top=1.1,
#                     wspace=0.4,
#                     hspace=0.2)

# # Set axis titles and ticks for each region (Figure 0)
# ax0[0].set_ylabel(ax0_ylabel)
# for n in range(0,len(ax0_titles),1):
#     if project_end==2050:
#         ax0[n].set_ylim([-0.7,ax0_ylim[0][n]])
#     else:
#         ax0[n].set_ylim([-0.7,ax0_ylim[1][n]])
#     ax0[n].set_xlim([1996,project_end+5])
#     ax0[n].set_title(ax0_titles[n], fontsize='medium')
#     ax0[n].xaxis.set_minor_locator(MultipleLocator(2))
#     if project_end==2050:
#         ax0[n].xaxis.set_major_locator(MultipleLocator(15))
#     else:
#         ax0[n].xaxis.set_major_locator(MultipleLocator(10))
#     ax0[n].yaxis.set_minor_locator(MultipleLocator(0.1))
#     ax0[n].yaxis.set_major_locator(MultipleLocator(0.5))
#     ax0[n].tick_params(direction="in", which="both", top=True, 
#                       bottom=True, left=True, right=True)
#     ax0[n].axhline(0, color='black', linewidth=0.8) 

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

fig6, ax6 = plt.subplots(1, 5, figsize=(17,3.5))
fig6.subplots_adjust(wspace=0.3)

ax2_ylim = [[4.0, 7.0, 1.9], [3.0, 4.0, 2.1]]
ax3_ylim = [[4.0, 5.0, 1.5], [3.0, 3.0, 1.5]]
ax4_ylim = [[4.0, 4.0, 1.1], [3.0, 2.0, 1.1]]
ax5_ylim = [[4.0, 2.5, 0.7], [3.0, 1.5, 0.7]]
 
ax2_nylim = [-0.2, -6, -8]
ax3_nylim = [-0.2, -6, -8]
ax4_nylim = [-0.2 -6, -8]
ax5_nylim = [-0.2, -6, -8]
ax2_xticklabels = ['SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
ssp_ticks = [0.5, 1.5, 2.5, 3.5]
ax2_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Arctic SAT Anomaly (°C)', 
             'Future Arctic Amplification']
ax3_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Regional SAT Anomaly (°C)', 
             'Future Regional Amplification']
ax4_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Regional SAT Anomaly (°C)', 
             'Future Regional Amplification']
ax5_ylabels = ['Future Global SAT Anomaly (°C)', 'Future Regional SAT Anomaly (°C)', 
             'Future Regional Amplification']
ax6_ylabels = ['Future Global SAT Anomaly (°C)', 'Arctic Amplification', 'NHML Amplification',
               'Tropical Amplification', 'SH Amplification']

region_names = ['Arctic', 'NH Midlatitudes', 'Tropical', 'SH Extratropical']
suptitle2 = ' SAT Anomalies vs Global Anomaly of AMAP Ensemble ('+str(project_start)+' - '+str(project_end)+')'
suptitle3 = 'AMAP Emulator Regional SAT Amplification ('+str(project_start)+' - '+str(project_end)+')'

fig2.suptitle(region_names[0]+suptitle2)
fig3.suptitle(region_names[1]+suptitle2)
fig4.suptitle(region_names[2]+suptitle2)
fig5.suptitle(region_names[3]+suptitle2)
fig6.suptitle(suptitle3)

fig2.subplots_adjust(wspace=0.3)
fig3.subplots_adjust(wspace=0.3)
fig4.subplots_adjust(wspace=0.3)
fig5.subplots_adjust(wspace=0.3)

for n in range(0,len(ax6_ylabels),1):
    # ax6[n].set_ylim([-0.5, ax6_ylim[0][n]])
    # if n!=0:     
    #     ax6[n].yaxis.set_minor_locator(MultipleLocator(0.01))
    #     ax6[n].yaxis.set_major_locator(MultipleLocator(0.025))
    if n==0:
        ax6[n].yaxis.set_minor_locator(MultipleLocator(0.1))
        ax6[n].yaxis.set_major_locator(MultipleLocator(0.5))
    # ax6[0].axhline(0, color='grey') 
    ax6[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
              rotation=20, minor=False)

    ax6[n].tick_params(direction="in", top=True, 
                  bottom=True, left=True, right=True)
    ax6[n].set_ylabel(ax6_ylabels[n])

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
            
            ax2[n].axhline(0, color='grey') 
            ax3[n].axhline(0, color='grey') 
            ax4[n].axhline(0, color='grey') 
            ax5[n].axhline(0, color='grey')
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
            ax2[n].axhline(0, color='grey') 
            ax3[n].axhline(0, color='grey') 
            ax4[n].axhline(0, color='grey') 
            ax5[n].axhline(0, color='grey')
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

    ax2[n].set_xlim([0, 4.0])
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
        
        # print(s, ' NH A ',
        #       "{:.2f}".format(tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]), 
        #       ' ± ',
        #       "{:.2f}".format(uncertainty) )
        
        
        x,y,d = xmaxes[r-1], ymaxes[r-1], delta[r-1]
        ax1[3].text(x=3.15, y=1.8-s*0.15, s=ssp_labels[s], color=colors0[s], 
                    fontsize='medium')
        
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
        slopes[r-1,s] = lsr[0]
        slopes_sigma[r-1,s] = lsr[2]
        # print('%.2f' % lsr[0] + ' ± %.2f' % lsr[2])
        
        # plot the global anomalies for the amplification plots
        ax2[0].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                       y=gm_tas_time_series_ensemble[:,project_end_i],
                       color=colors0[0], s=10)
        ax2[0].scatter(x=ssp_ticks[s],
                       y=gm_tas_time_series[project_end_i],
                       color='k', s=10)
        
        ax6[0].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                       y=gm_tas_time_series_ensemble[:,project_end_i],
                       color=colors0[0], s=10)
        ax6[0].scatter(x=ssp_ticks[s],
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
        
        a = curr_time_series[project_start_i:project_end_i]/gm_tas_time_series[project_start_i:project_end_i]
        # print(ax2_xticklabels[s], ' A ',
        #           "{:.2f}".format(np.mean(a)), 
        #           ' ± ',
        #           "{:.2f}".format(np.std(a)))
        
    
    # Mean amplifications for each year (across the ensemble)
    # curr_arc_amp = arc_tas_time_series/gm_tas_time_series
    # curr_nhml_amp = nhml_tas_time_series/gm_tas_time_series
    # curr_tropics_amp = tropics_tas_time_series/gm_tas_time_series
    # curr_sh_amp = sh_tas_time_series/gm_tas_time_series
    
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
    # curr_arc_amp_sigma = np.sqrt( np.square(arc_tas_time_series_sigma/arc_tas_time_series) +\
    #                               np.square(gm_tas_time_series_sigma/gm_tas_time_series) )
    # curr_nhml_amp_sigma = np.sqrt( np.square(nhml_tas_time_series_sigma/nhml_tas_time_series_sigma) +\
    #                               np.square(gm_tas_time_series_sigma/gm_tas_time_series) )
    # curr_tropics_amp_sigma = np.sqrt( np.square(tropics_tas_time_series_sigma/tropics_tas_time_series) +\
    #                               np.square(gm_tas_time_series_sigma/gm_tas_time_series) )
    # curr_sh_amp_sigma = np.sqrt( np.square(sh_tas_time_series_sigma/sh_tas_time_series) +\
    #                               np.square(gm_tas_time_series_sigma/gm_tas_time_series_sigma) )
    
    # # uncertainty of ensemble mean amplification    
    # amplifications_sigma[0,s] = curr_arc_amp_sigma[project_end_i]
    # amplifications_sigma[1,s] = curr_nhml_amp_sigma[project_end_i]
    # amplifications_sigma[2,s] = curr_tropics_amp_sigma[project_end_i]
    # amplifications_sigma[3,s] = curr_sh_amp_sigma[project_end_i]
        
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
    ax6[1].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                y=arc_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                color=colors0[s], s=10)
    ax6[1].scatter(x=ssp_ticks[s],
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
    ax6[2].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                y=nhml_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                color=colors0[s], s=10)
    ax6[2].scatter(x=ssp_ticks[s],
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
    ax6[3].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=tropics_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[s], s=10)
    ax6[3].scatter(x=ssp_ticks[s],
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
    ax6[4].scatter(x=ssp_ticks[s]*np.ones(n_ensemble),
                    y=sh_tas_time_series_ensemble[:,project_end_i]/gm_tas_time_series_ensemble[:,project_end_i],
                    color=colors0[s], s=10)
    ax6[4].scatter(x=ssp_ticks[s],
                    y=sh_tas_time_series[project_end_i]/gm_tas_time_series[project_end_i],
                    color='k', s=10) 


# Compute average ampliifcations for each region (across the SSPs)
avg_amp, amp_std = np.zeros(n_regions-1), np.zeros(n_regions-1)
amp_uncertainties = np.zeros(n_regions-1)

avg_m, avg_m_std = np.zeros(n_regions-1), np.zeros(n_regions-1)
avg_m_uncertainties = np.zeros(n_regions-1)
for r in range(0, n_regions-1, 1):
    avg_amp[r] = np.mean(amplifications[r,:]) # Mean amplification for the region
    amp_std[r] = np.std(amplifications[r,:])
    amp_uncertainties[r] = (1/n_ssps)*np.sqrt(
          np.sum( np.square(amplifications_std[r,:]) )
        ) 
    print('mean amplification = {:.2f}'.format(avg_amp[r]), 
          '\nstddev = {:.2f}'.format(amp_std[r]), 
          '\nuncertainty = {:.2f} \n'.format(amp_uncertainties[r]))  
  
    avg_m[r] = np.mean(slopes[r,:]) # Mean amplification for the region
    avg_m_std[r] = np.std(slopes[r,:])
    avg_m_uncertainties[r] = (1/n_ssps)*np.sqrt(
          np.sum( np.square(slopes_sigma[r,:]) ) # propagate uncertainty (uncertainty is the std of each ensemble mean amp)
        ) 
    # print('mean slope = {:.2f}'.format(avg_m[r]), 
    #       '\nstddev = {:.2f}'.format(avg_m_std[r]), 
    #       '\nuncertainty = {:.2f} \n'.format(avg_m_uncertainties[r]))
    
fig1.savefig('regional_sat_anomalies_vs_global_amap_ensemble_'+str(project_end)+'.png')
fig2.savefig('arctic_amplification_amap_ensemble_'+str(project_end)+'.png')
fig3.savefig('nhml_amplification_amap_ensemble_'+str(project_end)+'.png')
fig4.savefig('tropics_amplification_amap_ensemble_'+str(project_end)+'.png')
fig5.savefig('sh_amplification_amap_ensemble_'+str(project_end)+'.png')
fig6.savefig('amplification_amap_ensemble_'+str(project_end)+'.png')