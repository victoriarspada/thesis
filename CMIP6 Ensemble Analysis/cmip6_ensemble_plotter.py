# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

This plotter takes a CMIP6 ensemble for 4 SSPs and creates 6 figures.

A reduced major axis fit is performed for each of the SSPs 
regionl versus global SAT anomalies over 2015-2050, and the results are plotted
for each of the 4 regions.

The spread of Arctic, NHML, Tropical, and SH SAT anomalies and regional 
amplifications for each SSP are also plotted for the year 2050.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from netCDF4 import Dataset # http://unidata.github.io/netcdf4-python//
from regression import *
from rma import * 
from scipy.stats import pearsonr
import matplotlib.gridspec as gridspec

data_dir = 'cmip6_data/'
ssp_dirs = ['DATA_ssp126/', 'DATA_ssp245/', 'DATA_ssp370/', 'DATA_ssp585/']
n_ssps = len(ssp_dirs)

n_years = 451
project_start, project_end = 2015, 2050
ref_start, ref_end = 1995, 2014

ssp_labels = ['SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
colors0 = ['rebeccapurple', 'seagreen', 'darkgoldenrod', 'crimson', 'seagreen']


# Plot regional SAT anomalies vs the global SAT anomalies (Figure 1)
fig1 = plt.figure(figsize=(6,6))
gs = gridspec.GridSpec(2,2)
axs1 = fig1.add_subplot(gs[0, 0], )
axs2 = fig1.add_subplot(gs[0, 1])
axs3 = fig1.add_subplot(gs[1, 0])
axs4 = fig1.add_subplot(gs[1, 1])
gs.update(wspace=0.4, hspace=0.4)
ax1 = [axs1, axs2, axs3, axs4]

ax1_ylim = [[4.0, 2.3, 1.6, 1.4], [2.25, 1.25, 0.80, 0.75]]
ax1_xlabels = ['Arctic SAT Anomaly (°C)', 'NHML SAT Anomaly (°C)', 
             'Tropical SAT Anomaly (°C)', 'SH SAT Anomaly (°C)']
ax1_ylabel = 'Global SAT Anomaly (°C)'
suptitle1 = 'Regional SAT Anomalies vs Global Anomaly of CMIP6 Ensemble ('+str(project_start)+' - '+str(project_end)+')\n Relative to 1995-2014'
fig1.suptitle(suptitle1, y=0.98)
fig1.subplots_adjust(wspace=0.3)
    
# Set axis titles and ticks for each region (Figure 1)
for n in range(0,len(ax1_xlabels),1):
    ax1[n].set_xlabel(ax1_ylabel)
    ax1[n].set_ylabel(ax1_xlabels[n], fontsize='medium')
    if project_end==2050:
        ax1[n].set_xlim([0,1.74])
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
            ax1[n].yaxis.set_major_locator(MultipleLocator(0.05))
        if n==3:
            ax1[n].yaxis.set_minor_locator(MultipleLocator(0.005))
            ax1[n].yaxis.set_major_locator(MultipleLocator(0.01))
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

fig6 = plt.figure(figsize=(11,7), tight_layout=True)
gs = gridspec.GridSpec(2, 6)
gs.update(wspace=0.8, hspace=0.3)
axs1 = fig6.add_subplot(gs[0, :2], )
axs2 = fig6.add_subplot(gs[0, 2:4])
axs3 = fig6.add_subplot(gs[0, 4:])
axs4 = fig6.add_subplot(gs[1, 1:3])
axs5 = fig6.add_subplot(gs[1, 3:5])
ax6 = [axs1, axs2, axs3, axs4, axs5]

ax2_ylim = [[4.0, 8.0, 9.0], [3.0, 6.0, 9.0]] 
ax2_nylim = [-1, -5, -8]
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
suptitle2 = ' SAT Anomalies vs Global Anomaly of CMIP6 Ensemble ('+str(project_end)+')'
suptitle3 = 'CMIP6 Ensemble Regional SAT Amplification ('+str(project_end)+')'

fig2.suptitle(region_names[0]+suptitle2)
fig3.suptitle(region_names[1]+suptitle2)
fig4.suptitle(region_names[2]+suptitle2)
fig5.suptitle(region_names[3]+suptitle2)
fig6.suptitle(suptitle3, y=0.93, fontsize='large')
fig2.subplots_adjust(wspace=0.3)
fig3.subplots_adjust(wspace=0.3)
fig4.subplots_adjust(wspace=0.3)
fig5.subplots_adjust(wspace=0.3)

# Set axis titles and ticks for each subplot (Arctic Amplification) (Figure 2)
for n in range(0,len(ax2_ylabels),1):
    if project_end==2050:
        ax2[n].set_ylim([-1.0, ax2_ylim[0][n]])
        if n==2:
            ax2[n].set_ylim([-4.0, ax2_ylim[0][n]])
            ax6[1].set_ylim([-4.0, ax2_ylim[0][n]])
        ax2[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax2[n].yaxis.set_major_locator(MultipleLocator(1.0))
    else:
        ax2[n].set_ylim([ax2_nylim[n], ax2_ylim[1][n]])
        ax6[1].set_ylim([ax2_nylim[n], ax2_ylim[1][n]])
        ax2[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax2[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==1:
            ax2[n].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax2[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==2:
           ax2[n].yaxis.set_minor_locator(MultipleLocator(1.0))
           ax2[n].yaxis.set_major_locator(MultipleLocator(2.0))  

    ax2[n].set_xlim([0, 4.0])
    ax2[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                      rotation=20, minor=False)
    ax2[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax2[n].set_ylabel(ax2_ylabels[n])
    ax2[n].axhline(0, color='grey', linewidth=0.8)   
    
# Set axis titles and ticks for each subplot (NH Amplification) (Figure 3)
ax3_ylim = [[4.0, 4.0, 5.0], [3.0, 3.0, 8.0]] 
ax3_nylim = [-1, -1, -9]
for n in range(0,len(ax3_ylabels),1):
    if project_end==2050:
        ax3[n].set_ylim([-1.0, ax3_ylim[0][n]])
        if n==2:
            ax3[n].set_ylim([-1.0, ax3_ylim[0][n]])
            ax6[2].set_ylim([-1.0, ax3_ylim[0][n]])
        ax3[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax3[n].yaxis.set_major_locator(MultipleLocator(1.0))
    else:
        ax3[n].set_ylim([ax3_nylim[n], ax3_ylim[1][n]])
        ax6[2].set_ylim([ax3_nylim[n], ax3_ylim[1][n]])
        ax3[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax3[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==1:
            ax3[n].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax3[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==2:
            ax3[n].yaxis.set_minor_locator(MultipleLocator(1.0))
            ax3[n].yaxis.set_major_locator(MultipleLocator(2.0))  

    ax3[n].set_xlim([0, 4.0])
    ax3[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                      rotation=20, minor=False)
    ax3[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax3[n].set_ylabel(ax3_ylabels[n])
    ax3[n].axhline(0, color='grey', linewidth=0.8)      
    
# Set axis titles and ticks for each subplot (Tropical Amplification) (Figure 4)
ax4_ylim = [[4.0, 4.0, 3.0], [3.0, 3.0, 4.0]] 
ax4_nylim = [-1, -1, -3]
for n in range(0,len(ax4_ylabels),1):
    if project_end==2050:
        ax4[n].set_ylim([-1.0, ax4_ylim[0][n]])
        if n==2:
            ax4[n].set_ylim([-3.0, ax4_ylim[0][n]])
            ax6[3].set_ylim([-3.0, ax4_ylim[0][n]])
        ax4[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax4[n].yaxis.set_major_locator(MultipleLocator(1.0))
    else:
        ax4[n].set_ylim([ax4_nylim[n], ax4_ylim[1][n]])
        ax6[3].set_ylim([ax4_nylim[n], ax4_ylim[1][n]])
        ax4[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax4[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==1:
            ax4[n].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax4[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==2:
            ax4[n].yaxis.set_minor_locator(MultipleLocator(1.0))
            ax4[n].yaxis.set_major_locator(MultipleLocator(2.0))  

    ax4[n].set_xlim([0, 4.0])
    ax4[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                      rotation=20, minor=False)
    ax4[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax4[n].set_ylabel(ax4_ylabels[n])
    ax4[n].axhline(0, color='grey', linewidth=0.8)   

# Set axis titles and ticks for each subplot (SH Amplification) (Figure 5)
ax5_ylim = [[4.0, 4.5, 5.0], [3.0, 4.0, 6.0]] 
ax5_nylim = [-1, -2, -2]
for n in range(0,len(ax5_ylabels),1):
    if project_end==2050:
        ax5[n].set_ylim([-1.0, ax5_ylim[0][n]])
        if n==2 or n==1:
            ax5[n].set_ylim([-2, ax5_ylim[0][n]])
            ax6[4].set_ylim([-2, ax5_ylim[0][n]])
        ax5[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax5[n].yaxis.set_major_locator(MultipleLocator(1.0))
    else:
        ax5[n].set_ylim([ax5_nylim[n], ax5_ylim[1][n]])
        ax6[4].set_ylim([ax5_nylim[n], ax5_ylim[1][n]])
        ax5[n].yaxis.set_minor_locator(MultipleLocator(0.25))
        ax5[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==1:
            ax5[n].yaxis.set_minor_locator(MultipleLocator(0.5))
            ax5[n].yaxis.set_major_locator(MultipleLocator(1.0)) 
        if n==2:
            ax5[n].yaxis.set_minor_locator(MultipleLocator(1.0))
            ax5[n].yaxis.set_major_locator(MultipleLocator(2.0))  

    ax5[n].set_xlim([0, 4.0])
    ax5[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
                      rotation=20, minor=False)
    ax5[n].tick_params(direction="in", top=True, 
                      bottom=True, left=True, right=True)
    ax5[n].set_ylabel(ax5_ylabels[n])
    ax5[n].axhline(0, color='grey', linewidth=0.8)   
    
for n in range(0,len(ax6_ylabels),1):
    # ax6[n].set_ylim([-0.5, ax6_ylim[0][n]])
    if n==0:
        ax6[n].yaxis.set_minor_locator(MultipleLocator(0.1))
        ax6[n].yaxis.set_major_locator(MultipleLocator(0.5))
        if project_end==2030:
            ax6[n].set_ylim([-1,3])
        else:
            ax6[n].set_ylim([-0.5,3.5])
    ax6[n].set_xticks(ticks=ssp_ticks, labels=ax2_xticklabels, 
              rotation=20, minor=False)

    ax6[n].tick_params(direction="in", top=True, 
                  bottom=True, left=True, right=True)
    ax6[n].set_ylabel(ax6_ylabels[n])  
    ax6[n].grid(color='gainsboro')
    

region_prefixes = ['gma', 'arc', 'nhm', 'tro', 'she']     
n_regions = len(region_prefixes) # includes global

# Arrays for computing region median time series
years = np.arange(1850, 2301, 1)
li_years = list(years)
project_start_i, project_end_i = li_years.index(project_start), li_years.index(project_end) 
ref_start_i, ref_end_i = li_years.index(ref_start), li_years.index(ref_end)

gm_anomalies = np.zeros((n_ssps+1, n_years)) # mean anomalies for each ssp
arc_anomalies = np.zeros((n_ssps+1, n_years))
nh_anomalies = np.zeros((n_ssps+1, n_years))
trop_anomalies = np.zeros((n_ssps+1, n_years))
sh_anomalies = np.zeros((n_ssps+1, n_years))

# Average amplifications for each of the regions (for each of the SSPs)
amplifications = np.zeros((n_regions-1, n_ssps)) 
amplifications_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the amplification
amplifications_std = np.zeros((n_regions-1, n_ssps)) # std dev of amps in ensemble at year

slopes = np.zeros((n_regions-1, n_ssps)) 
slopes_sigma = np.zeros((n_regions-1, n_ssps)) # Propagated uncertainty for the slopes

# Plot the time series for each SSP
for s in range(0, n_ssps, 1):
    curr_data_dir = data_dir+ssp_dirs[s]
    curr_files = os.listdir(curr_data_dir)
    curr_files.remove('CanESM5') # folder with CanESM files copied
    curr_files.remove('Extras') # folders with duplicates
    # print(curr_files)
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
        # print(ds.__dict__) # Print metadata
        
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
    
    # Take the MEANS, each row is a region
    tas_time_series_mean = np.zeros((n_regions, n_years)) # Define mean temp. matrix
    tas_mean_anomaly = np.zeros((n_regions, n_years)) # Define mean temp. anomaly matrix
    sigma = np.zeros((n_regions, n_years))
    max_tas_anomaly = np.zeros((n_regions, n_years))
    min_tas_anomaly = np.zeros((n_regions, n_years))
    
    # Finished cycling through SSP files.
    mean_refs = np.zeros(n_regions) # Define array for regional mean temperatures of ref. period
    for r in range(0, n_regions, 1):
        tas_time_series_mean[r, :] = tas_time_series[r, :]/tas_count[r]
        # Compute the MEAN REFERENCE temperature for region r
        # Take mean value over reference period
        mean_refs[r] = np.mean(tas_time_series[r, ref_start_i:ref_end_i]/tas_count[r]) 
        ref = np.mean(tas_time_series[r, ref_start_i:ref_end_i]/tas_count[r]) 
        # Compute mean time series SAT anomaly for each region
        tas_mean_anomaly[r,:] = tas_time_series_mean[r,:] - ref
        # compute standard deviation of the anomaly for each region
        if r == 0:
            gm_ensemble_anomaly =  gm_tas_time_series - ref
            gm_anomalies[s,:] = tas_mean_anomaly[r,:]
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(gm_tas_time_series[:, y] - ref)
                max_tas_anomaly[r, y] = max(gm_tas_time_series[:, y] - ref)
                min_tas_anomaly[r, y] = min(gm_tas_time_series[:, y] - ref)
        if r == 1:
            arc_ensemble_anomaly =  arc_tas_time_series - ref
            arc_anomalies[s,:] = tas_mean_anomaly[r,:]
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(arc_tas_time_series[:, y] - ref)
                max_tas_anomaly[r, y] = max(arc_tas_time_series[:, y] - ref)
                min_tas_anomaly[r, y] = min(arc_tas_time_series[:, y] - ref)
        if r == 2:
            nh_ensemble_anomaly =  nh_tas_time_series - ref
            nh_anomalies[s,:] = tas_mean_anomaly[r,:]
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(nh_tas_time_series[:, y] - ref)
                max_tas_anomaly[r, y] = max(nh_tas_time_series[:, y] - ref)
                min_tas_anomaly[r, y] = min(nh_tas_time_series[:, y] - ref)
        if r == 3:
            trop_ensemble_anomaly =  trop_tas_time_series - ref
            trop_anomalies[s,:] = tas_mean_anomaly[r,:]
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(trop_tas_time_series[:, y] - ref)
                max_tas_anomaly[r, y] = max(trop_tas_time_series[:, y] - ref)
                min_tas_anomaly[r, y] = min(trop_tas_time_series[:, y] - ref)
        if r == 4:
            sh_ensemble_anomaly =  sh_tas_time_series - ref
            sh_anomalies[s,:] = tas_mean_anomaly[r,:]
            for y in range(0, n_years, 1): 
                sigma[r, y] = np.std(sh_tas_time_series[:, y] - ref)
                max_tas_anomaly[r, y] = max(sh_tas_time_series[:, y] - ref)
                min_tas_anomaly[r, y] = min(sh_tas_time_series[:, y] - ref)


    # Do the plotting now for the current SSP 
    for r in range(0, n_regions, 1): 
        # Plot the CMIP6 Regional vs Global temperature anomalies
        if r==1:
            ax1[r-1].plot(gm_anomalies[s, project_start_i:project_end_i],
                        arc_anomalies[s, project_start_i:project_end_i],
                        color=colors0[s],
                        linewidth=1)
            lsr = least_squares_fit(gm_anomalies[s, project_start_i:project_end_i],
                        arc_anomalies[s, project_start_i:project_end_i])
            
            m,b,sigma_m,sigma_b = reduced_major_axis(gm_anomalies[s, project_start_i:project_end_i],
                        arc_anomalies[s, project_start_i:project_end_i])
            corr, _ = pearsonr(gm_anomalies[s, project_start_i:project_end_i],
                                arc_anomalies[s, project_start_i:project_end_i])
            # print(s,'Arctic Pearsons correlation: %.2f' % corr)
            
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,     
                          x=0.12, y=3.6-s*0.29,
                          color=colors0[s],
                          fontsize=8)
            slopes[r-1,s] = m #lsr[0]
            slopes_sigma[r-1,s] = sigma_m #lsr[2]
        if r==2:
            ax1[r-1].plot(gm_anomalies[s, project_start_i:project_end_i],
                        nh_anomalies[s, project_start_i:project_end_i],
                        color=colors0[s],
                        linewidth=1)
            lsr = least_squares_fit(gm_anomalies[s, project_start_i:project_end_i],
                        nh_anomalies[s, project_start_i:project_end_i])
            m,b,sigma_m,sigma_b = reduced_major_axis(gm_anomalies[s, project_start_i:project_end_i],
                        nh_anomalies[s, project_start_i:project_end_i])
            
            corr, _ = pearsonr(gm_anomalies[s, project_start_i:project_end_i],
                                nh_anomalies[s, project_start_i:project_end_i])
            # print(s,'NH Pearsons correlation: %.2f' % corr)
            
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,
                          x=0.12, y=2.05-s*0.17,
                          color=colors0[s],
                          fontsize=8)
            slopes[r-1,s] = lsr[0]
            slopes_sigma[r-1,s] = lsr[2]
        if r==3:
            ax1[r-1].plot(gm_anomalies[s, project_start_i:project_end_i],
                        trop_anomalies[s, project_start_i:project_end_i],
                        color=colors0[s],
                        linewidth=1)
            lsr = least_squares_fit(gm_anomalies[s, project_start_i:project_end_i],
                        trop_anomalies[s, project_start_i:project_end_i])
            m,b,sigma_m,sigma_b = reduced_major_axis(gm_anomalies[s, project_start_i:project_end_i],
                        trop_anomalies[s, project_start_i:project_end_i])
            
            corr, _ = pearsonr(gm_anomalies[s, project_start_i:project_end_i],
                                trop_anomalies[s, project_start_i:project_end_i])
            # print(s,'Trop Pearsons correlation: %.2f' % corr)
            
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,
                          x=0.12, y=1.46-s*0.13,
                          color=colors0[s],
                          fontsize=8)
            slopes[r-1,s] = lsr[0]
            slopes_sigma[r-1,s] = lsr[2]
        if r==4:
            ax1[r-1].plot(gm_anomalies[s, project_start_i:project_end_i],
                        sh_anomalies[s, project_start_i:project_end_i],
                        color=colors0[s],
                        linewidth=1)
            lsr = least_squares_fit(gm_anomalies[s, project_start_i:project_end_i],
                        sh_anomalies[s, project_start_i:project_end_i])
            m,b,sigma_m,sigma_b = reduced_major_axis(gm_anomalies[s, project_start_i:project_end_i],
                        sh_anomalies[s, project_start_i:project_end_i])
            corr, _ = pearsonr(sh_anomalies[s, project_start_i:project_end_i],
                                nh_anomalies[s, project_start_i:project_end_i])
            # print(s,'SH Pearsons correlation: %.2f' % corr)
            
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,
                          x=0.12, y=1.26-s*0.102,
                          color=colors0[s],
                          fontsize=8)
            slopes[r-1,s] = lsr[0]
            slopes_sigma[r-1,s] = lsr[2]
        if r>0:
            qqq= True
            print(' %.2f' % lsr[0] + ' ± %.2f' % lsr[2])
        # Plots of regional amplification until 2050
        if r==0:
            # plot the spread for the SSP in the current region
            a = gm_tas_time_series[:, project_end_i].shape[0]
            aa = gm_tas_time_series[:, project_end_i]-mean_refs[r]
            # print(s, r, tas_mean_anomaly[r, project_end_i], np.mean(gm_tas_time_series[:, project_end_i]-mean_refs[r]))
            
            # scatter T_2050 - T_ref
            ax2[r].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=gm_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
            ax3[r].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=gm_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
            ax4[r].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=gm_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
            ax5[r].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=gm_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
            ax6[r].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=gm_tas_time_series[:,project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
        
            # plot the average for the SSP in the current region
            ax2[r].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')
            ax3[r].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')
            ax4[r].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')
            ax5[r].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')
            ax6[r].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k', s=10)
        if r==1: # Arctic amplification
            # plot the spread for the SSP in the current region
            a = gm_tas_time_series[:, project_end_i].shape[0]
            bb = arc_tas_time_series[:, project_end_i]-mean_refs[r]
            # print(s, r, tas_mean_anomaly[r, project_end_i], tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i])
            ax2[1].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=arc_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
        
            # plot the average for the SSP in the current region
            ax2[1].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')        
            
            # plot the spread for the SSP in the current region
            ax2[2].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           # y=(arc_tas_time_series[:, project_end_i]-mean_refs[r])/tas_mean_anomaly[r, project_end_i],
                           color=colors0[s], s=10)
            
            # plot the average for the SSP in the current region
            ax2[2].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k') 
            
            # Plotting for figure 6
            ax6[1].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            ax6[1].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k')
            
            mean_amp = tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]
            uncertainty = ( (sigma[r, project_end_i]/tas_mean_anomaly[r, project_end_i])**2 + \
                           (sigma[0, project_end_i]/tas_mean_anomaly[0, project_end_i])**2 )**0.5
            # print(ax2_xticklabels[s], ' ARC A ',
            #       "{:.2f}".format(tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]), 
            #       ' ± ',
            #       "{:.2f}".format(uncertainty) )
            
        if r==2: # NH amplification
            # plot the spread for the SSP in the current region
            a = gm_tas_time_series[:, project_end_i].shape[0]
            bb = nh_tas_time_series[:, project_end_i]-mean_refs[r]
            ax3[1].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=nh_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
        
            # plot the average for the SSP in the current region
            ax3[1].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')        
            
            # plot the spread for the SSP in the current region
            ax3[2].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            
            # plot the average for the SSP in the current region
            ax3[2].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k') 
            
            # Plotting for figure 6
            ax6[2].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            ax6[2].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k')
            
            mean_amp = tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]
            uncertainty = ( (sigma[r, project_end_i]/tas_mean_anomaly[r, project_end_i])**2 + \
                            (sigma[0, project_end_i]/tas_mean_anomaly[0, project_end_i])**2 )**0.5
            # print(ax2_xticklabels[s], ' NH A ',
            #       "{:.2f}".format(tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]), 
            #       ' ± ',
            #       "{:.2f}".format(uncertainty) )
            
        if r==3: # Tropical amplification
            # plot the spread for the SSP in the current region
            a = gm_tas_time_series[:, project_end_i].shape[0]
            bb = trop_tas_time_series[:, project_end_i]-mean_refs[r]
            ax4[1].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=trop_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
        
            # plot the average for the SSP in the current region
            ax4[1].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')        
            
            # plot the spread for the SSP in the current region
            ax4[2].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            
            # Plotting for figure 6
            ax6[3].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            ax6[3].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k')
            
            # plot the average for the SSP in the current region
            ax4[2].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k') 
            mean_amp = tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]
            uncertainty = ( (sigma[r, project_end_i]/tas_mean_anomaly[r, project_end_i])**2 + \
                            (sigma[0, project_end_i]/tas_mean_anomaly[0, project_end_i])**2 )**0.5
            # print(ax2_xticklabels[s], ' TROP A ',
            #       "{:.2f}".format(tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]), 
            #       ' ± ',
            #       "{:.2f}".format(uncertainty) )
            
            
        if r==4: # SH Extratropical amplification
            # plot the spread for the SSP in the current region
            a = gm_tas_time_series[:, project_end_i].shape[0]
            bb = sh_tas_time_series[:, project_end_i]-mean_refs[r]
            ax5[1].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=sh_tas_time_series[:, project_end_i]-mean_refs[r],
                           color=colors0[s], s=10)
        
            # plot the average for the SSP in the current region
            ax5[1].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i],
                           color='k')        
            
            # plot the spread for the SSP in the current region
            ax5[2].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            
            # Plotting for figure 6
            ax6[4].scatter(x=ssp_ticks[s]*np.ones(a),
                           y=np.divide(bb,aa),
                           color=colors0[s], s=10)
            ax6[4].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k')
            
            # plot the average for the SSP in the current region
            ax5[2].scatter(x=ssp_ticks[s],
                           y=tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i],
                           color='k') 
            mean_amp = tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]
            uncertainty = ( (sigma[r, project_end_i]/tas_mean_anomaly[r, project_end_i])**2 + \
                            (sigma[0, project_end_i]/tas_mean_anomaly[0, project_end_i])**2 )**0.5
            # print(ax2_xticklabels[s], ' SH A ',
            #       "{:.2f}".format(tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]), 
            #       ' ± ',
            #       "{:.2f}".format(uncertainty) + '\n' )
            
        if r>0:
            amplifications[r-1,s] = tas_mean_anomaly[r, project_end_i]/tas_mean_anomaly[0, project_end_i]
            amplifications_sigma[r-1,s] = uncertainty

    # Add SSP text labels
    if project_end==2050:
        ax1[1].text(x=1.8, y=2.1-s*0.22, s=ssp_labels[s], color=colors0[s], 
                    fontsize='medium') #, fontweight='bold')
    else:
        ax1[3].text(x=1.05, y=0.7-s*0.1, s=ssp_labels[s], color=colors0[s], fontsize='medium', fontweight='bold')

# Compute average ampliifcations for each region (across the SSPs)
avg_amp, amp_std = np.zeros(n_regions-1), np.zeros(n_regions-1)
amp_uncertainties = np.zeros(n_regions-1)

avg_m, avg_m_std = np.zeros(n_regions-1), np.zeros(n_regions-1)
avg_m_uncertainties = np.zeros(n_regions-1)
for r in range(0, n_regions-1, 1):
    avg_amp[r] = np.mean(amplifications[r,:]) # Mean amplification for the region
    amp_std[r] = np.std(amplifications[r,:])
    amp_uncertainties[r] = (1/n_ssps)*np.sqrt(
          np.sum( np.square(amplifications_sigma[r,:]) ) # propagate uncertainty (uncertainty is the std of each ensemble mean amp)
        ) 
    # print(r, 'mean amplification = {:.2f}'.format(avg_amp[r]), 
    #       '\nstddev = {:.2f}'.format(amp_std[r]), 
    #       '\nuncertainty = {:.2f} \n'.format(amp_uncertainties[r]))
    
    avg_m[r] = np.mean(slopes[r,:]) # Mean amplification for the region
    avg_m_std[r] = np.std(slopes[r,:])
    avg_m_uncertainties[r] = (1/n_ssps)*np.sqrt(
          np.sum( np.square(slopes_sigma[r,:]) ) # propagate uncertainty (uncertainty is the std of each ensemble mean amp)
        ) 
    print(r, 'mean slope = {:.2f}'.format(avg_m[r]), 
          '\nstddev = {:.2f}'.format(avg_m_std[r]), 
          '\nuncertainty = {:.2f} \n'.format(avg_m_uncertainties[r]))
    

fig1.savefig('./regional_sat_anomalies_vs_global_cmip6_ensemble_'+str(project_end)+'.png')
fig2.savefig('./arctic_amplification_cmip6_ensemble_'+str(project_end)+'.png')
fig3.savefig('./nh_amplification_cmip6_ensemble_'+str(project_end)+'.png')
fig4.savefig('./tropical_amplification_cmip6_ensemble_'+str(project_end)+'.png')
fig5.savefig('./sh_amplification_cmip6_ensemble_'+str(project_end)+'.png')   
fig6.savefig('amplification_cmip6_ensemble_'+str(project_end)+'.png') 
    
    
    