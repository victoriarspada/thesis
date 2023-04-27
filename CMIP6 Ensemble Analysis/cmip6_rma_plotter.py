# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

Plotting RMA fits of regional amplification for CMIP6 ensemble, for 4 SSPs.
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
ssp_labels = ['SSP1-2.6', 'SSP2-4.5', 'SSP3-7.0', 'SSP5-8.5']
colors0 = ['rebeccapurple', 'seagreen', 'darkgoldenrod', 'crimson', 'seagreen']

n_ssps = len(ssp_dirs)

n_years = 451
project_start, project_end = 1995, 2050
ref_start, ref_end = 1995, 2014

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
    
region_names = ['Arctic', 'NH Midlatitudes', 'Tropical', 'SH Extratropical']
suptitle2 = ' SAT Anomalies vs Global Anomaly of CMIP6 Ensemble ('+str(project_start)+' - '+str(project_end)+')'
suptitle3 = 'CMIP6 Emulator Regional SAT Amplification ('+str(project_end)+')'


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
    
    # Take the MEANS, each row is a region
    tas_time_series_mean = np.zeros((n_regions, n_years)) # Define mean temp. matrix
    tas_mean_anomaly = np.zeros((n_regions, n_years)) # Define mean temp. anomaly matrix
    sigma = np.zeros((n_regions, n_years))
    max_tas_anomaly = np.zeros((n_regions, n_years))
    min_tas_anomaly = np.zeros((n_regions, n_years))
    
    # Finished cycling through SSP files.
    # Now create the ensemble for each region.
    
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
            
            m,b,sigma_m,sigma_b = reduced_major_axis(
                        abs(gm_anomalies[s, project_start_i:project_end_i]),
                        abs(arc_anomalies[s, project_start_i:project_end_i]))
            corr, _ = pearsonr(gm_anomalies[s, project_start_i:project_end_i],
                                arc_anomalies[s, project_start_i:project_end_i])
            # print(s,'Arctic Pearsons correlation: %.2f' % corr)
            
            # ax1[r-1].text(s='m =  %.2f' % lsr[0] + ' ± %.2f' % lsr[2],
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,     
                          x=0.12, y=3.6-s*0.29,
                          #x=0.06, y=2.1-s*0.14,
                          color=colors0[s],
                          fontsize=8)
            slopes[r-1,s] = lsr[0]
            slopes_sigma[r-1,s] = lsr[2]
        if r==2:
            ax1[r-1].plot(gm_anomalies[s, project_start_i:project_end_i],
                        nh_anomalies[s, project_start_i:project_end_i],
                        color=colors0[s],
                        linewidth=1)
            lsr = least_squares_fit(gm_anomalies[s, project_start_i:project_end_i],
                        nh_anomalies[s, project_start_i:project_end_i])
            m,b,sigma_m,sigma_b = reduced_major_axis(
                        abs(gm_anomalies[s, project_start_i:project_end_i]),
                        abs(nh_anomalies[s, project_start_i:project_end_i]))
            
            corr, _ = pearsonr(gm_anomalies[s, project_start_i:project_end_i],
                                nh_anomalies[s, project_start_i:project_end_i])
            # print(s,'NH Pearsons correlation: %.2f' % corr)
            
            # ax1[r-1].text(s='m =  %.2f' % lsr[0] + ' ± %.2f' % lsr[2],
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,
                          x=0.12, y=2.05-s*0.17,
                          #x=0.06, y=1.15-s*0.08,
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
            m,b,sigma_m,sigma_b = reduced_major_axis(
                        abs(gm_anomalies[s, project_start_i:project_end_i]),
                        abs(trop_anomalies[s, project_start_i:project_end_i]))
            
            corr, _ = pearsonr(gm_anomalies[s, project_start_i:project_end_i],
                                trop_anomalies[s, project_start_i:project_end_i])
            # print(s,'Trop Pearsons correlation: %.2f' % corr)
            
            # ax1[r-1].text(s='m =  %.2f' % lsr[0] + ' ± %.2f' % lsr[2],
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,
                          x=0.12, y=1.46-s*0.13,
                          # x=0.06, y=0.74-s*0.05,
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
            m,b,sigma_m,sigma_b = reduced_major_axis(
                        abs(gm_anomalies[s, project_start_i:project_end_i]),
                        abs(sh_anomalies[s, project_start_i:project_end_i]))
            corr, _ = pearsonr(sh_anomalies[s, project_start_i:project_end_i],
                                nh_anomalies[s, project_start_i:project_end_i])
            # print(s,'SH Pearsons correlation: %.2f' % corr)
            
            # ax1[r-1].text(s='m =  %.2f' % lsr[0] + ' ± %.2f' % lsr[2],
            ax1[r-1].text(s='m =  %.2f' % m + ' ± %.2f' % sigma_m,
                          x=0.12, y=1.26-s*0.102,
                          # x=0.06, y=0.7-s*0.05,
                          color=colors0[s],
                          fontsize=8)
            slopes[r-1,s] = lsr[0]
            slopes_sigma[r-1,s] = lsr[2]
        if r>0:
            print(' %.2f' % lsr[0] + ' ± %.2f' % lsr[2])

    # Add SSP text labels
    if project_end==2050:
        ax1[1].text(x=1.8, y=2.1-s*0.22, s=ssp_labels[s], color=colors0[s], 
                    fontsize='medium') #, fontweight='bold')
    else:
        ax1[3].text(x=1.05, y=0.7-s*0.1, s=ssp_labels[s], color=colors0[s], fontsize='medium', fontweight='bold')

