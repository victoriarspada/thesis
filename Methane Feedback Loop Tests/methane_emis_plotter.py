"""
Created on Fri Apr 14 22:21:17 2023
@author: Victoria Spada

Plot the natural emissions resulting from AMAP emulator runs with different 
alpha values.
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator


fig0, ax0 = plt.subplots(2, 1, figsize=(9,7))
fig0.subplots_adjust(wspace=0.4,
                    hspace=0.5)
ax0[0].set_title('Natural Emissions of Methane')
ax0[1].set_title('Arctic Temperature Perturbations (°C)')

ax0[0].set_xlabel('Years')
ax0[1].set_xlabel('Years')
ax0[0].set_ylabel('Annual Pulse of Natural Methane Emissions (Tg)', fontsize='small')
ax0[1].set_ylabel('Arctic Temperature Perturbation (°C)', fontsize='small')

ax0[0].set_xlim([1750,2050])
ax0[1].set_xlim([1750,2050])

files = ['0_results.txt', '0.1_results.txt', '0.25_results.txt', '1_results.txt']
colors = ['firebrick', 'royalblue', 'darkgoldenrod', 'rebeccapurple']
labels = ['\u03B1=0.00', '\u03B1=0.10', '\u03B1=0.25', '\u03B1=1.00']
years = np.arange(1750,2050.1,1.0)
for f in range(0, len(files), 1):
    data = np.loadtxt(files[f])
    natemis = data[0,:][0:301] 
    dtemps = data[1,:] 
    ax0[0].plot(years, natemis, c=colors[f], label=labels[f])
    ax0[1].plot(years, dtemps, c=colors[f], label=labels[f])
    
ax0[0].legend()