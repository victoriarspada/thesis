# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 17:28:39 2023

@author: Victoria Spada
Print out the names of the model names in the CMIP6 Data set

Move the files representing ESMs that aren't included in all ESMs to an 
Extras directory
"""

import shutil
import os
# Plot the time series for each SSP

data_dir = 'cmip6_data/'
ssp_dirs = ['DATA_ssp126/', 'DATA_ssp245/', 'DATA_ssp370/', 'DATA_ssp585/']
n_ssps = len(ssp_dirs)

esms = []
esm_occurences = []
for s in range(0, n_ssps, 1): # For each of the SSPs
    curr_data_dir = data_dir+ssp_dirs[s]
    curr_files = os.listdir(curr_data_dir)
    curr_files.remove('CanESM5') # folder with CanESM files copied
    curr_files.remove('Extras') # folders with duplicates
    
    # print(curr_files, len(curr_files))
    curr_esms = []
    for i in curr_files: # gather all ESM names
        if i.split('_')[3] not in curr_esms:
            curr_esms += [i.split('_')[3]]
       
    for j in curr_esms:
        if j in esms: 
            esm_occurences[esms.index(j)] += 1
        else:
            esms += [j]
            esm_occurences += [1]
            
# print('all esms', len(esms), esms)

subset = []
for e in range(0, len(esms), 1):
    if esm_occurences[e] == 4:
        subset += [esms[e]]
        
# print('subset',len(subset), subset)


# Now print out a list of all the files that correspond to an ESM that we 
# don't have a realization of each SSP for
for s in range(0, n_ssps, 1): # For each of the SSPs
    curr_data_dir = data_dir+ssp_dirs[s]
    curr_files = os.listdir(curr_data_dir)
    curr_files.remove('CanESM5') # folder with CanESM files copied
    curr_files.remove('Extras') # folders with duplicates
    
    curr_esms = []
    for i in curr_files: # gather all ESM names
        if i.split('_')[3] not in subset: 
            source = './'+data_dir+ssp_dirs[s]+i
            destination = './'+data_dir+ssp_dirs[s]+'Extras/'+i
            print(source, destination)
            shutil.move(source, destination)
            
            
            
            