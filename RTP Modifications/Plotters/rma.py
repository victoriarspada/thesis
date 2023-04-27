# -*- coding: utf-8 -*-
"""
@author: Victoria Spada

Function for reduced major axis fitting of a 2-D uncertain dataset.
"""
import numpy as np
from scipy.stats import pearsonr

def reduced_major_axis(x, y):
    print(x.shape, y.shape)
    corr, _ = pearsonr(x,y)
    m = np.std(y)/np.std(x)
    b = (np.mean(y)-m*np.mean(x))/corr
    n = len(x)
    sigma_m = abs(m)*np.sqrt((1-corr**2)/n)
    sigma_b = np.std(y)*np.sqrt(((1-corr**2)/n)*(1+(np.mean(x)/np.std(x))**2))
    return m, b, sigma_m, sigma_b