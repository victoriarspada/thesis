# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 17:10:11 2023

@author: Victoria Spada 
RTP coefficient ratios and Mean Absolute Differences
"""
import numpy as np

pi = np.pi
# Latitude band boundaries dividing the latitude bands
latitudes = np.array([90, 60, 28, -28, -90])
theta = np.array([0, 30, 62, 118, 180])
theta = (pi/180)*theta

# Coefficients for relative area of Earth
A_arctic = (np.cos(theta[0]) - np.cos(theta[1])) 
A_nh = (np.cos(theta[1]) - np.cos(theta[2]))
A_sh = (np.cos(theta[2]) - np.cos(theta[3]))
A_ant = (np.cos(theta[3]) - np.cos(theta[4]))

A = np.array([A_arctic, A_nh, A_sh, A_ant]) # Array of relative areas
A = A/(2*np.pi)

# a = 11.17
# b = 120
# c = 150
# d = 200
# e = ((1/a) + (1/b) + (1/c) + (1/d))**(-1)
# f = 1/e

# scenarios = [SSP119, SSP126, SSP245, SSP370, SSP585, CLE, MFR, CFM, ]
# Load AMAP results from text file (2030)
data2030 = np.loadtxt('2030_results_amap.txt')

arctic2030 = data2030[0,:] #[1.379664654322452, 1.3996047047830464, 1.346670101811802, 1.553327686387983, 1.4568171844096274]
a3ae = data2030[1,:] #[0.21222776362769613,0.21175255240142854, 0.20183175964973366,0.2594324288018739, 0.22824551794848325]

global2030 = data2030[2,:] #[0.7755789205079199, 0.7689151555706294, 0.7353277683938022, 0.7353577501523292, 0.8028828818475411]
g3ae = data2030[3,:] #[0.10758207847201307, 0.10157770168944308, 0.09519806096469025, 0.09247995654807939, 0.11135195344936161]

# CMIP Results 
global2030_cmip = [0.63,0.57,0.63,0.62,0.66]
g3ce_min = [0.32,0.38,0.46,0.45,0.49]
g3ce_max = [0.91,0.78,0.80,0.78,0.89]

arctic2030_cmip = [1.92,1.30,1.53,1.36,1.46]
a3ce_min = [0.54,0.41,0.54,0.55,0.66]
a3ce_max = [2.87,2.08,2.20,2.12,2.40]

# Load AMAP results from text file (2050)
data2050 = np.loadtxt('2050_results_amap.txt')

arctic2050 = data2050[0,:] #[1.4924181556637672, 1.8575219859196235, 1.8762024364653482, 2.298714317872845, 2.507596706019554]
a5ae = data2050[1,:] #[0.1481944295431653,0.18845582285930462,0.17524494495351362,0.2703399914294953,0.30179501550914223]

global2050 = data2050[2,:] #[0.9999826382206163, 1.16866956754836, 1.2150631932216058, 1.2755085397038224, 1.529021431372502]
g5ae = data2050[3,:] #[0.11420359716195155,0.1218434724127263, 0.12409110347143448, 0.12258776325255605, 0.1792830812670383]


global2050_cmip = [0.67,0.80,1.14,1.22,1.40]
g5ce_min = [0.31,0.49,0.85,0.96,1.1]
g5ce_max = [1.03,1.20,1.42,1.56,1.79]

arctic2050_cmip = [2.04, 1.89, 2.55, 2.94, 3.20 ]
a5ce_min = [0.68,0.76,1.59,1.67,2.11]
a5ce_max = [2.73,3.22,3.64,3.95,4.41]

g3ce, g5ce, a3ce, a5ce = [], [], [], []
for i in range(0, len(global2030), 1):
    g3ce += [max( global2030_cmip[i] - g3ce_min[i], g3ce_max[i] - global2030_cmip[i])]
    g5ce += [max( global2050_cmip[i] - g5ce_min[i], g5ce_max[i] - global2050_cmip[i])]
    a3ce += [max( arctic2030_cmip[i] - a3ce_min[i], a3ce_max[i] - arctic2030_cmip[i])]
    a5ce += [max( arctic2050_cmip[i] - a5ce_min[i], a5ce_max[i] - arctic2050_cmip[i])]

def cv(a, b):
    array1, array2 = np.array(a), np.array(b)
    diff = array1 - array2 
    mean = (array1+array2)/2
    reldiff = diff/mean #a fraction ie, 10%
    print('reldiff',reldiff,'\n')
    mu = np.mean(reldiff) # mean relative difference
    sigma = np.std(reldiff) # std of difference
    return sigma/mean

def rel_diff(a,b):
    array1, array2 = np.array(a), np.array(b)
    diff = array1 - array2 
    mean = (array1+array2)/2
    N = len(a)
    c = np.divide(diff,mean)
    
    return (1/N)*np.sum(c)*100

def rel_diff_e(a,b,c,d):
    array1, array2 = np.array(a), np.array(b)
    diff = array1 - array2 
    mean = (array1+array2)/2
    N = len(a)
    j = np.divide(diff,mean)

    earray1, earray2 = np.array(c), np.array(d)
    k = np.sqrt(np.square(earray1) + np.square(earray2))
    l = 0.5*k
    m = np.multiply(j, np.sqrt(np.square( np.divide(l, mean)) + np.square(np.divide(k, diff))))
    n = np.square(m/N)
    o = np.sqrt( np.sum(n) )*100
    return o
    
print(rel_diff_e(arctic2030, arctic2030_cmip, a3ae, a3ce ))
print(rel_diff_e(arctic2050, arctic2050_cmip, a5ae, a5ce ))
print(rel_diff_e(global2030, global2030_cmip, g3ae, g3ce ))
print(rel_diff_e(global2050, global2050_cmip, g5ae, g5ce ))

print(rel_diff(arctic2030, arctic2030_cmip ))
print(rel_diff(arctic2050, arctic2050_cmip ))
print(rel_diff(global2030, global2030_cmip ))
print(rel_diff(global2050, global2050_cmip ))

# print(cv(arctic2030, arctic2030_cmip))
# print(cv(arctic2050, arctic2050_cmip))
# print(cv(global2030, global2030_cmip))
# print(cv(global2050, global2050_cmip))
