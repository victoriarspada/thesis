"""
Author: Victoria Spada

File to produce modified RTP coefficients with a specified D increase
to the area-weighted sum of the Arctic RTP coefficients.
"""

import numpy as np

# area fractions of latitude bands
area_frc_ar = 0.0669873                           # 60N-90N (Arctic) 
area_frc_ml = 0.1982769                           # 28N-60N (mid latitudes)
area_frc_tr = 0.4694716                           # 28N-28S (tropics)
area_frc_sh = 0.2652642                           # 28S-90S (southern hemisphere)

# regional temperature sensitivity (K (W/m^2)^-1, Shindell et al., 2010)
cs_ar_old = np.array([0.31, 0.17, 0.16, 0.06])    # 60N-90N (Arctic) 
cs_ml_old = np.array([0.06, 0.24, 0.17, 0.07])    # 28N-60N (mid latitudes)
cs_tr_old = np.array([0.02, 0.10, 0.24, 0.09])    # 28N-28S (tropics)
cs_sh_old = np.array([0.  , 0.02, 0.05, 0.19])    # 28S-90S (southern hemisphere)

# weightings for the RTP coefficient rescaling
D = 0.50 # percent increase for the area weighted sum of the Arctic RTP
delta = D*area_frc_ar*cs_ar_old
d1 = delta/(area_frc_ar*cs_ar_old)
d2 = -delta/(3*area_frc_ml*cs_ml_old)
d3 = -delta/(3*area_frc_tr*cs_tr_old)
d4 = -delta/(3*area_frc_sh*cs_sh_old)

cs_ar = list(np.array([0.31, 0.17, 0.16, 0.06])*(1+d1))    # 60N-90N (Arctic) 
cs_ml = list(np.array([0.06, 0.24, 0.17, 0.07])*(1+d2))    # 28N-60N (mid latitudes)
cs_tr = list(np.array([0.02, 0.10, 0.24, 0.09])*(1+d3))    # 28N-28S (tropics)
cs_sh = list(np.array([0.  , 0.02, 0.05, 0.19])*(1+d4))    # 28S-90S (southern hemisphere)

# account for cs_ch4_sh[0] = 0 (get inf.)
d4 = -delta[0]/(2*area_frc_ml*cs_ml_old[0])
d5 = -delta[0]/(2*area_frc_tr*cs_tr_old[0])

cs_ml[0] = cs_ml_old[0]*(1+d4)    # 28N-60N (mid latitudes)
cs_tr[0] = cs_tr_old[0]*(1+d5)
cs_sh[0] = 0

print(cs_ar)
print(cs_ml)
print(cs_tr)
print(cs_sh)

fname = str(D)+'_rtp_coeffs.txt'
np.savetxt(fname, [cs_ar,cs_ml,cs_tr,cs_sh], delimiter=',')

