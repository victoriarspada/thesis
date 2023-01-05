Contact:
--------

Knut von Salzen, ECCC, knut.vonsalzen@ec.gc.ca

Software licence:
-----------------

https://open.canada.ca/en/open-government-licence-canada

Instructions:
-------------

1) The emulator configuration is provided via the file emulator.csv
2) Various python scripts are available (use python3 to run). Separate emulator components for the impacts of different GHGs and SLCF on radiative forcings and temperature are provided (*_CO2.py, *_CH4.py, *_aerosol.py, *_O3.py). In addition, emulator_pm25.py provides PM2.5 projections.
3) Postprocessing and data analysis tools are provided in the data directory "postproc".
4) For verification of simulation results, compare with data in the data directory "verification".

Emissions:
----------

AMAP scenarios: ECLIPSE V6B (simulated by IIASA GAINS), combined with CMIP6 biomass burning emissions.

For ECLIPSE emissions, a distinction is made between emissions from upstream oil and gas production sources (“oil & gas flaring”) and combined fossil and bio fuel sources related to energy consumption (“fossil & bio fuel”). The former includes storage and distribution, including the intended venting and unintended leakage during extraction and transportation of oil and gas, the release of ventilation air CH4 during coal mining, and the flaring of excess gases and liquids. The latter includes land-based emissions from residential and commercial sources; agriculture and waste burning on fields; power plants, energy conversion and extraction; industrial combustion and processing; surface transportation; and waste processing. In addition to these upstream and downstream sources, emissions from international shipping sources (“shipping”) are considered, too.

Supporting files ("netcdf" folder):
-----------------------------------

grid2gains_amap_BC_2019.nc:

Grid2gains-region allocation mask file. Based on GAINS country/region allocation.

Regions:
AC members:
-       NORD, Nordic countries
-       CANA, Canada
-       USAM, USA
-       RUSS, Russia
Observers:
-       EURO, European observers
-       CHIN, China
-       INDI, India
-       ASIO, Other Asian observers
Other:
-       OEUR, Other European countries (non-member & non-observer)
-       ALL, containing all areas belonging to any of the above
-       ROW, rest of the world, all areas not included in the ALL
-       ROW_land, rest of the world, land areas only (included mainly due to technical reasons). 

grid2gains_amap_BC_2019_minimal.nc:

Minimal set of regions. Regions provided in grid2gains_amap_BC_2019.nc were grouped into larger geographic regions in order to come up with a smaller set of regions which could be used for modelling.

Regions:
Western AC members (AWEST):
-       CANA, Canada
-       USAM, USA
Eastern AC members (AEAST):
-       NORD, Nordic countries
-       RUSS, Russia
Rest of Europe (ROEUR):
-       EURO, European observers
-       OEUR, Other European countries (non-member & non-observer)
Asia (ASIA):
-       CHIN, China
-       INDI, India
-       ASIO, Other Asian observers
Rest of the world (ROW):
-       ROW, all remaining regions (identical to ROW in grid2gains_amap_BC_2019.nc)

arctic.nc:

Arctic region (ARCTIC, all latitudes >= 60 degree N). This mask can be used to account for Arctic vs. non-Arctic shipping. These two regions should be used in addition to the GAINS regions for simulations of shipping impacts.
