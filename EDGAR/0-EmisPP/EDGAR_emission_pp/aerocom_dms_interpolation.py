"""
Interpolate the aerocom DMS data for the year needed 
as the DMS data that the model uses has only data foe the first and last year (horrible)
"""

import netCDF4 as nc4
import numpy as np
import math
import sys
from mpl_toolkits.basemap import interp

file = '/exports/csce/datastore/geos/users/s1667168/EDGAR/aerocom_mam3_dms_surf_1849-2100_c111017.nc'

year_b=1944;year_e=2101;year=2010
slope = (year-year_b)/(year_e-year_b)

date_n = np.zeros((36));DMS_n = np.zeros((36,96,144))

nc_fid = nc4.Dataset(file,mode='r+')
date_o=nc_fid.variables['date'][:]
DMS_o=nc_fid.variables['DMS'][:]

date_n[0:12]=date_o[0:12];date_n[24:36]=date_o[12:24];
DMS_n[0:12,:,:]=DMS_o[0:12,:,:];DMS_n[24:36,:,:]=DMS_o[12:24,:,:];

for layer in range(12,24):
	date_n[layer]=date_o[layer]%10000+10000*2010
	DMS_n[layer,:,:] = DMS_o[layer-12,:,:]+slope*(DMS_o[layer,:,:]-DMS_o[layer-12,:,:])

nc_fid.variables['date'][:]=date_n
nc_fid.variables['DMS'][:]=DMS_n

nc_fid.close()