"""
This is to plot the ozone concentration difference between 2010 and 1970 
the 1970 data is from ozone_1.9x2.5_L26_1850-2005_c090803.nc
the 2010 data is from ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202.nc
"""

import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib.pyplot as plt


def get_ozone():
	input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/ozone/'
	ozone1970 = input_path+'ozone_1.9x2.5_L26_1850-2005_c090803.nc'
	nc_fid = nc4.Dataset(ozone1970,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	lev = nc_fid.variables['lev'][:]
	o1970 = nc_fid.variables['O3'][156:168,:,:,:]
	nc_fid.close()
	ozone2010 = input_path+'ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202.nc'
	nc_fid = nc4.Dataset(ozone2010,mode='r')
	o2010 = nc_fid.variables['O3'][48:60,:,:,:]
	nc_fid.close()
	ozone_change = o2010-o1970

	return lat,lev,o1970,ozone_change
	

lat,lev,o1970,ozone_change = get_ozone();

## Here replace the  Trop ozone in 2010 file with that from 1970
## the way to determine trop ozone is that under 150ppb in 1970
input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/ozone/'
ozone1970 = input_path+'ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202_1970_Trop_ozoe.nc'
nc_fid = nc4.Dataset(ozone1970,mode='r+')
o2010 = nc_fid.variables['O3'][48:60,:,:,:]
o2010[o1970*10**(9)<150] =  o1970[o1970*10**(9)<150]
nc_fid.variables['O3'][48:60,:,:,:] = o2010

nc_fid.close()



