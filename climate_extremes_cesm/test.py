import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import stats
import math
# from scipy.interpolate import interp2d
from lib import *

oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

variable_name ='AODVIS'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/cloud_physics/'
file_name = 'ensumble_mean_'+variable_name+'_200602_210101.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
time = nc_fid.variables['time'][:]
lons = nc_fid.variables['lon'][:]
lats = nc_fid.variables['lat'][:]

#RCP8.5

att_rcp85 = nc_fid.variables['rcp85'][:]
att_rcp85 = np.multiply(att_rcp85,oceanmask)
att_rcp85 = stats.nanmean(att_rcp85,axis = 0)
att_rcp85 = np.ma.masked_where(np.isnan(att_rcp85),att_rcp85)
#FixA
rcp85_fixA = nc_fid.variables['rcp85_fixA'][:]
rcp85_fixA = np.multiply(rcp85_fixA,oceanmask)
lons,lats,rcp85_fixA = range_clip(40,65,18,30,lons,lats,rcp85_fixA)
rcp85_fixA = stats.nanmean(rcp85_fixA[4:7,:,:],axis = 0)
rcp85_fixA = np.ma.masked_where(np.isnan(rcp85_fixA),rcp85_fixA)
print stats.nanmean(rcp85_fixA)
# for value in rcp85_fixA:
	# print value
print np.shape(rcp85_fixA)
cmap = cm.jet
cmap.set_bad('white',1.)

plt.imshow(att_rcp85,origin='lower',cmap=cmap)
plt.show()

plt.imshow(rcp85_fixA,origin='lower',cmap=cmap)
plt.show()