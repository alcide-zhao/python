# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to convert the world ocean_land map mask from 1440*720 into needed resolution
@author: Alcide.Zhao
"""

import numpy as np
import scipy.io as spio
# from scipy.interpolate import interp2d
from scipy.interpolate import RectBivariateSpline as RectBivariateSpline

import netCDF4 as nc4
import matplotlib.pyplot as plt


lat_res = 0.5;lon_res = 0.5;
lat_grid= int(180/lat_res);lon_grid= int(360/lon_res);
oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720')['landocean']
y = np.arange(0,lat_grid,lat_grid/720.0)
x = np.arange(0,lon_grid,lon_grid/1440.0)
# xi, yi = np.meshgrid(x, y)
y_new = range(lat_grid)
x_new = range(lon_grid)
print np.shape(y)
print np.shape(x)
f = RectBivariateSpline(y, x, oceanmask)

ocean_mask = np.empty((lat_grid,lon_grid))
ocean_mask = np.flipud(f(y_new,x_new))
cache = np.empty((lat_grid,lon_grid))

# make the longitudes from 0 to 360
cache[:,0:int(lon_grid/2)] = ocean_mask[:,int(lon_grid/2):lon_grid]
cache[:,int(lon_grid/2):lon_grid] = ocean_mask[:,0:int(lon_grid/2)]
ocean_mask = cache
# make the land to be 1 and ocean to be 0
ocean_mask[ocean_mask<125] =0
ocean_mask[ocean_mask>=125] =1
ocean_mask[0,:]=1


plt.imshow(ocean_mask, origin = 'lower')
plt.show()
file_name='/home/s1667168/coding/python/climate_extremes_cesm/external_data/\
landoceanmask_'+str(lat_grid)+'_'+str(lon_grid)+'.mat'
mdict={'landoceanmask':ocean_mask}
spio.savemat(file_name, mdict)