# -*- coding: utf-8 -*-
"""
Plot the annual total of NO,CO and NH3 moclecules 
"""

import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import math
import site
import scipy.io as sio

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

#######################################################
# 0.2 functions and subroutines                       #
#######################################################
def movingaverage (values, window=3):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result
	
def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth area in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

def MonthWeight(annual_data,area):
	"""
	weight the days to each month and sum the mass up throughout the year
	"""
	calendar_day = {'1':31,'2':28,'3':31,'4':30,'5':31,'6':30,'7':31,'8':31,'9':30,'10':31,'11':30,'12':31}
	month_sum= np.empty((np.shape(annual_data)))
	for month in range(0,12):
		day_no = calendar_day[str(month+1)]
		month_sum[month,:,:] =annual_data[month,:,:]*day_no*24*60*60/(10**(9))
	month_total = np.nansum(np.nansum(np.multiply(month_sum,area),axis=2),axis=1)
	return month_total

def get_RCP85_emissions(filename):
	"""
	get the sum of emissions over 2010 from all the three sectors
	"""
	rcp85_path = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
	rcp85 = rcp85_path+filename
	print rcp85
	nc_fid = nc4.Dataset(rcp85,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	emission = nc_fid.variables['anthro'][:]  	
	lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
	lons,lats = np.meshgrid (lon,lat)  #meshgrid
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	AreaWeighMass =np.nansum(np.nansum(np.multiply(emission,area),axis=2),axis=1)[:]
	print np.shape(AreaWeighMass)
	return AreaWeighMass
	
CO = get_RCP85_emissions('rcp85_finn_2002-2015_CO_woBiog_1.9x2.5_mol_c20160128.nc');
NH3 = get_RCP85_emissions('rcp85_finn_2002-2015_NH3_1.9x2.5_mol_c20160128.nc');
NO = get_RCP85_emissions('rcp85_finn_2002-2015_NO_1.9x2.5_mol_c20160128.nc');
print np.nanmean(CO[97:109]),np.nanmean(CO[107:119])
print np.nanmean(NH3[97:109]),np.nanmean(NH3[107:119])
print np.nanmean(NO[97:109]),np.nanmean(NO[107:119])	
time_series=range(0,170)
fig1 = plt.figure(facecolor='White',figsize=[5,3]);plot_setup();
ax1 = plt.subplot(1,1,1);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
# lins1 = ax1.plot(time_series,CO,'-',color="k",label='CO',linewidth=1.5)
lins2 = ax1.plot(time_series,NH3,'-',color="b",label='NH3',linewidth=1.5)
# lins3 = ax1.plot(time_series,NO,'-',color="g",label='NO',linewidth=1.5)

plt.show()
