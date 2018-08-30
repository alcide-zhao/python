# -*- coding: utf-8 -*-
'''
This is to plot the gradients of HW and CS magnitude against one degree of warming from 
	hiS (1920-2005) GHG(RCp8.5_FixA) GHG+AEROSOL (RCP8.5) and AEROSOL (RCP8.5-FixA)
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

def movingaverage (values, window=10):
	boudary = int(math.floor(window/2))
	result = np.empty((np.shape(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index,:,:] = np.nanmean(values[index:index+window,:,:],axis=0)
	for index in range(boudary,np.shape(values)[0]):
		result[index,:,:] = np.nanmean(values[index-window:index,:,:],axis=0)
	# for index in range(-1*boudary-1,0):
		# result[index,:,:] = np.nanmean(values[index:-1,:,:],axis=0)
	# for index in range(boudary,len(values)-boudary):
		# result[index,:,:] = np.nanmean(values[index-boudary:index+boudary+1,:,:],axis=0)
	return result

"""
def get_annual_monthly_oscilation():
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_TS_192002_200512.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	lon=nc_fid.variables['lon'];
	lat=nc_fid.variables['lat'];
	his=nc_fid.variables['TS'][912:1032,:,:];
	nc_fid.close()
	
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_TS_200602_210101.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	lon=nc_fid.variables['lon'];
	lat=nc_fid.variables['lat'];
	R85=np.concatenate([his,nc_fid.variables['rcp85']],axis=0);
	FiX=np.concatenate([his,nc_fid.variables['rcp85_fixA']],axis=0);
	rcp85_T=np.empty((105,192,288)); rcp85_T[:]=np.nan
	FixA_T=np.empty((105,192,288)); FixA_T[:]=np.nan
	rcp85_O=np.empty((13,95,192,288)); rcp85_T[:]=np.nan
	FixA_O=np.empty((13,95,192,288)); FixA_T[:]=np.nan	
	for iyear in range(0,105):
		rcp85_T[iyear,:,:] = stats.nanmean(R85[iyear*12:(iyear+1)*12,:,:],axis=0)
		FixA_T[iyear,:,:] = stats.nanmean(FiX[iyear*12:(iyear+1)*12,:,:],axis=0)
	rcp85_O[0,:,:,:] = movingaverage(rcp85_T)[10:105,:,:];
	FixA_O[0,:,:,:] = movingaverage(FixA_T)[10:105,:,:];

	for imonth in range(1,13):
		rcp85_T[:] = R85[range(imonth-1,1260,12),:,:]
		FixA_T[:] = FiX[range(imonth-1,1260,12),:,:]
		rcp85_O[imonth,:,:,:] = movingaverage(rcp85_T)[10:105,:,:];
		FixA_O[imonth,:,:,:] = movingaverage(FixA_T)[10:105,:,:];
	return lon,lat,rcp85_O,FixA_O
	
lon,lat,rcp85_o,fixa_o = get_annual_monthly_oscilation()

file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/TS_oscilation_rcp85_fixa_2006_2100.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4') 
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('year', 95)
f.createDimension('scale', 13)

lats = f.createVariable('lat',np.float32, ('lat'))
lons = f.createVariable('lon',np.float32, ('lon'))
years = f.createVariable('year',np.float32, ('year'))
scales = f.createVariable('scale',np.float32, ('scale'))

rcp85_os = f.createVariable('rcp85',np.float32,('scale','year','lat','lon'))
fixa_os = f.createVariable('fixa',np.float32,('scale','year','lat','lon'))
his_6190_means = f.createVariable('his_6190_mean',np.float32,('lat','lon'))

lats[:] = lat
lons[:] = lon
years[:] = range(2006,2101)
scales[:] = range(0,13)
rcp85_os[:] = rcp85_o
fixa_os[:] = fixa_o

scales.long_name='0-annual, 1-12 for each month'
f.description = 'TS oscilaitons in the future under RCP8.5 and RCP8.5_FIXA, decadal smoothed (5 years either side)'
f.close()

	
"""

def get_annual_monthly_oscilation():
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_TS_192002_200512.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	lon=nc_fid.variables['lon'];
	lat=nc_fid.variables['lat'];
	TS=nc_fid.variables['TS'][372:1032,:,:];
	print np.shape(TS)
	his_T=np.empty((55,192,288)); his_T[:]=np.nan
	his_O=np.empty((13,45,192,288)); his_O[:]=np.nan
	for iyear in range(0,55):
		his_T[iyear,:,:] = stats.nanmean(TS[iyear*12:(iyear+1)*12,:,:],axis=0)
	his_O[0,:,:,:] = movingaverage(his_T)[10:55,:,:];
	TS_6190_A =  stats.nanmean(his_T[10:40,:,:],axis=0)
	
	for imonth in range(1,13):
		his_T[:] = TS[range(imonth-1,660,12),:,:]
		his_O[imonth,:,:,:] = movingaverage(his_T)[10:55,:,:];
	
	return lon,lat,his_O,TS_6190_A
	nc_fid.close();
	
lon,lat,his_o,TS_6190_A = get_annual_monthly_oscilation()

file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/TS_oscilation_rhis_1960_2005.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4') 
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('year', 45)
f.createDimension('scale', 13)

lats = f.createVariable('lat',np.float32, ('lat'))
lons = f.createVariable('lon',np.float32, ('lon'))
years = f.createVariable('year',np.float32, ('year'))
scales = f.createVariable('scale',np.float32, ('scale'))
his_os = f.createVariable('his',np.float32,('scale','year','lat','lon'))
TS_6190_As = f.createVariable('TS_6190_A',np.float32,('lat','lon'))

lats[:] = lat
lons[:] = lon
years[:] = range(1961,2006)
scales[:] = range(0,13)
his_os[:] = his_o
TS_6190_As[:] = TS_6190_A
scales.long_name='0-annual, 1-12 for each month'
TS_6190_As.long_name="Inter-annual mean (1961-1990) of the annually averaged TS "
f.description = 'TS oscilaitons 1961-2005, decadal smoothed (5 years either side)'
f.close()

	
	
	
	
	
	

