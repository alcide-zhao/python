# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This script is t preprocess the surface temperature, prcipitation, 
and 850hPa winds into 1971-2000 JJA mean for model evaluation
The precipitation and temperature mean and std are derived from CESM model products
The precipitation mean are derived from the Aphridate observations over only AM region
The temperature and 850hPa winds are derived from NCEP reanalysis
"""

import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as sio
import os 
import site
# from scipy.interpolate import interp2d

lib_path = os.path.join(
	os.path.realpath(
		os.path.dirname(__file__)
	), 
	os.path.pardir, 
	os.path.pardir,
    # 'lib'
)

site.addsitedir(lib_path)

from lib import *
#leapyear judgement
def leapyear(iyear):
	leapyear = False
	if (iyear % 400 == 0):
		leapyear = True
	else:
		if (iyear % 100 == 0):
			leapyear = False
		else:
			if (iyear % 4 == 0):
				leapyear = True
			else:
				leapyear = False

def JJA_mean_CESM_variables(variable):
	file_name =  '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/ensumble_mean_'+variable+'_192002_200512.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	variable_cache = nc_fid.variables[variable][:] 
	variable_MV = nc_fid.variables[variable].missing_value
	if (np.rank(variable_cache) == 4):
		variable =variable_cache[:,23,:,:] # 850hpa winds for U and V 
	else:
		variable = variable_cache
	variable[variable == variable_MV] = np.nan
	# units =  nc_fid.variables[variable].units
	nc_fid.close()
	size = np.shape(variable)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2005]
	###JJA mean
	
	variable_JJAmean = np.empty((86,size[1],size[2])); variable_JJAmean[:] =np.nan
	
	layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
	# print year_series
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = variable[layer_b:layer_e+1,:,:]
		variable_JJAmean[iyear-1920,:,:] = stats.nanmean(cache,axis=0)
	variable_71_00 = variable_JJAmean[66:85,:,:]
	year_series = year_series[66:85]
	return year_series,lon,lat,variable_71_00

def JJA_mean_NCEP_variables(variable):
	file_name =  '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/NCEP_'+variable+'.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	variable_cache = nc_fid.variables[variable][:]  # cache 
	variable_MV = nc_fid.variables[variable].missing_value
	variable_valid_range = nc_fid.variables[variable].valid_range
	if (np.rank(variable_cache) == 4):
		variable =variable_cache[:,2,:,:] # 850hpa winds for U and V 
	else:
		variable = variable_cache
	variable[variable == variable_MV ] = np.nan;
	variable[variable>variable_valid_range[1] ] = np.nan;variable[variable<variable_valid_range[0] ] = np.nan;
	# units =  nc_fid.variables[variable].units
	nc_fid.close()
	size = np.shape(variable)
	NCEP_time = datetime_second_datetime_integer( time*60*60,reference = '1800-01-01 00:00:00')
	year_series= np.unique(map(int,NCEP_time/10000))
	###JJA mean
	variable_JJAmean = np.empty((70,size[1],size[2])); variable_JJAmean[:] =np.nan # 1948-2017
	layer_e = -5  # CESMN MONTHLY DATA BEGINS FROM Jan OF 1948
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = variable[layer_b:layer_e+1,:,:]
		variable_JJAmean[iyear-1948,:,:] = stats.nanmean(cache,axis=0)
	variable_86_05 = variable_JJAmean[38:57,:,:]
	year_series = year_series[38:57]

	return year_series,lon,lat,variable_86_05
	

# CESM precipitation mean and std
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat_CESM = nc_fid.variables['lat'][:]
lon_CESM = nc_fid.variables['lon'][:]
time_CESM = nc_fid.variables['time'][66:85]
prep_mean_CESM = nc_fid.variables['mean_precip'][66:85,:,:]
nc_fid.close()

CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/ensumble_mean_PREPSTD_192001-200512.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
prep_std_CESM = nc_fid.variables['PREPSTD'][66:85,:,:]
nc_fid.close()

# CESM SURFACE TEMPERATURE mean and std	
variable = 'TS'; _,_,_,TS_JJA_mean_CESM=JJA_mean_CESM_variables(variable);
variable = 'TSSTD'; _,_,_,TS_JJA_std_CESM=JJA_mean_CESM_variables(variable);	
# CESM 850hPa winds			
variable = 'U'; _,_,_,U_JJA_mean_CESM=JJA_mean_CESM_variables(variable);	
variable = 'V'; _,_,_,V_JJA_mean_CESM=JJA_mean_CESM_variables(variable);	

# APHRODITE precipitation mean
APHRO_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(APHRO_file,mode='r')
lat_APHRO = nc_fid.variables['lat'][:]
lon_APHRO = nc_fid.variables['lon'][:]
date_APHRO = nc_fid.variables['time'][20:50]
prep_mean_APHRO = nc_fid.variables['mean_precip'][20:50,:,:]


# APHRODITE TEMPERATURE 
APHRO_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/APHRO_T_MA_050deg_V1101_1961_2006.nc'

nc_fid = nc4.Dataset(APHRO_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
temp = nc_fid.variables['temp'][:]
time = nc_fid.variables['time'][:]
missing_value = nc_fid.variables['temp'].missing_value
temp[temp == missing_value] = np.nan
print missing_value
print temp

temp_mean_APHRO	= np.empty((20,140,180));temp_mean_APHRO[:]=np.nan;
for iyear in range(1986,2005):
	if (leapyear(iyear)):
		layer_b = [l for l in range(len(time)) if time[l] == iyear*1000+152 ][0]
		layer_e = layer_b + 123
	else:
		layer_b = [l for l in range(len(time)) if time[l] == iyear*1000+151 ][0]
		layer_e = layer_b + 123
	temp_cache = temp[layer_b:layer_e,:,:]
	temp_mean_APHRO[iyear-1986,:,:] = stats.nanmean(temp_cache,axis=0)


# NCEP surface tamperature
variable = 'TS'; year_series,lon_NCEP,lat_NCEP,TS_JJA_mean_NCEP=JJA_mean_NCEP_variables(variable);
# NCEP 850hPa Winds		
variable = 'U'; _,_,_,U_JJA_mean_NCEP=JJA_mean_NCEP_variables(variable);	
variable = 'V'; _,_,_,V_JJA_mean_NCEP=JJA_mean_NCEP_variables(variable);	

# output the data intp .mat file
output_PATH ='/exports/csce/datastore/geos/users/s1667168/PP/'
sio.savemat(output_PATH+'TS_P_850UV_model_eval.mat', {'time':year_series,\
'lon_CESM':lon_CESM,'lat_CESM':lat_CESM,'prep_mean_CESM':prep_mean_CESM,'prep_std_CESM':prep_std_CESM,\
'TS_JJA_mean_CESM':TS_JJA_mean_CESM,'TS_JJA_std_CESM':TS_JJA_std_CESM,'U_JJA_mean_CESM':U_JJA_mean_CESM,'V_JJA_mean_CESM':V_JJA_mean_CESM,\
'lon_APHRO':lon_APHRO,'lat_APHRO':lat_APHRO,'prep_mean_APHRO':prep_mean_APHRO,'temp_mean_APHRO':temp_mean_APHRO,\
'lon_NCEP':lon_NCEP,'lat_NCEP':lat_NCEP,'TS_JJA_mean_NCEP':TS_JJA_mean_NCEP,'U_JJA_mean_NCEP':U_JJA_mean_NCEP,'V_JJA_mean_NCEP':V_JJA_mean_NCEP})

# data = sio.loadmat(output_PATH+'Aerosol_Burden_AOD_JJA.mat')
# print sio.whosmat(output_PATH+"Aerosol_Burden_AOD_JJA.mat")