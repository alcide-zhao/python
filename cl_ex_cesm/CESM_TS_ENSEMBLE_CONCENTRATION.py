# -*- coding: utf-8 -*-
"""
This script is to preprocess the CESM ensembles into JJA mean and std
The variables which are of interests here arre the surfacr temperature, the JJA total precipitation and the 
derived extreme precipitation indixes : RX5DAY, R99P, CWD, and SDII
the year period are 2006-2100
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


# def JJA_mean_CESM_variables(,variable,file_name):
	# nc_fid = nc4.Dataset(file_name,mode='r')
	# time = nc_fid.variables['time'][:]
	# lat = nc_fid.variables['lat'][:]
	# lon = nc_fid.variables['lon'][:]
	# rcp85_cache = nc_fid.variables['rcp85'][:] 
	# fixa_cache = nc_fid.variables['rcp85_fixA'][:] 
	# variable_MV = nc_fid.variables['rcp85'].missing_value
	# if (np.rank(rcp85_cache) == 4):
		# rcp85 =rcp85_cache[:,23,:,:] # 850hpa winds for U and V 
		# fixa =fixa_cache[:,23,:,:] # 850hpa winds for U and V 
	# else:
		# rcp85 = rcp85_cache;
		# fixa = fixa_cache
	# rcp85[rcp85 == variable_MV] = np.nan
	# fixa[fixa == variable_MV] = np.nan
	# # units =  nc_fid.variables[variable].units
	# nc_fid.close()
	# size = np.shape(rcp85)
	# # year= [ iyear/10000 for iyear in map(int,time)]
	# # year_series = [value for value in np.unique(year) if value <=2100]
	# # print year_series
	# year_series = range(2006,2100)
	# ###JJA mean
	# rcp85_JJAmean = np.empty((95,size[1],size[2])); rcp85_JJAmean[:] =np.nan
	# fixa_JJAmean = np.empty((95,size[1],size[2])); fixa_JJAmean[:] =np.nan
	# layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
	# # print year_series
	# for iyear in year_series:
		# layer_b = layer_e + 10
		# layer_e = layer_b + 2
		# cache = rcp85[layer_b:layer_e+1,:,:]
		# rcp85_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
		# cache = fixa[layer_b:layer_e+1,:,:]
		# fixa_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)		
	# return year_series,lon,lat,rcp85_JJAmean, fixa_JJAmean

	
def JJA_mean_CESM_variables(file_name,variable):
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	variable_cache = nc_fid.variables[variable][:] 
	if (np.rank(variable_cache) == 4):
		variable =variable_cache[:,23,:,:] # 850hpa winds for U and V 
	else:
		variable = variable_cache
	# units =  nc_fid.variables[variable].units
	nc_fid.close()
	size = np.shape(variable)
	year_series = range(2006,2101)
	###JJA mean
	
	variable_JJAmean = np.empty((95,size[1],size[2])); variable_JJAmean[:] =np.nan
	
	layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
	# print year_series
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = variable[layer_b:layer_e+1,:,:]
		variable_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	return year_series,lon,lat,variable_JJAmean


	
# CESM SURFACE TEMPERATURE mean and std	
import glob
variable= 'TS'
data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'+variable+'/rcp85/*.nc' 
files=sorted(glob.glob(data_path))
TS_rcp85 = np.empty((30,95,192,288))
ensemble = 0
for file in files:
	time,lon,lat,TS_rcp85[ensemble,:,:,:] = JJA_mean_CESM_variables(file,variable)
	ensemble  = ensemble+1
	
data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'+variable+'/fixa/*.nc' 
files=sorted(glob.glob(data_path))
TS_fixa = np.empty((15,95,192,288))
ensemble = 0
for file in files:
	time,lon,lat,TS_fixa[ensemble,:,:,:] = JJA_mean_CESM_variables(file,variable)
	ensemble  = ensemble+1

file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/ensenmble_TS_200602_210101.nc'	
f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write

f.createDimension('time', 95)
f.createDimension('ensumble_r', 30)
f.createDimension('ensumble_f', 15)
f.createDimension('lon', len(lon))
f.createDimension('lat', len(lat))

ensumble_rs = f.createVariable('ensumble_r',np.float64, ('ensumble_r'))
ensumble_fs = f.createVariable('ensumble_f',np.float64, ('ensumble_f'))
times = f.createVariable('time',np.float64, ('time'))
lons = f.createVariable('lon',np.float64, ('lon'))
lats = f.createVariable('lat',np.float64, ('lat'))
TS_rcp85s = f.createVariable('TS_rcp85',np.float32,('ensumble_r','time','lat','lon'))
TS_fixas = f.createVariable('TS_fixa',np.float32,('ensumble_f','time','lat','lon'))

ensumble_rs[:] = range(1,31)
ensumble_fs[:] = range(1,16)
times[:] = time
lons[:] =lon;
lats[:] =lat;
TS_rcp85s[:] = TS_rcp85
TS_fixas[:] = TS_fixa

f.description = 'JJA mean surface temperature'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()