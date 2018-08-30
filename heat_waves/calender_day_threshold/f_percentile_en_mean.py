'''
The calenday day 95/90 percentiles of TX and TN were derived from all the 
'''

import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import os
import glob



#############################################
# 2 ensumble mean 4D data
#############################################

def get_en_mean_var(variable):
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/threshold/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + 'file_list.txt', "r")
	text_content = text_file.readlines()
	en_store = np.empty((10,365,192,288))
	for en_no in range(0,10):
		nc_f = text_content[en_no][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		en_store[en_no,:,:,:] = nc_fid.variables[variable][:]
	var_mean = stats.nanmean(en_store,axis=0)
	del en_store
	return lat,lon,var_mean
	

lat,lon,TX95P = get_en_mean_var('TX95P')
lat,lon,TN95P = get_en_mean_var('TN95P')
lat,lon,TX5P = get_en_mean_var('TX5P')
lat,lon,TN5P = get_en_mean_var('TN5P')

output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/threshold/'
file_name_out = output_path+'TmTn_percentile_calender_enMean.nc'
f = nc4.Dataset(file_name_out,'w', format='NETCDF4')

f.createDimension('day', 365)
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))


days = f.createVariable('year',np.float32, ('day'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))

days[:] = range(1,366)
latitudes[:] = lat
longitudes[:] = lon

TX95Ps= f.createVariable('TX95P',np.float32,('day','lat','lon'))
TN95Ps= f.createVariable('TN95P',np.float32,('day','lat','lon'))
TX5Ps= f.createVariable('TX5P',np.float32,('day','lat','lon'))
TN5Ps= f.createVariable('TN5P',np.float32,('day','lat','lon'))

TX95Ps[:] = TX95P
TN95Ps[:] = TN95P
TX5Ps[:] =TX5P
TN5Ps[:]= TN5P


f.description = 'degree kelvin'
f.close()


