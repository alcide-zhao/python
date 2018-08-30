
# -*- coding: utf-8 -*-
"""
This is to caculate the 90 and 95 threshold of the global temp and precipi
for computing the extremes. the baseline fot the threshold is 1961-1950
The results are going to be stored into nc files
"""

import site
import os
import time as clock
import numpy as np
import netCDF4 as nc4
from scipy import stats

	
########################################################
#0. setting variables
########################################################
# linuc path
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/historical/'


########################################################
#1. find all .nc files that need to be extracted to be read
########################################################
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()


p90 = np.empty((30,192,288)); p90[:] = np.nan   
p95 = np.empty((30,192,288)); p95[:] = np.nan 
for ensumble_member in range(0,28): 
	nc_f = text_content[ensumble_member][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	time = nc_fid.variables['date']
	TX = nc_fid.variables['TX'][:]
	TN = nc_fid.variables['TN'][:]
	
	size_data = np.shape(TX)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if (value <=1950 and value >=1961)]
	
	#  varibale initialization 
	cache = np.empty((92*30,size_data[1],size_data[2])); cache[:] = np.nan


	layer_output = 0 # the time dimenssion of the output variable
	
	for iyear in range(1961,1951):	
		"""
		#for precip
		layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+601][0]
		layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+831][0]
		cache[(iyear-1961)*92:(iyear+1-1961)*92,:,:] = precip[layer_b:layer_e+1,:,:]*24*60*60*1000 
		"""
	p90[ensumble_member,:,:] = np.percentile(cache,90,axis = 0)	
	p95[ensumble_member,:,:] = np.percentile(cache,95,axis = 0)	

r90pt = stats.nanmean(p90,axis=0)
r95pt = stats.nanmean(p95,axis=0)

file_name = input_path+'precip_9095_thre.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	

f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
r90pts = f.createVariable('r90pt',np.float32,('lat','lon'))
r95pts = f.createVariable('r95pt',np.float32,('lat','lon'))

latitudes[:] = lat
longitudes[:] = lon
r90pts[:] = r90pt
r95pts[:] = r95pt


f.description = '90 and 95 precipitation threshold 1961-1950 of 30 CESM ensumble mean'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'

latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'
longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

r90pts.standard_name = 'Daily Maximum 90 percentile'
r90pts.units = 'mm'
r95pts.standard_name = 'Daily Maximum 90 percentile'
r95pts.units = 'mm'

f.close()