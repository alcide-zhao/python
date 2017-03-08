
# -*- coding: utf-8 -*-
"""
This is to caculate the 95 and 99 threshold of the global temp and precipi
for computing the extremes. the baseline fot the threshold is 1961-1990
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
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/PRECT/'


########################################################
#1. find all .nc files that need to be extracted to be read
########################################################
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()


p95 = np.empty((30,192,288)); p95[:] = np.nan   
p99 = np.empty((30,192,288)); p99[:] = np.nan 
for ensumble_member in range(0,29): 
	nc_f = text_content[ensumble_member][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	time = nc_fid.variables['date']
	precip = nc_fid.variables['PRECT'][:,:,:]
	precip[precip <0]=np.nan
	
	#######################
	# 2.2 calculation 
	#######################
	size_data = np.shape(precip)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if (value <=1990 and value >=1961)]
	
	#  varibale initialization 
	cache = np.empty((92*30,size_data[1],size_data[2])); cache[:] = np.nan


	layer_output = 0 # the time dimenssion of the output variable
	
	for iyear in range(1961,1991):	
		layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+601][0]
		layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+831][0]
		cache[(iyear-1961)*92:(iyear+1-1961)*92,:,:] = precip[layer_b:layer_e+1,:,:]*24*60*60*1000 
		# print iyear
	p95[ensumble_member,:,:] = np.percentile(cache,95,axis = 0)	
	p99[ensumble_member,:,:] = np.percentile(cache,99,axis = 0)	

r95pt = stats.nanmean(p95,axis=0)
r99pt = stats.nanmean(p99,axis=0)

file_name = input_path+'precip_9599_thre.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	

f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
r95pts = f.createVariable('r95pt',np.float32,('lat','lon'))
r99pts = f.createVariable('r99pt',np.float32,('lat','lon'))

latitudes[:] = lat
longitudes[:] = lon
r95pts[:] = r95pt
r99pts[:] = r99pt


f.description = '95 and 99 precipitation threshold 1961-1990 of 30 CESM ensumble mean'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'

latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'
longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

r95pts.standard_name = 'precipitation 95 percentile'
r95pts.units = 'mm'
r99pts.standard_name = 'precipitation 99 percentile'
r99pts.units = 'mm'

f.close()