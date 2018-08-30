import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob
import os  


def JJAS_daily_CESM_TmTn(date,TREFHTMN,TREFHTMX):
	# print date
	year= [ iyear/10000 for iyear in map(int,date)]
	year_series = [value for value in np.unique(year) if value <=2100]
	TREFHTMN_JJAS = np.empty((len(year_series),365,np.shape(TREFHTMN)[1],np.shape(TREFHTMN)[2]))
	TREFHTMX_JJAS = np.empty((len(year_series),365,np.shape(TREFHTMN)[1],np.shape(TREFHTMN)[2]))
	for iyear in year_series:
		# print iyear
		layer_b = [layer for layer in range(len(date)) if date[layer] == iyear*10000+101][0]
		layer_e = layer_b+365
		TREFHTMN_JJAS[(iyear-year_series[0]),:,:,:] = TREFHTMN[layer_b:layer_e,:,:]
		TREFHTMX_JJAS[(iyear-year_series[0]),:,:,:] = TREFHTMX[layer_b:layer_e,:,:]
		iyear = iyear+1
	return year_series,TREFHTMN_JJAS,TREFHTMX_JJAS	


#####################
##     Main        ##
#####################	
input_path_mn = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/rcp45/TREFHTMN/'
input_path_mx = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/rcp45/TREFHTMX/'

os.system('find ' + os.path.abspath(input_path_mn) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path_mn + '/file_list.txt')
text_file_mn = open(input_path_mn + '/file_list.txt', "r")
text_content_mn = text_file_mn.readlines()

os.system('find ' + os.path.abspath(input_path_mx) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path_mx + '/file_list.txt')
text_file_mx = open(input_path_mx + '/file_list.txt', "r")
text_content_mx = text_file_mx.readlines()

for en_no in range(0,(len(text_content_mn))):
	# TREFHTMN
	nc_f_mn = text_content_mn[en_no][:-1]
	nc_fid_mn = nc4.Dataset(nc_f_mn,mode='r')
	lat = nc_fid_mn.variables['lat'][:]
	lon = nc_fid_mn.variables['lon'][:]
	date = nc_fid_mn.variables['date'][:]
	TREFHTMN = nc_fid_mn.variables['TREFHTMN'][:]-273.15
	TREFHTMN_ln= nc_fid_mn.variables['TREFHTMN'].long_name
	nc_fid_mn.close()
	# TREFHTMX
	nc_f_mx = text_content_mx[en_no][:-1]
	nc_fid_mx = nc4.Dataset(nc_f_mx,mode='r')	
	TREFHTMX = nc_fid_mx.variables['TREFHTMX'][:]-273.15
	TREFHTMX_ln= nc_fid_mx.variables['TREFHTMX'].long_name
	year_series,TREFHTMN_JJAS,TREFHTMX_JJAS = JJAS_daily_CESM_TmTn(date,TREFHTMN,TREFHTMX)
	nc_fid_mx.close()
	en_no = en_no+1
	#writting the results into nc files
	output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/rcp45/'
	file_name_out = output_path+'TmTn_rcp45_'+str(en_no)+'.nc'
	print file_name_out
	f = nc4.Dataset(file_name_out,'w', format='NETCDF4') #'w' stands for write
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	f.createDimension('year', len(year_series))
	f.createDimension('day', 365)

	years = f.createVariable('year',np.float64, ('year'))
	days = f.createVariable('day',np.float32, ('day'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	TXs = f.createVariable('TX',np.float32,('year','day','lat','lon'))
	TNs = f.createVariable('TN',np.float32,('year','day','lat','lon'))
		
	years[:] = year_series
	days[:] = range(1,366)
	latitudes[:] = lat
	longitudes[:] = lon
	TXs[:]= TREFHTMX_JJAS
	TNs[:]= TREFHTMN_JJAS

	latitudes.long_name = 'latitude'
	latitudes.units = 'degree_north'
	longitudes.long_name = 'longitude'
	longitudes.units = 'degree_east'
	TXs.units='degree Celsius'
	TXs.long_name =TREFHTMX_ln
	TNs.units='degree Celsius'
	TNs.long_name =TREFHTMN_ln
	f.description = 'CESM daily TX and TN'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	f.close()
	os.remove(nc_f_mn); os.remove(nc_f_mx)
	










