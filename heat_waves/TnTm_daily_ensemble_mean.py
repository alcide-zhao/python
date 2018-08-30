import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  

#############################################
# 1 ensumble mean 3D data
#############################################
def JJAS_daily_extract(data,date):
	year= [ iyear/10000 for iyear in map(int,date)]
	year_series = [value for value in np.unique(year) if value <=2100]
	data_JJAS = np.empty((len(year_series)*122,np.shape(data)[1],np.shape(data)[2]))
	time_JJAS = np.empty((len(year_series)*122))
	for iyear in year_series:
		layer_b = [layer for layer in range(len(date)) if date[layer] == iyear*10000+601][0]
		layer_e = layer_b+122
		# print iyear
		# print date[layer_b]
		# print date[layer_e]
		data_JJAS[(iyear-year_series[0])*122:(iyear-year_series[0]+1)*122,:,:] = data[layer_b:layer_e,:,:]
		time_JJAS[(iyear-year_series[0])*122:(iyear-year_series[0]+1)*122] = date[layer_b:layer_e]
		iyear = iyear+1
	return data_JJAS, time_JJAS

def ensumble_mean_of_CESM_T(files,variable):
	att_ensumble_mean = np.zeros((34675,192,288)); 
	ensumble_no = 0
	for file in files:
		# print file
		print file
		nc_fid = nc4.Dataset(file,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		date = nc_fid.variables['date'][:]
		value = nc_fid.variables[variable][0:34675,:,:]
		units = nc_fid.variables[variable].units
		long_name = nc_fid.variables[variable].long_name
		ensumble_no=ensumble_no+1
		# print ensumble_no
		att_ensumble_mean = att_ensumble_mean+value
		nc_fid.close()
	att_ensumble_mean = att_ensumble_mean/ensumble_no
	data_JJAS, time_JJAS = JJAS_daily_extract(att_ensumble_mean,date)
	return lon,lat,time_JJAS,data_JJAS,units,long_name

variable = 	'TREFHTMN'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/rcp85/'+variable+'/*.nc'
files=sorted(glob.glob(input_path))
lon,lat,time_JJAS,TN_JJAS,units,long_name = ensumble_mean_of_CESM_T(files,variable)
variable = 	'TREFHTMX'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/rcp85/'+variable+'/*.nc'
files=sorted(glob.glob(input_path))
lon,lat,time_JJAS,TX_JJAS,units,long_name = ensumble_mean_of_CESM_T(files,variable)

######################################
#1.1 writting the results into nc files
######################################
output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/'
file_name_out = output_path+'ensumble_mean_daily_TmTn_2006_2100_JJAS_rcp85.nc'

f = nc4.Dataset(file_name_out,'w', format='NETCDF4') #'w' stands for write
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('time', len(time_JJAS))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
TXs = f.createVariable('TX',np.float32,('time','lat','lon'))
TNs = f.createVariable('TN',np.float32,('time','lat','lon'))
	
times[:] = time_JJAS
latitudes[:] = lat
longitudes[:] = lon
TXs[:]= TX_JJAS
TNs[:]= TN_JJAS

times.long_name = 'Month'
times.units = 'Month'
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'
longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'
TXs.units=units
TXs.long_name =long_name
TNs.units=units
TNs.long_name =long_name
f.description = 'Ensemble mean of CESM daily TX and TN'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()