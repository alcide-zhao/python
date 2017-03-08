"""
This is ro concentrate the APHRO_MA V1101 daily precipitation with 0.50deg 
grids 1951_2007 daily precipitationd data into a single nc file
"""
import numpy as np
import netCDF4 as nc4
import time as clock

################################
# variables and functions      #
################################

input_path = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/APHRO/'

value_dic ={'lon':[],'lat':[]}
time_list =[]
precip_list = []

# leapyear judgement
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
	return leapyear

data_len = 0
for iyear in range(1951,2008):
	if (leapyear(iyear)):
		data_len += 366
	else:
		data_len += 365

year_series = np.empty((data_len)); precip_series = np.empty((data_len,140,180))
layer_e = 0
for iyear in range(1951,2008):
	file= input_path + 'APHRO_MA_050deg_V1101.'+str(iyear)+'.nc'
	nc_fid = nc4.Dataset(file,mode='r')
	lat = nc_fid.variables['latitude'][:]
	lon = nc_fid.variables['longitude'][:]
	precip = nc_fid.variables['precip'][:]
	units = nc_fid.variables['precip'].units
	long_name = nc_fid.variables['precip'].long_name
	missing_value = nc_fid.variables['precip'].missing_value
	time = nc_fid.variables['time'][:]+1; time=iyear*1000+time
	if (leapyear(iyear)):
		layer_b = layer_e
		layer_e = layer_b + 366
	else:
		layer_b = layer_e
		layer_e = layer_b + 365
	year_series[layer_b:layer_e] = time
	precip_series[layer_b:layer_e,:,:] = precip

file_name_out = input_path+'APHRO_MA_050deg_V1101_1951_2007.nc'


f = nc4.Dataset(file_name_out,'w', format='NETCDF4') #'w' stands for write
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('time',len(year_series))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
precips = f.createVariable('precip',np.float32,('time','lat','lon'))
	
times[:] = year_series
latitudes[:] = lat
longitudes[:] = lon

precips[:]= precip_series
precips.units=units
precips.missing_value=missing_value
precips.long_name =long_name

f.description = 'APHRO_MA V1101 daily precipitation with 0.50deg grids 1951_2007'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()

	
