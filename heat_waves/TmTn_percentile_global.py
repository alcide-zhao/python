"""
This is to caculate the 5th, 10th, 90th and 95th  percentile threshold of the TX and TN
the baseline fot the threshold is 1961-1995
Ensemble mean are output in the TmTn_percentile_365.nc file
"""

import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import matplotlib.pyplot as plt

# linuc path
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/his/'
#os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()

for en_no in range(12,23): 
	cache_TX= np.empty((20,395,192,288));cache_TX[:] = np.nan 
	cache_TN= np.empty((20,395,192,288));cache_TN[:] = np.nan
	nc_f = text_content[en_no][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	TX = nc_fid.variables['TX'][66:86,:,:,:];
	cache_TX[:,0:15,:,:] = TX[:,350:365,:,:];
	cache_TX[:,15:380,:,:] = TX[:];
	cache_TX[:,380:395,:,:]= TX[:,0:15,:,:];
	del TX
	TN = nc_fid.variables['TN'][66:86,:,:,:]
	cache_TN[:,0:15,:,:] = TN[:,350:365,:,:];
	cache_TN[:,15:380,:,:] = TN[:];
	cache_TN[:,380:395,:,:]= TN[:,0:15,:,:];
	del TN
	TX95P = np.empty((365,192,288));TX95P[:] = np.nan   
	TN95P = np.empty((365,192,288));TN95P[:] = np.nan
	TX5P = np.empty((365,192,288));TX5P[:] = np.nan   
	TN5P = np.empty((365,192,288));TN5P[:] = np.nan 
	for iday in range(15,380):
		for lat_index in range(len(lat)):
			for lon_index in range(len(lon)):
				data_cache =  np.reshape(cache_TX[:,iday-15:iday+16,lat_index,lon_index],20*31)
				TX95P[iday-15,lat_index,lon_index] = np.nanpercentile(data_cache,95,axis = 0)
				TX5P[iday-15,lat_index,lon_index] = np.nanpercentile(data_cache,5,axis = 0)			
				data_cache =  np.reshape(cache_TN[:,iday-15:iday+16,lat_index,lon_index],20*31)
				TN95P[iday-15,lat_index,lon_index] = np.nanpercentile(data_cache,95,axis = 0)
				TN5P[iday-15,lat_index,lon_index] = np.nanpercentile(data_cache,5,axis = 0)
				del data_cache
	del cache_TN,cache_TX
	## write into .nc file
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/threshold/TmTn_percentile_calender'+str(['_0'+str(en_no+1) if en_no<9 else '_'+str(en_no+1)])[2:5]+'.nc'
	f = nc4.Dataset(file_name,'w', format='NETCDF4') 
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	f.createDimension('day', 365)

	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	days = f.createVariable('day',np.float32, ('day'))
	TX95ps = f.createVariable('TX95P',np.float32,('day','lat','lon'))
	TN95ps = f.createVariable('TN95P',np.float32,('day','lat','lon'))
	TX5ps = f.createVariable('TX5P',np.float32,('day','lat','lon'))
	TN5ps = f.createVariable('TN5P',np.float32,('day','lat','lon'))

	latitudes[:] = lat
	longitudes[:] = lon
	days[:]=range(1,366)
	TX95ps[:] = TX95P
	TN95ps[:] = TN95P
	TX5ps[:] = TX5P
	TN5ps[:] = TN5P

	f.description = 'TxTn 5, 10, 90, 95 percentiles for 1961-1950 of 30 CESM ensumble mean'
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	f.close()
	nc_fid.close()
	print en_no+' finished'
"""

##########
#NCEP  #
##########

# linuc path
input_path = '/exports/csce/datastore/geos/users/s1667168/obs/NCEP/'

nc_f = input_path+'NCEP_TXTN.2m.gauss.1948_2017.nc'
nc_fid = nc4.Dataset(nc_f,mode='r')
TX = nc_fid.variables['TX'][13:43,:,:,:]
nc_fid.close()

nc_f = input_path+'NCEP_TXTN.2m.gauss.1948_2017.nc'
nc_fid = nc4.Dataset(nc_f,mode='r')
lat = nc_fid.variables['lat']
lon = nc_fid.variables['lon']
TN = nc_fid.variables['TN'][13:43,:,:,:]

cache = np.reshape(TX,(30*365,94,192))
TX95P = np.nanpercentile(cache,95,axis = 0)
TX90P = np.nanpercentile(cache,90,axis = 0)
TX10P = np.nanpercentile(cache,10,axis = 0)
TX5P  = np.nanpercentile(cache,5,axis = 0)
cache = np.reshape(TN,(30*365,94,192))
TN95P = np.nanpercentile(cache,95,axis = 0)	
TN90P = np.nanpercentile(cache,90,axis = 0)
TN10P = np.nanpercentile(cache,10,axis = 0)
TN5P  = np.nanpercentile(cache,5,axis = 0)

## write into .nc file
file_name = input_path+'NCEP_TmTn_percentile_365_nan.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4') 
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
TX95ps = f.createVariable('TX95P',np.float32,('lat','lon'))
TN95ps = f.createVariable('TN95P',np.float32,('lat','lon'))
TX90ps = f.createVariable('TX90P',np.float32,('lat','lon'))
TN90ps = f.createVariable('TN90P',np.float32,('lat','lon'))
TX10ps = f.createVariable('TX10P',np.float32,('lat','lon'))
TN10ps = f.createVariable('TN10P',np.float32,('lat','lon'))
TX5ps = f.createVariable('TX5P',np.float32,('lat','lon'))
TN5ps = f.createVariable('TN5P',np.float32,('lat','lon'))

latitudes[:] = lat
longitudes[:] = lon
TX95ps[:] = TX95P
TN95ps[:] = TN95P
TX90ps[:] = TX90P
TN90ps[:] = TN90P
TX10ps[:] = TX10P
TN10ps[:] = TN10P
TX5ps[:] = TX5P
TN5ps[:] = TN5P

f.description = 'TxTn 5, 10, 90, 95 percentiles for 1961-1990 NCEP'
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()
"""
