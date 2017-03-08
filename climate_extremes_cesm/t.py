import netCDF4 as nc4
import numpy as np
import math
import time as clock

variable='FSDS'
output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/'
file_name = output_path+'ensumble_mean_'+variable+'_200602_210101.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time_series = nc_fid.variables['time'][:]

rcp85 = nc_fid.variables['rcp85'][:]
rcp85_MV = nc_fid.variables['rcp85'].missing_value
rcp85_fixA = nc_fid.variables['rcp85_fixA'][:]
RCP85_fixA_MV = nc_fid.variables['rcp85_fixA'].missing_value
units =  nc_fid.variables['rcp85_fixA'].units
long_name =  'Atmospheric Incident Solar Radiation'
# rcp85_MV = np.nan
# units = '#/m2'
nc_fid.close()


f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('time',len(time_series))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))

rcp85s = f.createVariable('rcp85',np.float32,('time','lat','lon'))
rcp85_fixAs = f.createVariable('rcp85_fixA',np.float32,('time','lat','lon'))
	
times[:] = time_series
latitudes[:] = lat
longitudes[:] = lon
rcp85s[:]= rcp85
rcp85_fixAs[:]= rcp85_fixA

times.long_name = 'Month'
times.units = 'Month'
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'
longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

rcp85s.units=units
rcp85s.missing_value=rcp85_MV
rcp85s.long_name =long_name
rcp85_fixAs.units=units
rcp85_fixAs.missing_value=rcp85_MV
rcp85_fixAs.long_name =long_name

f.description = 'Time series of nsumble mean of '+variable+' calculated using the CESM model outputs'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()

