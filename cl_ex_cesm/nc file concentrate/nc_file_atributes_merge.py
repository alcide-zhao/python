"""
Created on tuesday 22 2016
This is to Merge different varibales in different files to one nc file
"""
import numpy as np
import netCDF4 as nc4
import time as clock

data_path = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/' 
file = data_path + 'NCEP_tmax.2m.gauss.1971_2001.nc'
nc_fid = nc4.Dataset(file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
year_series = nc_fid.variables['time'][:]

TX = nc_fid.variables['tmax'][:]
# cwd = nc_fid.variables['cwd'][:]
# precptot = nc_fid.variables['precptot'][:]
# r10 = nc_fid.variables['r10'][:]
# r20 = nc_fid.variables['r20'][:]
# rnm = nc_fid.variables['rnm'][:]
# r95p = nc_fid.variables['r95p'][:]
# r99p = nc_fid.variables['r99p'][:]
# rx1day = nc_fid.variables['rx1day'][:]
# rx5day = nc_fid.variables['rx5day'][:]
# sdii = nc_fid.variables['sdii'][:]
nc_fid.close()

file = data_path + 'NCEP_tmin.2m.gauss.1971_2001.nc'
nc_fid = nc4.Dataset(file,mode='r')
TN = nc_fid.variables['tmin'][:]
# mean_precip = nc_fid.variables['mean_precip'][:]
# std_precip = nc_fid.variables['std_precip'][:]
# total_precip = nc_fid.variables['total_precip'][:]
nc_fid.close()

file_output = data_path+'TmTn_NCEP_2m_1971_2001.nc'

f = nc4.Dataset(file_output,'w', format='NETCDF4') #'w' stands for write

f.createDimension('time', len(year_series))
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
rx5days = f.createVariable('rx5day',np.float32,('time','lat','lon'))
sdiis = f.createVariable('sdii',np.float32,('time','lat','lon'))
# cwds = f.createVariable('cwd',np.float32,('time','lat','lon'))
# r99ps = f.createVariable('r99p',np.float32,('time','lat','lon'))
# total_precips = f.createVariable('total_precip',np.float32,('time','lat','lon'))
# rx1days = f.createVariable('rx1day',np.float32,('time','lat','lon'))
# r10s = f.createVariable('r10',np.float32,('time','lat','lon'))
# r20s = f.createVariable('r20',np.float32,('time','lat','lon'))
# rnms = f.createVariable('rnm',np.float32,('time','lat','lon'))
# cdds = f.createVariable('cdd',np.float32,('time','lat','lon'))
# r95ps = f.createVariable('r95p',np.float32,('time','lat','lon'))
# prcptots = f.createVariable('precptot',np.float32,('time','lat','lon'))
# r95prs = f.createVariable('r95pr',np.float32,('time','lat','lon'))
# mean_precips = f.createVariable('mean_precip',np.float32,('time','lat','lon'))
# std_precips = f.createVariable('std_precip',np.float32,('time','lat','lon'))
	
times[:] = year_series
latitudes[:] = lat
longitudes[:] = lon
rx5days[:] = rx5day
sdiis[:] =sdii
cwds[:] = cwd
r99ps[:] = r99p
total_precips[:] = total_precip
r95ps[:] = r95p
prcptots[:] = precptot
rx1days[:] = rx1day
r10s[:] = r10
r20s[:] = r20
rnms[:] = rnm
cdds[:] =cdd
r95prs[:] = np.divide(r95p,total_precip-r95p)
mean_precips[:] = mean_precip
std_precips[:] = std_precip

f.description = 'Precipitation extrme indeces calculated using the CESM model outputs'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()
