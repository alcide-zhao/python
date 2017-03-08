"""
Created on tuesday 22 2016
This is to produce figures relating the indices
inputs  : extreme indices nc files
outputs : saved figures
"""

import site
import os
import numpy as np
from scipy import stats
import scipy.io as spio
from scipy.interpolate import interp2d
import scipy

import netCDF4 as nc4
import glob  
import time as clock
import matplotlib.pyplot as plt
from lib import *




def time_seeries_of_spatial_mean(time_s, time_e, time, data):
	"""
	for the input array, values of the interesed geographical domain is to be averaged
	the time seris of the averaged value is then produced 
	"""
	time_series = np.empty((time_e-time_s+1))
	
	layer_s = [layer for layer in range(len(time)) if time[layer] == time_s]
	layer_e = [layer for layer in range(len(time)) if time[layer] == time_e]
	layer_s = layer_s[0]
	layer_e = layer_e[0]
	for layer in range(layer_s , layer_e+1):
		# print data[layer,:,:]
		time_series[layer-layer_s] = stats.nanmean(stats.nanmean(data[layer,:,:]))
	return time_series

def spatial_mean_of_time_period(time_s, time_e, time, data):
	"""
	for each point inside the interested dommain of the input data, the mean of the 
	specified time period is calculated, the map of the mean of the time period is produced
	"""	
	layer_s = [layer for layer in range(len(time)) if time[layer] == time_s][0]
	layer_e = [layer for layer in range(len(time)) if time[layer] == time_e][0]
	spatial_mean = stats.nanmean(data[layer_s:layer_e,:,:], axis = 0)
	return spatial_mean
	

########################################################
# 0. setting variables
########################################################

# World ocean avs land map cover

# oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/normal_data/world.oceanmask.1440x720')['landocean']
# y = np.arange(0,192,0.266666666666666666)
# x = np.arange(0,288,0.2)
# xi, yi = np.meshgrid(x, y)
# y_new = range(192)
# x_new = range(288)
# x_newi,y_newi =np.meshgrid(y_new, x_new)
# f = scipy.interpolate.RectBivariateSpline(y, x, oceanmask)

# ocean_mask = np.empty((192,288))
# print np.shape(ocean_mask)
# ocean_mask = np.flipud(f(y_new,x_new))
# temperal = np.empty((192,288))
# print np.shape(temperal)
# temperal[:,0:144] = ocean_mask[:,144:288]
# temperal[:,144:288] = ocean_mask[:,0:144]
# ocean_mask = temperal
# ocean_mask[ocean_mask<125] =0
# ocean_mask[ocean_mask>=125] =1
# ocean_mask[0,:]=1


# file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/ensumble_mean_Precip_extremes_global_2006_2080.nc'
# nc_fid = nc4.Dataset(file_name,mode='r')
# lats = nc_fid.variables['lat'][:]
# lons = nc_fid.variables['lon'][:]
# _,_,ocean_mask=range_clip(60,150,0,55,lons,lats,ocean_mask)
# print np.shape(ocean_mask)
# plt.figure(figsize = (20,10))

# plt.imshow(ocean_mask, origin = 'lower')
# plt.show()
# file_name='/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/landoceanmask_asia.mat'
# mdict={'landoceanmask':ocean_mask}
# spio.savemat(file_name, mdict)

# linuc path
# input_path = file_path + file_name
time_s = 2006
time_e = 2100

oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']

#########Precip
output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
date_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/fixa/*.nc'   
file_name_t = output_path+'time_series_Precip_extremes_asia_'+str(time_s)+'_'+str(time_e)+'.nc'
file_name_s = output_path+'ensumble_mean_PEI_global_'+str(time_s)+'_'+str(time_e)+'_fixa.nc'	
print file_name_s
#########Precip
"""
#########TEMP
output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
date_path = '/exports/csce/datastore/geos/users/s1667168/CESM/prep/8000/*.nc'   
file_name_t = output_path+'time_series_Prep_extremes_global_'+str(time_s)+'_'+str(time_e)+'.nc'
file_name_s = output_path+'Spatial_ensumble_mean_Prep_extremes_global_'+str(time_s)+'_'+str(time_e)+'.nc'	
#########TEMP
"""

#######################################################
# 1. data preprocessing                               #
#######################################################
files=sorted(glob.glob(date_path))
# print files


#########Precip
att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[],\
'precptot':[],'total_precip':[]} #'mean_precip':[],'std_precip':[]
att_list = ['rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot','total_precip'] #'mean_precip','std_precip'
#########Precip

#########Temp
# att_dic = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_list = ['txx', 'txn', 'tnx', 'tnn', 'dtr', 'fd0', 'su25','id0', 'tr20', 'tn10p', 'tn90p','tx10p','tx90p']

#########Temp
"""
#######################################################
# 1.1 TIME SERIES OF 30 ENSUMBLES                     #
#######################################################
for att_name in att_list:
	ensumble_series = 0
	time_series = np.empty([len(files),time_e-time_s+1])
	for file in files: 
		nc_fid = nc4.Dataset(file,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		year_series = nc_fid.variables['time'][:]
		att_value = nc_fid.variables[att_name][:]
		att_value=np.multiply(att_value,oceanmask)
		nc_fid.close()
		lons,lats,att_clipped = range_clip(0,360,-90,90,lon,lat,att_value) # Global
		# lons,lats,att_clipped = range_clip(60,150,0,55,lon,lat,att_value)
		# lon_clipped,lat_clipped,att_clipped = range_clip(70,130,20,55,lon,lat,att_value) # China
		time_series[ensumble_series,:] = time_seeries_of_spatial_mean(time_s, time_e, year_series, att_clipped)
		# print ensumble_series
		# print file
		ensumble_series = ensumble_series+1		
	att_dic[att_name] = time_series
	# print att_name+'_finished'

year_series = range(time_s,time_e+1)
f = nc4.Dataset(file_name_t,'w', format='NETCDF4') #'w' stands for write

f.createDimension('time', time_e-time_s+1)
f.createDimension('ensumble_no', len(files))
ensumble_nos = f.createVariable('ensumble_no',np.float64, ('ensumble_no'))
times = f.createVariable('time',np.float64, ('time'))


#########Precip
rx1days = f.createVariable('rx1day',np.float32,('ensumble_no','time'))
rx5days = f.createVariable('rx5day',np.float32,('ensumble_no','time'))
sdiis = f.createVariable('sdii',np.float32,('ensumble_no','time'))
r10s = f.createVariable('r10',np.float32,('ensumble_no','time'))
r20s = f.createVariable('r20',np.float32,('ensumble_no','time'))
rnms = f.createVariable('rnm',np.float32,('ensumble_no','time'))
cdds = f.createVariable('cdd',np.float32,('ensumble_no','time'))
cwds = f.createVariable('cwd',np.float32,('ensumble_no','time'))
r95ps = f.createVariable('r95p',np.float32,('ensumble_no','time'))
r99ps = f.createVariable('r99p',np.float32,('ensumble_no','time'))
prcptots = f.createVariable('precptot',np.float32,('ensumble_no','time'))

times[:] = year_series
ensumble_nos[:] = range(len(files))
rx1days[:] = att_dic['rx1day']
rx5days[:] = att_dic['rx5day']
sdiis[:] = att_dic['sdii']
r10s[:] = att_dic['r10']
r20s[:] = att_dic['r20']
rnms[:] = att_dic['rnm']
cdds[:] = att_dic['cdd']
cwds[:] = att_dic['cwd']
r95ps[:] = att_dic['r95p']
r99ps[:] =  att_dic['r99p']
prcptots[:] =  att_dic['precptot']


#########temp

txxs = f.createVariable('txx',np.float32,('ensumble_no','time'))
txns = f.createVariable('txn',np.float32,('ensumble_no','time'))
tnxs = f.createVariable('tnx',np.float32,('ensumble_no','time'))
tnns = f.createVariable('tnn',np.float32,('ensumble_no','time'))
dtrs = f.createVariable('dtr',np.float32,('ensumble_no','time'))
fd0s = f.createVariable('fd0',np.float32,('ensumble_no','time'))
su25s = f.createVariable('su25',np.float32,('ensumble_no','time'))
id0s = f.createVariable('id0',np.float32,('ensumble_no','time'))
tr20s = f.createVariable('tr20',np.float32,('ensumble_no','time'))
tn10ps = f.createVariable('tn10p',np.float32,('ensumble_no','time'))
tn90ps = f.createVariable('tn90p',np.float32,('ensumble_no','time'))
tx10ps = f.createVariable('tx10p',np.float32,('ensumble_no','time'))
tx90ps = f.createVariable('tx90p',np.float32,('ensumble_no','time'))

times[:] = year_series
ensumble_nos[:] = range(len(files))
txxs[:] = att_dic['txx']
txns[:] = att_dic['txn']
tnxs[:] = att_dic['tnx']
tnns[:] = att_dic['tnn']
dtrs[:] = att_dic['dtr']
fd0s[:] = att_dic['fd0']
su25s[:] = att_dic['su25']
id0s[:] = att_dic['id0']
tr20s[:] = att_dic['tr20']
tn10ps[:] =  att_dic['tn10p']
tn90ps[:] =  att_dic['tn90p']
tx10ps[:] =  att_dic['tx10p']
tx90ps[:] =  att_dic['tx90p']

f.description = 'time series of temperature extrme indeces calculated using the CESM model outputs'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()
"""


######################
#ENSUMBLE MEAN
######################

for att_name in att_list:	
	att_ensumble_mean = np.empty((time_e-time_s+1,192,288))
	for iyear in range(time_s,time_e+1):
		att_value = np.empty((len(files),192,288))
		ensumble_no = 0
		for file in files:
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			year_series = nc_fid.variables['time'][:]
			att_value[ensumble_no,:,:] = nc_fid.variables[att_name][iyear-time_s]
			ensumble_no =ensumble_no+1
			nc_fid.close()
		att_ensumble_mean[iyear-time_s,:,:] = stats.nanmean(att_value,axis = 0)
		# print stats.nanmean(att_value,axis = 0)
	att_dic[att_name] = att_ensumble_mean
	
	
	
f = nc4.Dataset(file_name_s,'w', format='NETCDF4') #'w' stands for write

f.createDimension('time', len(year_series))
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))

#########Precip
rx1days = f.createVariable('rx1day',np.float32,('time','lat','lon'))
rx5days = f.createVariable('rx5day',np.float32,('time','lat','lon'))
sdiis = f.createVariable('sdii',np.float32,('time','lat','lon'))
r10s = f.createVariable('r10',np.float32,('time','lat','lon'))
r20s = f.createVariable('r20',np.float32,('time','lat','lon'))
rnms = f.createVariable('rnm',np.float32,('time','lat','lon'))
cdds = f.createVariable('cdd',np.float32,('time','lat','lon'))
cwds = f.createVariable('cwd',np.float32,('time','lat','lon'))
r95ps = f.createVariable('r95p',np.float32,('time','lat','lon'))
r99ps = f.createVariable('r99p',np.float32,('time','lat','lon'))
prcptots = f.createVariable('precptot',np.float32,('time','lat','lon'))
mean_precips = f.createVariable('mean_precip',np.float32,('time','lat','lon'))
# mean_precips = f.createVariable('mean_precip',np.float32,('time','lat','lon'))
# std_precips = f.createVariable('std_precip',np.float32,('time','lat','lon'))
	
times[:] = year_series
latitudes[:] = lat
longitudes[:] = lon
rx1days[:] = att_dic['rx1day']
rx5days[:] = att_dic['rx5day']
sdiis[:] = att_dic['sdii']
r10s[:] = att_dic['r10']
r20s[:] = att_dic['r20']
rnms[:] = att_dic['rnm']
cdds[:] = att_dic['cdd']
cwds[:] = att_dic['cwd']
r95ps[:] = att_dic['r95p']
r99ps[:] =  att_dic['r99p']
prcptots[:] =  att_dic['precptot']
# mean_precips[:] =  att_dic['mean_precip']
# std_precips[:] =  att_dic['std_precip']
mean_precips[:] =  att_dic['total_precip']/91.0   # 91days of JJA
"""
#########
temp
#########

txxs = f.createVariable('txx',np.float32,('time','lat','lon'))
txns = f.createVariable('txn',np.float32,('time','lat','lon'))
tnxs = f.createVariable('tnx',np.float32,('time','lat','lon'))
tnns = f.createVariable('tnn',np.float32,('time','lat','lon'))
dtrs = f.createVariable('dtr',np.float32,('time','lat','lon'))
fd0s = f.createVariable('fd0',np.float32,('time','lat','lon'))
su25s = f.createVariable('su25',np.float32,('time','lat','lon'))
id0s = f.createVariable('id0',np.float32,('time','lat','lon'))
tr20s = f.createVariable('tr20',np.float32,('time','lat','lon'))
tn10ps = f.createVariable('tn10p',np.float32,('time','lat','lon'))
tn90ps = f.createVariable('tn90p',np.float32,('time','lat','lon'))
tx10ps = f.createVariable('tx10p',np.float32,('time','lat','lon'))
tx90ps = f.createVariable('tx90p',np.float32,('time','lat','lon'))

times[:] = year_series
latitudes[:] = lat
longitudes[:] = lon
txxs[:] = att_dic['txx']
txns[:] = att_dic['txn']
tnxs[:] = att_dic['tnx']
tnns[:] = att_dic['tnn']
dtrs[:] = att_dic['dtr']
fd0s[:] = att_dic['fd0']
su25s[:] = att_dic['su25']
id0s[:] = att_dic['id0']
tr20s[:] = att_dic['tr20']
tn10ps[:] =  att_dic['tn10p']
tn90ps[:] =  att_dic['tn90p']
tx10ps[:] =  att_dic['tx10p']
tx90ps[:] =  att_dic['tx90p']
"""

f.description = 'Precipitation extrme indeces calculated using the CESM model outputs'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()


