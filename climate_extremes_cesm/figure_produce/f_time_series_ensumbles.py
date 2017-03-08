
import site
import os
import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt

from scipy import stats
from matplotlib import rcParams
from scipy.interpolate import interp2d
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
	
# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/'

oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
time_s = 2006
time_e = 2100

precipitation extremes
att_ref = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
att_list=['rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot']
unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

# temperature extremes
# att_ref stores the reference value which is used to remove climatology
# att_ref = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_dic = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_list = ['txx', 'txn', 'tnx', 'tnn', 'dtr', 'fd0', 'su25','id0','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
# unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}


time_s = 2006
time_e = 2100
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'India':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}

region = rergion_dic['SCI'][:]

	
#######################################################
# 1. time series of ensumble mean                     #
#######################################################
# spst = [2,3]
# input_path = '/home/s1667168/scratch/extremes_indices/temp/'
# plt.figure(1, facecolor='White',figsize=[20,16])

# file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100_Rcp8.5.nc'
# nc_fid = nc4.Dataset(input_path+file_name,mode='r')
# sub_plot=1

# 
# region = region_dic['GLOBE'][:]
# for att_name in att_list:
	# time = nc_fid.variables['time'][:]
	# lons = nc_fid.variables['lon'][:]
	# lats = nc_fid.variables['lat'][:]
	# att_value = nc_fid.variables[att_name][:]
	# if (att_name == 'tn10p'or att_name == 'tn90p' or att_name == 'tx10p' or att_name == 'tx90p'):
		# att_value = att_value*100  # 92 is the length of the days of JJA. 
	# # print np.shape(att_value)
	# att_value = np.multiply(att_value,oceanmask)
	# lons,lats,att_value = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	# time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value)
	
	# att_ref[att_name] = time_series[0]
	# # time_series = time_series -att_ref[att_name]
	# # print time_series
	# ax1=plt.subplot(spst[0],spst[1],sub_plot)
	# # print np.shape(time)
	# # print np.shape(time_series)
	# ax1.plot(time,time_series,'-',color="red",linewidth=2,label='RCP8.5')
	# sub_plot = sub_plot+1
# nc_fid.close()


# file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100_fixA.nc'
# nc_fid = nc4.Dataset(input_path+file_name,mode='r')
# plt.figure(1, facecolor='White',figsize=[20,16])	
# sub_plot=1	
# for att_name in att_list:
	# time = nc_fid.variables['time'][:]
	# lons = nc_fid.variables['lon'][:]
	# lats = nc_fid.variables['lat'][:]
	# att_value = nc_fid.variables[att_name][:]
	# if (att_name == 'tn10p'or att_name == 'tn90p' or att_name == 'tx10p' or att_name == 'tx90p'):
		# att_value = att_value*100
	# att_value = np.multiply(att_value,oceanmask)
	# lons,lats,att_value = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	# time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value)
	# att_ref[att_name] = time_series[0]
	# # time_series = time_series -att_ref[att_name]
	
	# ax1=plt.subplot(spst[0],spst[1],sub_plot)
	# ax1.plot(time,time_series,'-',color="blue",linewidth=2,label='RCP8.5_fixA')
	# ax1.set_xlabel('Year',fontsize=30)
	# ax1.set_xlim([2000,2100])
	# # ax1.set_ylim([4,7.5])
	# ax1.set_ylabel(unit_dic[att_name].title(),fontsize=30)
	# # ax1.set(aspect=10)
	# ax1.set_title(att_name.title(),fontsize=30)
	# ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	# ax1.tick_params(labelsize=25,axis='both', which='major', pad=30)
	# sub_plot = sub_plot+1
# nc_fid.close()

# legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(0.95, 0.8),fontsize=20)	
# # plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/time_series_global.tif')
# plt.tight_layout()

# plt.show()

######################################################
# time series of the ensumble range
######################################################
input_path = '/home/s1667168/scratch/extremes_indices/temp/'
plt.figure(1, facecolor='White',figsize=[10,5])
file_name = 'time_series_Temp_extremes_global_2006_2100_Rcp8.5.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
sub_plot=1
att_list = ['txx','tr20']
for att_name in att_list:
	time_series_range = np.empty([3,time_e-time_s+1])
	time = nc_fid.variables['time'][:]
	att_value = nc_fid.variables[att_name][:]
	# att_value = np.multiply(att_value,oceanmask)
	time_series_range[1,:] = stats.nanmean(att_value,axis = 0)
	time_series_range[0,:] = np.nanmax(att_value,axis = 0)
	time_series_range[2,:] = np.nanmin(att_value,axis = 0)	
	att_ref[att_name] = time_series_range[1,0]
	time_series_range[1,:]=time_series_range[1,:]-att_ref[att_name]
	time_series_range[0,:]=time_series_range[0,:]-att_ref[att_name]
	time_series_range[2,:] = time_series_range[2,:] -att_ref[att_name]
	att_dic[att_name] = time_series_range
	
	ax1=plt.subplot(1,2,sub_plot)
	ax1.plot(time,time_series_range[1,:],'-',color="red",linewidth=2,label='RCP8.5')
	ax1.fill_between(time,time_series_range[0,:],time_series_range[2,:], where=time_series_range[0,:]>=time_series_range[2,:], color="red", alpha = 0.2) 
	sub_plot = sub_plot+1

nc_fid.close()


file_name = 'time_series_Temp_extremes_global_2006_2100_fixA.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

sub_plot=1

nc_fid = nc4.Dataset(input_path+file_name,mode='r')

plt.figure(1, facecolor='White',figsize=[20,16])


for att_name in att_list:
	time_series_range = np.empty([3,time_e-time_s+1])
	time = nc_fid.variables['time'][:]
	att_value = nc_fid.variables[att_name][:]
	# att_value = np.multiply(att_value,oceanmask)
	time_series_range[1,:] = stats.nanmean(att_value,axis = 0)
	time_series_range[0,:] = np.nanmax(att_value,axis = 0)
	time_series_range[2,:] = np.nanmin(att_value,axis = 0)	
	att_ref[att_name] = time_series_range[1,0]
	time_series_range[1,:]=time_series_range[1,:]-att_ref[att_name]
	time_series_range[0,:]=time_series_range[0,:]-att_ref[att_name]
	time_series_range[2,:] = time_series_range[2,:] -att_ref[att_name]
	att_dic[att_name] = time_series_range
	
	ax1=plt.subplot(1,2,sub_plot)
	ax1.plot(time,time_series_range[1,:],'-',color="blue",linewidth=2,label='RCP8.5_fixA')
	ax1.fill_between(time,time_series_range[0,:],time_series_range[2,:], where=time_series_range[0,:]>=time_series_range[2,:], color="blue", alpha = 0.2) 
	ax1.set_title(att_name.title(),fontsize=30)
	ax1.set_xlabel('Year',fontsize=30)
	ax1.set_xlim([2000,2100])
	# ax1.set_ylim([4,7.5])
	ax1.set_ylabel(unit_dic[att_name].title(),fontsize=30)
	# ax1.set(aspect=10)
	sub_plot = sub_plot+1
	ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	ax1.tick_params(labelsize=25,axis='both', which='major', pad=30)
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.95),fontsize=20)	

plt.tight_layout()

plt.show()