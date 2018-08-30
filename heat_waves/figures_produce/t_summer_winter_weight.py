# -*- coding: utf-8 -*-
'''
This is to plot the boxplots of HW and CS characteristics over 7 different countries (regions)
His(1986-2005), rcp(2081-2100) and FiXA (2081-2100) are compared
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *



def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# plt.imshow(area);plt.show()
	return area

def mask_match(country_key,lon,lat):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask= ocean_mask[country_key][:]
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	# print np.unique(mask)
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
	mask = f(lon, lat);
	mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid (lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
	mask=np.multiply(mask,area);  #crop the interested region
	mask=np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
	return mask
	

def summer_winter_weight_print(country_key):
	def annual_read(scenario):
		file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/' 
		nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_'+scenario+'.nc',mode='r')
		lon = nc_fid.variables['lon'][:]
		lat = nc_fid.variables['lat'][:]
		year = nc_fid.variables['year'][:]
		mask = mask_match(country_key,lon,lat)
		HWDC =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWDC'][0,-51:-1,:,:],axis=0),mask),axis=1),axis=0)
		HWDuM =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWDuM'][0,-51:-1,:,:],axis=0),mask),axis=1),axis=0)
		HWI_NO =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWI_NO'][0,-51:-1,:,:],axis=0),mask),axis=1),axis=0)
		HWInM =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWInM'][0,-51:-1,:,:],axis=0),mask),axis=1),axis=0)
		result = np.empty((2,4));result[0,:]=np.array([HWInM,HWDuM,HWI_NO,HWDC]);result[1,:]=np.array([HWInM,HWDuM,HWI_NO,HWDC]);
		# print result
		return result
	def summer_winter_read(scenario):
		def summer_winter_hemisphere_swap(variable,country_key):
			summer_id='HW_MJJASO_'+variable;winter_id='HW_NDJFMA_'+variable;
			summer = np.empty((50,192,288));winter = np.empty((50,192,288));
			file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr_season/'  
			nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_'+scenario+'.nc',mode='r')
			lon = nc_fid.variables['lon'][:]
			lat = nc_fid.variables['lat'][:];#print lat
			summer[:,0:96,:] = nc_fid.variables[winter_id][1,-51:-1,0:96,:]
			summer[:,96:192,:] = nc_fid.variables[summer_id][1,-51:-1,96:192,:]
			winter[:,0:96,:] = nc_fid.variables[summer_id][1,-51:-1,0:96,:]
			winter[:,96:192,:] = nc_fid.variables[winter_id][1,-51:-1,96:192,:]
			# print np.nanmean(np.nanmean(np.nanmean(winter-summer,axis=0),axis=1),axis=0)
			# plt.imshow(a);plt.show(a)
			mask = mask_match(country_key,lon,lat)
			summer_key = np.nansum(np.nansum(np.multiply(stats.nanmean(summer,axis=0),mask),axis=1),axis=0)
			winter_key = np.nansum(np.nansum(np.multiply(stats.nanmean(winter,axis=0),mask),axis=1),axis=0)
			return summer_key,winter_key
		HWDC_S,HWDC_W = summer_winter_hemisphere_swap('DC',country_key)
		HWDu_S,HWDu_W = summer_winter_hemisphere_swap('DuM',country_key)
		HWNO_S,HWNO_W = summer_winter_hemisphere_swap('NO',country_key)
		HWIn_S,HWIn_W = summer_winter_hemisphere_swap('InM',country_key)		

		summer = np.array([HWIn_S,HWDu_S,HWNO_S,HWDC_S]);#print np.shape(summer)
		winter = np.array([HWIn_W,HWDu_W,HWNO_W,HWDC_W]);#print np.shape(winter)
		# print summer-winter
		# result = np.concatenate([summer,winter],axis=0)
		result = np.empty((2,4));result[0,:]=summer;result[1,:]=winter;
		return result
		
	print country_key
	annual_rcp = annual_read('rcp85');season_rcp = summer_winter_read('rcp85');
	rcp_ratio =season_rcp[0,:]-season_rcp[1,:]; print np.round(rcp_ratio,2)
	rcp_ratio =100*np.divide(season_rcp[0,:]-season_rcp[1,:],season_rcp[1,:]); print np.round(rcp_ratio,2)
	
	# annual_fix = annual_read('fixa');season_fix = summer_winter_read('fixa');
	# fix_ratio =season_fix[0,:]-season_fix[1,:]; print np.round(fix_ratio,2)
	# fix_ratio =100*np.divide(season_fix[0,:]-season_fix[1,:],season_fix[1,:]); print np.round(fix_ratio,2)
	
	
summer_winter_weight_print('Globe')
summer_winter_weight_print('Australia')
summer_winter_weight_print('Brazil')
summer_winter_weight_print('China')
summer_winter_weight_print('Europe')
summer_winter_weight_print('India')
summer_winter_weight_print('Southern African')
summer_winter_weight_print('USA')





