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
import matplotlib.ticker as ticker


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

"""
land_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
lon=land_mask['lon'][:];lat=land_mask['lat'][:]
# print lon, lat
# GLO=np.empty((360,720));GLO[:]=1
GLO=land_mask['Globe'][:];Russia=land_mask['Russia'][:];
AUS=land_mask['Australia'][:];EUR=land_mask['Europe'][:];CHA=land_mask['China'][:];
USA=land_mask['USA'][:];IND=land_mask['India'][:];BRZ=land_mask['Brazil'][:]; SAF=land_mask['Southern African'][:]; 

CHA[np.isnan(CHA)]=0;	CHA[CHA>0]=1;
import scipy.io as sio
data_path = './'
sio.savemat(data_path+'Euro_USA_AUS_BRICS_STA_720_360.mat',{'lon':lon,'lat':lat,\
'Globe':GLO,'Europe':EUR,'USA':USA,'India':IND,'Brazil':BRZ,'China':CHA,'Russia':Russia,'Australia':AUS,'Southern African':SAF})   

"""

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
	

def HWI_CSI_boxplot_pp(scenario,country_key,directory):
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp'+directory  #_InterOsi
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_'+scenario+'.nc',mode='r')
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	year = nc_fid.variables['year'][:]
	mask = mask_match(country_key,lon,lat)
	HWDC =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWDC'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWDuM =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWDuX'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWI_NO =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWI_NO'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWInM =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['HWInX'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSDC =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['CSDC'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSDuM =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['CSDuX'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSI_NO =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['CSI_NO'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSInM =  np.nansum(np.nansum(np.multiply(stats.nanmean(nc_fid.variables['CSInX'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	return HWDC,HWDuM,HWI_NO,HWInM,CSDC,CSDuM,CSI_NO,CSInM
	
	
ratio_store=np.empty((2,8,12));	
	

#globes
HWDC_GLO_h,HWDuM_GLO_h,HWI_NO_GLO_h,HWInM_GLO_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Globe','/CalDayThr/')
HWDC_GLO_r,HWDuM_GLO_r,HWI_NO_GLO_r,HWInM_GLO_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Globe','/CalDayThr/')
HWDC_GLO_f,HWDuM_GLO_f,HWI_NO_GLO_f,HWInM_GLO_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Globe','/CalDayThr/')
print "GLO"
ratio_store[0,0,:] = np.array([HWInM_GLO_h[1],HWDuM_GLO_h[1],HWI_NO_GLO_h[1],HWDC_GLO_h[1],\
					HWInM_GLO_r[1],HWDuM_GLO_r[1],HWI_NO_GLO_r[1],HWDC_GLO_r[1],\
					HWInM_GLO_f[1],HWDuM_GLO_f[1],HWI_NO_GLO_f[1],HWDC_GLO_f[1]])
#Australia
HWDC_AUS_h,HWDuM_AUS_h,HWI_NO_AUS_h,HWInM_AUS_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Australia','/CalDayThr/')
HWDC_AUS_r,HWDuM_AUS_r,HWI_NO_AUS_r,HWInM_AUS_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Australia','/CalDayThr/')
HWDC_AUS_f,HWDuM_AUS_f,HWI_NO_AUS_f,HWInM_AUS_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Australia','/CalDayThr/')
print "AUS"
ratio_store[0,1,:] = np.array([HWInM_AUS_h[1],HWDuM_AUS_h[1],HWI_NO_AUS_h[1],HWDC_AUS_h[1],\
					HWInM_AUS_r[1],HWDuM_AUS_r[1],HWI_NO_AUS_r[1],HWDC_AUS_r[1],\
					HWInM_AUS_f[1],HWDuM_AUS_f[1],HWI_NO_AUS_f[1],HWDC_AUS_f[1]])
					
HWDC_BRA_h,HWDuM_BRA_h,HWI_NO_BRA_h,HWInM_BRA_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Brazil','/CalDayThr/')
HWDC_BRA_r,HWDuM_BRA_r,HWI_NO_BRA_r,HWInM_BRA_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Brazil','/CalDayThr/')
HWDC_BRA_f,HWDuM_BRA_f,HWI_NO_BRA_f,HWInM_BRA_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Brazil','/CalDayThr/')
print "BRA"
ratio_store[0,2,:] = np.array([HWInM_BRA_h[1],HWDuM_BRA_h[1],HWI_NO_BRA_h[1],HWDC_BRA_h[1],\
					HWInM_BRA_r[1],HWDuM_BRA_r[1],HWI_NO_BRA_r[1],HWDC_BRA_r[1],\
					HWInM_BRA_f[1],HWDuM_BRA_f[1],HWI_NO_BRA_f[1],HWDC_BRA_f[1]])
					
#China
HWDC_CHA_h,HWDuM_CHA_h,HWI_NO_CHA_h,HWInM_CHA_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','China','/CalDayThr/')
HWDC_CHA_r,HWDuM_CHA_r,HWI_NO_CHA_r,HWInM_CHA_r, _,_,_,_ = HWI_CSI_boxplot_pp('rcp85','China','/CalDayThr/')
HWDC_CHA_f,HWDuM_CHA_f,HWI_NO_CHA_f,HWInM_CHA_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','China','/CalDayThr/')
print "CHA"
ratio_store[0,3,:] = np.array([HWInM_CHA_h[1],HWDuM_CHA_h[1],HWI_NO_CHA_h[1],HWDC_CHA_h[1],\
					HWInM_CHA_r[1],HWDuM_CHA_r[1],HWI_NO_CHA_r[1],HWDC_CHA_r[1],\
					HWInM_CHA_f[1],HWDuM_CHA_f[1],HWI_NO_CHA_f[1],HWDC_CHA_f[1]])

#Europe
HWDC_EUR_h,HWDuM_EUR_h,HWI_NO_EUR_h,HWInM_EUR_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Europe','/CalDayThr/')
HWDC_EUR_r,HWDuM_EUR_r,HWI_NO_EUR_r,HWInM_EUR_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Europe','/CalDayThr/')
HWDC_EUR_f,HWDuM_EUR_f,HWI_NO_EUR_f,HWInM_EUR_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Europe','/CalDayThr/')
print "EUR"
ratio_store[0,4,:] = np.array([HWInM_EUR_h[1],HWDuM_EUR_h[1],HWI_NO_EUR_h[1],HWDC_EUR_h[1],\
					HWInM_EUR_r[1],HWDuM_EUR_r[1],HWI_NO_EUR_r[1],HWDC_EUR_r[1],\
					HWInM_EUR_f[1],HWDuM_EUR_f[1],HWI_NO_EUR_f[1],HWDC_EUR_f[1]])

#India
HWDC_IND_h,HWDuM_IND_h,HWI_NO_IND_h,HWInM_IND_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','India','/CalDayThr/')
HWDC_IND_r,HWDuM_IND_r,HWI_NO_IND_r,HWInM_IND_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','India','/CalDayThr/')
HWDC_IND_f,HWDuM_IND_f,HWI_NO_IND_f,HWInM_IND_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','India','/CalDayThr/')
print "IND"
ratio_store[0,5,:] = np.array([HWInM_IND_h[1],HWDuM_IND_h[1],HWI_NO_IND_h[1],HWDC_IND_h[1],\
					HWInM_IND_r[1],HWDuM_IND_r[1],HWI_NO_IND_r[1],HWDC_IND_r[1],\
					HWInM_IND_f[1],HWDuM_IND_f[1],HWI_NO_IND_f[1],HWDC_IND_f[1]])

HWDC_SAF_h,HWDuM_SAF_h,HWI_NO_SAF_h,HWInM_SAF_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Southern African','/CalDayThr/')
HWDC_SAF_r,HWDuM_SAF_r,HWI_NO_SAF_r,HWInM_SAF_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Southern African','/CalDayThr/')
HWDC_SAF_f,HWDuM_SAF_f,HWI_NO_SAF_f,HWInM_SAF_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Southern African','/CalDayThr/')
print "SAF"

ratio_store[0,6,:] = np.array([HWInM_SAF_h[1],HWDuM_SAF_h[1],HWI_NO_SAF_h[1],HWDC_SAF_h[1],\
					HWInM_SAF_r[1],HWDuM_SAF_r[1],HWI_NO_SAF_r[1],HWDC_SAF_r[1],\
					HWInM_SAF_f[1],HWDuM_SAF_f[1],HWI_NO_SAF_f[1],HWDC_SAF_f[1]])
					
#USA
HWDC_USA_h,HWDuM_USA_h,HWI_NO_USA_h,HWInM_USA_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','USA','/CalDayThr/')
HWDC_USA_r,HWDuM_USA_r,HWI_NO_USA_r,HWInM_USA_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','USA','/CalDayThr/')
HWDC_USA_f,HWDuM_USA_f,HWI_NO_USA_f,HWInM_USA_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','USA','/CalDayThr/')
print "USA"
ratio_store[0,7,:] = np.array([HWInM_USA_h[1],HWDuM_USA_h[1],HWI_NO_USA_h[1],HWDC_USA_h[1],\
					HWInM_USA_r[1],HWDuM_USA_r[1],HWI_NO_USA_r[1],HWDC_USA_r[1],\
					HWInM_USA_f[1],HWDuM_USA_f[1],HWI_NO_USA_f[1],HWDC_USA_f[1]])
					
#globes
HWDC_GLO_h,HWDuM_GLO_h,HWI_NO_GLO_h,HWInM_GLO_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Globe','/CalDayThr_InterOsi/')
HWDC_GLO_r,HWDuM_GLO_r,HWI_NO_GLO_r,HWInM_GLO_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Globe','/CalDayThr_InterOsi/')
HWDC_GLO_f,HWDuM_GLO_f,HWI_NO_GLO_f,HWInM_GLO_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Globe','/CalDayThr_InterOsi/')
print "GLO"
ratio_store[1,0,:] = np.array([HWInM_GLO_h[1],HWDuM_GLO_h[1],HWI_NO_GLO_h[1],HWDC_GLO_h[1],\
					HWInM_GLO_r[1],HWDuM_GLO_r[1],HWI_NO_GLO_r[1],HWDC_GLO_r[1],\
					HWInM_GLO_f[1],HWDuM_GLO_f[1],HWI_NO_GLO_f[1],HWDC_GLO_f[1]])
#Australia
HWDC_AUS_h,HWDuM_AUS_h,HWI_NO_AUS_h,HWInM_AUS_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Australia','/CalDayThr_InterOsi/')
HWDC_AUS_r,HWDuM_AUS_r,HWI_NO_AUS_r,HWInM_AUS_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Australia','/CalDayThr_InterOsi/')
HWDC_AUS_f,HWDuM_AUS_f,HWI_NO_AUS_f,HWInM_AUS_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Australia','/CalDayThr_InterOsi/')
print "AUS"
ratio_store[1,1,:] = np.array([HWInM_AUS_h[1],HWDuM_AUS_h[1],HWI_NO_AUS_h[1],HWDC_AUS_h[1],\
					HWInM_AUS_r[1],HWDuM_AUS_r[1],HWI_NO_AUS_r[1],HWDC_AUS_r[1],\
					HWInM_AUS_f[1],HWDuM_AUS_f[1],HWI_NO_AUS_f[1],HWDC_AUS_f[1]])
					
HWDC_BRA_h,HWDuM_BRA_h,HWI_NO_BRA_h,HWInM_BRA_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Brazil','/CalDayThr_InterOsi/')
HWDC_BRA_r,HWDuM_BRA_r,HWI_NO_BRA_r,HWInM_BRA_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Brazil','/CalDayThr_InterOsi/')
HWDC_BRA_f,HWDuM_BRA_f,HWI_NO_BRA_f,HWInM_BRA_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Brazil','/CalDayThr_InterOsi/')
print "BRA"
ratio_store[1,2,:] = np.array([HWInM_BRA_h[1],HWDuM_BRA_h[1],HWI_NO_BRA_h[1],HWDC_BRA_h[1],\
					HWInM_BRA_r[1],HWDuM_BRA_r[1],HWI_NO_BRA_r[1],HWDC_BRA_r[1],\
					HWInM_BRA_f[1],HWDuM_BRA_f[1],HWI_NO_BRA_f[1],HWDC_BRA_f[1]])
					
#China
HWDC_CHA_h,HWDuM_CHA_h,HWI_NO_CHA_h,HWInM_CHA_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','China','/CalDayThr_InterOsi/')
HWDC_CHA_r,HWDuM_CHA_r,HWI_NO_CHA_r,HWInM_CHA_r, _,_,_,_ = HWI_CSI_boxplot_pp('rcp85','China','/CalDayThr_InterOsi/')
HWDC_CHA_f,HWDuM_CHA_f,HWI_NO_CHA_f,HWInM_CHA_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','China','/CalDayThr_InterOsi/')
print "CHA"
ratio_store[1,3,:] = np.array([HWInM_CHA_h[1],HWDuM_CHA_h[1],HWI_NO_CHA_h[1],HWDC_CHA_h[1],\
					HWInM_CHA_r[1],HWDuM_CHA_r[1],HWI_NO_CHA_r[1],HWDC_CHA_r[1],\
					HWInM_CHA_f[1],HWDuM_CHA_f[1],HWI_NO_CHA_f[1],HWDC_CHA_f[1]])

#Europe
HWDC_EUR_h,HWDuM_EUR_h,HWI_NO_EUR_h,HWInM_EUR_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Europe','/CalDayThr_InterOsi/')
HWDC_EUR_r,HWDuM_EUR_r,HWI_NO_EUR_r,HWInM_EUR_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Europe','/CalDayThr_InterOsi/')
HWDC_EUR_f,HWDuM_EUR_f,HWI_NO_EUR_f,HWInM_EUR_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Europe','/CalDayThr_InterOsi/')
print "EUR"
ratio_store[1,4,:] = np.array([HWInM_EUR_h[1],HWDuM_EUR_h[1],HWI_NO_EUR_h[1],HWDC_EUR_h[1],\
					HWInM_EUR_r[1],HWDuM_EUR_r[1],HWI_NO_EUR_r[1],HWDC_EUR_r[1],\
					HWInM_EUR_f[1],HWDuM_EUR_f[1],HWI_NO_EUR_f[1],HWDC_EUR_f[1]])

#India
HWDC_IND_h,HWDuM_IND_h,HWI_NO_IND_h,HWInM_IND_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','India','/CalDayThr_InterOsi/')
HWDC_IND_r,HWDuM_IND_r,HWI_NO_IND_r,HWInM_IND_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','India','/CalDayThr_InterOsi/')
HWDC_IND_f,HWDuM_IND_f,HWI_NO_IND_f,HWInM_IND_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','India','/CalDayThr_InterOsi/')
print "IND"
ratio_store[1,5,:] = np.array([HWInM_IND_h[1],HWDuM_IND_h[1],HWI_NO_IND_h[1],HWDC_IND_h[1],\
					HWInM_IND_r[1],HWDuM_IND_r[1],HWI_NO_IND_r[1],HWDC_IND_r[1],\
					HWInM_IND_f[1],HWDuM_IND_f[1],HWI_NO_IND_f[1],HWDC_IND_f[1]])

HWDC_SAF_h,HWDuM_SAF_h,HWI_NO_SAF_h,HWInM_SAF_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','Southern African','/CalDayThr_InterOsi/')
HWDC_SAF_r,HWDuM_SAF_r,HWI_NO_SAF_r,HWInM_SAF_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','Southern African','/CalDayThr_InterOsi/')
HWDC_SAF_f,HWDuM_SAF_f,HWI_NO_SAF_f,HWInM_SAF_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','Southern African','/CalDayThr_InterOsi/')
print "SAF"

ratio_store[1,6,:] = np.array([HWInM_SAF_h[1],HWDuM_SAF_h[1],HWI_NO_SAF_h[1],HWDC_SAF_h[1],\
					HWInM_SAF_r[1],HWDuM_SAF_r[1],HWI_NO_SAF_r[1],HWDC_SAF_r[1],\
					HWInM_SAF_f[1],HWDuM_SAF_f[1],HWI_NO_SAF_f[1],HWDC_SAF_f[1]])
					
#USA
HWDC_USA_h,HWDuM_USA_h,HWI_NO_USA_h,HWInM_USA_h,_,_,_,_  = HWI_CSI_boxplot_pp('his','USA','/CalDayThr_InterOsi/')
HWDC_USA_r,HWDuM_USA_r,HWI_NO_USA_r,HWInM_USA_r,_,_,_,_  = HWI_CSI_boxplot_pp('rcp85','USA','/CalDayThr_InterOsi/')
HWDC_USA_f,HWDuM_USA_f,HWI_NO_USA_f,HWInM_USA_f,_,_,_,_  = HWI_CSI_boxplot_pp('fixa','USA','/CalDayThr_InterOsi/')
print "USA"
ratio_store[1,7,:] = np.array([HWInM_USA_h[1],HWDuM_USA_h[1],HWI_NO_USA_h[1],HWDC_USA_h[1],\
					HWInM_USA_r[1],HWDuM_USA_r[1],HWI_NO_USA_r[1],HWDC_USA_r[1],\
					HWInM_USA_f[1],HWDuM_USA_f[1],HWI_NO_USA_f[1],HWDC_USA_f[1]])

def pcolor_fraction(ax,data,x_label=False):
	y, x = np.meshgrid(range(0,5),range(0,9))
	colormap='RdYlBu';colormap = reverse_colourmap(colormap);cmap = discrete_cmap(20,colormap)
	cmap.set_over('k');
	colormesh=ax.pcolor(x, y, data,cmap=cmap, vmin=0, vmax=25)
	# ax.grid(True,linewidth=2,linestyle='-')
	ax.set_xlim([0,8]);ax.set_ylim([0,4]);
	ax.set_xticks(np.arange(0.5,8.5,1));ax.set_yticks(np.arange(0.5,4.5,1))
	ax.set_yticklabels(('Intensity','Duration','frequency','Total days'),rotation=0)
	if x_label:
		ax.set_xticklabels(('GLO','AUS','BRA','CHA','EUR','IND','SAF','USA'),rotation = 90);
	else:
		ax.set_xticklabels([])
	return colormesh
	
a =ratio_store[0,:,4]-ratio_store[0,:,0]
b =ratio_store[1,:,4]-ratio_store[1,:,0]
print b/a*1.0

baseline_a = ratio_store[0,:,0:4]
ratio_store[0,:,4:8] =ratio_store[0,:,4:8]- baseline_a
ratio_store[0,:,8:12] =ratio_store[0,:,8:12]- baseline_a
baseline_v = ratio_store[1,:,0:4]
ratio_store[1,:,4:8] =ratio_store[1,:,4:8]- baseline_v
ratio_store[1,:,8:12] =ratio_store[1,:,8:12]- baseline_v

ratio = np.divide(ratio_store[1,:,:],ratio_store[0,:,:])*100
# ratio_his =ratio[:,0:4]; print ratio_his
ratio_rcp =ratio[:,4:8]; 
ratio_fix =ratio[:,8:12]; 
	
	
fig = plt.figure(facecolor='White',figsize=[5,5]);plot_setup();pad= 5;
# ax = plt.subplot(3,1,1);
# ax.annotate('(a) His (1986-2005)',xy=(0.02,1.02), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)		
# pcolor_fraction(ax,ratio_his,x_label=False)
ax = plt.subplot(2,1,1);
ax.annotate('(a) RCP8.5 (2081-2100)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
pcolor_fraction(ax,ratio_rcp,x_label=False)

ax = plt.subplot(2,1,2);
ax.annotate('(b) RCP8.5_FixA (2081-2100)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

colormesh=pcolor_fraction(ax,ratio_fix,x_label=True)
cbar_ax = fig.add_axes([0.90, 0.16, 0.02, 0.7])
char = fig.colorbar(colormesh,orientation='vertical',extend='max',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.2)*25,2))
char.ax.set_ylabel('%', rotation=90)

plt.subplots_adjust(left=0.15, bottom=0.1, right=0.85, top=0.92, wspace=0.2, hspace=0.2); 
plt.savefig('figureS2_HW_TempVar_fraction.png', format='png', dpi=1000)


