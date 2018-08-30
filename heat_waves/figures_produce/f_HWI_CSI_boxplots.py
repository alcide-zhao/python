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
	
	
def boxplot_plot(his,rcp,fix):
	# ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.set_xlim([-1,20]);ax.set_xticks(np.arange(0.5,20,2.5));ax.set_xticklabels(('GLO','AUS','BRA','CHA','EUR','IND','SAF','USA'));
	medians=his[:,1];maxes=his[:,2];mins=his[:,3];p75=his[:,4];p25=his[:,5];
	ax.errorbar(np.arange(0,20,2.5), medians, [medians - p25, p75 - medians], color='g',fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7)
	ax.errorbar(np.arange(0,20,2.5), medians, [medians - mins, maxes - medians],fmt='+w', capsize=4,ecolor='g', lw=1,markeredgecolor='w',ms=7)
	medians=rcp[:,1];maxes=rcp[:,2];mins=rcp[:,3];p75=rcp[:,4];p25=rcp[:,5];
	ax.errorbar(np.arange(0.5,20,2.5), medians, [medians - p25, p75 - medians], color='r',fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7)
	ax.errorbar(np.arange(0.5,20,2.5), medians, [medians - mins, maxes - medians],fmt='+w', capsize=4,ecolor='r', lw=1,markeredgecolor='w',ms=7)
	medians=fix[:,1];maxes=fix[:,2];mins=fix[:,3];p75=fix[:,4];p25=fix[:,5];
	ax.errorbar(np.arange(1,20,2.5), medians, [medians - p25, p75 - medians], color='b',fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7)
	ax.errorbar(np.arange(1,20,2.5), medians, [medians - mins, maxes - medians],fmt='+w', capsize=4,ecolor='b', lw=1,markeredgecolor='w',ms=7)
	
#globes
HWDC_GLO_h,HWDuM_GLO_h,HWI_NO_GLO_h,HWInM_GLO_h,CSDC_GLO_h,CSDuM_GLO_h,CSI_NO_GLO_h,CSInM_GLO_h = HWI_CSI_boxplot_pp('his','Globe','/CalDayThr/')
HWDC_GLO_r,HWDuM_GLO_r,HWI_NO_GLO_r,HWInM_GLO_r,CSDC_GLO_r,CSDuM_GLO_r,CSI_NO_GLO_r,CSInM_GLO_r = HWI_CSI_boxplot_pp('rcp85','Globe','/CalDayThr/')
HWDC_GLO_f,HWDuM_GLO_f,HWI_NO_GLO_f,HWInM_GLO_f,CSDC_GLO_f,CSDuM_GLO_f,CSI_NO_GLO_f,CSInM_GLO_f = HWI_CSI_boxplot_pp('fixa','Globe','/CalDayThr/')
print "GLO"
print HWInM_GLO_h[1],HWInM_GLO_h[4],HWInM_GLO_h[5]
print HWInM_GLO_r[1],HWInM_GLO_r[4],HWInM_GLO_r[5]
print HWInM_GLO_f[1],HWInM_GLO_f[4],HWInM_GLO_f[5]

#China
HWDC_CHA_h,HWDuM_CHA_h,HWI_NO_CHA_h,HWInM_CHA_h,CSDC_CHA_h,CSDuM_CHA_h,CSI_NO_CHA_h,CSInM_CHA_h = HWI_CSI_boxplot_pp('his','China','/CalDayThr/')
HWDC_CHA_r,HWDuM_CHA_r,HWI_NO_CHA_r,HWInM_CHA_r,CSDC_CHA_r,CSDuM_CHA_r,CSI_NO_CHA_r,CSInM_CHA_r = HWI_CSI_boxplot_pp('rcp85','China','/CalDayThr/')
HWDC_CHA_f,HWDuM_CHA_f,HWI_NO_CHA_f,HWInM_CHA_f,CSDC_CHA_f,CSDuM_CHA_f,CSI_NO_CHA_f,CSInM_CHA_f = HWI_CSI_boxplot_pp('fixa','China','/CalDayThr/')
print "CHA"
print HWInM_CHA_h[1]
print HWInM_CHA_r[1]
print HWInM_CHA_f[1]

#Australia
HWDC_AUS_h,HWDuM_AUS_h,HWI_NO_AUS_h,HWInM_AUS_h,CSDC_AUS_h,CSDuM_AUS_h,CSI_NO_AUS_h,CSInM_AUS_h = HWI_CSI_boxplot_pp('his','Australia','/CalDayThr/')
HWDC_AUS_r,HWDuM_AUS_r,HWI_NO_AUS_r,HWInM_AUS_r,CSDC_AUS_r,CSDuM_AUS_r,CSI_NO_AUS_r,CSInM_AUS_r = HWI_CSI_boxplot_pp('rcp85','Australia','/CalDayThr/')
HWDC_AUS_f,HWDuM_AUS_f,HWI_NO_AUS_f,HWInM_AUS_f,CSDC_AUS_f,CSDuM_AUS_f,CSI_NO_AUS_f,CSInM_AUS_f = HWI_CSI_boxplot_pp('fixa','Australia','/CalDayThr/')
print "AUS"
print HWInM_AUS_h[1],HWInM_AUS_h[4],HWInM_AUS_h[5]
print HWInM_AUS_r[1],HWInM_AUS_r[4],HWInM_AUS_r[5]
print HWInM_AUS_f[1],HWInM_AUS_f[4],HWInM_AUS_f[5]
#Europe
HWDC_EUR_h,HWDuM_EUR_h,HWI_NO_EUR_h,HWInM_EUR_h,CSDC_EUR_h,CSDuM_EUR_h,CSI_NO_EUR_h,CSInM_EUR_h = HWI_CSI_boxplot_pp('his','Europe','/CalDayThr/')
HWDC_EUR_r,HWDuM_EUR_r,HWI_NO_EUR_r,HWInM_EUR_r,CSDC_EUR_r,CSDuM_EUR_r,CSI_NO_EUR_r,CSInM_EUR_r = HWI_CSI_boxplot_pp('rcp85','Europe','/CalDayThr/')
HWDC_EUR_f,HWDuM_EUR_f,HWI_NO_EUR_f,HWInM_EUR_f,CSDC_EUR_f,CSDuM_EUR_f,CSI_NO_EUR_f,CSInM_EUR_f = HWI_CSI_boxplot_pp('fixa','Europe','/CalDayThr/')
print "EUR"
print HWInM_EUR_h[1],HWInM_EUR_h[4],HWInM_EUR_h[5]
print HWInM_EUR_r[1],HWInM_EUR_r[4],HWInM_EUR_r[5]
print HWInM_EUR_f[1],HWInM_EUR_f[4],HWDuM_EUR_f[5]

#CHA
HWDC_USA_h,HWDuM_USA_h,HWI_NO_USA_h,HWInM_USA_h,CSDC_USA_h,CSDuM_USA_h,CSI_NO_USA_h,CSInM_USA_h = HWI_CSI_boxplot_pp('his','USA','/CalDayThr/')
HWDC_USA_r,HWDuM_USA_r,HWI_NO_USA_r,HWInM_USA_r,CSDC_USA_r,CSDuM_USA_r,CSI_NO_USA_r,CSInM_USA_r = HWI_CSI_boxplot_pp('rcp85','USA','/CalDayThr/')
HWDC_USA_f,HWDuM_USA_f,HWI_NO_USA_f,HWInM_USA_f,CSDC_USA_f,CSDuM_USA_f,CSI_NO_USA_f,CSInM_USA_f = HWI_CSI_boxplot_pp('fixa','USA','/CalDayThr/')
print "USA"
print HWInM_USA_h[1]
print HWInM_USA_r[1]
print HWInM_USA_f[1]
#India
HWDC_IND_h,HWDuM_IND_h,HWI_NO_IND_h,HWInM_IND_h,CSDC_IND_h,CSDuM_IND_h,CSI_NO_IND_h,CSInM_IND_h = HWI_CSI_boxplot_pp('his','India','/CalDayThr/')
HWDC_IND_r,HWDuM_IND_r,HWI_NO_IND_r,HWInM_IND_r,CSDC_IND_r,CSDuM_IND_r,CSI_NO_IND_r,CSInM_IND_r = HWI_CSI_boxplot_pp('rcp85','India','/CalDayThr/')
HWDC_IND_f,HWDuM_IND_f,HWI_NO_IND_f,HWInM_IND_f,CSDC_IND_f,CSDuM_IND_f,CSI_NO_IND_f,CSInM_IND_f = HWI_CSI_boxplot_pp('fixa','India','/CalDayThr/')
print "IND"
print HWInM_IND_h[1]
print HWInM_IND_r[1]
print HWInM_IND_f[1]
#India
HWDC_BRA_h,HWDuM_BRA_h,HWI_NO_BRA_h,HWInM_BRA_h,CSDC_BRA_h,CSDuM_BRA_h,CSI_NO_BRA_h,CSInM_BRA_h = HWI_CSI_boxplot_pp('his','Brazil','/CalDayThr/')
HWDC_BRA_r,HWDuM_BRA_r,HWI_NO_BRA_r,HWInM_BRA_r,CSDC_BRA_r,CSDuM_BRA_r,CSI_NO_BRA_r,CSInM_BRA_r = HWI_CSI_boxplot_pp('rcp85','Brazil','/CalDayThr/')
HWDC_BRA_f,HWDuM_BRA_f,HWI_NO_BRA_f,HWInM_BRA_f,CSDC_BRA_f,CSDuM_BRA_f,CSI_NO_BRA_f,CSInM_BRA_f = HWI_CSI_boxplot_pp('fixa','Brazil','/CalDayThr/')
print "BRA"
print HWInM_BRA_h[1],HWInM_BRA_h[4],HWInM_EUR_h[5]
print HWInM_BRA_r[1],HWInM_BRA_r[4],HWInM_BRA_r[5]
print HWInM_BRA_f[1],HWInM_BRA_f[4],HWInM_BRA_f[5]

HWDC_SAF_h,HWDuM_SAF_h,HWI_NO_SAF_h,HWInM_SAF_h,CSDC_SAF_h,CSDuM_SAF_h,CSI_NO_SAF_h,CSInM_SAF_h = HWI_CSI_boxplot_pp('his','Southern African','/CalDayThr/')
HWDC_SAF_r,HWDuM_SAF_r,HWI_NO_SAF_r,HWInM_SAF_r,CSDC_SAF_r,CSDuM_SAF_r,CSI_NO_SAF_r,CSInM_SAF_r = HWI_CSI_boxplot_pp('rcp85','Southern African','/CalDayThr/')
HWDC_SAF_f,HWDuM_SAF_f,HWI_NO_SAF_f,HWInM_SAF_f,CSDC_SAF_f,CSDuM_SAF_f,CSI_NO_SAF_f,CSInM_SAF_f = HWI_CSI_boxplot_pp('fixa','Southern African','/CalDayThr/')
print "SAF"
print HWInM_SAF_h[1]
print HWInM_SAF_r[1]
print HWInM_SAF_f[1]


fig = plt.figure(facecolor='White',figsize=[10,13]);plot_setup();pad= 5;
ax = plt.subplot(4,2,1);
ax.annotate('Annual peak intensity\n($^\circ$C)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
# ax.annotate('Heat Waves',xy=(0.5, 1.05), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',rotation='horizontal')
ax.annotate('(a)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
			
his = np.array([HWInM_GLO_h,HWInM_AUS_h,HWInM_BRA_h,HWInM_CHA_h,HWInM_EUR_h,HWInM_IND_h,HWInM_SAF_h,HWInM_USA_h])
rcp = np.array([HWInM_GLO_r,HWInM_AUS_r,HWInM_BRA_r,HWInM_CHA_r,HWInM_EUR_r,HWInM_IND_r,HWInM_SAF_r,HWInM_USA_r])
fix = np.array([HWInM_GLO_f,HWInM_AUS_f,HWInM_BRA_f,HWInM_CHA_f,HWInM_EUR_f,HWInM_IND_f,HWInM_SAF_f,HWInM_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,8]);ax.set_yticks(np.arange(0,8.1,2));

ax = plt.subplot(4,2,3);
ax.annotate('Anuual maximum duration\n'+r'($\mathrm{day\/event^{-1}}$)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(b)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
his = np.array([HWDuM_GLO_h,HWDuM_AUS_h,HWDuM_BRA_h,HWDuM_CHA_h,HWDuM_EUR_h,HWDuM_IND_h,HWDuM_SAF_h,HWDuM_USA_h])
rcp = np.array([HWDuM_GLO_r,HWDuM_AUS_r,HWDuM_BRA_r,HWDuM_CHA_r,HWDuM_EUR_r,HWDuM_IND_r,HWDuM_SAF_r,HWDuM_USA_r])
fix = np.array([HWDuM_GLO_f,HWDuM_AUS_f,HWDuM_BRA_f,HWDuM_CHA_f,HWDuM_EUR_f,HWDuM_IND_f,HWDuM_SAF_f,HWDuM_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,120]);ax.set_yticks(np.arange(0,120.1,30));

	
ax = plt.subplot(4,2,5);
ax.annotate('Anuual frequency\n'+r'($\mathrm{events\/year^{-1}}$)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
his = np.array([HWI_NO_GLO_h,HWI_NO_AUS_h,HWI_NO_BRA_h,HWI_NO_CHA_h,HWI_NO_EUR_h,HWI_NO_IND_h,HWI_NO_SAF_h,HWI_NO_USA_h])
rcp = np.array([HWI_NO_GLO_r,HWI_NO_AUS_r,HWI_NO_BRA_r,HWI_NO_CHA_r,HWI_NO_EUR_r,HWI_NO_IND_r,HWI_NO_SAF_r,HWI_NO_USA_r])
fix = np.array([HWI_NO_GLO_f,HWI_NO_AUS_f,HWI_NO_BRA_f,HWI_NO_CHA_f,HWI_NO_EUR_f,HWI_NO_IND_f,HWI_NO_SAF_f,HWI_NO_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,30]);ax.set_yticks(np.arange(0,30.1,10));

		
ax = plt.subplot(4,2,7);
ax.annotate('Total hot days\n'+r'($\mathrm{days\/year^{-1}}$)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(d)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
				
his = np.array([HWDC_GLO_h,HWDC_AUS_h,HWDC_BRA_h,HWDC_CHA_h,HWDC_EUR_h,HWDC_IND_h,HWDC_SAF_h,HWDC_USA_h])
rcp = np.array([HWDC_GLO_r,HWDC_AUS_r,HWDC_BRA_r,HWDC_CHA_r,HWDC_EUR_r,HWDC_IND_r,HWDC_SAF_r,HWDC_USA_r])
fix = np.array([HWDC_GLO_f,HWDC_AUS_f,HWDC_BRA_f,HWDC_CHA_f,HWDC_EUR_f,HWDC_IND_f,HWDC_SAF_f,HWDC_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,360]);ax.set_yticks(np.arange(0,360.1,90));



#China
HWDC_CHA_h,HWDuM_CHA_h,HWI_NO_CHA_h,HWInM_CHA_h,CSDC_CHA_h,CSDuM_CHA_h,CSI_NO_CHA_h,CSInM_CHA_h = HWI_CSI_boxplot_pp('his','China','/CalDayThr_InterOsi/')
HWDC_CHA_r,HWDuM_CHA_r,HWI_NO_CHA_r,HWInM_CHA_r,CSDC_CHA_r,CSDuM_CHA_r,CSI_NO_CHA_r,CSInM_CHA_r = HWI_CSI_boxplot_pp('rcp85','China','/CalDayThr_InterOsi/')
HWDC_CHA_f,HWDuM_CHA_f,HWI_NO_CHA_f,HWInM_CHA_f,CSDC_CHA_f,CSDuM_CHA_f,CSI_NO_CHA_f,CSInM_CHA_f = HWI_CSI_boxplot_pp('fixa','China','/CalDayThr_InterOsi/')
print "CHA"
print HWDuM_CHA_h[1]
print HWDuM_CHA_r[1]
print HWDuM_CHA_f[1]

#globes
HWDC_GLO_h,HWDuM_GLO_h,HWI_NO_GLO_h,HWInM_GLO_h,CSDC_GLO_h,CSDuM_GLO_h,CSI_NO_GLO_h,CSInM_GLO_h = HWI_CSI_boxplot_pp('his','Globe','/CalDayThr_InterOsi/')
HWDC_GLO_r,HWDuM_GLO_r,HWI_NO_GLO_r,HWInM_GLO_r,CSDC_GLO_r,CSDuM_GLO_r,CSI_NO_GLO_r,CSInM_GLO_r = HWI_CSI_boxplot_pp('rcp85','Globe','/CalDayThr_InterOsi/')
HWDC_GLO_f,HWDuM_GLO_f,HWI_NO_GLO_f,HWInM_GLO_f,CSDC_GLO_f,CSDuM_GLO_f,CSI_NO_GLO_f,CSInM_GLO_f = HWI_CSI_boxplot_pp('fixa','Globe','/CalDayThr_InterOsi/')
print "GLO"
print HWDuM_GLO_h[1]
print HWDuM_GLO_r[1]
print HWDuM_GLO_f[1]
#Australia
HWDC_AUS_h,HWDuM_AUS_h,HWI_NO_AUS_h,HWInM_AUS_h,CSDC_AUS_h,CSDuM_AUS_h,CSI_NO_AUS_h,CSInM_AUS_h = HWI_CSI_boxplot_pp('his','Australia','/CalDayThr_InterOsi/')
HWDC_AUS_r,HWDuM_AUS_r,HWI_NO_AUS_r,HWInM_AUS_r,CSDC_AUS_r,CSDuM_AUS_r,CSI_NO_AUS_r,CSInM_AUS_r = HWI_CSI_boxplot_pp('rcp85','Australia','/CalDayThr_InterOsi/')
HWDC_AUS_f,HWDuM_AUS_f,HWI_NO_AUS_f,HWInM_AUS_f,CSDC_AUS_f,CSDuM_AUS_f,CSI_NO_AUS_f,CSInM_AUS_f = HWI_CSI_boxplot_pp('fixa','Australia','/CalDayThr_InterOsi/')
print "AUS"
print HWDuM_AUS_h[1]
print HWDuM_AUS_r[1]
print HWDuM_AUS_f[1]
#Europe
HWDC_EUR_h,HWDuM_EUR_h,HWI_NO_EUR_h,HWInM_EUR_h,CSDC_EUR_h,CSDuM_EUR_h,CSI_NO_EUR_h,CSInM_EUR_h = HWI_CSI_boxplot_pp('his','Europe','/CalDayThr_InterOsi/')
HWDC_EUR_r,HWDuM_EUR_r,HWI_NO_EUR_r,HWInM_EUR_r,CSDC_EUR_r,CSDuM_EUR_r,CSI_NO_EUR_r,CSInM_EUR_r = HWI_CSI_boxplot_pp('rcp85','Europe','/CalDayThr_InterOsi/')
HWDC_EUR_f,HWDuM_EUR_f,HWI_NO_EUR_f,HWInM_EUR_f,CSDC_EUR_f,CSDuM_EUR_f,CSI_NO_EUR_f,CSInM_EUR_f = HWI_CSI_boxplot_pp('fixa','Europe','/CalDayThr_InterOsi/')
print "EUR"
print HWDuM_EUR_h[1]
print HWDuM_EUR_r[1]
print HWDuM_EUR_f[1]

#CHA
HWDC_USA_h,HWDuM_USA_h,HWI_NO_USA_h,HWInM_USA_h,CSDC_USA_h,CSDuM_USA_h,CSI_NO_USA_h,CSInM_USA_h = HWI_CSI_boxplot_pp('his','USA','/CalDayThr_InterOsi/')
HWDC_USA_r,HWDuM_USA_r,HWI_NO_USA_r,HWInM_USA_r,CSDC_USA_r,CSDuM_USA_r,CSI_NO_USA_r,CSInM_USA_r = HWI_CSI_boxplot_pp('rcp85','USA','/CalDayThr_InterOsi/')
HWDC_USA_f,HWDuM_USA_f,HWI_NO_USA_f,HWInM_USA_f,CSDC_USA_f,CSDuM_USA_f,CSI_NO_USA_f,CSInM_USA_f = HWI_CSI_boxplot_pp('fixa','USA','/CalDayThr_InterOsi/')
print "USA"
print HWDuM_USA_h[1]
print HWDuM_USA_r[1]
print HWDuM_USA_f[1]
#India
HWDC_IND_h,HWDuM_IND_h,HWI_NO_IND_h,HWInM_IND_h,CSDC_IND_h,CSDuM_IND_h,CSI_NO_IND_h,CSInM_IND_h = HWI_CSI_boxplot_pp('his','India','/CalDayThr_InterOsi/')
HWDC_IND_r,HWDuM_IND_r,HWI_NO_IND_r,HWInM_IND_r,CSDC_IND_r,CSDuM_IND_r,CSI_NO_IND_r,CSInM_IND_r = HWI_CSI_boxplot_pp('rcp85','India','/CalDayThr_InterOsi/')
HWDC_IND_f,HWDuM_IND_f,HWI_NO_IND_f,HWInM_IND_f,CSDC_IND_f,CSDuM_IND_f,CSI_NO_IND_f,CSInM_IND_f = HWI_CSI_boxplot_pp('fixa','India','/CalDayThr_InterOsi/')
print "IND"
print HWDuM_IND_h[1]
print HWDuM_IND_r[1]
print HWDuM_IND_f[1]
#India
HWDC_BRA_h,HWDuM_BRA_h,HWI_NO_BRA_h,HWInM_BRA_h,CSDC_BRA_h,CSDuM_BRA_h,CSI_NO_BRA_h,CSInM_BRA_h = HWI_CSI_boxplot_pp('his','Brazil','/CalDayThr_InterOsi/')
HWDC_BRA_r,HWDuM_BRA_r,HWI_NO_BRA_r,HWInM_BRA_r,CSDC_BRA_r,CSDuM_BRA_r,CSI_NO_BRA_r,CSInM_BRA_r = HWI_CSI_boxplot_pp('rcp85','Brazil','/CalDayThr_InterOsi/')
HWDC_BRA_f,HWDuM_BRA_f,HWI_NO_BRA_f,HWInM_BRA_f,CSDC_BRA_f,CSDuM_BRA_f,CSI_NO_BRA_f,CSInM_BRA_f = HWI_CSI_boxplot_pp('fixa','Brazil','/CalDayThr_InterOsi/')
print "BRA"
print HWDuM_BRA_h[1]
print HWDuM_BRA_r[1]
print HWDuM_BRA_f[1]

HWDC_SAF_h,HWDuM_SAF_h,HWI_NO_SAF_h,HWInM_SAF_h,CSDC_SAF_h,CSDuM_SAF_h,CSI_NO_SAF_h,CSInM_SAF_h = HWI_CSI_boxplot_pp('his','Southern African','/CalDayThr_InterOsi/')
HWDC_SAF_r,HWDuM_SAF_r,HWI_NO_SAF_r,HWInM_SAF_r,CSDC_SAF_r,CSDuM_SAF_r,CSI_NO_SAF_r,CSInM_SAF_r = HWI_CSI_boxplot_pp('rcp85','Southern African','/CalDayThr_InterOsi/')
HWDC_SAF_f,HWDuM_SAF_f,HWI_NO_SAF_f,HWInM_SAF_f,CSDC_SAF_f,CSDuM_SAF_f,CSI_NO_SAF_f,CSInM_SAF_f = HWI_CSI_boxplot_pp('fixa','Southern African','/CalDayThr_InterOsi/')
print "SAF"
print HWDuM_SAF_h[1]
print HWDuM_SAF_r[1]
print HWDuM_SAF_f[1]


ax = plt.subplot(4,2,2);
ax.annotate('(e)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
			
his = np.array([HWInM_GLO_h,HWInM_AUS_h,HWInM_BRA_h,HWInM_CHA_h,HWInM_EUR_h,HWInM_IND_h,HWInM_SAF_h,HWInM_USA_h])
rcp = np.array([HWInM_GLO_r,HWInM_AUS_r,HWInM_BRA_r,HWInM_CHA_r,HWInM_EUR_r,HWInM_IND_r,HWInM_SAF_r,HWInM_USA_r])
fix = np.array([HWInM_GLO_f,HWInM_AUS_f,HWInM_BRA_f,HWInM_CHA_f,HWInM_EUR_f,HWInM_IND_f,HWInM_SAF_f,HWInM_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,8]);ax.set_yticks(np.arange(0,8.1,2));

ax = plt.subplot(4,2,4);
ax.annotate('(f)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
his = np.array([HWDuM_GLO_h,HWDuM_AUS_h,HWDuM_BRA_h,HWDuM_CHA_h,HWDuM_EUR_h,HWDuM_IND_h,HWDuM_SAF_h,HWDuM_USA_h])
rcp = np.array([HWDuM_GLO_r,HWDuM_AUS_r,HWDuM_BRA_r,HWDuM_CHA_r,HWDuM_EUR_r,HWDuM_IND_r,HWDuM_SAF_r,HWDuM_USA_r])
fix = np.array([HWDuM_GLO_f,HWDuM_AUS_f,HWDuM_BRA_f,HWDuM_CHA_f,HWDuM_EUR_f,HWDuM_IND_f,HWDuM_SAF_f,HWDuM_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,120]);ax.set_yticks(np.arange(0,120.1,30));
	
ax = plt.subplot(4,2,6);
ax.annotate('(g)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
his = np.array([HWI_NO_GLO_h,HWI_NO_AUS_h,HWI_NO_BRA_h,HWI_NO_CHA_h,HWI_NO_EUR_h,HWI_NO_IND_h,HWI_NO_SAF_h,HWI_NO_USA_h])
rcp = np.array([HWI_NO_GLO_r,HWI_NO_AUS_r,HWI_NO_BRA_r,HWI_NO_CHA_r,HWI_NO_EUR_r,HWI_NO_IND_r,HWI_NO_SAF_r,HWI_NO_USA_r])
fix = np.array([HWI_NO_GLO_f,HWI_NO_AUS_f,HWI_NO_BRA_f,HWI_NO_CHA_f,HWI_NO_EUR_f,HWI_NO_IND_f,HWI_NO_SAF_f,HWI_NO_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,30]);ax.set_yticks(np.arange(0,30.1,10));

			
ax = plt.subplot(4,2,8);
ax.annotate('h)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
				
his = np.array([HWDC_GLO_h,HWDC_AUS_h,HWDC_BRA_h,HWDC_CHA_h,HWDC_EUR_h,HWDC_IND_h,HWDC_SAF_h,HWDC_USA_h])
rcp = np.array([HWDC_GLO_r,HWDC_AUS_r,HWDC_BRA_r,HWDC_CHA_r,HWDC_EUR_r,HWDC_IND_r,HWDC_SAF_r,HWDC_USA_r])
fix = np.array([HWDC_GLO_f,HWDC_AUS_f,HWDC_BRA_f,HWDC_CHA_f,HWDC_EUR_f,HWDC_IND_f,HWDC_SAF_f,HWDC_USA_f])
boxplot_plot(his,rcp,fix)
ax.set_ylim([0,360]);ax.set_yticks(np.arange(0,360.1,90));


plt.subplots_adjust(left=0.1, bottom=0.05, right=0.98, top=0.95, wspace=0.2, hspace=0.2); 
plt.savefig('HW_boxplot.png', format='png', dpi=1000)


"""
ax = plt.subplot(4,2,2);
ax.annotate('Cold Spells',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(b)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
			
his = np.array([CSInM_GLO_h,CSInM_AUS_h,CSInM_BRA_h,CSInM_CHA_h,CSInM_EUR_h,CSInM_IND_h,CSInM_SAF_h,CSInM_USA_h])
rcp = np.array([CSInM_GLO_r,CSInM_AUS_r,CSInM_BRA_r,CSInM_CHA_r,CSInM_EUR_r,CSInM_IND_r,CSInM_SAF_r,CSInM_USA_r])
fix = np.array([CSInM_GLO_f,CSInM_AUS_f,CSInM_BRA_f,CSInM_CHA_f,CSInM_EUR_f,CSInM_IND_f,CSInM_SAF_f,CSInM_USA_f])
boxplot_plot(his,rcp,fix)


ax = plt.subplot(4,2,4);
ax.annotate('(d)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
his = np.array([CSDuM_GLO_h,CSDuM_AUS_h,CSDuM_BRA_h,CSDuM_CHA_h,CSDuM_EUR_h,CSDuM_IND_h,CSDuM_SAF_h,CSDuM_USA_h])
rcp = np.array([CSDuM_GLO_r,CSDuM_AUS_r,CSDuM_BRA_r,CSDuM_CHA_r,CSDuM_EUR_r,CSDuM_IND_r,CSDuM_SAF_r,CSDuM_USA_r])
fix = np.array([CSDuM_GLO_f,CSDuM_AUS_f,CSDuM_BRA_f,CSDuM_CHA_f,CSDuM_EUR_f,CSDuM_IND_f,CSDuM_SAF_f,CSDuM_USA_f])

boxplot_plot(his,rcp,fix)
	
ax = plt.subplot(4,2,6);
ax.annotate('(f)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
his = np.array([CSI_NO_GLO_h,CSI_NO_AUS_h,CSI_NO_BRA_h,CSI_NO_CHA_h,CSI_NO_EUR_h,CSI_NO_IND_h,CSI_NO_SAF_h,CSI_NO_USA_h])
rcp = np.array([CSI_NO_GLO_r,CSI_NO_AUS_r,CSI_NO_BRA_r,CSI_NO_CHA_r,CSI_NO_EUR_r,CSI_NO_IND_r,CSI_NO_SAF_r,CSI_NO_USA_r])
fix = np.array([CSI_NO_GLO_f,CSI_NO_AUS_f,CSI_NO_BRA_f,CSI_NO_CHA_f,CSI_NO_EUR_f,CSI_NO_IND_f,CSI_NO_SAF_f,CSI_NO_USA_f])

boxplot_plot(his,rcp,fix)
				
ax = plt.subplot(4,2,8);

ax.annotate('(h)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
				
his = np.array([CSDC_GLO_h,CSDC_AUS_h,CSDC_BRA_h,CSDC_CHA_h,CSDC_EUR_h,CSDC_IND_h,CSDC_SAF_h,CSDC_USA_h])
rcp = np.array([CSDC_GLO_r,CSDC_AUS_r,CSDC_BRA_r,CSDC_CHA_r,CSDC_EUR_r,CSDC_IND_r,CSDC_SAF_r,CSDC_USA_r])
fix = np.array([CSDC_GLO_f,CSDC_AUS_f,CSDC_BRA_f,CSDC_CHA_f,CSDC_EUR_f,CSDC_IND_f,CSDC_SAF_f,CSDC_USA_f])
boxplot_plot(his,rcp,fix)
"""
