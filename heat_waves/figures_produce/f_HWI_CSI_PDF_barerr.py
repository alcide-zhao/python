# -*- coding: utf-8 -*-
'''
This is to plot the PDF shifts of HWI and CSI magnitudes over 7 different regions
using bar chart with error bars (25-75 ensemble percentile) on them
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as sio
import math
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


ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_720_360.mat')
lon_mask = ocean_mask['lon'][0,:];


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
	

def HWI_CSI_boxplot_pp(scenario,country_key):
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_'+scenario+'.nc',mode='r')
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	year = nc_fid.variables['year'][:]
	mask = mask_match(country_key,lon,lat)
	"""
	Take the sum of no of HW events over 20 years over each grid-box and then spatially averaged into the interested region
	return the means, 75th and 25th percentiles
	"""

	HWCM =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCM'][:,-21:-1,:,:]+nc_fid.variables['HWCN'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWCD =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCD'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWCS =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCS'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWCE =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCE'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWCV =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCV'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWCP =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCP'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	HWCU =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['HWCU'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSCM =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['CSCM'][:,-21:-1,:,:]+nc_fid.variables['CSCN'][:,-21:-1,:,:]+nc_fid.variables['CSCD'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSCS =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['CSCS'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSCE =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['CSCE'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSCV =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['CSCV'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSCP =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['CSCP'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)
	CSCU =  np.nansum(np.nansum(np.multiply(np.nansum(nc_fid.variables['CSCU'][:,-21:-1,:,:],axis=1),mask),axis=2),axis=1)	
	cache = np.array([CSCU,CSCP,CSCV,CSCE,CSCS,CSCM,HWCM,HWCD,HWCS,HWCE,HWCV,HWCP,HWCU])
	# result = np.array([cache[:,0],cache[:,4],cache[:,5]])
	# print cache[:,5]
	return cache
	
	
def bar_err_plot(his,rcp,fix):
	# ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
	ax.yaxis.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
	ax.axvline(x=0,color='k',lw=2)
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.tick_params(axis='x', which='major', pad=0.5);
	# ax.set_xlim([-6,7]);ax.set_xticks(np.arange(-5.5,6.6,1));
	ax.set_xlim([0,7]);ax.set_xticks(np.arange(0.5,7.1,1));
	# ax.set_ylim([0,90]);ax.set_yticks(np.arange(0,91,30));
	ax.set_ylim([0,150]);ax.set_yticks(np.arange(0,151,25));
	ax.set_xticklabels(('(-inf,0]','(0,2]','(2,4]','(4,8]','(8,16]','(16,32]','(32,inf)'),rotation = 90);
	# ax.set_xticklabels(('[16,inf)','[8,16]','[4,8]','[2,4]','[0,2]','(-inf,0]','(-inf,0]','[0,2]','[2,4]','[4,8]','[8,16]','[16,32]','[32,inf)'),rotation = 90);
	height = rcp[:,0];yerr=np.array([rcp[:,1]-rcp[:,5],rcp[:,4]- rcp[:,1]]); 
	# rects21 = ax.bar(np.arange(-6,0), height[0:6], width=1,yerr=yerr[:,0:6],color='r',ecolor='r',edgecolor='none',lw=0,capsize=3,alpha=0.6)
	rects22 = ax.bar(np.arange(0,7), height[6:13], width=1,yerr=yerr[:,6:13],color='r',ecolor='r',edgecolor='none',lw=0,capsize=3,alpha=0.6)
	height = fix[:,0];yerr=np.array([fix[:,1]-fix[:,5],fix[:,4]- fix[:,1]]); 
	# rects31 = ax.bar(np.arange(-6,0), height[0:6], width=1,yerr=yerr[:,0:6],color='b',ecolor='b',edgecolor='none',lw=0,capsize=4,alpha=0.4)
	rects32 = ax.bar(np.arange(0,7), height[6:13], width=1,yerr=yerr[:,6:13],color='b',ecolor='b',edgecolor='none',lw=0,capsize=4,alpha=0.4)
	height = his[:,0];yerr=np.array([his[:,1]-his[:,5],his[:,4]- his[:,1]]); 
	# rects11 = ax.bar(np.arange(-6,0), height[0:6], width=1,yerr=yerr[:,0:6],color='g',ecolor='g',edgecolor='none',lw=0,capsize=2,alpha=0.6)
	rects12 = ax.bar(np.arange(0,7), height[6:13], width=1,yerr=yerr[:,6:13],color='g',ecolor='g',edgecolor='none',lw=0,capsize=2,alpha=0.6)
	print height,yerr
	
	return rects12,rects22,rects32
	
#globes

MagHis_GLO_h = HWI_CSI_boxplot_pp('his','Globe');MagHis_GLO_r = HWI_CSI_boxplot_pp('rcp85','Globe');MagHis_GLO_f = HWI_CSI_boxplot_pp('fixa','Globe');
#Australia
MagHis_AUS_h = HWI_CSI_boxplot_pp('his','Australia');MagHis_AUS_r = HWI_CSI_boxplot_pp('rcp85','Australia');MagHis_AUS_f = HWI_CSI_boxplot_pp('fixa','Australia');
#Europe
MagHis_EUR_h = HWI_CSI_boxplot_pp('his','Europe');MagHis_EUR_r = HWI_CSI_boxplot_pp('rcp85','Europe');MagHis_EUR_f = HWI_CSI_boxplot_pp('fixa','Europe');
#China
MagHis_CHA_h = HWI_CSI_boxplot_pp('his','China');MagHis_CHA_r = HWI_CSI_boxplot_pp('rcp85','China');MagHis_CHA_f = HWI_CSI_boxplot_pp('fixa','China');
#USA
MagHis_USA_h = HWI_CSI_boxplot_pp('his','USA');MagHis_USA_r = HWI_CSI_boxplot_pp('rcp85','USA');MagHis_USA_f = HWI_CSI_boxplot_pp('fixa','USA');
#India
MagHis_IND_h = HWI_CSI_boxplot_pp('his','India');MagHis_IND_r = HWI_CSI_boxplot_pp('rcp85','India');MagHis_IND_f = HWI_CSI_boxplot_pp('fixa','India');
#Brazil
MagHis_BRA_h = HWI_CSI_boxplot_pp('his','Brazil');MagHis_BRA_r = HWI_CSI_boxplot_pp('rcp85','Brazil');MagHis_BRA_f = HWI_CSI_boxplot_pp('fixa','Brazil');
#Southern African
MagHis_SAF_h = HWI_CSI_boxplot_pp('his','Southern African');MagHis_SAF_r = HWI_CSI_boxplot_pp('rcp85','Southern African');MagHis_SAF_f = HWI_CSI_boxplot_pp('fixa','Southern African');


fig = plt.figure(facecolor='White',figsize=[8,13]);plot_setup();pad= 5;
ax = plt.subplot(4,2,1);
ax.annotate('(a) GLO',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
rects11,rects21,rects31=bar_err_plot(MagHis_GLO_h,MagHis_GLO_r,MagHis_GLO_f)
legend=ax.legend((rects11[0], rects21[0],rects31[0]),('His (1986-2005)','RCP8.5 (2081-2100)','RCP8.5_FixA (2081-2100)'),ncol=1,labelspacing=1,markerscale =15,bbox_to_anchor=(1.0, 1.0))
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)

print "------------------------"
ax = plt.subplot(4,2,2);
ax.annotate('(b) AUS',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
bar_err_plot(MagHis_AUS_h,MagHis_AUS_r,MagHis_AUS_f)

ax = plt.subplot(4,2,3);
ax.annotate('(c) BRA',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
bar_err_plot(MagHis_BRA_h,MagHis_BRA_r,MagHis_BRA_f)

ax = plt.subplot(4,2,4);
ax.annotate('(d) CHA',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
bar_err_plot(MagHis_CHA_h,MagHis_CHA_r,MagHis_CHA_f)

ax = plt.subplot(4,2,5);
ax.annotate('(e) EUR',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
bar_err_plot(MagHis_EUR_h,MagHis_EUR_r,MagHis_EUR_f)
				
ax = plt.subplot(4,2,6);
ax.annotate('(f) IND',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
bar_err_plot(MagHis_IND_h,MagHis_IND_r,MagHis_IND_f)

ax = plt.subplot(4,2,7);
ax.annotate('(g) SAF',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
bar_err_plot(MagHis_SAF_h,MagHis_SAF_r,MagHis_SAF_f)

				
ax = plt.subplot(4,2,8);
ax.annotate('(h) USA',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
bar_err_plot(MagHis_USA_h,MagHis_USA_r,MagHis_USA_f)		

plt.subplots_adjust(left=0.1, bottom=0.03, right=0.98, top=0.98, wspace=0.2, hspace=0.3); 
plt.savefig('Fig4_HW_PDF_Bar_err.png', format='png', dpi=1000)  #_CalDayThr_InterOsi

