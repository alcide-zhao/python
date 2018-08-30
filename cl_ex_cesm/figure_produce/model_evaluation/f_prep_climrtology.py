# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This is to plot the climatology mean of precipitation and precipitation extremes
The observational precipitation data are from the APHRO spans from 1951 to 2007 with a resolution of 
The reference monthly RX5day precipitation indices are derived from the AROPH from 1901 to 2010
The model data are from the CESM Ensemble  means
@author: Alcide.Zhao
"""
import scipy

import glob
import netCDF4 as nc4
import numpy as np
import scipy.io as spio
from scipy import arange
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d
import os
from scipy import stats
import math
import site
# from scipy.interpolate import interp2d

lib_path = os.path.join(
	os.path.realpath(
		os.path.dirname(__file__)
	), 
	os.path.pardir, 
	os.path.pardir,
    # 'lib'
)

site.addsitedir(lib_path)

from lib import *

########################################
# 0. Vaeiables and functions
########################################
#leapyear judgement
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

rergion_dic={'ASIA':[60,150,-5,55],'EA':[110,135,22.5,45],'SA':[70,105,10,30]}
region = rergion_dic['ASIA'][:]

# ocean land mask
oceanmask=spio.loadmat('/home/s1667168/coding/python/cl_ex_cesm/external_data/world.oceanmask.1440x720')['landocean']
ocean_mask = np.flipud(oceanmask)
cache = np.empty((720,1440))
cache[:,0:720] = ocean_mask[:,720:1440]; cache[:,720:1440] = ocean_mask[:,0:720]
ocean_mask_720_1440 = cache

#######################################
# DATA INPUT and preprocessing 1971-2000
#######################################

# CESM Ensemble  mean AND rx5day 0.9375*1.25

CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat_CESM = nc_fid.variables['lat'][:]
lon_CESM = nc_fid.variables['lon'][:]
time_CESM = nc_fid.variables['time'][66:85]
mean_prep_CESM = nc_fid.variables['mean_precip'][66:85,:,:]
std_prep_CESM = nc_fid.variables['std_precip'][66:85,:,:]
RX5DAY_CESM = nc_fid.variables['rx5day'][66:85,:,:]

# interpolate the ocean_mask to CESM resolution
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440,kind='linear')
ocean_mask_192_288 = f(lon_CESM, lat_CESM)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1
lon_CESMs,lat_CESMs,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,ocean_mask_192_288)		

# 20 years mean
mean_CESMs = stats.nanmean(mean_prep_CESM,axis =0)
std_CESMs = stats.nanmean(std_prep_CESM,axis =0)
# mask off the ocean
mean_prep_CESM_masked = np.multiply(mean_CESMs,ocean_mask_192_288)
std_prep_CESM_masked = np.multiply(std_CESMs,ocean_mask_192_288)
RX5DAY_3D_CESM_masked = np.multiply(RX5DAY_CESM,ocean_mask_192_288)
lon_CESMs,lat_CESMs,mean_prep_CESM_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,mean_prep_CESM_masked)
_,_,std_prep_CESM_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,std_prep_CESM_masked)
_,_,RX5DAY_3D_CESM_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,RX5DAY_3D_CESM_masked)
RX5DAY_CESM_masked_clipped = stats.nanmean(RX5DAY_3D_CESM_masked_clipped,axis = 0)

"""
#The APHRO observational precipitatin data 0.5*0.5 degree 
"""
APHRO_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/APHRO_P_MA_050deg_V1101_1951_2007.nc'
nc_fid = nc4.Dataset(APHRO_file,mode='r')
lat_APHRO = nc_fid.variables['lat'][:]
lon_APHRO = nc_fid.variables['lon'][:]
date_APHRO = nc_fid.variables['time'][:]
prep_APHRO = nc_fid.variables['precip'][:]
# missing_value = nc_fid.variables['PREC'].missing_value
# fail_value = nc_fid.variables['PREC']._FillValue
prep_APHRO[prep_APHRO<0] =0

mean_prep_APHRO = np.empty((20,len(lat_APHRO),len(lon_APHRO)));mean_prep_APHRO[:] = np.nan
std_prep_APHRO = np.empty((20,len(lat_APHRO),len(lon_APHRO)));std_prep_APHRO[:] = np.nan
for iyear in range(1986,2005):
	# print iyear
	if leapyear(iyear):
		layer_b = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+153][0]
		layer_e = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+243][0]
	else:
		layer_b = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+152][0]
		layer_e = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+242][0]
	precp_year_data = prep_APHRO[layer_b:layer_e+1,:,:] 
	mean_prep_APHRO[iyear-1986,:,:] = stats.nanmean(precp_year_data,axis = 0)
	std_prep_APHRO[iyear-1986,:,:] = stats.nanstd(precp_year_data,axis = 0)
	
# 20 years mean	
mean_APHROs = stats.nanmean(mean_prep_APHRO,axis =0)
std_APHROs = stats.nanmean(std_prep_APHRO,axis =0)
# interpolate the APHRO resolution mean and std
f = interp2d(lon_APHRO, lat_APHRO, mean_APHROs,kind='linear')
mean_prep_APHRO_interp = f(lon_CESMs, lat_CESMs)
f = interp2d(lon_APHRO, lat_APHRO, std_APHROs,kind='linear')
std_prep_APHRO_interp = f(lon_CESMs, lat_CESMs)
# mask off the ocean: IT IS NOT NECESSSARYily true FOR THE aphro DATA
# mean_prep_APHRO_masked = mean_prep_APHRO_interp
# std_prep_APHRO_masked = std_prep_APHRO_interp
mean_prep_APHRO_masked = np.multiply(mean_prep_APHRO_interp,ocean_mask_192_288s)
std_prep_APHRO_masked = np.multiply(std_prep_APHRO_interp,ocean_mask_192_288s)	
# clip SCI region for plot
_,_,mean_prep_APHRO_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESMs,lat_CESMs,mean_prep_APHRO_masked)
_,_,std_prep_APHRO_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESMs,lat_CESMs,std_prep_APHRO_masked)
mean_prep_APHRO_masked_clipped[mean_prep_APHRO_masked_clipped<0]=np.nan
std_prep_APHRO_masked_clipped[std_prep_APHRO_masked_clipped<0]=np.nan


"""
# AROPH RX5DAY  0.5*0.5 degree
"""
AROPH_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(AROPH_file,mode='r')
lat_AROPH = nc_fid.variables['lat'][:]
lon_AROPH = nc_fid.variables['lon'][:]
time_AROPH = nc_fid.variables['time'][20:50]
RX5DAY_AROPH = nc_fid.variables['rx5day'][20:50]


# interpolate the  AROPH RX5DAY to CESM MODEL resolution 
size_RX5DAY_AROPH = np.shape(RX5DAY_AROPH)
RX5DAY_AROPH_3D_interp = np.empty((size_RX5DAY_AROPH[0],len(lat_CESMs),len(lon_CESMs)))
for ilayer in range(size_RX5DAY_AROPH[0]):
	cache = RX5DAY_AROPH[ilayer,:,:]
	f = interp2d(lon_AROPH, lat_AROPH, cache,kind='linear')
	RX5DAY_AROPH_3D_interp[ilayer,:,:] = f(lon_CESMs, lat_CESMs)

RX5DAY_AROPH_3D_interp_masked = np.multiply(RX5DAY_AROPH_3D_interp,ocean_mask_192_288s)
# RX5DAY_AROPH_3D_interp_masked = RX5DAY_AROPH;
lon_AROPHs,lat_AROPHs,RX5DAY_AROPH_3D_interp_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESMs,lat_CESMs,RX5DAY_AROPH_3D_interp_masked)
RX5DAY_AROPH_interp_masked_clipped = stats.nanmean(RX5DAY_AROPH_3D_interp_masked_clipped, axis = 0)
# print RX5DAY_AROPH_interp_masked_clipped



######################################
##          plotting            ######
######################################
# climatology mean

fig = plt.figure(facecolor='White',figsize=[6.5,4]);plot_setup();
p_value = np.zeros((np.shape(mean_prep_APHRO_masked_clipped)))
ax = plt.subplot(2,3,1)
colormap ='YlGnBu';pad= 5;cb_min = 0;cb_max =15
spatial_figure(ax,mean_prep_CESM_masked_clipped,lon_CESMs,lat_CESMs,colormap,cb_min,cb_max,p_value, tb_lef=True,tb_bot=False )
ax.annotate('Mean',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='vertical')
ax.annotate('CESM1 ensemble mean',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax = plt.subplot(2,3,2);
spatial_figure(ax,mean_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,colormap,cb_min,cb_max,p_value, tb_lef=False,tb_bot=False)
ax.annotate('APHRODITE',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
ax = plt.subplot(2,3,4);
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('Std',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='vertical')
spatial_figure(ax,std_prep_CESM_masked_clipped,lon_CESMs,lat_CESMs,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=True)

				
ax = plt.subplot(2,3,5)
colormesh1 =spatial_figure(ax,std_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
cbar_ax = fig.add_axes([0.10, 0.10, 0.55, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,15.1,2.5))
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$')

ax = plt.subplot(2,3,3);colormap ='BrBG'; cb_min = -5;cb_max =5
spatial_figure_norm(ax,mean_prep_CESM_masked_clipped-mean_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('CESM1 - APHRODITE',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
ax = plt.subplot(2,3,6);
colormesh2 =spatial_figure_norm(ax,std_prep_CESM_masked_clipped-std_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.70, 0.10, 0.28, 0.02])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both',orientation='horizontal',ticks=np.arange(-5.0,5.1,2.0))
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$')

plt.subplots_adjust(left=0.08, bottom=0.20, right=0.98, top=0.90, wspace=0.05, hspace=0.01); 
plt.savefig('Fig2.png', format='png', dpi=1000)
plt.show()
######################################
##       statistics         ##########
######################################

print 'ASIA'
print ' MEAN' 
region =rergion_dic['ASIA'][:]
_,_,mean_ccesm = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,mean_prep_CESM_masked_clipped)
_,_,mean_aphro = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,mean_prep_APHRO_masked_clipped)
print stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print (stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0))/stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print ' std' 
_,_,mean_ccesm = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,std_prep_CESM_masked_clipped)
_,_,mean_aphro = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,std_prep_APHRO_masked_clipped) 
print stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print (stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0))/stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)


print 'EA'
print ' MEAN' 
region =rergion_dic['EA'][:]
_,_,mean_ccesm = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,mean_prep_CESM_masked_clipped)
_,_,mean_aphro = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,mean_prep_APHRO_masked_clipped)
print stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print (stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0))/stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print ' std' 
_,_,mean_ccesm = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,std_prep_CESM_masked_clipped)
_,_,mean_aphro = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,std_prep_APHRO_masked_clipped) 
print stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print (stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0))/stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)


print 'SA'
print ' MEAN' 
region =rergion_dic['SA'][:]
_,_,mean_ccesm = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,mean_prep_CESM_masked_clipped)
_,_,mean_aphro = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,mean_prep_APHRO_masked_clipped)
print stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print (stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0))/stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print ' std' 
_,_,mean_ccesm = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,std_prep_CESM_masked_clipped)
_,_,mean_aphro = range_clip(region[0],region[1],region[2],region[3],lon_AROPHs,lat_AROPHs,std_prep_APHRO_masked_clipped) 
print stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)
print (stats.nanmean(stats.nanmean(mean_ccesm,axis=1),axis=0) - stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0))/stats.nanmean(stats.nanmean(mean_aphro,axis=1),axis=0)


"""

#RX5DAY :  mean and pdf
plt.figure(2, facecolor='White',figsize=[15,10])
spst=[2,3]
plt.subplot(spst[0],spst[1],1)
colormap ='jet'; figure_title = 'AROPH mean'
spatial_figure(RX5DAY_AROPH_interp_masked_clipped,lon_AROPHs,lat_AROPHs,figure_title,colormap,0.1,200,'mm')
plt.subplot(spst[0],spst[1],2)
colormap ='jet'; figure_title = 'CESM mean'
spatial_figure(RX5DAY_CESM_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,0.1,200,'mm')

####PDF
from scipy.stats.kde import gaussian_kde
from numpy import linspace
ax3=plt.subplot(spst[0],spst[1],3)
size =np.shape(RX5DAY_3D_CESM_masked_clipped)
RX5DAY_AROPH_m_c_reshape = np.reshape(RX5DAY_3D_CESM_masked_clipped,(size[0]*size[1]*size[2],1))
RX5DAY_AROPH_m_c_reshape [RX5DAY_AROPH_m_c_reshape <=0]= np.nan
mask = ~np.isnan(RX5DAY_AROPH_m_c_reshape); 
RX5DAY_AROPH_m_c_reshape = RX5DAY_AROPH_m_c_reshape[mask]
kde = gaussian_kde( RX5DAY_AROPH_m_c_reshape )
dist_space = linspace( np.nanmin(RX5DAY_AROPH_m_c_reshape), np.nanmax(RX5DAY_AROPH_m_c_reshape), 100 )
ax3.plot( dist_space, kde(dist_space),'k-',alpha=1, label = 'CESM' )
# rv = norm()
# plt.plot(RX5DAY_AROPH_m_c_reshape, rv.pdf(RX5DAY_AROPH_m_c_reshape), 'k-', lw=2,label = 'AROPH')

size =np.shape(RX5DAY_AROPH_3D_interp_masked_clipped)
RX5DAY_CESM_m_c_reshape = np.reshape(RX5DAY_AROPH_3D_interp_masked_clipped,(size[0]*size[1]*size[2],1))
RX5DAY_CESM_m_c_reshape [RX5DAY_CESM_m_c_reshape <=0]= np.nan
mask = ~np.isnan(RX5DAY_CESM_m_c_reshape); RX5DAY_CESM_m_c_reshape = RX5DAY_CESM_m_c_reshape[mask]
kde = gaussian_kde( RX5DAY_CESM_m_c_reshape )
dist_space = linspace( np.nanmin(RX5DAY_CESM_m_c_reshape), np.nanmax(RX5DAY_CESM_m_c_reshape), 100 )
ax3.plot( dist_space, kde(dist_space),'r-',alpha=1,label = 'AROPH' )
ax3.set_xlim([-20,400])ax = plt.subplot(3,1,2)
legend = ax3.legend(loc=10, shadow=True,fontsize=15)	
# rv = norm()
# plt.plot(RX5DAY_CESM_m_c_reshape, rv.pdf(RX5DAY_CESM_m_c_reshape), 'r-', lw=2,label = 'CESM')

#### Trend of the rx5day 
size = np.shape(RX5DAY_3D_CESM_masked_clipped)
slope_CESM =np.empty((size[1],size[2]));
for row in range(size[1]):
	for column in range(size[2]):
		cache = RX5DAY_3D_CESM_masked_clipped[:,row,column]
		slope_CESM[row,column], _, _, _, _ = stats.linregress(time_CESM,cache)
slope_CESM_masked =np.multiply(slope_CESM,ocean_mask_192_288s)
slope_CESM_masked[slope_CESM_masked==0] = np.nan

size = np.shape(RX5DAY_AROPH_3D_interp_masked_clipped)
slope_AROPH =np.empty((size[1],size[2]))
for row in range(size[1]):
	for column in range(size[2]):
		cache = RX5DAY_AROPH_3D_interp_masked_clipped[:,row,column]
		slope_AROPH[row,column], _, _, _, _ = stats.linregress(time_AROPH,cache)

slope_AROPH_masked =np.multiply(slope_AROPH,ocean_mask_192_288s)
slope_AROPH_masked[slope_AROPH_masked==0] = np.nan
		
plt.subplot(spst[0],spst[1],4)
colormap ='bwr'; figure_title = 'AROPH TREND'
spatial_figure(slope_AROPH_masked,lon_AROPHs,lat_AROPHs,figure_title,colormap,-0.5,0.8,'mm/year')
plt.subplot(spst[0],spst[1],5)
colormap ='bwr'; figure_title = 'CESM TREND'
spatial_figure(slope_CESM_masked,lon_CESMs,lat_CESMs,figure_title,colormap,-0.5,0.8,'mm/year')
plt.show()
"""
"""

# CESM MODEL Ensemble  RANGE VS APHPR MEAN FROM 1951-2005
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical/*.nc'

files=sorted(glob.glob(input_path))
Ensemble _no=0
mean_prep = np.empty((30,55)); mean_prep[:]=np.nan
for file in files:
	nc_fid = nc4.Dataset(file,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	time = nc_fid.variables['lon'][31:86]
	value = nc_fid.variables['mean_precip'][31:86,:,:]
	value_masked =np.multiply(value,ocean_mask_192_288)
	lons,lats,value_m_clipped =  range_clip(region[0],region[1],region[2],region[3],lon,lat,value_masked)
	mean_prep[Ensemble _no,:] = stats.nanmean(stats.nanmean(value_m_clipped,axis=2),axis=1)
	Ensemble _no =Ensemble _no+1

# APHRO mean prep evelution
APHRO_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(APHRO_file,mode='r')
lat_APHRO = nc_fid.variables['lat'][:]
lon_APHRO = nc_fid.variables['lon'][:]
date_APHRO = nc_fid.variables['time'][0:55]
prep_APHRO = nc_fid.variables['mean_precip'][0:55,:,:]

prep_APHRO[prep_APHRO<0] =np.nan
APHRO_series = stats.nanmean(stats.nanmean(prep_APHRO,axis=2),axis=1)
nc_fid.close()

# print mean_prep_APHRO_series
time_series_range = np.empty([3,55])
# att_value = np.multiply(att_value,oceanmask)
time_series_range[1,:] = stats.nanmean(mean_prep,axis = 0)
time_series_range[0,:] = np.nanmax(mean_prep,axis = 0)
time_series_range[2,:] = np.nanmin(mean_prep,axis = 0)	

time_series_range[1,:]=time_series_range[1,:]
time_series_range[0,:]=time_series_range[0,:]
time_series_range[2,:] = time_series_range[2,:]

plt.figure(1, facecolor='White',figsize=[10,10])
sub_plot=1
ax1=plt.subplot(1,1,sub_plot)
ax1.plot(date_APHRO,time_series_range[1,:],'-',color="red",linewidth=2,label='CESM')
ax1.fill_between(date_APHRO,time_series_range[0,:],time_series_range[2,:], where=time_series_range[0,:]>=time_series_range[2,:], color="red", alpha = 0.2) 

ax1.plot(date_APHRO,APHRO_series,'-',color="black",linewidth=2,label='APHRO')
"""


