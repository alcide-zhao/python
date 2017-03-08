# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This is to plot the climatology mean of precipitation and precipitation extremes
The observational precipitation data are from the APHRO spans from 1951 to 2007 with a resolution of 
The reference monthly RX5day precipitation indices are derived from the AROPH from 1901 to 2010
The model data are from the CESM ensumble means
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
from scipy.stats import norm
from mpl_toolkits.basemap import Basemap
import math
import site
# from scipy.interpolate import interp2d

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)

site.addsitedir(lib_path)

from lib import *

########################################
# 0. Vaeiables and functions
########################################

prep_att_ref = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_dic = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_list=['rx5day','r99p','sdii','cwd']
prep_unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}


region = rergion_dic['ASIA'][:]

# ocean land mask
oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720')['landocean']
ocean_mask = np.flipud(oceanmask)
cache = np.empty((720,1440))
cache[:,0:720] = ocean_mask[:,720:1440]; cache[:,720:1440] = ocean_mask[:,0:720]
ocean_mask_720_1440 = cache


def spatial_figure(data,lons,lats,plot_title,colormap,colorbar_min,colorbar_max,colorbar_title):
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	
	# calculate the origin of the map
	lon_0 = lons.mean()
	lat_0 = lats.mean()
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = int((lon_e -lon_b)/3)
	lat_bin = int((lat_e -lat_b)/3)
	
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, data)
	
	# Add Grid Lines
	map.drawparallels(np.arange(round(lat_b,1), round(lat_e,1), lat_bin), labels=[1,0,0,0],linewidth=0.0, fontsize=10)
	map.drawmeridians(np.arange(round(lon_b,1), round(lon_e,1), lon_bin), labels=[0,0,0,1],linewidth=0.0, fontsize=10)
	map.drawcoastlines()
	cmap = discrete_cmap(20,colormap)
	# cmap= reverse_colourmap(cmap)
	cmap.set_bad('w',alpha = 1.0)
	cmap.set_over('k')
	cmap.set_under('w')
	# cmap.set_over((1., 0., 0.))
	# cmap.set_under((0., 0., 1.))
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)
	cbar = map.colorbar(colormesh, location='bottom', pad="20%",extend='both')
	cbar.set_label(colorbar_title,fontsize=10)
	# Add Title
	plt.title(plot_title)
	
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

#######################################
# DATA INPUT and preprocessing 1956-2005
#######################################

# CESM ensumble mean AND rx5day 0.9375*1.25

CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat_CESM = nc_fid.variables['lat'][:]
lon_CESM = nc_fid.variables['lon'][:]
time_CESM = nc_fid.variables['time'][51:81]
mean_prep_CESM = nc_fid.variables['mean_precip'][51:81,:,:]
std_prep_CESM = nc_fid.variables['std_precip'][51:81,:,:]
RX5DAY_CESM = nc_fid.variables['rx5day'][51:81,:,:]

# interpolate the ocean_mask to CESM resolution
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440,kind='linear')
ocean_mask_192_288 = f(lon_CESM, lat_CESM)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1
lon_CESMs,lat_CESMs,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,ocean_mask_192_288)		

# 30 years mean
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
APHRO_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/APHRO_MA_050deg_V1101_1951_2007.nc'
nc_fid = nc4.Dataset(APHRO_file,mode='r')
lat_APHRO = nc_fid.variables['lat'][:]
lon_APHRO = nc_fid.variables['lon'][:]
date_APHRO = nc_fid.variables['time'][:]
prep_APHRO = nc_fid.variables['precip'][:]
# missing_value = nc_fid.variables['PREC'].missing_value
# fail_value = nc_fid.variables['PREC']._FillValue
prep_APHRO[prep_APHRO<0] =0

mean_prep_APHRO = np.empty((30,len(lat_APHRO),len(lon_APHRO)));mean_prep_APHRO[:] = np.nan
std_prep_APHRO = np.empty((30,len(lat_APHRO),len(lon_APHRO)));std_prep_APHRO[:] = np.nan
for iyear in range(1971,2001):
	# print iyear
	if leapyear(iyear):
		layer_b = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+153][0]
		layer_e = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+243][0]
	else:
		layer_b = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+152][0]
		layer_e = [layer for layer in range(len(date_APHRO)) if date_APHRO[layer] == iyear*1000+242][0]
	precp_year_data = prep_APHRO[layer_b:layer_e+1,:,:] 
	mean_prep_APHRO[iyear-1971,:,:] = stats.nanmean(precp_year_data,axis = 0)
	std_prep_APHRO[iyear-1971,:,:] = stats.nanstd(precp_year_data,axis = 0)
	
# 30 years mean	
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



# climatology mean
plt.figure(1, facecolor='White',figsize=[15,10])
spst=[2,3]
plt.subplot(spst[0],spst[1],1)
colormap ='jet'; figure_title = 'APHRO mean'
spatial_figure(mean_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,0.1,14,'mm/day')
plt.subplot(spst[0],spst[1],2)
colormap ='jet'; figure_title = 'CESM mean'
spatial_figure(mean_prep_CESM_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,0.1,14,'mm/day')
plt.subplot(spst[0],spst[1],3)
colormap ='PiYG'; figure_title = 'CESM-APHRO'
spatial_figure(mean_prep_CESM_masked_clipped-mean_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,-3,3,'mm/day')
plt.subplot(spst[0],spst[1],4)
colormap ='jet'; figure_title = 'APHRO std'
spatial_figure(std_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,0.1,14,'mm/day')
plt.subplot(spst[0],spst[1],5)
colormap ='jet'; figure_title = 'CESM std'
spatial_figure(std_prep_CESM_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,0.1,14,'mm/day')
plt.subplot(spst[0],spst[1],6)
colormap ='PiYG'; figure_title = 'CESM-APHRO'
spatial_figure(std_prep_CESM_masked_clipped-std_prep_APHRO_masked_clipped,lon_CESMs,lat_CESMs,figure_title,colormap,-3.5,3.5,'mm/day')

# plt.show()

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
ax3.set_xlim([-20,400])
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

# CESM MODEL ENSUMBLE RANGE VS APHPR MEAN FROM 1951-2005
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical/*.nc'

files=sorted(glob.glob(input_path))
ensumble_no=0
mean_prep = np.empty((30,55)); mean_prep[:]=np.nan
for file in files:
	nc_fid = nc4.Dataset(file,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	time = nc_fid.variables['lon'][31:86]
	value = nc_fid.variables['mean_precip'][31:86,:,:]
	value_masked =np.multiply(value,ocean_mask_192_288)
	lons,lats,value_m_clipped =  range_clip(region[0],region[1],region[2],region[3],lon,lat,value_masked)
	mean_prep[ensumble_no,:] = stats.nanmean(stats.nanmean(value_m_clipped,axis=2),axis=1)
	ensumble_no =ensumble_no+1

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

plt.show()
"""
