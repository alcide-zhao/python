# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This is to plot the precipitation extreme climatologies
"""
import scipy
import numpy as np
import netCDF4 as nc4

import math
from scipy import stats
from scipy import arange
from scipy.interpolate import interp2d  as interp2d

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from scipy.stats.kde import gaussian_kde
from numpy import linspace

# from scipy.interpolate import interp2d
import os;import site
lib_path = os.path.join(
	os.path.realpath(
		os.path.dirname(__file__)
	), 
	os.path.pardir, 
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


		
def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,p_value,tb_lef=True,tb_bot=True ): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	# lons[lons>180]-=360; 
	# calculate the origin of the map
	lon_0 = lons.mean(); 
	lat_0 = lats.mean(); 
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	# lon_bin = int((lon_e -lon_b)/5)
	# lat_bin = int((lat_e -lat_b)/5)
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, data)
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='.',color='k',zorder=10)
	
	# Add Grid Lines
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines();# map.drawstates(); map.drawcountries()
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); cmap.set_under('b')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min-0.01, vmax=colorbar_max) #,latlon=True
	return colormesh


########################################
# 0. Vaeiables and functions
########################################
rergion_dic={'ASIA':[60,155,-5,55],'EA':[110,145,20,50],'SA':[65,100,5,28]}
region = rergion_dic['ASIA'][:]

# ocean land mask
import scipy.io as spio
oceanmask=spio.loadmat('/home/s1667168/coding/python/external_data/world.oceanmask.1440x720')['landocean']
ocean_mask = np.flipud(oceanmask)
cache = np.empty((720,1440))
cache[:,0:720] = ocean_mask[:,720:1440]; cache[:,720:1440] = ocean_mask[:,0:720]
ocean_mask_720_1440 = cache


#######################################
# DATA INPUT and preprocessing 1996-2005
#######################################

"""
# CESM ensumble mean AND RX5DAY 0.9375*1.25
"""
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][66:86]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:86,:,:]-nc_fid.variables['rx5day'][30:50,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][66:86,:,:]-nc_fid.variables['r95p'][30:50,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][66:86,:,:] -nc_fid.variables['cwd'][30:50,:,:] ,axis=0);
r10 = stats.nanmean(nc_fid.variables['r10'][66:86,:,:]-nc_fid.variables['r10'][30:50,:,:],axis=0);
ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:] - nc_fid.variables['total_precip'][30:50,:,:] ,axis=0);
# r95p= 100*np.divide(r95p,ptot)

# cdd = nc_fid.variables['cdd'][51:81,:,:]
# r10 = nc_fid.variables['r10'][51:81,:,:]

# interpolate the ocean_mask to required resolution
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440, kind='linear')
ocean_mask_192_288 = f(lon, lat)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1
RX5DAY_3D_masked = np.multiply(RX5DAY,ocean_mask_192_288)
lons,lats,RX5DAY_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,RX5DAY_3D_masked)
r10_3D_masked = np.multiply(r10,ocean_mask_192_288)
lons,lats,r10_masked_clipped= range_clip(region[0],region[1],region[2],region[3],lon,lat,r10_3D_masked)
cwd_3D_masked = np.multiply(cwd,ocean_mask_192_288)
lons,lats,cwd_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cwd_3D_masked)
r95p_3D_masked = np.multiply(r95p,ocean_mask_192_288)
lons,lats,r95p_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r95p_3D_masked)

"""
# AROPH RX5DAY  0.5*0.5 degree
"""
_,_,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon,lat,ocean_mask_192_288)	
AROPH_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(AROPH_file,mode='r')
lat_AROPH = nc_fid.variables['lat'][:]
lon_AROPH = nc_fid.variables['lon'][:]
time_AROPH = nc_fid.variables['time'][35:55]
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][35:55]-nc_fid.variables['rx5day'][0:20],axis=0);
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][35:55]-nc_fid.variables['cwd'][0:20],axis=0)
r10_AROPH = stats.nanmean(nc_fid.variables['r10'][35:55]-nc_fid.variables['r10'][0:20],axis=0); r10_AROPH[np.isnan(r10_AROPH) ]=0.0
r95p_AROPH = stats.nanmean(nc_fid.variables['r95p'][35:55]-nc_fid.variables['r95p'][0:20],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0
ptot_AROPH = stats.nanmean(nc_fid.variables['total_precip'][35:55]-nc_fid.variables['total_precip'][0:20],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0
# r95p_AROPH =100*np.divide(r95p_AROPH,ptot_AROPH);r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0


# CDD_AROPH = stats.nanmean(nc_fid.variables['cdd'][20:50],axis=0)
# R10_AROPH = stats.nanmean(nc_fid.variables['r10'][20:50],axis=0)

f1 = interp2d(lon_AROPH, lat_AROPH, RX5DAY_AROPH,kind='quintic')
RX5DAY_AROPH_interp = f1(lons, lats)
f2 = interp2d(lon_AROPH, lat_AROPH, r10_AROPH,kind='quintic')
r10_AROPH_interp = f2(lons, lats)
f3 = interp2d(lon_AROPH, lat_AROPH, CWD_AROPH,kind='quintic')
CWD_AROPH_interp = f3(lons, lats)
f4 = interp2d(lon_AROPH, lat_AROPH, r95p_AROPH,kind='quintic')
r95p_AROPH_interp = f4(lons, lats)

RX5DAY_AROPH_interp_masked = np.multiply(RX5DAY_AROPH_interp,ocean_mask_192_288s)
cwd_AROPH_interp_masked = np.multiply(CWD_AROPH_interp,ocean_mask_192_288s)
r10_AROPH_interp_masked = np.multiply(r10_AROPH_interp,ocean_mask_192_288s)
r95p_AROPH_interp_masked = np.multiply(r95p_AROPH_interp,ocean_mask_192_288s)



####################################
###plotting spatial mapss
####################################
# fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(10, 7), facecolor='White');plot_setup();
fig = plt.figure(facecolor='White',figsize=[5,7]);plot_setup();
p_value = np.zeros((np.shape(RX5DAY_masked_clipped)))
#############################
#RX5DAY
#############################
ax = plt.subplot(4,2,1);colormap ='BrBG';pad= 5;cb_min = -15;cb_max =15
colormesh1 =spatial_figure(ax,RX5DAY_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate('(a)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('CESM1',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('RX5DAY (mm)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
		
ax = plt.subplot(4,2,2);
spatial_figure(ax,RX5DAY_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('(b)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('APHRODITE',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
cbar_ax = fig.add_axes([0.15, 0.745, 0.80, 0.01])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)

#############################
# r10
#############################
ax = plt.subplot(4,2,3);colormap ='BrBG';pad= 5;cb_min = -8;cb_max =8
colormesh4 =spatial_figure(ax,r10_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate('R10 (days)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
# char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$',x=0.62, y=0.515)
ax = plt.subplot(4,2,4);
spatial_figure(ax,r10_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value, tb_lef=False,tb_bot=False)
ax.annotate('(d)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.15, 0.51, 0.80, 0.01])
char = fig.colorbar(colormesh4,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)
#############################
# CWD
#############################
ax = plt.subplot(4,2,5);colormap ='BrBG';pad= 5;cb_min = -8;cb_max =8
colormesh3 =spatial_figure(ax,cwd_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate('CWD (days)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(e)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

# char.set_label('day')
ax = plt.subplot(4,2,6);
spatial_figure(ax,cwd_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('(f)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.15, 0.275, 0.80, 0.01])
char = fig.colorbar(colormesh3,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)
#############################
# r95p
#############################
ax = plt.subplot(4,2,7);colormap ='BrBG';pad= 5;cb_min = -30;cb_max =30
colormesh2 =spatial_figure(ax,r95p_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=True)
ax.annotate('(g)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('R95P (mm)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax = plt.subplot(4,2,8);
spatial_figure(ax,r95p_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
ax.annotate('(h)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.15, 0.025, 0.80, 0.01])
char = fig.colorbar(colormesh2,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)


					
plt.subplots_adjust(left=0.15, bottom=0.06, right=0.98, top=0.95, wspace=0.15, hspace=0.22); 
plt.savefig('historical_changes.png', format='png', dpi=800)
plt.show()
