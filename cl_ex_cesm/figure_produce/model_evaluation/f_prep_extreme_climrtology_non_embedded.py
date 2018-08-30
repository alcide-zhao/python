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
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:86,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][66:86,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][66:86,:,:],axis=0);
r10 = stats.nanmean(nc_fid.variables['r10'][66:86,:,:],axis=0);
ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
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
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][35:55],axis=0);
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][35:55],axis=0)
r10_AROPH = stats.nanmean(nc_fid.variables['r10'][35:55],axis=0); r10_AROPH[np.isnan(r10_AROPH) ]=0.0
r95p_AROPH = stats.nanmean(nc_fid.variables['r95p'][35:55],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0
ptot_AROPH = stats.nanmean(nc_fid.variables['total_precip'][35:55],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0
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
fig = plt.figure(facecolor='White',figsize=[10,7]);plot_setup();
p_value = np.zeros((np.shape(RX5DAY_masked_clipped)))
#############################
#RX5DAY
#############################
ax = plt.subplot(4,4,1);colormap ='YlGnBu';pad= 5;cb_min = 0.000000000001;cb_max =250
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
		
ax = plt.subplot(4,4,2);
spatial_figure(ax,RX5DAY_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('(b)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('APHRODITE',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
cbar_ax = fig.add_axes([0.12, 0.745, 0.4, 0.01])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)

#############################
# r10
#############################
ax = plt.subplot(4,4,5);colormap ='YlGnBu';pad= 5;cb_min = 0;cb_max =40
colormesh4 =spatial_figure(ax,r10_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate('R10 (days)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(e)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
# char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$',x=0.62, y=0.515)
ax = plt.subplot(4,4,6);
spatial_figure(ax,r10_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value, tb_lef=False,tb_bot=False)
ax.annotate('(f)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.12, 0.51, 0.4, 0.01])
char = fig.colorbar(colormesh4,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)
#############################
# CWD
#############################
ax = plt.subplot(4,4,9);colormap ='YlGnBu';pad= 5;cb_min = 0;cb_max =80
colormesh3 =spatial_figure(ax,cwd_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate('CWD (days)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(i)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

# char.set_label('day')
ax = plt.subplot(4,4,10);
spatial_figure(ax,cwd_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('(j)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.12, 0.275, 0.4, 0.01])
char = fig.colorbar(colormesh3,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)
#############################
# r95p
#############################
ax = plt.subplot(4,4,13);colormap ='YlGnBu';pad= 5;cb_min = 0;cb_max =350
colormesh2 =spatial_figure(ax,r95p_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=True)
ax.annotate('(m)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('R95P (mm)',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax = plt.subplot(4,4,14);
spatial_figure(ax,r95p_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
ax.annotate('(n)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.12, 0.025, 0.4, 0.01])
char = fig.colorbar(colormesh2,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.1,0.2)*(cb_max-cb_min)+cb_min)

region = rergion_dic['EA'][:]
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][66:86]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:86,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][66:86,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][66:86,:,:],axis=0);
r10 = stats.nanmean(nc_fid.variables['r10'][66:86,:,:],axis=0);
ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
# r95p= 100*np.divide(r95p,ptot)
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
time_AROPH = nc_fid.variables['time'][35:54]
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][35:54],axis=0);
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][35:54],axis=0)
r10_AROPH = stats.nanmean(nc_fid.variables['r10'][35:54],axis=0); r10_AROPH[np.isnan(r10_AROPH) ]=0.0
r95p_AROPH = stats.nanmean(nc_fid.variables['r95p'][35:54],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0
ptot_AROPH = stats.nanmean(nc_fid.variables['total_precip'][35:55],axis=0); ptot_AROPH[np.isnan(ptot_AROPH) ]=0.0
# r95p_AROPH =100*np.divide(r95p_AROPH,ptot_AROPH);r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0

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


#############################
#RX5DAY
#############################
ax = plt.subplot(4,4,4);


size =np.shape(RX5DAY_AROPH_interp_masked)
hist_AROPH = np.reshape(RX5DAY_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(RX5DAY_masked_clipped)
hist_CESM = np.reshape(RX5DAY_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( -20, 250, 200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,250.1,50.0));ax.set_yticks(np.arange(0.0,1.21,0.4))
ax.set_xlim([0,250]);ax.set_ylim([0,1.2])
ax.annotate('(d)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('East Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')

print 'EA','RX5',np.nanmean(np.nanmean(RX5DAY_masked_clipped-RX5DAY_AROPH_interp_masked,axis=1),axis=0)			
#############################
# r10
#############################
size =np.shape(r10_AROPH_interp_masked)
hist_AROPH = np.reshape(r10_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(r10_masked_clipped)
hist_CESM = np.reshape(r10_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
	
ax = plt.subplot(4,4,8)						
dist_space = linspace( 0, 40, 200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,40.1,8));ax.set_yticks(np.arange(0.0,7.51,2.5))
ax.set_xlim([0,40.1]);ax.set_ylim([0,7.5])
ax.annotate('(h)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
print 'EA','R10',np.nanmean(np.nanmean(r10_masked_clipped-r10_AROPH_interp_masked,axis=1),axis=0)
#############################
# CWD
#############################
ax = plt.subplot(4,4,12);
# Histogram embedded in the subplot
size =np.shape(cwd_AROPH_interp_masked)
hist_AROPH = np.reshape(cwd_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 80, 200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0.0,80.1,16.0));ax.set_yticks(np.arange(0.0,16.1,5.0))
ax.set_xlim([0,80.1]);ax.set_ylim([0,15.01])
ax.annotate('(l)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
print 'EA','CWD',np.nanmean(np.nanmean(cwd_masked_clipped-cwd_AROPH_interp_masked,axis=1),axis=0)
				
#############################
# r95p
#############################
ax = plt.subplot(4,4,16);
# Histogram embedded in the subplot
size =np.shape(r95p_AROPH_interp_masked)
hist_AROPH = np.reshape(r95p_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 350,200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xlim([0,350]);ax.set_xticks(np.arange(0.0,350.1,70))
ax.set_ylim([0,0.9]);ax.set_yticks(np.arange(0.0,0.91,0.3))
ax.annotate('(p)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
print 'EA','95',np.nanmean(np.nanmean(r95p_masked_clipped-r95p_AROPH_interp_masked,axis=1),axis=0)

									
region = rergion_dic['SA'][:]
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][66:85]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:85,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][66:85,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][66:85,:,:],axis=0);
r10 = stats.nanmean(nc_fid.variables['r10'][66:85,:,:],axis=0);
ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
# r95p= 100*np.divide(r95p,ptot)
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
time_AROPH = nc_fid.variables['time'][35:54]
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][35:54],axis=0);
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][35:54],axis=0)
r10_AROPH = stats.nanmean(nc_fid.variables['r10'][35:54],axis=0); r10_AROPH[np.isnan(r10_AROPH) ]=0.0
r95p_AROPH = stats.nanmean(nc_fid.variables['r95p'][35:54],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0
ptot_AROPH = stats.nanmean(nc_fid.variables['total_precip'][35:54],axis=0); ptot_AROPH[np.isnan(ptot_AROPH) ]=0.0
# r95p_AROPH =100*np.divide(r95p_AROPH,ptot_AROPH);r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0

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

#############################
#RX5DAY
#############################
ax = plt.subplot(4,4,3);
from scipy.stats.kde import gaussian_kde
from numpy import linspace

size =np.shape(RX5DAY_AROPH_interp_masked)
hist_AROPH = np.reshape(RX5DAY_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(RX5DAY_masked_clipped)
hist_CESM = np.reshape(RX5DAY_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( -20, 250, 200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,250.1,50.0));ax.set_yticks(np.arange(0.0,0.91,0.3))
ax.set_xlim([0,250]);ax.set_ylim([0,0.9])

ax.annotate('(c)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('South Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
legend=plt.legend(loc=8,shadow=False);
legend.get_frame().set_facecolor('None');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)
print 'SA','RX5',np.nanmean(np.nanmean(RX5DAY_masked_clipped-RX5DAY_AROPH_interp_masked,axis=1),axis=0)


#############################
# r10
#############################
size =np.shape(r10_AROPH_interp_masked)
hist_AROPH = np.reshape(r10_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(r10_masked_clipped)
hist_CESM = np.reshape(r10_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
	
ax = plt.subplot(4,4,7)						
dist_space = linspace( 0, 40, 200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,40.1,8));ax.set_yticks(np.arange(0.0,4.21,1.4))
ax.set_xlim([0,40]);ax.set_ylim([0,4.2])
ax.annotate('(g)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
print 'SA','R10',np.nanmean(np.nanmean(r10_masked_clipped-r10_AROPH_interp_masked,axis=1),axis=0)
				

#############################
# CWD
#############################
ax = plt.subplot(4,4,11);
# Histogram embedded in the subplot
size =np.shape(cwd_AROPH_interp_masked)
hist_AROPH = np.reshape(cwd_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 80, 200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0.0,80.1,16.0));ax.set_yticks(np.arange(0.0,3.1,1))
ax.set_xlim([0,80]);ax.set_ylim([0,3])
ax.annotate('(k)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
print 'SA','CWD',np.nanmean(np.nanmean(cwd_masked_clipped-cwd_AROPH_interp_masked,axis=1),axis=0)
					
#############################
# r95p
#############################
ax = plt.subplot(4,4,15);
# Histogram embedded in the subplot
size =np.shape(r95p_AROPH_interp_masked)
hist_AROPH = np.reshape(r95p_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 350,200)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'CESM1' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xlim([0,350]);ax.set_xticks(np.arange(0.0,350.1,70));
ax.set_ylim([0,0.75]);ax.set_yticks(np.arange(0.0,0.751,0.25))
ax.annotate('(o)',xy=(0.90, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
print 'SA','95',np.nanmean(np.nanmean(r95p_masked_clipped-r95p_AROPH_interp_masked,axis=1),axis=0)


#################################################################
######PDFs
#################################################################
region = rergion_dic['EA'][:]
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/result2.0/his2/'
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()
_,_,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon,lat,ocean_mask_192_288)
for ensumble_member in range(0,25): 
	nc_f = text_content[ensumble_member][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	time = nc_fid.variables['time'][66:86]
	RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:86,:,:],axis=0);
	r95p =stats.nanmean( nc_fid.variables['r95p'][66:86,:,:],axis=0);
	cwd = stats.nanmean(nc_fid.variables['cwd'][66:86,:,:],axis=0);
	r10 = stats.nanmean(nc_fid.variables['r10'][66:86,:,:],axis=0);
	ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
	# r95p= 100*np.divide(r95p,ptot)
	# interpolate the ocean_mask to required resolution
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
	
	ax = plt.subplot(4,4,4);
	size =np.shape(RX5DAY_masked_clipped)
	hist_CESM = np.reshape(RX5DAY_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
	dist_space = linspace( 0, 250, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )
	
	ax = plt.subplot(4,4,8)
	size =np.shape(r10_masked_clipped)
	hist_CESM = np.reshape(r10_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	dist_space = linspace( 0, 40, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )
	
	ax = plt.subplot(4,4,12);
	size =np.shape(cwd_masked_clipped)
	hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]				
	dist_space = linspace( 0, 80, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )
	
	ax = plt.subplot(4,4,16);
	size =np.shape(r95p_masked_clipped)
	hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]			
	dist_space = linspace( 0, 350,100)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )


region = rergion_dic['SA'][:]
input_path = '//exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/result2.0/his2/'
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()
_,_,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon,lat,ocean_mask_192_288)
for ensumble_member in range(0,25): 
	nc_f = text_content[ensumble_member][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	time = nc_fid.variables['time'][66:85]
	RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:85,:,:],axis=0);
	r95p =stats.nanmean( nc_fid.variables['r95p'][66:85,:,:],axis=0);
	cwd = stats.nanmean(nc_fid.variables['cwd'][66:85,:,:],axis=0);
	r10 = stats.nanmean(nc_fid.variables['r10'][66:85,:,:],axis=0);
	ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
	# r95p= 100*np.divide(r95p,ptot)
	# interpolate the ocean_mask to required resolution
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
		
	ax = plt.subplot(4,4,3);
	size =np.shape(RX5DAY_masked_clipped)
	hist_CESM = np.reshape(RX5DAY_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
	dist_space = linspace( 0, 250, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )
	
	ax = plt.subplot(4,4,7)
	size =np.shape(r10_masked_clipped)
	hist_CESM = np.reshape(r10_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	dist_space = linspace( 0, 40, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )
	
	ax = plt.subplot(4,4,11);
	size =np.shape(cwd_masked_clipped)
	hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]				
	dist_space = linspace( 0, 80, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )
	
	ax = plt.subplot(4,4,15);
	size =np.shape(r95p_masked_clipped)
	hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
	hist_CESM [hist_CESM <=0]= np.nan
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]			
	dist_space = linspace( 0, 350, 200)
	kde = gaussian_kde( hist_CESM )
	plt.plot(dist_space, 100*kde(dist_space),'k-',alpha=0.1, label = 'CESM1' )	
					
plt.subplots_adjust(left=0.10, bottom=0.06, right=0.98, top=0.95, wspace=0.15, hspace=0.22); 
plt.savefig('Fig3.png', format='png', dpi=800)
plt.show()
