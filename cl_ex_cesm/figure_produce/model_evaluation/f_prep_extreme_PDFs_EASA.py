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
rergion_dic={'ASIA':[60,150,-15,55],'EA':[105,135,22.5,45],'SA':[70,105,10,28]}
region = rergion_dic['EA'][:]

# ocean land mask
import scipy.io as spio
oceanmask=spio.loadmat('/home/s1667168/coding/python/cl_ex_cesm/external_data/world.oceanmask.1440x720')['landocean']
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
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/10/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][51:81]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][51:81,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][51:81,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][51:81,:,:],axis=0);
sdii = stats.nanmean(nc_fid.variables['sdii'][51:81,:,:],axis=0);

# interpolate the ocean_mask to required resolution
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440, kind='linear')
ocean_mask_192_288 = f(lon, lat)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1
RX5DAY_3D_masked = np.multiply(RX5DAY,ocean_mask_192_288)
lons,lats,RX5DAY_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,RX5DAY_3D_masked)
sdii_3D_masked = np.multiply(sdii,ocean_mask_192_288)
lons,lats,sdii_masked_clipped= range_clip(region[0],region[1],region[2],region[3],lon,lat,sdii_3D_masked)
cwd_3D_masked = np.multiply(cwd,ocean_mask_192_288)
lons,lats,cwd_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cwd_3D_masked)
r95p_3D_masked = np.multiply(r95p,ocean_mask_192_288)
lons,lats,r95p_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r95p_3D_masked)

"""
# AROPH RX5DAY  0.5*0.5 degree
"""
_,_,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon,lat,ocean_mask_192_288)	
AROPH_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(AROPH_file,mode='r')
lat_AROPH = nc_fid.variables['lat'][:]
lon_AROPH = nc_fid.variables['lon'][:]
time_AROPH = nc_fid.variables['time'][20:50]
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][20:50],axis=0);
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][20:50],axis=0)
SDII_AROPH = stats.nanmean(nc_fid.variables['sdii'][20:50],axis=0); SDII_AROPH[np.isnan(SDII_AROPH) ]=0.0
r95p_AROPH = stats.nanmean(nc_fid.variables['r95p'][20:50],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0

# CDD_AROPH = stats.nanmean(nc_fid.variables['cdd'][20:50],axis=0)
# R10_AROPH = stats.nanmean(nc_fid.variables['r10'][20:50],axis=0)

f1 = interp2d(lon_AROPH, lat_AROPH, RX5DAY_AROPH,kind='quintic')
RX5DAY_AROPH_interp = f1(lons, lats)
f2 = interp2d(lon_AROPH, lat_AROPH, SDII_AROPH,kind='quintic')
SDII_AROPH_interp = f2(lons, lats)
f3 = interp2d(lon_AROPH, lat_AROPH, CWD_AROPH,kind='quintic')
CWD_AROPH_interp = f3(lons, lats)
f4 = interp2d(lon_AROPH, lat_AROPH, r95p_AROPH,kind='quintic')
r95p_AROPH_interp = f4(lons, lats)

RX5DAY_AROPH_interp_masked = np.multiply(RX5DAY_AROPH_interp,ocean_mask_192_288s)
cwd_AROPH_interp_masked = np.multiply(CWD_AROPH_interp,ocean_mask_192_288s)
sdii_AROPH_interp_masked = np.multiply(SDII_AROPH_interp,ocean_mask_192_288s)
r95p_AROPH_interp_masked = np.multiply(r95p_AROPH_interp,ocean_mask_192_288s)


fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(7, 7.5), facecolor='White');plot_setup();
pad= 5
#############################
#RX5DAY
#############################
ax = plt.subplot(4,2,2);
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
dist_space = linspace( -20, 200, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,210,40.0));ax.set_yticks(np.arange(0.0,1.8,0.4))
ax.set_xlim([0,200]);ax.set_ylim([0,1.8])
ax.annotate('(e)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('East Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
#############################
# SDII
#############################
size =np.shape(sdii_AROPH_interp_masked)
hist_AROPH = np.reshape(sdii_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(sdii_masked_clipped)
hist_CESM = np.reshape(sdii_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
	
ax = plt.subplot(4,2,4)						
dist_space = linspace( 0, 19, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,21,4.0));ax.set_yticks(np.arange(0.0,25.0,5.0))
ax.set_xlim([0,20]);ax.set_ylim([0,23])
ax.annotate('(f)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
#############################
# CWD
#############################
ax = plt.subplot(4,2,6);
# Histogram embedded in the subplot
size =np.shape(cwd_AROPH_interp_masked)
hist_AROPH = np.reshape(cwd_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 60, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0.0,41.0,10.0));ax.set_yticks(np.arange(0.0,18.0,3.0))
ax.set_xlim([0,40]);ax.set_ylim([0,16])
ax.annotate('(g)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
#############################
# r95p
#############################
ax = plt.subplot(4,2,8);
# Histogram embedded in the subplot
size =np.shape(r95p_AROPH_interp_masked)
hist_AROPH = np.reshape(r95p_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 250, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0.0,260,50));ax.set_yticks(np.arange(0.0,1.2,0.2))
ax.set_xlim([0,250]);ax.set_ylim([0,1.2])
ax.annotate('(h)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)


#############################################################
#############################################################
region = rergion_dic['SA'][:]

# ocean land mask
import scipy.io as spio
oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720')['landocean']
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
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/10/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][51:81]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][51:81,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][51:81,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][51:81,:,:],axis=0);
sdii = stats.nanmean(nc_fid.variables['sdii'][51:81,:,:],axis=0);

# interpolate the ocean_mask to required resolution
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440, kind='linear')
ocean_mask_192_288 = f(lon, lat)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1
RX5DAY_3D_masked = np.multiply(RX5DAY,ocean_mask_192_288)
lons,lats,RX5DAY_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,RX5DAY_3D_masked)
sdii_3D_masked = np.multiply(sdii,ocean_mask_192_288)
lons,lats,sdii_masked_clipped= range_clip(region[0],region[1],region[2],region[3],lon,lat,sdii_3D_masked)
cwd_3D_masked = np.multiply(cwd,ocean_mask_192_288)
lons,lats,cwd_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cwd_3D_masked)
r95p_3D_masked = np.multiply(r95p,ocean_mask_192_288)
lons,lats,r95p_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r95p_3D_masked)

"""
# AROPH RX5DAY  0.5*0.5 degree
"""
_,_,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon,lat,ocean_mask_192_288)	
AROPH_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(AROPH_file,mode='r')
lat_AROPH = nc_fid.variables['lat'][:]
lon_AROPH = nc_fid.variables['lon'][:]
time_AROPH = nc_fid.variables['time'][20:50]
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][20:50],axis=0);
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][20:50],axis=0)
SDII_AROPH = stats.nanmean(nc_fid.variables['sdii'][20:50],axis=0); SDII_AROPH[np.isnan(SDII_AROPH) ]=0.0
r95p_AROPH = stats.nanmean(nc_fid.variables['r95p'][20:50],axis=0); r95p_AROPH[np.isnan(r95p_AROPH) ]=0.0

# CDD_AROPH = stats.nanmean(nc_fid.variables['cdd'][20:50],axis=0)
# R10_AROPH = stats.nanmean(nc_fid.variables['r10'][20:50],axis=0)

f1 = interp2d(lon_AROPH, lat_AROPH, RX5DAY_AROPH,kind='quintic')
RX5DAY_AROPH_interp = f1(lons, lats)
f2 = interp2d(lon_AROPH, lat_AROPH, SDII_AROPH,kind='quintic')
SDII_AROPH_interp = f2(lons, lats)
f3 = interp2d(lon_AROPH, lat_AROPH, CWD_AROPH,kind='quintic')
CWD_AROPH_interp = f3(lons, lats)
f4 = interp2d(lon_AROPH, lat_AROPH, r95p_AROPH,kind='quintic')
r95p_AROPH_interp = f4(lons, lats)

RX5DAY_AROPH_interp_masked = np.multiply(RX5DAY_AROPH_interp,ocean_mask_192_288s)
cwd_AROPH_interp_masked = np.multiply(CWD_AROPH_interp,ocean_mask_192_288s)
sdii_AROPH_interp_masked = np.multiply(SDII_AROPH_interp,ocean_mask_192_288s)
r95p_AROPH_interp_masked = np.multiply(r95p_AROPH_interp,ocean_mask_192_288s)

#############################
#RX5DAY
#############################
ax = plt.subplot(4,2,1);
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
dist_space = linspace( -20, 200, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,250,50.0));ax.set_yticks(np.arange(0.0,1.0,0.3))
ax.set_xlim([0,200]);ax.set_ylim([0,0.9])

ax.annotate('(a)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('South Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('RX5DAY (mm)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
legend=plt.legend(loc=1,shadow=False);
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)

#############################
# SDII
#############################
size =np.shape(sdii_AROPH_interp_masked)
hist_AROPH = np.reshape(sdii_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(sdii_masked_clipped)
hist_CESM = np.reshape(sdii_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
	
ax = plt.subplot(4,2,3)						
dist_space = linspace( 0, 19, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0,21,4.0));ax.set_yticks(np.arange(0.0,22.0,5.0))
ax.set_xlim([0,20]);ax.set_ylim([0,16])
ax.annotate(r'SDII ($\mathrm{\mathsf{mm\/day^{-1}}}$)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(b)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
#############################
# CWD
#############################
ax = plt.subplot(4,2,5);
# Histogram embedded in the subplot
size =np.shape(cwd_AROPH_interp_masked)
hist_AROPH = np.reshape(cwd_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 60, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0.0,61.0,12.0));ax.set_yticks(np.arange(0.0,4.0,1.0))
ax.set_xlim([0,61]);ax.set_ylim([0,4])
ax.annotate('CWD (days)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
#############################
# r95p
#############################
ax = plt.subplot(4,2,7);
# Histogram embedded in the subplot
size =np.shape(r95p_AROPH_interp_masked)
hist_AROPH = np.reshape(r95p_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]
					
dist_space = linspace( 0, 250, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax.set_xticks(np.arange(0.0,300,50));ax.set_yticks(np.arange(0.0,1.0,0.2))
ax.set_xlim([0,250]);ax.set_ylim([0,0.75])
ax.annotate('(d)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('R95P',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')

plt.subplots_adjust(left=0.10, bottom=0.06, right=0.98, top=0.95, wspace=0.12, hspace=0.22); 
plt.savefig('A1_PDFs.png', format='png', dpi=800)
plt.show()
