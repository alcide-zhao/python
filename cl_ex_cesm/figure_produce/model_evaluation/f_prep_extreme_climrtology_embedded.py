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
prep_att_list=['RX5DAY','R99P','sdii','cwd']
prep_unit_dic = {'rx1day':'Mm', 'RX5DAY':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'R99P':'Mm', 'R99P':'Mm', 'precptot':'Mm'}

rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}


region = rergion_dic['ASIA'][:]

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
# cdd = nc_fid.variables['cdd'][51:81,:,:]
# r10 = nc_fid.variables['r10'][51:81,:,:]

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
# f = interp2d(lon_AROPH, lat_AROPH, CDD_AROPH,kind='quintic')
# CDD_AROPH_interp = f(lons, lats)
# f = interp2d(lon_AROPH, lat_AROPH, R10_AROPH,kind='quintic')
# R10_AROPH_interp = f(lons, lats)

RX5DAY_AROPH_interp_masked = np.multiply(RX5DAY_AROPH_interp,ocean_mask_192_288s)
cwd_AROPH_interp_masked = np.multiply(CWD_AROPH_interp,ocean_mask_192_288s)
sdii_AROPH_interp_masked = np.multiply(SDII_AROPH_interp,ocean_mask_192_288s)
r95p_AROPH_interp_masked = np.multiply(r95p_AROPH_interp,ocean_mask_192_288s)
# CDD_AROPH_interp_masked = np.multiply(CDD_AROPH_interp,ocean_mask_192_288s)
# R10_AROPH_interp_masked = np.multiply(R10_AROPH_interp,ocean_mask_192_288s)

# SDII_AROPH_interp_masked = SDII_AROPH
# r95p_AROPH_interp_masked = r95p_AROPH
# RX5DAY_AROPH_interp_masked = RX5DAY_AROPH
# CWD_AROPH_interp_masked = CWD_AROPH
# CDD_AROPH_interp_masked = CDD_AROPH
# R10_AROPH_interp_masked = R10_AROPH


# climatology mean


fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(12.5, 5.5), facecolor='White');plot_setup();

p_value = np.zeros((np.shape(RX5DAY_masked_clipped)))
#############################
#RX5DAY
#############################
ax = plt.subplot(2,4,1);
colormap ='GnBu';pad= 5;cb_min = 0.000000000001;cb_max =150

spatial_figure(ax,RX5DAY_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate('RX5DAY',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('APHRODITE',xy=(-0.2, 0.6), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='vertical')
ax = plt.subplot(2,4,5)
colormesh1 =spatial_figure(ax,RX5DAY_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=True)
ax.annotate('CESM ensembles',xy=(-0.2, 0.7), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='vertical')
cbar_ax = fig.add_axes([0.1, 0.08, 0.18, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.25,0.25)*(cb_max-cb_min)+cb_min)
char.set_label('mm')

# Histogram embedded in the subplot
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

ax_his = inset_axes(ax, 
                    width="37%", # width = 30% of parent_bbox
                    height='26%', # height : 27% of parent_bbox
                    loc=3)		
					
dist_space = linspace( -20, 200, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =0.05,which= 'major')
ax_his.set_xticks(np.arange(0,200,50.0));ax_his.set_yticks(np.arange(0.0,1.2,0.4))
ax_his.set_xlim([-10,200]);ax_his.set_ylim([0,1.4])
ax_his.xaxis.tick_top();ax_his.yaxis.tick_right()

legend=plt.legend(loc=1,fontsize=5,shadow=False);
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)

#############################
# r95p
#############################
ax = plt.subplot(2,4,4);
colormap ='YlOrRd';pad= 5;cb_min = 0;cb_max =200
spatial_figure(ax,r95p_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('R95P',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax = plt.subplot(2,4,8);
colormesh2 =spatial_figure(ax,r95p_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
cbar_ax = fig.add_axes([0.78, 0.08, 0.18, 0.02])
char = fig.colorbar(colormesh2,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.25,0.25)*(cb_max-cb_min)+cb_min)

# Histogram embedded in the subplot
size =np.shape(r95p_AROPH_interp_masked)
hist_AROPH = np.reshape(r95p_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]

ax_his = inset_axes(ax, 
                    width="37%", # width = 30% of parent_bbox
                    height='26%', # height : 27% of parent_bbox
                    loc=3)		
					
dist_space = linspace( 0, 250, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =0.05,which= 'major')
ax_his.set_xticks(np.arange(0.0,250,120));ax_his.set_yticks(np.arange(0.0,1.0,0.5))
ax_his.set_xlim([-5,280]);ax_his.set_ylim([0,1.0])
ax_his.xaxis.tick_top();ax_his.yaxis.tick_right()

legend=plt.legend(loc=1,fontsize=5,shadow=False);
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)
#############################
# CWD
#############################
ax = plt.subplot(2,4,3);
colormap ='YlOrBr';pad= 5;cb_min = 0;cb_max =40
spatial_figure(ax,cwd_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('CWD',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax = plt.subplot(2,4,7);
colormesh3 =spatial_figure(ax,cwd_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
cbar_ax = fig.add_axes([0.55, 0.08, 0.18, 0.02])
char = fig.colorbar(colormesh3,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.25,0.25)*(cb_max-cb_min)+cb_min)
char.set_label('day')

# Histogram embedded in the subplot
size =np.shape(cwd_AROPH_interp_masked)
hist_AROPH = np.reshape(cwd_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]

ax_his = inset_axes(ax, 
                    width="37%", # width = 30% of parent_bbox
                    height='26%', # height : 27% of parent_bbox
                    loc=3)		
					
dist_space = linspace( 0, 60, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =0.05,which= 'major')
ax_his.set_xticks(np.arange(0.0,80.0,30.0));ax_his.set_yticks(np.arange(0.0,9.0,4.0))
ax_his.set_xlim([-1,60]);ax_his.set_ylim([0,9.0])
ax_his.xaxis.tick_top();ax_his.yaxis.tick_right()

legend=plt.legend(loc=1,fontsize=5,shadow=False);
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)
#############################
# SDII
#############################
ax = plt.subplot(2,4,2);
colormap ='YlGn';pad= 5;cb_min = 0;cb_max =12
spatial_figure(ax,sdii_AROPH_interp_masked,lons,lats,colormap,cb_min,cb_max,p_value, tb_lef=False,tb_bot=False)
ax.annotate('SDII',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax = plt.subplot(2,4,6);
colormesh4 =spatial_figure(ax,sdii_masked_clipped,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
cbar_ax = fig.add_axes([0.325, 0.08, 0.18, 0.02])
char = fig.colorbar(colormesh4,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(0,1.25,0.25)*(cb_max-cb_min)+cb_min)
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$')

# Histogram embedded in the subplot
size =np.shape(sdii_AROPH_interp_masked)
hist_AROPH = np.reshape(sdii_AROPH_interp_masked,(size[0]*size[1],1))
hist_AROPH [hist_AROPH <=0]= np.nan;mask = ~np.isnan(hist_AROPH); 
hist_AROPH = hist_AROPH[mask]

size =np.shape(sdii_masked_clipped)
hist_CESM = np.reshape(sdii_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]

ax_his = inset_axes(ax, 
                    width="37%", # width = 30% of parent_bbox
                    height='26%', # height : 27% of parent_bbox
                    loc=3)		
					
dist_space = linspace( 0, 19, 30)
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, 100*kde(dist_space),'r-',alpha=1, label = 'CESM' )
kde = gaussian_kde(hist_AROPH)
plt.plot( dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'APHRODITE' )
plt.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =0.05,which= 'major')
ax_his.set_xticks(np.arange(0.0,19,5.0));ax_his.set_yticks(np.arange(0.0,24.0,8.0))
ax_his.set_xlim([1,19]);ax_his.set_ylim([0,24])
ax_his.xaxis.tick_top();ax_his.yaxis.tick_right()
legend=plt.legend(loc=1,fontsize=5,shadow=False);
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)

plt.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.95, wspace=0.02, hspace=0.002); 
plt.savefig('Fig3_Precipitation_extremes_CESM_observation.png', format='png', dpi=800)
plt.show()
