# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This is to plot the precipitation extreme climatologies
"""
import scipy

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

#######################################
# DATA INPUT and preprocessing 1996-2005
#######################################

"""
# CESM ensumble mean AND rx5day 0.9375*1.25
"""
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][51:81]
rx5day = nc_fid.variables['rx5day'][51:81,:,:]
sdii = nc_fid.variables['sdii'][51:81,:,:]
cwd = nc_fid.variables['cwd'][51:81,:,:]
r99p = nc_fid.variables['r99p'][51:81,:,:]
cdd = nc_fid.variables['cdd'][51:81,:,:]
r10 = nc_fid.variables['r10'][51:81,:,:]

# interpolate the ocean_mask to required resolution
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440, kind='linear')
ocean_mask_192_288 = f(lon, lat)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1


rx5day_3D_masked = np.multiply(rx5day,ocean_mask_192_288)
lons,lats,rx5day_3D_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rx5day_3D_masked)
rx5day_masked_clipped = stats.nanmean(rx5day_3D_masked_clipped,axis = 0)

sdii_3D_masked = np.multiply(sdii,ocean_mask_192_288)
lons,lats,sdii_3D_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,sdii_3D_masked)
sdii_masked_clipped = stats.nanmean(sdii_3D_masked_clipped,axis = 0)

cwd_3D_masked = np.multiply(cwd,ocean_mask_192_288)
lons,lats,cwd_3D_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cwd_3D_masked)
cwd_masked_clipped = stats.nanmean(cwd_3D_masked_clipped,axis = 0)

r99p_3D_masked = np.multiply(r99p,ocean_mask_192_288)
lons,lats,r99p_3D_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r99p_3D_masked)
r99p_masked_clipped = stats.nanmean(r99p_3D_masked_clipped,axis = 0)

r10_3D_masked = np.multiply(r10,ocean_mask_192_288)
lons,lats,r10_3D_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r10_3D_masked)
r10_masked_clipped = stats.nanmean(r10_3D_masked_clipped,axis = 0)

cdd_3D_masked = np.multiply(cdd,ocean_mask_192_288)
lons,lats,cdd_3D_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cdd_3D_masked)
cdd_masked_clipped = stats.nanmean(cdd_3D_masked_clipped,axis = 0)

"""
# AROPH RX5DAY  0.5*0.5 degree
"""
_,_,ocean_mask_192_288s = range_clip(region[0],region[1],region[2],region[3],lon,lat,ocean_mask_192_288)	
AROPH_file = '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/ARHPO_PrepEI_1951_2007_JJA.nc'
nc_fid = nc4.Dataset(AROPH_file,mode='r')
lat_AROPH = nc_fid.variables['lat'][:]
lon_AROPH = nc_fid.variables['lon'][:]
time_AROPH = nc_fid.variables['time'][20:50]
RX5DAY_AROPH = stats.nanmean(nc_fid.variables['rx5day'][20:50],axis=0)
CWD_AROPH = stats.nanmean(nc_fid.variables['cwd'][20:50],axis=0)
SDII_AROPH = stats.nanmean(nc_fid.variables['sdii'][20:50],axis=0)
R99P_AROPH = stats.nanmean(nc_fid.variables['r99p'][20:50],axis=0)
CDD_AROPH = stats.nanmean(nc_fid.variables['cdd'][20:50],axis=0)
R10_AROPH = stats.nanmean(nc_fid.variables['r10'][20:50],axis=0)

f = interp2d(lon_AROPH, lat_AROPH, RX5DAY_AROPH,kind='quintic')
RX5DAY_AROPH_interp = f(lons, lats)
f = interp2d(lon_AROPH, lat_AROPH, SDII_AROPH,kind='quintic')
SDII_AROPH_interp = f(lons, lats)
f = interp2d(lon_AROPH, lat_AROPH, CWD_AROPH,kind='quintic')
CWD_AROPH_interp = f(lons, lats)

f = interp2d(lon_AROPH, lat_AROPH, R99P_AROPH,kind='quintic')
R99P_AROPH_interp = f(lons, lats)
f = interp2d(lon_AROPH, lat_AROPH, CDD_AROPH,kind='quintic')
CDD_AROPH_interp = f(lons, lats)
f = interp2d(lon_AROPH, lat_AROPH, R10_AROPH,kind='quintic')
R10_AROPH_interp = f(lons, lats)

# RX5DAY_AROPH_interp_masked = np.multiply(RX5DAY_AROPH_interp,ocean_mask_192_288s)
# CWD_AROPH_interp_masked = np.multiply(CWD_AROPH_interp,ocean_mask_192_288s)
# SDII_AROPH_interp_masked = np.multiply(SDII_AROPH_interp,ocean_mask_192_288s)
# R99P_AROPH_interp_masked = np.multiply(R99P_AROPH_interp,ocean_mask_192_288s)
# CDD_AROPH_interp_masked = np.multiply(CDD_AROPH_interp,ocean_mask_192_288s)
# R10_AROPH_interp_masked = np.multiply(R10_AROPH_interp,ocean_mask_192_288s)

SDII_AROPH_interp_masked = SDII_AROPH
R99P_AROPH_interp_masked = R99P_AROPH
RX5DAY_AROPH_interp_masked = RX5DAY_AROPH
CWD_AROPH_interp_masked = CWD_AROPH
CDD_AROPH_interp_masked = CDD_AROPH
R10_AROPH_interp_masked = R10_AROPH


# climatology mean

plt.figure(1, facecolor='White',figsize=[18,6])
spst=[2,6]
plt.subplot(spst[0],spst[1],1)
colormap ='gist_rainbow'; figure_title = 'CESM RX5DAY'
spatial_figure(rx5day_masked_clipped,lons,lats,figure_title,colormap,0.1,150,'mm')

plt.subplot(spst[0],spst[1],7)
colormap ='gist_rainbow'; figure_title = 'AROPH RX5DAY'
spatial_figure(RX5DAY_AROPH_interp_masked,lon_AROPH,lat_AROPH,figure_title,colormap,0.1,150,'mm')


plt.subplot(spst[0],spst[1],2)
colormap ='gist_rainbow'; figure_title = 'CESM R99P'
spatial_figure(r99p_masked_clipped,lons,lats,figure_title,colormap,0.1,80,'mm')

plt.subplot(spst[0],spst[1],8)
colormap ='gist_rainbow'; figure_title = 'AROPH R99P'
spatial_figure(R99P_AROPH_interp_masked,lon_AROPH,lat_AROPH,figure_title,colormap,0.1,80,'mm')

plt.subplot(spst[0],spst[1],3)
colormap ='gist_rainbow'; figure_title = 'CESM R10'
spatial_figure(r10_masked_clipped,lons,lats,figure_title,colormap,0.1,30,'dayS')

plt.subplot(spst[0],spst[1],9)
colormap ='gist_rainbow'; figure_title = 'AROPH R10'
spatial_figure(R10_AROPH_interp_masked,lon_AROPH,lat_AROPH,figure_title,colormap,0.1,30,'days')

plt.subplot(spst[0],spst[1],4)
colormap ='gist_rainbow'; figure_title = 'CESM SDII'
spatial_figure(sdii_masked_clipped,lons,lats,figure_title,colormap,3,10,'mm/day')

plt.subplot(spst[0],spst[1],10)
colormap ='gist_rainbow'; figure_title = 'AROPH SDII'
spatial_figure(SDII_AROPH_interp_masked,lon_AROPH,lat_AROPH,figure_title,colormap,3,10,'mm/day')

plt.subplot(spst[0],spst[1],5)
colormap ='gist_rainbow'; figure_title = 'CESM CWD'
spatial_figure(cwd_masked_clipped,lons,lats,figure_title,colormap,0.1,45,'days')

plt.subplot(spst[0],spst[1],11)
colormap ='gist_rainbow'; figure_title = 'AROPH CWD'
spatial_figure(CWD_AROPH_interp_masked,lon_AROPH,lat_AROPH,figure_title,colormap,0.1,45,'days')

plt.subplot(spst[0],spst[1],6)
colormap ='gist_rainbow'; figure_title = 'CESM CDD'
spatial_figure(cdd_masked_clipped,lons,lats,figure_title,colormap,0.1,50,'days')

plt.subplot(spst[0],spst[1],12)
colormap ='gist_rainbow'; figure_title = 'AROPH CDD'
spatial_figure(CDD_AROPH_interp_masked,lon_AROPH,lat_AROPH,figure_title,colormap,0.1,50,'days')



plt.show()
