# -*- coding: utf-8 -*-
"""
This is to show the variation of Aerosol burdens and AOD variations 
in 2036-2045 and 2090-2100 with respct to current day (2006-2015) levels
@author: Alcide.Zhao
"""

import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import os
from scipy import stats
import site
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
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


time_s = 2006
time_e = 2100
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}

region = rergion_dic['ASIA'][:]



def Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name):
	'''
	This function is to process thw diagnostic field into JJA mean from the CESM
	monthly data
	The input is the file which contains both rcp85 and fixA simulations
	The output is should be normaly lons, lats, time_series and range clipped data
	the lev information willl also be returned for 4D field
	'''
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	rcp85 = nc_fid.variables['rcp85'][:] 
	RCP85_MV = nc_fid.variables['rcp85'].missing_value
	rcp85[rcp85 == RCP85_MV] = np.nan
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:] 
	RCP85_fixA_MV = nc_fid.variables['rcp85_fixA'].missing_value
	rcp85_fixA[rcp85_fixA == RCP85_fixA_MV] = np.nan
	units =  nc_fid.variables['rcp85_fixA'].units
	# long_name =  nc_fid.variables['rcp85_fixA'].long_name
	long_name= 'Surface Temperature (radiative)'
	diff = rcp85-rcp85_fixA
	lons,lats,diff_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,diff)
	
	size = np.shape(diff_clipped)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	###JJA mean
	if (np.rank(diff_clipped) == 4):
		levs = nc_fid.variables['lev'][:]
		variable_JJAmean = np.empty((95,30,size[2],size[3])); variable_JJAmean[:] =np.nan
		
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = diff_clipped[layer_b:layer_e+1,:,:,:]
			variable_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
	else:
		levs=np.nan
		variable_JJAmean = np.empty((95,size[1],size[2])); variable_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = diff_clipped[layer_b:layer_e+1,:,:]
			variable_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	nc_fid.close()
	return lons, lats, levs, year_series, variable_JJAmean,units,long_name

#BC	
variable = 'BURDENBC'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, BC_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)


# oc
variable = 'BURDENSOA'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, SOA_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

variable = 'BURDENPOM'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, POM_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
OC_JJAmean = SOA_JJAmean + POM_JJAmean

variable = 'BURDENSO4'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, SO4_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

variable = 'AODVIS'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, AOD_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)


BC_0615_D = stats.nanmean(BC_JJAmean[0:10,:,:],axis=0)*10**5
OC_0615_D = stats.nanmean(OC_JJAmean[0:10,:,:],axis=0)*10**5
SO4_0615_D = stats.nanmean(SO4_JJAmean[0:10,:,:],axis=0)*10**5
AOD_0615_D = stats.nanmean(AOD_JJAmean[0:10,:,:],axis=0)
BC_3645_D = stats.nanmean(BC_JJAmean[30:40,:,:],axis=0)*10**5
OC_3645_D = stats.nanmean(OC_JJAmean[30:40,:,:],axis=0)*10**5
SO4_3645_D = stats.nanmean(SO4_JJAmean[30:40,:,:],axis=0)*10**5
AOD_3645_D = stats.nanmean(AOD_JJAmean[30:40,:,:],axis=0)
BC_9100_D = stats.nanmean(BC_JJAmean[85:95,:,:],axis=0)*10**5
OC_9100_D = stats.nanmean(OC_JJAmean[85:95,:,:],axis=0)*10**5
SO4_9100_D = stats.nanmean(SO4_JJAmean[85:95,:,:],axis=0)*10**5
AOD_9100_D = stats.nanmean(AOD_JJAmean[85:95,:,:],axis=0)


# now plotting
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(6, 8), facecolor='White')
fig.tight_layout()
# tight_layout doesn't take these labels into account. We'll need 
# to make some room. These numbers are are manually tweaked. 
# You could automatically calculate them, but it's a pain.
fig.subplots_adjust(left=0.15, top=0.9)


colormap ='PRGn'; p_value = np.zeros((np.shape(BC_0615_D))); 
pad= 5 
# BC Burden
colorbar_min=-0.30; colorbar_max = 0.05;
ax = plt.subplot(4,3,1)
spatial_figure_norm(BC_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('2006_2015',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax.annotate('BC',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(4,3,2)
spatial_figure_norm(BC_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('2036_2045',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(4,3,3)
colormesh1=spatial_figure_norm(BC_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('2090_2010',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
cbar_ax = fig.add_axes([0.93, 0.72, 0.01, 0.18])
char = fig.colorbar(colormesh1, cax=cbar_ax,extend='both')
char.set_label('10^-5*kg/m2',fontsize=15)



# OC Burden		
colorbar_min=-0.7; colorbar_max = 0.40;		
ax=plt.subplot(4,3,4)
spatial_figure_norm(OC_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('OC',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(4,3,5)
spatial_figure_norm(OC_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
plt.subplot(4,3,6)
colormesh2 =spatial_figure_norm(OC_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
cbar_ax = fig.add_axes([0.93, 0.50, 0.01, 0.18])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both')
char.set_label('10^-5*kg/m2',fontsize=15)

# SO4
colorbar_min=-3.0; colorbar_max = 0.7;	
ax=plt.subplot(4,3,7)
spatial_figure_norm(SO4_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('SO4',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(4,3,8)
spatial_figure_norm(SO4_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
plt.subplot(4,3,9)
colormesh3 =spatial_figure_norm(SO4_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
cbar_ax = fig.add_axes([0.93, 0.27, 0.01, 0.18])
char = fig.colorbar(colormesh3, cax=cbar_ax,extend='both')
char.set_label('10^-5*kg/m2',fontsize=15)

# AOD550
colormap= 'RdBu'
colorbar_min=-0.26; colorbar_max = 0.07;	
ax=plt.subplot(4,3,10)
spatial_figure_norm(AOD_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('AOD',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(4,3,11)
spatial_figure_norm(AOD_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
plt.subplot(4,3,12)
colormesh4=spatial_figure_norm(AOD_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)		
cbar_ax = fig.add_axes([0.93, 0.05, 0.01, 0.18])
char = fig.colorbar(colormesh4, cax=cbar_ax,extend='both')
char.set_label('',fontsize=20)

plt.show()
