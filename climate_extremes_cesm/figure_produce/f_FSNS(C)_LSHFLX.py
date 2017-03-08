# -*- coding: utf-8 -*-
"""
This is to show the variation in  shortwave solar flux at surface from both clear-sky and all-sky 
and the surface latent and sensible flux
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
	
	lons,lats,rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85)
	lons,lats,fixA_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85_fixA)	
	size = np.shape(rcp85_clipped)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	###JJA mean
	if (np.rank(rcp85_clipped) == 4):
		levs = nc_fid.variables['lev'][:]
		rcp85_JJAmean = np.empty((95,30,size[2],size[3])); rcp85_JJAmean[:] =np.nan
		fixa_JJAmean = np.empty((95,30,size[2],size[3])); fixa_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = rcp85_clipped[layer_b:layer_e+1,:,:,:]
			rcp85_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
			cache = fixA_clipped[layer_b:layer_e+1,:,:,:]
			fixa_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
	else:
		levs=np.nan
		rcp85_JJAmean = np.empty((95,size[1],size[2])); rcp85_JJAmean[:] =np.nan
		fixa_JJAmean = np.empty((95,size[1],size[2])); fixa_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = rcp85_clipped[layer_b:layer_e+1,:,:]
			rcp85_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
			cache = fixA_clipped[layer_b:layer_e+1,:,:]
			fixa_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	nc_fid.close()
	return lons, lats, levs, year_series,rcp85_JJAmean,fixa_JJAmean,units,long_name

#FSNS
variable = 'FSNS'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, FSNS_rcp85_JJAmean,FSNS_fixa_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

FSNS_0615_D = stats.nanmean((FSNS_rcp85_JJAmean-FSNS_fixa_JJAmean)[0:10,:,:],axis=0)
_, FSNS_0615_pvalue =scipy.stats.ttest_1samp(FSNS_rcp85_JJAmean[0:10,:,:], stats.nanmean(FSNS_fixa_JJAmean[0:10,:,:],axis=0), axis=0)
FSNS_0615_pvalue[FSNS_0615_pvalue>0.05] =np.nan; FSNS_0615_pvalue[FSNS_0615_pvalue<=0.05] =1

FSNS_3645_D = stats.nanmean((FSNS_rcp85_JJAmean-FSNS_fixa_JJAmean)[30:40,:,:],axis=0)
_, FSNS_3645_pvalue =scipy.stats.ttest_1samp(FSNS_rcp85_JJAmean[30:40,:,:], stats.nanmean(FSNS_fixa_JJAmean[30:40,:,:],axis=0), axis=0)
FSNS_3645_pvalue[FSNS_3645_pvalue>0.05] =np.nan; FSNS_3645_pvalue[FSNS_3645_pvalue<=0.05] =1

FSNS_9100_D = stats.nanmean((FSNS_rcp85_JJAmean-FSNS_fixa_JJAmean)[85:95,:,:],axis=0)
_, FSNS_9100_pvalue =scipy.stats.ttest_1samp(FSNS_rcp85_JJAmean[85:95,:,:], stats.nanmean(FSNS_fixa_JJAmean[85:95,:,:],axis=0), axis=0)
FSNS_9100_pvalue[FSNS_9100_pvalue>0.05] =np.nan; FSNS_9100_pvalue[FSNS_9100_pvalue<=0.05] =1

#FSNSC
variable = 'FSNSC'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, FSNSC_rcp85_JJAmean,FSNSC_fixa_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

FSNSC_0615_D = stats.nanmean((FSNSC_rcp85_JJAmean-FSNSC_fixa_JJAmean)[0:10,:,:],axis=0)
_, FSNSC_0615_pvalue =scipy.stats.ttest_1samp(FSNSC_rcp85_JJAmean[0:10,:,:], stats.nanmean(FSNSC_fixa_JJAmean[0:10,:,:],axis=0), axis=0)
FSNSC_0615_pvalue[FSNSC_0615_pvalue>0.05] =np.nan; FSNSC_0615_pvalue[FSNSC_0615_pvalue<=0.05] =1

FSNSC_3645_D = stats.nanmean((FSNSC_rcp85_JJAmean-FSNSC_fixa_JJAmean)[30:40,:,:],axis=0)
_, FSNSC_3645_pvalue =scipy.stats.ttest_1samp(FSNSC_rcp85_JJAmean[30:40,:,:], stats.nanmean(FSNSC_fixa_JJAmean[30:40,:,:],axis=0), axis=0)
FSNSC_3645_pvalue[FSNSC_3645_pvalue>0.05] =np.nan; FSNSC_3645_pvalue[FSNSC_3645_pvalue<=0.05] =1

FSNSC_9100_D = stats.nanmean((FSNSC_rcp85_JJAmean-FSNSC_fixa_JJAmean)[85:95,:,:],axis=0)
_, FSNSC_9100_pvalue =scipy.stats.ttest_1samp(FSNSC_rcp85_JJAmean[85:95,:,:], stats.nanmean(FSNSC_fixa_JJAmean[85:95,:,:],axis=0), axis=0)
FSNSC_9100_pvalue[FSNSC_9100_pvalue>0.05] =np.nan; FSNSC_9100_pvalue[FSNSC_9100_pvalue<=0.05] =1


#LHFLUX
variable = 'LHFLX'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, LHFLUX_rcp85_JJAmean,LHFLUX_fixa_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

LHFLUX_0615_D = stats.nanmean((LHFLUX_rcp85_JJAmean-LHFLUX_fixa_JJAmean)[0:10,:,:],axis=0)
_, LHFLUX_0615_pvalue =scipy.stats.ttest_1samp(LHFLUX_rcp85_JJAmean[0:10,:,:], stats.nanmean(LHFLUX_fixa_JJAmean[0:10,:,:],axis=0), axis=0)
LHFLUX_0615_pvalue[LHFLUX_0615_pvalue>0.05] =np.nan; LHFLUX_0615_pvalue[LHFLUX_0615_pvalue<=0.05] =1

LHFLUX_3645_D = stats.nanmean((LHFLUX_rcp85_JJAmean-LHFLUX_fixa_JJAmean)[30:40,:,:],axis=0)
_, LHFLUX_3645_pvalue =scipy.stats.ttest_1samp(LHFLUX_rcp85_JJAmean[30:40,:,:], stats.nanmean(LHFLUX_fixa_JJAmean[30:40,:,:],axis=0), axis=0)
LHFLUX_3645_pvalue[LHFLUX_3645_pvalue>0.05] =np.nan; LHFLUX_3645_pvalue[LHFLUX_3645_pvalue<=0.05] =1

LHFLUX_9100_D = stats.nanmean((LHFLUX_rcp85_JJAmean-LHFLUX_fixa_JJAmean)[85:95,:,:],axis=0)
_, LHFLUX_9100_pvalue =scipy.stats.ttest_1samp(LHFLUX_rcp85_JJAmean[85:95,:,:], stats.nanmean(LHFLUX_fixa_JJAmean[85:95,:,:],axis=0), axis=0)
LHFLUX_9100_pvalue[LHFLUX_9100_pvalue>0.05] =np.nan; LHFLUX_9100_pvalue[LHFLUX_9100_pvalue<=0.05] =1


#SHFLUX
variable = 'SHFLX'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, SHFLUX_rcp85_JJAmean,SHFLUX_fixa_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

SHFLUX_0615_D = stats.nanmean((SHFLUX_rcp85_JJAmean-SHFLUX_fixa_JJAmean)[0:10,:,:],axis=0)
_, SHFLUX_0615_pvalue =scipy.stats.ttest_1samp(SHFLUX_rcp85_JJAmean[0:10,:,:], stats.nanmean(SHFLUX_fixa_JJAmean[0:10,:,:],axis=0), axis=0)
SHFLUX_0615_pvalue[SHFLUX_0615_pvalue>0.05] =np.nan; SHFLUX_0615_pvalue[SHFLUX_0615_pvalue<=0.05] =1

SHFLUX_3645_D = stats.nanmean((SHFLUX_rcp85_JJAmean-SHFLUX_fixa_JJAmean)[30:40,:,:],axis=0)
_, SHFLUX_3645_pvalue =scipy.stats.ttest_1samp(SHFLUX_rcp85_JJAmean[30:40,:,:], stats.nanmean(SHFLUX_fixa_JJAmean[30:40,:,:],axis=0), axis=0)
SHFLUX_3645_pvalue[SHFLUX_3645_pvalue>0.05] =np.nan; SHFLUX_3645_pvalue[SHFLUX_3645_pvalue<=0.05] =1

SHFLUX_9100_D = stats.nanmean((SHFLUX_rcp85_JJAmean-SHFLUX_fixa_JJAmean)[85:95,:,:],axis=0)
_, SHFLUX_9100_pvalue =scipy.stats.ttest_1samp(SHFLUX_rcp85_JJAmean[85:95,:,:], stats.nanmean(SHFLUX_fixa_JJAmean[85:95,:,:],axis=0), axis=0)
SHFLUX_9100_pvalue[SHFLUX_9100_pvalue>0.05] =np.nan; SHFLUX_9100_pvalue[SHFLUX_9100_pvalue<=0.05] =1


# now plotting
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(6, 10), facecolor='White')
fig.tight_layout()
# tight_layout doesn't take these labels into account. We'll need 
# to make some room. These numbers are are manually tweaked. 
# You could automatically calculate them, but it's a pain.
fig.subplots_adjust(left=0.15, top=0.9)


colormap ='seismic'; 
pad= 5 
#FSNS
colorbar_min=-5; colorbar_max = 25;
ax = plt.subplot(4,3,1)
spatial_figure_norm(FSNS_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNS_0615_pvalue)
ax.annotate('2006_2015',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax.annotate('FSNS',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(4,3,2)
spatial_figure_norm(FSNS_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNS_3645_pvalue)
ax.annotate('2036_2045',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(4,3,3)
colormesh1=spatial_figure_norm(FSNS_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNS_9100_pvalue)
ax.annotate('2090_2010',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
cbar_ax = fig.add_axes([0.93, 0.74, 0.01, 0.18])
char = fig.colorbar(colormesh1, cax=cbar_ax,extend='both')
char.set_label('W/m2',fontsize=15)

#FSNSC	
colorbar_min=-5; colorbar_max = 25;		
ax=plt.subplot(4,3,4)
spatial_figure_norm(FSNSC_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSC_0615_pvalue)
ax.annotate('FSNSC',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(4,3,5)
spatial_figure_norm(FSNSC_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSC_3645_pvalue)
plt.subplot(4,3,6)
colormesh2 =spatial_figure_norm(FSNSC_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSC_9100_pvalue)
cbar_ax = fig.add_axes([0.93, 0.52, 0.01, 0.18])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both')
char.set_label('W/m2',fontsize=15)

colormap= 'seismic'
# SHFLUX
colorbar_min=-2.0; colorbar_max = 12;	
ax=plt.subplot(4,3,7)
spatial_figure_norm(SHFLUX_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,SHFLUX_0615_pvalue)
ax.annotate('SHFLUX',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(4,3,8)
spatial_figure_norm(SHFLUX_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,SHFLUX_3645_pvalue)
plt.subplot(4,3,9)
colormesh3 =spatial_figure_norm(SHFLUX_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,SHFLUX_9100_pvalue)
cbar_ax = fig.add_axes([0.93, 0.29, 0.01, 0.18])
char = fig.colorbar(colormesh3, cax=cbar_ax,extend='both')
char.set_label('W/m2',fontsize=15)

# LHFLUX

colorbar_min=-2; colorbar_max = 12;	
ax=plt.subplot(4,3,10)
spatial_figure_norm(LHFLUX_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,LHFLUX_0615_pvalue)
ax.annotate('LHFLUX',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(4,3,11)
spatial_figure_norm(LHFLUX_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,LHFLUX_3645_pvalue)
plt.subplot(4,3,12)
colormesh4=spatial_figure_norm(LHFLUX_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,LHFLUX_9100_pvalue)		
cbar_ax = fig.add_axes([0.93, 0.07, 0.01, 0.18])
char = fig.colorbar(colormesh4, cax=cbar_ax,extend='both')
char.set_label('W/m2',fontsize=15)

plt.show()
