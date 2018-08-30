# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to show the correlation between extremes and aerosol burden and aerosol emissions 
@author: Alcide.Zhao
"""

import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
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
    # 'lib'
)
site.addsitedir(lib_path)
from lib import *

#######################################
# 0.0 data input
#######################################
oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
#######################################
# 0.1 varibale definition
#######################################
time_s = 2006
time_e = 2100

# precipitation extremes
prep_att_ref = {'cwd':[],'r99p':[],'sdii':[],'rx5day':[]}
prep_att_dic = {'cwd':[],'r99p':[],'sdii':[],'rx5day':[]}
prep_att_list=['cwd','r99p','sdii','rx5day']
prep_unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

# temperature extremes
# temp_att_ref stores the reference value which is used to remove climatology
temp_att_ref = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_dic = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_list = ['txn','su25','dtr','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
temp_unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.],\
'ASIA':[0,360,0,90]}

region = rergion_dic['GLOBE'][:]
att_name = 'rx5day'

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'

file_name = 'Spatial_ensumble_mean_Prep_extremes_global_2006_2100_Rcp8.5.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
time = nc_fid.variables['time'][:]
lon = nc_fid.variables['lon'][:]
lat = nc_fid.variables['lat'][:]
att_value = nc_fid.variables[att_name][:]
att_value = np.multiply(att_value,oceanmask)
lons,lats,rcp85_value_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)

file_name = 'Spatial_ensumble_mean_Prep_extremes_global_2006_2100_fixA.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
time = nc_fid.variables['time'][:]
lon = nc_fid.variables['lon'][:]
lat = nc_fid.variables['lat'][:]
att_value = nc_fid.variables[att_name][:]
att_value = np.multiply(att_value,oceanmask)
lons,lats,fixa_value_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)

exreme_variation = rcp85_value_clipped - fixa_value_clipped
exreme_variation = stats.nanmean(stats.nanmean(exreme_variation,axis = 2),axis=1)


'''
Time series of AOD at 550nm
'''
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/cloud_physics/'
def time_series_of_aerosol_physics(variable_name):
	file_name = 'ensumble_mean_'+variable_name+'_200602_210101.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]

	#RCP8.5
	att_rcp85 = nc_fid.variables['rcp85'][:]
	att_rcp85 = np.multiply(att_rcp85,oceanmask)
	lons,lats,att_rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_rcp85)
	# time_series_att_rcp85 = stats.nanmean(stats.nanmean(att_rcp85_clipped,axis=2),axis=1)
	#FixA
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:]
	rcp85_fixA = np.multiply(rcp85_fixA,oceanmask)
	lons,lats,rcp85_fixA_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85_fixA)
	# time_series_rcp85_fixA = stats.nanmean(stats.nanmean(rcp85_fixA_clipped,axis=2),axis=1)
	nc_fid.close()
	
	###JJA mean
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	Aerosol_variation = np.empty((95,192,288))
	Aerosol_variation[:] =np.nan
	layer_e = -6   # In order to let the first year behins from the 152 day 274-213=151
	# print year_series
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = att_rcp85_clipped[layer_b:layer_e+1,:,:]-rcp85_fixA_clipped[layer_b:layer_e+1,:,:]
		# print np.shape(cache)
		Aerosol_variation[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	return stats.nanmean(stats.nanmean(Aerosol_variation,axis = 2),axis=1)


AOD550= time_series_of_aerosol_physics(variable_name='AODVIS')
BURDENBC= time_series_of_aerosol_physics(variable_name='BURDENBC')
# print BURDENBC
BURDENDUST = time_series_of_aerosol_physics('BURDENDUST')
BURDENPOM = time_series_of_aerosol_physics('BURDENPOM')
BURDENSEASALT = time_series_of_aerosol_physics('BURDENSEASALT')
BURDENSO4= time_series_of_aerosol_physics('BURDENSO4')
BURDENSOA = time_series_of_aerosol_physics('BURDENSOA')
BURDEN_TOTAL = BURDENBC+BURDENDUST+BURDENPOM+BURDENSEASALT+BURDENSO4+BURDENSOA

################################################
# RCP emissions
################################################

oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask_360_720.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
# emission_3d_dic = {'BC':[],'OC':[],'SO4':[]}
emission_TiSe_dic = {'BC':np.zeros((96)),'OC':np.zeros((96)),'SO2':np.zeros((96))}
input_path = '/exports/csce/datastore/geos/users/s1667168/RCP/'
for key in emission_TiSe_dic.keys():
	file_name ='accmip_interpolated_emissions_RCP85_'+key+'_2005_2100_0.5x0.5.nc'
	file = input_path+file_name
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	ships = nc_fid.variables['ships'][:]
	anthropogenic = nc_fid.variables['anthropogenic'][:]
	biomass_burning = nc_fid.variables['biomass_burning'][:]
	value = np.multiply(ships+anthropogenic+biomass_burning,oceanmask)
	_,_,value_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,value)
	value_clipped_TiSe = stats.nanmean(stats.nanmean(value_clipped,axis=2),axis=1) # time series
	
	#JJA mean
	layer_e = -5   # RCP hae the record for 2005 and begins from Jan
	# print year_series
	# yeaar_se_cache = np.empty((96))
	for iyear in range(2005,2101):
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = value_clipped_TiSe[layer_b:layer_e+1]
		# print stats.nanmean(cache,axis=0)
		emission_TiSe_dic[key][iyear-2005] = stats.nanmean(cache,axis=0)

emission_sum =emission_TiSe_dic['BC']+emission_TiSe_dic['OC']+emission_TiSe_dic['SO2']



Variable_dic ={'EMISSION_BC':emission_TiSe_dic['BC'][1:]-emission_TiSe_dic['BC'][0],\
'EMISSION_OC':emission_TiSe_dic['OC'][1:]-emission_TiSe_dic['OC'][0],\
'EMISSION_SO2':emission_TiSe_dic['SO2'][1:]-emission_TiSe_dic['SO2'][0],\
'BURDEN_BC':BURDENBC,'BURDEN_OC':BURDENSOA+BURDENPOM,'BURDEN_SO4':BURDENSO4,'BURDENSE_OTHERS':BURDENSEASALT+BURDENDUST,\
'EMISSION_TOTAL':emission_sum[1:]-emission_sum[0],'BURDEN_TOTAL':BURDEN_TOTAL,'AOD550':AOD550}

spst = [3,4]
fig_loc_dic ={'EMISSION_BC':1,'EMISSION_OC':2,'EMISSION_SO2':3,'BURDEN_BC':5,'BURDEN_OC':6,'BURDEN_SO4':7,'BURDENSE_OTHERS':8,\
'EMISSION_TOTAL':9,'BURDEN_TOTAL':10,'AOD550':11}
window =0
for key in Variable_dic.keys():
	# print Variable_dic[Aerosol_name]
	ax=plt.subplot(spst[0],spst[1],fig_loc_dic[key])
	# print np.shape(movingaverage(time_series,5))
	ax.scatter(exreme_variation, Variable_dic[key]*10**12)
	ax.set_xlabel(att_name.title(),fontsize=15)
	# ax1.set_ylim([4,7.5])
	ax.set_ylabel('*10^-12',fontsize=15)
	# ax1.set(aspect=10)
	ax.set_title(key,fontsize=15)
	window = window+1
plt.show()
	
# plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/ectreme_aerosol_correlation_scatterplot.png')

