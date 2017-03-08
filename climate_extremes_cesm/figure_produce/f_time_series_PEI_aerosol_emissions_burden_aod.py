# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to show the time evelution of selected ectreme indices 
				the evelution of aerosol emissions, burden and AOD 
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
prep_att_ref = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_dic = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_list=['rx5day','r99p','sdii','cwd']
prep_unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

# temperature extremes
# temp_att_ref stores the reference value which is used to remove climatology
temp_att_ref = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_dic = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_list = ['txn','su25','dtr','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
temp_unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],'EA':[100,145,20,50],'SA':[65,100,5,30],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}


region = rergion_dic['EA'][:]

#######################################################
# 0.2 functions and subroutines                       #
#######################################################
def time_seeries_of_spatial_mean(time_s, time_e, time, data):
	"""
	for the input array, values of the interesed geographical domain is to be averaged
	the time seris of the averaged value is then produced 
	"""
	time_series = np.empty((np.shape(data)[0]))
	
	layer_s = [layer for layer in range(len(time)) if time[layer] == time_s]
	layer_e = [layer for layer in range(len(time)) if time[layer] == time_e]
	layer_s = layer_s[0]
	layer_e = layer_e[0]
	for layer in range(layer_s , layer_e+1):
		# print data[layer,:,:]
		time_series[layer-layer_s] = stats.nanmean(stats.nanmean(data[layer,:,:]))
	return time_series

def movingaverage (values, window=3):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result
	
'''
#######################################################
# 1. time evolution                                   #
#######################################################

spst = [2,2]   # subplot style
"""
Time series of temperature extreme indecies 
"""

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/'
plt.figure(1, facecolor='White',figsize=[120,30])

file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100_Rcp8.5.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

sub_plot=1
for att_name in temp_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	temp_att_ref[att_name] = time_series[0]
	time_series = time_series -temp_att_ref[att_name]
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	# print np.shape(movingaverage(time_series,5))
	ax1.plot(time,movingaverage(time_series),'-',color="red",linewidth=2,label='RCP8.5')
	sub_plot = sub_plot+1
nc_fid.close()


file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100_fixA.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	
sub_plot=1	
for att_name in temp_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	# temp_att_ref[att_name] = time_series[0]
	time_series = time_series -temp_att_ref[att_name]
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	ax1.plot(time,movingaverage(time_series),'-',color="blue",linewidth=2,label='RCP8.5_fixA')
	ax1.set_xlabel('Year',fontsize=10)
	ax1.set_xlim([2000,2100])
	# ax1.set_ylim([4,7.5])
	ax1.set_ylabel(temp_unit_dic[att_name].title(),fontsize=15)
	# ax1.set(aspect=10)
	ax1.set_title(att_name.title(),fontsize=15)
	ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
	sub_plot = sub_plot+1
nc_fid.close()

"""
Time series of precipitation ectreme indecies 
"""
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'


file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
spst = [2,4]
sub_plot=1
for att_name in prep_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	prep_att_ref[att_name] = time_series[0]
	time_series = time_series -prep_att_ref[att_name]
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	ax1.plot(time,movingaverage(time_series),'-',color="red",linewidth=2,label='RCP8.5')
	sub_plot = sub_plot+1
nc_fid.close()


file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

sub_plot=1	
for att_name in prep_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	# prep_att_ref[att_name] = time_series[0]
	time_series = time_series -prep_att_ref[att_name]
	
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	ax1.plot(time,movingaverage(time_series),'-',color="blue",linewidth=2,label='RCP8.5_fixA')
	ax1.set_xlabel('Year',fontsize=10)
	ax1.set_xlim([2000,2100])
	# ax1.set_ylim([4,7.5])
	ax1.set_ylabel(prep_unit_dic[att_name].title(),fontsize=15)
	# ax1.set(aspect=10)
	ax1.set_title(att_name.title(),fontsize=15)
	ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
	sub_plot = sub_plot+1
nc_fid.close()

# legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(0.95, 0.8),fontsize=20)	
# plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/time_series_global.tif')

'''

'''
# Time series of AOD 550nm, aerosol burden and Emissions
'''


################################################
# Aerosol burden and AOD from CESM model outputs
################################################
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
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
	time_series_att_rcp85 = stats.nanmean(stats.nanmean(att_rcp85_clipped,axis=2),axis=1)
	#FixA
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:]
	rcp85_fixA = np.multiply(rcp85_fixA,oceanmask)
	lons,lats,rcp85_fixA_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85_fixA)
	time_series_rcp85_fixA = stats.nanmean(stats.nanmean(rcp85_fixA_clipped,axis=2),axis=1)
	nc_fid.close()
	
	###JJA mean
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	rcp85_jja_time_series = np.empty((95))
	rcp85_jja_time_series[:] =np.nan
	fixA_jja_time_series = np.empty((95))
	fixA_jja_time_series[:] =np.nan
	layer_e = -6  
	# print year_series
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = time_series_att_rcp85[layer_b:layer_e+1]
		rcp85_jja_time_series[iyear-2006] = stats.nanmean(cache,axis=0)
		cache = time_series_rcp85_fixA[layer_b:layer_e+1]
		fixA_jja_time_series[iyear-2006] = stats.nanmean(cache,axis=0)
	return year_series,fixA_jja_time_series,rcp85_jja_time_series
	
year_series,AODVIS_fixA_jja_time_series,AODVIS_rcp85_jja_time_series = time_series_of_aerosol_physics(variable_name='AODVIS')
year_series,BURDENBC_fixA_jja_time_series,BURDENBC_rcp85_jja_time_series = time_series_of_aerosol_physics(variable_name='BURDENBC')
#year_series,BURDENDUST_fixA_jja_time_series,BURDENDUST_rcp85_jja_time_series = time_series_of_aerosol_physics('BURDENDUST')
year_series,BURDENPOM_fixA_jja_time_series,BURDENPOM_rcp85_jja_time_series = time_series_of_aerosol_physics('BURDENPOM')
#year_series,BURDENSEASALT_fixA_jja_time_series,BURDENSEASALT_rcp85_jja_time_series = time_series_of_aerosol_physics('BURDENSEASALT')
year_series,BURDENSO4_fixA_jja_time_series,BURDENSO4_rcp85_jja_time_series = time_series_of_aerosol_physics('BURDENSO4')
year_series,BURDENSOA_fixA_jja_time_series,BURDENSOA_rcp85_jja_time_series = time_series_of_aerosol_physics('BURDENSOA')

# burden_sum_fixA = BURDENBC_fixA_jja_time_series+BURDENDUST_fixA_jja_time_series+BURDENSEASALT_fixA_jja_time_series+\
# BURDENPOM_fixA_jja_time_series+BURDENSO4_fixA_jja_time_series+BURDENSOA_fixA_jja_time_series

# burden_sum_rcp85 = BURDENBC_rcp85_jja_time_series+BURDENDUST_rcp85_jja_time_series+BURDENSEASALT_rcp85_jja_time_series+\
# BURDENPOM_rcp85_jja_time_series+BURDENSO4_rcp85_jja_time_series+BURDENSOA_rcp85_jja_time_series

burden_sum_fixA = BURDENBC_fixA_jja_time_series+\
BURDENPOM_fixA_jja_time_series+BURDENSO4_fixA_jja_time_series+BURDENSOA_fixA_jja_time_series

burden_sum_rcp85 = BURDENBC_rcp85_jja_time_series+\
BURDENPOM_rcp85_jja_time_series+BURDENSO4_rcp85_jja_time_series+BURDENSOA_rcp85_jja_time_series

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
	value_clipped_TiSe = stats.nanmean(stats.nanmean(value_clipped,axis=2),axis=1)
		###JJA mean
	
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
#################################
# plotting the figures
################################
# emissions_sum
fig, ax1 = plt.subplots(1, facecolor='White',figsize=[5,5])
lins1 = ax1.plot(year_series,10**11*movingaverage(emission_sum[1:]-emission_sum[0]),'-',color="r",linewidth=5,label='Total')
lins2 = ax1.plot(year_series,10**11*movingaverage(emission_TiSe_dic['OC'][1:]-emission_TiSe_dic['OC'][0]),'-',color="g",linewidth=5,label='OC')
lins3 = ax1.plot(year_series,10**11*movingaverage(emission_TiSe_dic['SO2'][1:]-emission_TiSe_dic['SO2'][0]),'-',color="y",linewidth=5,label='SO2')
lins4 = ax1.plot(year_series,10**12*movingaverage(emission_TiSe_dic['BC'][1:]-emission_TiSe_dic['BC'][0]),'-',color="b",linewidth=5,label='BC*10^-1')
ax1.set_ylabel('Emissions (10^-11*kg/M2/$)', color='k')
ax1.set_xlabel('Year',fontsize=20);ax1.set_xlim([2000,2100]);
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

legend1 = ax1.legend(loc=10, shadow=True, bbox_to_anchor=(0.1, 0.1),fontsize=20)	

ax2 = ax1.twinx()
lins5 = ax2.plot(year_series,10**5*movingaverage(burden_sum_rcp85-burden_sum_fixA),'-.',color="r",linewidth=5,label='B_Total')
lins6 = ax2.plot(year_series,10**5*movingaverage(BURDENSOA_rcp85_jja_time_series+BURDENPOM_rcp85_jja_time_series-\
BURDENSOA_fixA_jja_time_series-BURDENPOM_fixA_jja_time_series),'-.',color="g",linewidth=5,label='B_OC')
lins7 = ax2.plot(year_series,10**5*movingaverage(BURDENSO4_rcp85_jja_time_series-BURDENSO4_fixA_jja_time_series),'-.',color="y",linewidth=5,label='B_SO2')
lins8 = ax2.plot(year_series,10**6*movingaverage(BURDENBC_rcp85_jja_time_series-BURDENBC_fixA_jja_time_series),'-.',color="b",linewidth=5,label='B_BC*10^-1')
lins9 = ax2.plot(range(2000,2101),0*np.ones((101)),'-',color="k",linewidth=3,label='B_BC*10^-1')
ax2.set_ylabel('Burden (10^-5*kg/M2)', color='k')
ax2.set_title('Aerosol Emission and Burden over East Asia (RCP8.5-FixA)',fontsize=20)

align_yaxis(ax1, 0, ax2, 0)

"""
# emissions_sum
sub_plot=5
spst = [2,4]
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_sum[1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_sum[0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('Emissions',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
# aerosol_burden sum
sub_plot=6

ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(burden_sum_fixA)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)

ax1.plot(year_series,movingaverage(burden_sum_rcp85)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)


#aerosol optical depth
sub_plot=7
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(AODVIS_fixA_jja_time_series),'-',color="blue",linewidth=2,label='RCP8.5_fixA')
ax1.plot(year_series,movingaverage(AODVIS_rcp85_jja_time_series),'-',color="red",linewidth=2,label='RCP8.5')
# ax1.hold(True)
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
# ax1.set_ylabel('AOD550'.title(),fontsize=30)
# ax1.set(aspect=10)
ax1.set_title('AOD550',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(1.5, 0.8),fontsize=15)	


spst =[2,4]
plt.figure(2, facecolor='White',figsize=[10,10])

# BC emissions
sub_plot=1
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_TiSe_dic['BC'][1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_TiSe_dic['BC'][0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('BC',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# OC emissions
sub_plot=2
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_TiSe_dic['OC'][1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_TiSe_dic['OC'][0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('OC',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# BC emissions
sub_plot=3
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_TiSe_dic['SO2'][1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_TiSe_dic['SO2'][0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('SO2',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(1.5, 0.8),fontsize=15)	

# BC BURDEN
sub_plot=5
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENBC_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENBC_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('BC Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# POM+SOA BURDEN
sub_plot = 6
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENPOM_fixA_jja_time_series+BURDENSOA_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENPOM_rcp85_jja_time_series+BURDENSOA_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('POA+SOA Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(1.5, 0.9),fontsize=15)	


# SO4 BURDEN
sub_plot = 7
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENSO4_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENSO4_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('SO4 Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# SEASALT  BURDEN
# sub_plot = 8
# ax1=plt.subplot(spst[0],spst[1],sub_plot)
# ax1.plot(year_series,movingaverage(BURDENSEASALT_fixA_jja_time_series+BURDENDUST_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
# ax1.hold(True)
# ax1=plt.subplot(spst[0],spst[1],sub_plot)
# ax1.plot(year_series,movingaverage(BURDENSEASALT_rcp85_jja_time_series+BURDENDUST_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
# ax1.set_xlabel('Year',fontsize=10)
# ax1.set_xlim([2000,2100])
# # ax1.set_ylim([4,7.5])
# ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# # ax1.set(aspect=10)
# ax1.set_title('SEASALT Aerosol Burden',fontsize=15)
# ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
# ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/time_series_global.tif')
"""
plt.show()
