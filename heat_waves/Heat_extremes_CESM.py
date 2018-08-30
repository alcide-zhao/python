# -*- coding: utf-8 -*-
'''
This is to produce the heat wave related indices over the interested regions over JJAS
the input
	(a) The reference threshold (from 1961-1990) which are used to define the extreme heat and tropical night
	(b) The daily TX and Tn 
The definitions:
	(a) Extreme heat: Days when the TX > the percentiles of the TX
	(b) Tropical night: Nights when the TN> the percentiles o f the referring TN
	(c) Heat wave: when extreme heat coincide with tropical night
The output:
	The intensity,duration adn frequency of the heat waves
		(a) frequency: the proportion of days of heat waves to the total days(1222days)
		(b) durations: in Three categories based on the normal distribution of the duration 
			namely : < the 25th, between the 25th and 75th, and >75th percentiles
		(c) Intensity: here consider maximum and mean 
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


######################################################
# calculate HW for the year interested
######################################################
def zeros_lookups(data):
    # Create an array that is 1 where data is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(data, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges
	
def mask_match(country_key,lon,lat):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	from scipy.interpolate import interp2d  as interp2d
	import scipy.io as sio
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/countries_mask_lon_lat_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask= ocean_mask[country_key][:]
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
	mask = f(lon, lat); 
	mask[mask > 0] = 1;
	mask[mask < 1] = np.nan
	return mask
	
def EH_TN_HW(TD,TN,TD_ll,TD_ul,TN_ll,TN_ul):
	def percentile_selection (data):
		Imin= np.nanmin(data);
		Imea = stats.nanmean(data);
		Imax= np.nanmax(data)
		return Imin,Imea,Imax
	
	TD_index = np.ones(len(TD));TN_index = np.ones(len(TD));HW_index = np.ones(len(TD));
	if (np.isnan(TD_ul) or np.isnan(TN_ul)):
		tag = [item for item in range(len(TD)) if TD[item]>=TD_ll];TD_index[tag] = 0;
		tag = [item for item in range(len(TN)) if TN[item]>=TN_ll]; TN_index[tag] = 0;
		tag = [item for item in range(len(TD)) if (TD[item]>=TD_ll and TN[item]>=TN_ll)];HW_index [tag] = 0
	else:
		tag = [item for item in range(len(TD)) if (TD[item]>=TD_ll and TD[item]<TD_ul)];TD_index[tag] = 0;
		tag = [item for item in range(len(TN)) if (TN[item]>=TN_ll and TN[item]<TN_ul)]; TN_index[tag] = 0;
		tag = [item for item in range(len(TD)) if ((TD[item]>=TD_ll and TD[item]<TD_ul) and (TN[item]>=TN_ll and TN[item]<TN_ul))];HW_index [tag] = 0
	
	# heat extremes
	count_TD = len(TD)- np.sum(TD_index);TD_intensity = TD -TD_ll;
	ranges = zeros_lookups(TD_index); durations = ranges[:,1]-ranges[:,0];
	if count_TD ==0:
		TD_stats = np.empty((5));TD_stats[:]=np.nan
	else:
		Ievents = np.empty((np.shape(ranges)[0])); Ievents[:]=np.nan
		for ino in range(np.shape(ranges)[0]):
			Ievents[ino] = np.sum(TD_intensity[ranges[ino,0]:ranges[ino,1]])
		_,Imea,Imax = percentile_selection(Ievents)
		_,Dmea,Dmax = percentile_selection (durations)
		TD_stats = np.array([count_TD,Dmea,Dmax,Imea,Imax])
	# Tropical nights
	count_TN = len(TN)- np.sum(TN_index);TN_intensity = TN -TN_ll;
	ranges = zeros_lookups(TN_index); durations = ranges[:,1]-ranges[:,0];
	if count_TN ==0:
		TN_stats = np.empty((5));TN_stats[:]=np.nan
	else:
		Ievents = np.empty((np.shape(ranges)[0])); Ievents[:]=np.nan
		for ino in range(np.shape(ranges)[0]):
			Ievents[ino] = np.sum(TN_intensity[ranges[ino,0]:ranges[ino,1]])
		_,Imea,Imax = percentile_selection(Ievents)
		_,Dmea,Dmax = percentile_selection (durations)
		TN_stats = np.array([count_TN,Dmea,Dmax,Imea,Imax])
	# heat waves
	#1. screen out duration less than 3 days
	ranges = zeros_lookups(HW_index); durations = ranges[:,1]-ranges[:,0];
	for ino in range(np.shape(ranges)[0]):
		if (ranges[ino,1]-ranges[ino,0]  <  3):
				HW_index[ranges[ino,0]:ranges[ino,1]] = 1			
	Tdays_HW = len(TN)- np.sum(HW_index);
	#2. count the duratin, intensity and countsba
	if Tdays_HW ==0:
		HW_stats = np.empty((5));HW_stats[:]=np.nan
	else:
		ranges = zeros_lookups(HW_index); durations = ranges[:,1]-ranges[:,0];
		frequeency = np.shape(ranges)[0]
		_,Dmea,Dmax = percentile_selection (durations)
		Ievents = np.empty((np.shape(ranges)[0])); Ievents[:]=np.nan
		HW_intensity = TD+TN-TN_ll-TD_ll;
		for ino in range(np.shape(ranges)[0]):
			Ievents[ino] = np.sum(HW_intensity[ranges[ino,0]:ranges[ino,1]])
		_,Imea,Imax = percentile_selection(Ievents)
		HW_stats = np.array([frequeency,Tdays_HW,Dmax,Imea,Imax])	
	return TD_stats,TN_stats,HW_stats

def EH_TN_HW_TiSe(file_name,country_key,region_box,TX90,TN90,TX99,TN99):
	nc_fid = nc4.Dataset(file_name,mode='r')
	year = nc_fid.variables['year'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	mask = mask_match(country_key,lon,lat)
	TX = np.multiply(nc_fid.variables['TX'][:],mask)
	TN = np.multiply(nc_fid.variables['TN'][:],mask)
	region = region_dic[region_box][:]
	_,_,TX_region = range_clip(region[0],region[1],region[2],region[3],lon,lat,TX)
	_,_,TN_region = range_clip(region[0],region[1],region[2],region[3],lon,lat,TN)
	TD_TiSe = stats.nanmean(stats.nanmean(TX_region,axis=3),axis=2)
	TN_TiSe = stats.nanmean(stats.nanmean(TN_region,axis=3),axis=2)
	TD_stats = np.empty((len(year),2,5));TN_stats = np.empty((len(year),2,5));HW_stats = np.empty((len(year),2,5));
	year_base =year[0]
	for iyear in year:
		TD = TD_TiSe[iyear-year_base,:];TN = TN_TiSe[iyear-year_base,:];
		TD_stats_moderate,TN_stats_moderate,HW_stats_moderate = EH_TN_HW(TD,TN,TD_ll=TX90,TD_ul=np.nan,TN_ll=TN90,TN_ul=np.nan)
		TD_stats_extreme,TN_stats_extreme,HW_stats_extreme = EH_TN_HW(TD,TN,TD_ll=TX99,TD_ul=np.nan,TN_ll=TN99,TN_ul=np.nan)
		TD_stats[iyear-year_base,:,:]= np.array([TD_stats_moderate,TD_stats_extreme])
		TN_stats[iyear-year_base,:,:]= np.array([TN_stats_moderate,TN_stats_extreme])
		HW_stats[iyear-year_base,:,:]= np.array([HW_stats_moderate,HW_stats_extreme])
	return year, TD_stats,TN_stats,HW_stats

	
##################
####   Main   ####
##################	
	
region_dic = {'NorthIndia':[70,90,22,33],'SouthIndia':[70,90,8,22],'SouthEastChina':[105,125,20,30],'CentralEastChina':[105,125,30,40]}
TnTm_90_99 = {'NorthIndia':[23.28,32.83,23.76,34.40],'SouthIndia':[23.92,30.41,24.25,31.91],'SouthEastChina':[23.21,30.32,23.47,30.76],'CentralEastChina':[21.14,28.95,21.61,29.62]}

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/rcp85/'
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()	
LEN = len(text_content) # NUMBER OF ENSEMBLE WHIHC CAN BE OBTAINED

TD_CEC = np.zeros((LEN,95,2,5));TN_CEC = np.zeros((LEN,95,2,5));HW_CEC = np.zeros((LEN,95,2,5));
TD_SEC = np.zeros((LEN,95,2,5));TN_SEC = np.zeros((LEN,95,2,5));HW_SEC = np.zeros((LEN,95,2,5));
TD_NID = np.zeros((LEN,95,2,5));TN_NID = np.zeros((LEN,95,2,5));HW_NID = np.zeros((LEN,95,2,5));
TD_SID = np.zeros((LEN,95,2,5));TN_SID = np.zeros((LEN,95,2,5));HW_SID = np.zeros((LEN,95,2,5));

for en_no in range(LEN):
	print en_no
	file_name = text_content[en_no][:-1]
	region_box = 'CentralEastChina';country_key='China';
	TN90=TnTm_90_99[region_box][0];TX90=TnTm_90_99[region_box][1];TN99=TnTm_90_99[region_box][2];TX99=TnTm_90_99[region_box][3]
	year_series, TD_CEC[en_no,:,:],TN_CEC[en_no,:,:],HW_CEC[en_no,:,:] = EH_TN_HW_TiSe(file_name,country_key,region_box,TX90,TN90,TX99,TN99)
	region_box = 'SouthEastChina';country_key='China'
	TN90=TnTm_90_99[region_box][0];TX90=TnTm_90_99[region_box][1];TN99=TnTm_90_99[region_box][2];TX99=TnTm_90_99[region_box][3]
	_, TD_SEC[en_no,:,:],TN_SEC[en_no,:,:],HW_SEC[en_no,:,:] = EH_TN_HW_TiSe(file_name,country_key,region_box,TX90,TN90,TX99,TN99)
	region_box = 'NorthIndia';country_key='India'
	TN90=TnTm_90_99[region_box][0];TX90=TnTm_90_99[region_box][1];TN99=TnTm_90_99[region_box][2];TX99=TnTm_90_99[region_box][3]
	_, TD_NID[en_no,:,:],TN_NID[en_no,:,:],HW_NID[en_no,:,:] = EH_TN_HW_TiSe(file_name,country_key,region_box,TX90,TN90,TX99,TN99)
	region_box = 'SouthIndia';country_key='India'
	TN90=TnTm_90_99[region_box][0];TX90=TnTm_90_99[region_box][1];TN99=TnTm_90_99[region_box][2];TX99=TnTm_90_99[region_box][3]
	_, TD_SID[en_no,:,:],TN_SID[en_no,:,:],HW_SID[en_no,:,:] = EH_TN_HW_TiSe(file_name,country_key,region_box,TX90,TN90,TX99,TN99)


output_name = 'Heat_TD_TN_HW_EastChina_India_rcp85.nc'
f = nc4.Dataset(output_name,'w', format='NETCDF4') #'w' stands for write
f.createDimension('ensemble', LEN)
f.createDimension('year', len(year_series))
f.createDimension('category', 2)
f.createDimension('statistics', 5)

ensembles =  f.createVariable('ensemble',np.float64, ('ensemble'))
years = f.createVariable('year',np.float64, ('year'))
categorys = f.createVariable('category',np.float64, ('category'))
statisticss = f.createVariable('statistics',np.float64, ('statistics'))

TD_CECs = f.createVariable('TD_CEC',np.float32,('ensemble','year','category','statistics'))
TN_CECs = f.createVariable('TN_CEC',np.float32,('ensemble','year','category','statistics'))
HW_CECs = f.createVariable('HW_CEC',np.float32,('ensemble','year','category','statistics'))
TD_SECs = f.createVariable('TD_SEC',np.float32,('ensemble','year','category','statistics'))
TN_SECs = f.createVariable('TN_SEC',np.float32,('ensemble','year','category','statistics'))
HW_SECs = f.createVariable('HW_SEC',np.float32,('ensemble','year','category','statistics'))

TD_SIDs = f.createVariable('TD_SID',np.float32,('ensemble','year','category','statistics'))
TN_SIDs = f.createVariable('TN_SID',np.float32,('ensemble','year','category','statistics'))
HW_SIDs = f.createVariable('HW_SID',np.float32,('ensemble','year','category','statistics'))
TD_NIDs = f.createVariable('TD_NID',np.float32,('ensemble','year','category','statistics'))
TN_NIDs = f.createVariable('TN_NID',np.float32,('ensemble','year','category','statistics'))
HW_NIDs = f.createVariable('HW_NID',np.float32,('ensemble','year','category','statistics'))

ensembles = range(en_no)
years[:] = year_series
categorys[:] = np.array([0,1])
statisticss[:] = np.array([0,1,2,3,4])
TD_CECs[:]= TD_CEC;TN_CECs[:]= TN_CEC;HW_CECs[:]= HW_CEC
TD_SECs[:]= TD_SEC;TN_SECs[:]= TN_SEC;HW_SECs[:]= HW_SEC
TD_SIDs[:]= TD_SID;TN_SIDs[:]= TN_SID;HW_SIDs[:]= HW_SID
TD_NIDs[:]= TD_NID;TN_NIDs[:]= TN_NID;HW_NIDs[:]= HW_NID

import time as clock
f.description = 'Time series of Heat events over East China and India'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()
