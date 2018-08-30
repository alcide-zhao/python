# -*- coding: utf-8 -*-
'''
This scipt is to compute the baseline heat wave intensity (HWI) which is defined as
 the accumulated temperature of (TX+Tn-TX90P-TN90P)
The baseline is difined as the 30 ensemble mean of the mean of HWI over the 30 yrs 
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import time as clock
import math as math

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio

# functions
def zeros_lookups(data):
    # Create an array that is 1 where data is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(data, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where abs diff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def get_all_baseline_data(scenario,scale,directory):
	"""
	reading all the necessary date including the thresholds, 
	the 1961-1990 baselines and the inter annual/seasonal osculations upon requet
	"""
	def get_oscilations(scale):
		"""
		get the TS osculations
		"""
		if scenario=='his':
			osci_data =data_path+'/TS_oscilation_rhis_1960_2005.nc'
		else:
			osci_data = data_path+'/TS_oscilation_rcp85_fixa_2006_2100.nc'
		nc_fid = nc4.Dataset(osci_data)
		if scale=='interannual':
			osci= nc_fid.variables[scenario][0,:,:,:]
		elif scale== 'interseasonal':
			osci= nc_fid.variables[scenario][1:13,:,:,:]
		return osci
		nc_fid.close()

	## data read ins 
	ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
	ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;ocean_mask_CESM[0:27,:]=np.nan
	data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp'
	# TmTn 95percentiles for 1961-1995
	threshold_data =data_path+'/Temp_pp/TmTn_percentile_calender_enMean.nc'
	nc_fid = nc4.Dataset(threshold_data,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	TX95P = np.multiply(nc_fid.variables['TX95P'],ocean_mask_CESM)
	TN95P = np.multiply(nc_fid.variables['TN95P'],ocean_mask_CESM)
	TX5P = np.multiply(nc_fid.variables['TX5P'],ocean_mask_CESM)
	TN5P = np.multiply(nc_fid.variables['TN5P'],ocean_mask_CESM)
	nc_fid.close()
	
	# baseline data
	baseline_data =data_path+'/Temp_pp/'+directory+'/HWI_CSI_1961_1990_0595_'+directory+'.nc' #_InterOsi
	nc_fid = nc4.Dataset(baseline_data,mode='r')
	HWIM= stats.nanmean(nc_fid.variables['HWIM'],axis=0);
	HWIS= stats.nanmean(nc_fid.variables['HWIS'],axis=0);
	CSIM= stats.nanmean(nc_fid.variables['CSIM'],axis=0);
	CSIS= stats.nanmean(nc_fid.variables['CSIS'],axis=0);
	nc_fid.close()
	if bool(scale):
		osci = get_oscilations(scale)
		osci_data =data_path+'/TS_oscilation_rhis_1960_2005.nc'
		nc_fid = nc4.Dataset(osci_data)
		TS_6190_A= nc_fid.variables['TS_6190_A'][:]
		TX95P=TX95P-TS_6190_A;
		TN95P=TN95P-TS_6190_A;
		TX5P=TX5P-TS_6190_A;
		TN5P=TN5P-TS_6190_A;
	else:
		osci=0
	return TX95P,TN95P,TX5P,TN5P,HWIM,HWIS,CSIM,CSIS,osci

################
#     Main     #
################
scenario = 'fixa';scale = '';directory='CalDayThr' # _InterOsi

scenario_dic ={'his':[23,41,86],'rcp85':[22,0,95],'fixa':[14,0,95],'rcp45':[15,0,75]}
ensemble = scenario_dic[scenario][0]; layer_s = scenario_dic[scenario][1]; layer_e = scenario_dic[scenario][2]
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/'+scenario
# os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()

##get all data input
TX95P,TN95P,TX5P,TN5P,HWIM,HWIS,CSIM,CSIS,osci = get_all_baseline_data(scenario,scale,directory) 
for en_no in range(9,ensemble): 
	print en_no
	nc_f = text_content[en_no][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	year = nc_fid.variables['year'][layer_s:layer_e]
	TX = nc_fid.variables['TX'][layer_s:layer_e,:,:,:];
	TN = nc_fid.variables['TN'][layer_s:layer_e,:,:,:];  
	for iyear in range(len(year)):
		if bool(scale):
			oscilation = osci[iyear,:,:]
		else:
			oscilation = np.zeros((192,288))
		TX[iyear,:,:,:] = TX[iyear,:,:,:]-oscilation;
		TN[iyear,:,:,:] = TN[iyear,:,:,:]-oscilation;
	TX_HW = TX-TX95P; TX_CS = TX5P-TX; 
	TN_HW = TN-TN95P; TN_CS = TN5P-TN;
	print 'interannual oscilation removed'
	
	HWI_MJJASO = np.zeros((len(year),192,288,8));HWI_MJJASO[:]=np.nan
	HWI_NDJFMA = np.zeros((len(year),192,288,8));HWI_NDJFMA[:]=np.nan
	CSI_MJJASO = np.zeros((len(year),192,288,8));CSI_MJJASO[:]=np.nan
	CSI_NDJFMA = np.zeros((len(year),192,288,8));CSI_NDJFMA[:]=np.nan

	for lat_index in range(192):
		for lon_index in range(288):
			if (~np.isnan(TX95P[0,lat_index,lon_index])):
				for iyear in range(len(year)):
					##########HEAT WAVES
					HW_events = np.array([]);HW_intensity = np.array([]);event_no =0;HW_index = np.ones(184);
					TX_cache = TX_HW[iyear,120:304,lat_index,lon_index]; TN_cache = TN_HW[iyear,120:304,lat_index,lon_index];
					TXTN = (TX_cache+TN_cache)/2
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];HW_index[tag] = 0;
					HWDC = 184-np.sum(HW_index)  #Total days of heatwaves
					ranges = zeros_lookups(HW_index);
					for ino in range(np.shape(ranges)[0]): 
						# exclude events which last only two days and mark them with 1
						if (ranges[ino,1]-ranges[ino,0] < 3): 
							HW_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN)- np.sum(HW_index)==0):  # no days meet the conditions
						HWI_M=np.nan;HWI_X=np.nan; HWI_NO=np.nan
						Du_X=np.nan; Du_M=np.nan;In_X=np.nan; In_M=np.nan;
					else:
						ranges = zeros_lookups(HW_index);
						duraton= ranges[:,1]-ranges[:,0]; #print duraton
						Du_X=np.nanmax(duraton); Du_M=stats.nanmean(duraton);
						for ino in range(np.shape(ranges)[0]):
							event_no=event_no+1
							HW_events=np.append(HW_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
							HW_intensity=np.append(HW_intensity,stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
						In_X=np.nanmax(HW_intensity); In_M=stats.nanmean(HW_intensity);	
						HWI_NO = event_no
						## HW magnitude
						HW_scaled=np.divide(HW_events-HWIM[lat_index,lon_index],HWIS[lat_index,lon_index])
						HWI_M = stats.nanmean(HW_scaled);HWI_X = np.nanmax(HW_scaled);
						# tag = [item for item in range(len(HW_scaled)) if ( HW_scaled[item]<0)];HW_scaled[tag]=0;
					HWI_MJJASO[iyear,lat_index,lon_index,:] = np.array([HWI_M,HWI_X,HWI_NO,HWDC,Du_X,Du_M,In_X,In_M])					
					
					
					HW_events = np.array([]);HW_intensity = np.array([]);event_no =0;HW_index = np.ones(181);
					TX_cache =  np.concatenate((TX_HW[iyear,304:365,lat_index,lon_index],TX_HW[iyear,0:120,lat_index,lon_index]),axis=0);
					TN_cache =  np.concatenate((TN_HW[iyear,304:365,lat_index,lon_index],TN_HW[iyear,0:120,lat_index,lon_index]),axis=0);
					TXTN = (TX_cache+TN_cache)/2
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];HW_index[tag] = 0;
					HWDC = 181-np.sum(HW_index)  #Total days of heatwaves
					ranges = zeros_lookups(HW_index);
					for ino in range(np.shape(ranges)[0]): 
						# exclude events which last only two days and mark them with 1
						if (ranges[ino,1]-ranges[ino,0] < 3): 
							HW_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN)- np.sum(HW_index)==0):  # no cays meet the conditions
						HWI_M=np.nan;HWI_X=np.nan; HWI_NO=np.nan
						Du_X=np.nan; Du_M=np.nan;In_X=np.nan; In_M=np.nan;
					else:
						ranges = zeros_lookups(HW_index);
						duraton= ranges[:,1]-ranges[:,0]; Du_X=np.nanmax(duraton); Du_M=stats.nanmean(duraton);
						for ino in range(np.shape(ranges)[0]):
							event_no=event_no+1
							HW_events=np.append(HW_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
							HW_intensity=np.append(HW_intensity,stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
						In_X=np.nanmax(HW_intensity); In_M=stats.nanmean(HW_intensity);	
						HWI_NO = event_no
						## HW magnitude
						HW_scaled=np.divide(HW_events-HWIM[lat_index,lon_index],HWIS[lat_index,lon_index])
						HWI_M = stats.nanmean(HW_scaled);HWI_X = np.nanmax(HW_scaled);
					HWI_NDJFMA[iyear,lat_index,lon_index,:] = np.array([HWI_M,HWI_X,HWI_NO,HWDC,Du_X,Du_M,In_X,In_M])					
					
					##########COLD SPELLS
					CS_events = np.array([]);CS_intensity = np.array([]);event_no =0;CS_index = np.ones(184);
					TX_cache = TX_CS[iyear,120:304,lat_index,lon_index]; TN_cache = TN_CS[iyear,120:304,lat_index,lon_index];
					TXTN = (TX_cache+TN_cache)/2
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];CS_index[tag] = 0;
					CSDC = 184- np.sum(CS_index)
					ranges = zeros_lookups(CS_index);
					for ino in range(np.shape(ranges)[0]):
						if (ranges[ino,1]-ranges[ino,0] < 3):
							CS_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN) - np.sum(CS_index)==0):
						CSI_M=np.nan;CSI_X=np.nan;CSI_NO=np.nan;Du_X=np.nan;Du_M=np.nan;In_X=np.nan;In_M=np.nan;
					else:
						ranges = zeros_lookups(CS_index);
						duraton= ranges[:,1]-ranges[:,0]; Du_X=np.nanmax(duraton); Du_M=stats.nanmean(duraton);
						for ino in range(np.shape(ranges)[0]):
							event_no=event_no+1
							CS_events=np.append(CS_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
							CS_intensity=np.append(CS_intensity,-1*stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
						CSI_NO = event_no;
						In_X=np.nanmax(CS_intensity); In_M=stats.nanmean(CS_intensity);					
						## CS Magnitude
						CS_scaled=np.divide(CS_events-CSIM[lat_index,lon_index],CSIS[lat_index,lon_index])
						CSI_M = -1*stats.nanmean(CS_scaled);CSI_X = -1*np.nanmax(CS_scaled);
					CSI_MJJASO[iyear,lat_index,lon_index,:] = np.array([CSI_M,CSI_X,CSI_NO,CSDC,Du_X,Du_M,In_X,In_M])

					CS_events = np.array([]);CS_intensity = np.array([]);event_no =0;CS_index = np.ones(181);
					TX_cache =  np.concatenate((TX_CS[iyear,304:365,lat_index,lon_index],TX_CS[iyear,0:120,lat_index,lon_index]),axis=0);
					TN_cache =  np.concatenate((TN_CS[iyear,304:365,lat_index,lon_index],TN_CS[iyear,0:120,lat_index,lon_index]),axis=0);
					TXTN = (TX_cache+TN_cache)/2
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];CS_index[tag] = 0;
					CSDC = 181- np.sum(CS_index)
					ranges = zeros_lookups(CS_index);
					for ino in range(np.shape(ranges)[0]):
						if (ranges[ino,1]-ranges[ino,0] < 3):
							CS_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN) - np.sum(CS_index)==0):
						CSI_M=np.nan;CSI_X=np.nan;CSI_NO=np.nan;Du_X=np.nan;Du_M=np.nan;In_X=np.nan;In_M=np.nan;
					else:
						ranges = zeros_lookups(CS_index);
						duraton= ranges[:,1]-ranges[:,0]; Du_X=np.nanmax(duraton); Du_M=stats.nanmean(duraton);
						for ino in range(np.shape(ranges)[0]):
							event_no=event_no+1
							CS_events=np.append(CS_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
							CS_intensity=np.append(CS_intensity,-1*stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
						CSI_NO = event_no;
						In_X=np.nanmax(CS_intensity); In_M=stats.nanmean(CS_intensity);					
						## CS Magnitude
						CS_scaled=np.divide(CS_events-CSIM[lat_index,lon_index],CSIS[lat_index,lon_index])
						CSI_M = -1*stats.nanmean(CS_scaled);CSI_X = -1*np.nanmax(CS_scaled);
					CSI_NDJFMA[iyear,lat_index,lon_index,:] = np.array([CSI_M,CSI_X,CSI_NO,CSDC,Du_X,Du_M,In_X,In_M])
	year = nc_fid.variables['year'][layer_s:layer_e]
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	
	# writting each ensemble results into a .nc file
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr_season/'+scenario+'/'+'HWI_CSI'+str(['_0'+str(en_no+1) if en_no<9 else '_'+str(en_no+1)])[2:5]+'_'+scenario+'.nc'
	f = nc4.Dataset(file_name,'w', format='NETCDF4')
	f.createDimension('time', len(year))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	f.createDimension('stat', 8)
	
	times = f.createVariable('time',np.float32, ('time'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	statss = f.createVariable('stat',np.float32, ('stat'))
	HWI_MJJASOs = f.createVariable('HWI_MJJASO',np.float32,('time','lat','lon','stat'))
	HWI_NDJFMAs = f.createVariable('HWI_NDJFMA',np.float32,('time','lat','lon','stat'))	
	CSI_MJJASOs = f.createVariable('CSI_MJJASO',np.float32,('time','lat','lon','stat'))
	CSI_NDJFMAs = f.createVariable('CSI_NDJFMA',np.float32,('time','lat','lon','stat'))

	times[:] = year
	latitudes[:] = lat
	longitudes[:] = lon

	HWI_MJJASOs[:] = HWI_MJJASO
	HWI_NDJFMAs[:] = HWI_NDJFMA
	CSI_MJJASOs[:] = CSI_MJJASO
	CSI_NDJFMAs[:] = CSI_NDJFMA
	statss[:]=range(8)
	f.close()
