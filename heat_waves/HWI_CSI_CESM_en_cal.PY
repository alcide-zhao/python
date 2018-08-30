# -*- coding: utf-8 -*-
'''
This scipt is to compute the baseline heat wave intensity (HWI) which is defined as the accumulated temperature of (TX+Tn-TX90P-TN90P)
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
scenario = 'rcp85';scale = 'interannual';directory='CalDayThr_InterOsi' # _InterOsi

scenario_dic ={'his':[23,41,86],'rcp85':[22,0,95],'fixa':[14,0,95],'rcp45':[15,0,75]}
ensemble = scenario_dic[scenario][0]; layer_s = scenario_dic[scenario][1]; layer_e = scenario_dic[scenario][2]
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/'+scenario+'/'
# os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()

##get all data input
TX95P,TN95P,TX5P,TN5P,HWIM,HWIS,CSIM,CSIS,osci = get_all_baseline_data(scenario,scale,directory) 
for en_no in range(3,ensemble): 
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
	
	HWI = np.zeros((len(year),192,288,21));HWI[:]=np.nan
	CSI = np.zeros((len(year),192,288,21));CSI[:]=np.nan
	for lat_index in range(192):
		for lon_index in range(288):
			if (~np.isnan(TX95P[0,lat_index,lon_index])):
				for iyear in range(len(year)):
					##########HEAT WAVES
					HW_events = np.array([]);RF_events = np.array([]);HW_intensity = np.array([]);event_no =0;HW_index = np.ones(365);
					TX_cache = TX_HW[iyear,:,lat_index,lon_index]; TN_cache = TN_HW[iyear,:,lat_index,lon_index];
					TXTN = (TX_cache+TN_cache)/2
					## TX-TN statstics
					rf = TX[iyear,:,lat_index,lon_index]-TN[iyear,:,lat_index,lon_index]
					RFSM = stats.nanmean(rf[151:273],axis=0);
					rf_DJFM =  np.concatenate((rf[0:120],rf[334:-1]),axis=0);RFWT=stats.nanmean(rf_DJFM)
					## criteria for HW
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];HW_index[tag] = 0;
					HWDC = 365-np.sum(HW_index)  #Total days of heatwaves
					ranges = zeros_lookups(HW_index);
					for ino in range(np.shape(ranges)[0]): 
						# exclude events which last only two days and mark them with 1
						if (ranges[ino,1]-ranges[ino,0] < 3): 
							HW_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN)- np.sum(HW_index)==0):  # no cays meet the conditions
						HWI_M=np.nan;HWI_N=np.nan;HWI_X=np.nan; HWI_NO=np.nan;
						HWCM=np.nan;HWCN=np.nan;HWCD=np.nan;HWCS=np.nan; HWCE=np.nan;
						HWCV=np.nan;HWCP=np.nan;HWCU=np.nan;RF_M=np.nan; RF_N=np.nan; RF_X=np.nan;
						Du_X=np.nan; Du_M=np.nan;In_X=np.nan; In_M=np.nan;
					else:
						ranges = zeros_lookups(HW_index);
						duraton= ranges[:,1]-ranges[:,0]; Du_X=np.nanmax(duraton); Du_M=stats.nanmean(duraton);
						for ino in range(np.shape(ranges)[0]):
							event_no=event_no+1
							HW_events=np.append(HW_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
							RF_events=np.append(RF_events,stats.nanmean(rf[ranges[ino,0]:ranges[ino,1]]))
							HW_intensity=np.append(HW_intensity,stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
						In_X=np.nanmax(HW_intensity); In_M=stats.nanmean(HW_intensity);	
						HWI_NO = event_no
						## HW magnitude
						HW_scaled=np.divide(HW_events-HWIM[lat_index,lon_index],HWIS[lat_index,lon_index])
						HWI_M = stats.nanmean(HW_scaled);HWI_X = np.nanmax(HW_scaled);HWI_N = np.nanmin(HW_scaled);
						RF_M = stats.nanmean(RF_events);RF_N =np.nanmin(RF_events);RF_X =np.nanmax(RF_events)
						# tag = [item for item in range(len(HW_scaled)) if ( HW_scaled[item]<0)];HW_scaled[tag]=0;
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]<=-1)];HWCM= len(tag);
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]>-1 and HW_scaled[item]<=0)];HWCN= len(tag);
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]>0 and HW_scaled[item]<=2)];HWCD= len(tag); 
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]>2 and HW_scaled[item]<=4)];HWCS= len(tag);
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]>4 and HW_scaled[item]<=8)];HWCE= len(tag);
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]>8 and HW_scaled[item]<=16)];HWCV= len(tag);
						tag = [item for item in range(len(HW_scaled)) if (HW_scaled[item]>16 and HW_scaled[item]<=32)];HWCP= len(tag);
						tag = [item for item in range(len(HW_scaled)) if HW_scaled[item]>32];HWCU= len(tag);
					HWI[iyear,lat_index,lon_index,:] = np.array([HWI_M,HWI_N,HWI_X,HWI_NO,HWCM,HWCN,HWCD,HWCS,HWCE,HWCV,HWCP,HWCU,RF_M,RF_N,RF_X,HWDC,RFSM, Du_X,Du_M,In_X,In_M])					
					##########COLD SPELLS
					CS_events = np.array([]);RF_events = np.array([]);CS_intensity = np.array([]);event_no =0;CS_index = np.ones(365);
					TX_cache = TX_CS[iyear,:,lat_index,lon_index]; TN_cache = TN_CS[iyear,:,lat_index,lon_index];
					TXTN = (TX_cache+TN_cache)/2
					recovery_factor = TX[iyear,:,lat_index,lon_index]-TN[iyear,:,lat_index,lon_index]
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];CS_index[tag] = 0;
					CSDC = 365- np.sum(CS_index)
					ranges = zeros_lookups(CS_index);
					for ino in range(np.shape(ranges)[0]):
						if (ranges[ino,1]-ranges[ino,0] < 3):
							CS_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN) - np.sum(CS_index)==0):
						CSI_M=np.nan;CSI_N=np.nan;CSI_X=np.nan;CSI_NO=np.nan;CSCM=np.nan;CSCN=np.nan;
						CSCD=np.nan;CSCS=np.nan;CSCE=np.nan;CSCV=np.nan;CSCP=np.nan;CSCU=np.nan;RF_M=np.nan;
						RF_N=np.nan;RF_X=np.nan;Du_X=np.nan;Du_M=np.nan;In_X=np.nan;In_M=np.nan;
					else:
						ranges = zeros_lookups(CS_index);
						duraton= ranges[:,1]-ranges[:,0]; Du_X=np.nanmax(duraton); Du_M=stats.nanmean(duraton);
						for ino in range(np.shape(ranges)[0]):
							event_no=event_no+1
							CS_events=np.append(CS_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
							RF_events=np.append(RF_events,stats.nanmean(recovery_factor[ranges[ino,0]:ranges[ino,1]]))
							CS_intensity=np.append(CS_intensity,-1*stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
						CSI_NO = event_no;In_X=np.nanmax(CS_intensity); In_M=stats.nanmean(CS_intensity);					
						## CS Magnitude
						CS_scaled=np.divide(CS_events-CSIM[lat_index,lon_index],CSIS[lat_index,lon_index])
						CSI_M = -1*stats.nanmean(CS_scaled);CSI_N = -1*np.nanmin(CS_scaled);CSI_X = -1*np.nanmax(CS_scaled);
						RF_M = stats.nanmean(RF_events);RF_N =np.nanmin(RF_events);RF_X =np.nanmax(RF_events)
						# tag = [item for item in range(len(CS_scaled)) if ( CS_scaled[item]<0)];CS_scaled[tag]=0;
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]<=-2)];CSCM= len(tag);
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]>-2 and CS_scaled[item]<=-1)];CSCN= len(tag);
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]>-1 and CS_scaled[item]<=0)];CSCD= len(tag); 
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]>0 and CS_scaled[item]<=2)];CSCS= len(tag);
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]>2 and CS_scaled[item]<=4)];CSCE= len(tag);
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]>4 and CS_scaled[item]<=8)];CSCV= len(tag);
						tag = [item for item in range(len(CS_scaled)) if (CS_scaled[item]>8 and CS_scaled[item]<=16)];CSCP= len(tag);
						tag = [item for item in range(len(CS_scaled)) if CS_scaled[item]>16];CSCU= len(tag);
					CSI[iyear,lat_index,lon_index,:] = np.array([CSI_M,CSI_N,CSI_X,CSI_NO,CSCM,CSCN,CSCD,CSCS,CSCE,CSCV,CSCP,CSCU,RF_M,RF_N,RF_X,CSDC,RFWT,Du_X,Du_M,In_X,In_M])

	year = nc_fid.variables['year'][layer_s:layer_e]
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	
	# writting each ensemble results into a .nc file
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'+directory+'/'+scenario+'/'+'HWI_CSI'+str(['_0'+str(en_no+1) if en_no<9 else '_'+str(en_no+1)])[2:5]+'_'+scenario+'.nc'
	f = nc4.Dataset(file_name,'w', format='NETCDF4')
	f.createDimension('time', len(year))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	f.createDimension('stat', 21)
	
	times = f.createVariable('time',np.float32, ('time'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	statss = f.createVariable('stat',np.float32, ('stat'))
	HWIs = f.createVariable('HWI',np.float32,('time','lat','lon','stat'))
	CSIs = f.createVariable('CSI',np.float32,('time','lat','lon','stat'))
	times[:] = year
	latitudes[:] = lat
	longitudes[:] = lon
	statss[:]=range(21)
	HWIs[:] = HWI
	CSIs[:] = CSI
	
	f.description = 'Heat wave indexs'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	HWIs.standard_name = 'Statistics for HWI'
	HWIs.long_name = '0-annual mean; 1- annual minimum; 2-annual maximum; 3-annual counts;\
	4-no.of mild(<=1); 5-no.of nornal(1-2); 6-no.of moderate(2-4); 7-no. of severe(4-8);\
	8-no.of V. extreme(8-16);  9-no. of V. extreme(16-32);  10-no. of S. extreme(32-64);  11-no. of U. extreme(>64);\
	12-RF mean; 13- RF minimum; 14-RF maximum 15-day count of HW, 16-Summer mean RF;\
	17- maximum duration; 18- mean duration 19- maximum Intensity,20-mean Intensity'
	CSIs.standard_name = 'Statistics for CSI'
	CSIs.long_name = '0-annual mean; 1- annual minimum; 2-annual maximum; 3-annual counts;\
	4-no.of mild(<=1); 5-no.of nornal(1-2); 6-no.of moderate(2-4); 7-no. of severe(4-8);\
	8-no.of V. extreme(8-16);  9-no. of V. extreme(16-32);  10-no. of S. extreme(32-64);  11-no. of U. extreme(>64);\
	12-RF mean; 13- RF minimum; 14-RF maximum 15-day count of CS, 16-winter mean RF;\
	17- maximum duration; 18- mean duration 19- maximum Intensity,20-mean Intensity'
	nc_fid.close()
	f.close()
