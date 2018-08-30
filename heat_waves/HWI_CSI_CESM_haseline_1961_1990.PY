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
## data readin 
ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;ocean_mask_CESM[0:27,:]=np.nan

def zeros_lookups(data):
    # Create an array that is 1 where data is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(data, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def get_all_baseline_data(scenario,scale):
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
			osci= nc_fid.variables[scenario][0,0:30,:,:]
		elif scale== 'interseasonal':
			osci= nc_fid.variables[scenario][1:13,0:30,:,:]
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
	return TX95P,TN95P,TX5P,TN5P,osci

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/his/'
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()


scenario = 'his';scale = ''
TX95P,TN95P,TX5P,TN5P,osci = get_all_baseline_data(scenario,scale)
 
ensembles =23
HWIM= np.empty((ensembles,192,288));HWIM[:]=np.nan;
HWIS= np.empty((ensembles,192,288));HWIS[:]=np.nan;
CSIM= np.empty((ensembles,192,288));CSIM[:]=np.nan;
CSIS= np.empty((ensembles,192,288));CSIS[:]=np.nan;

for en_no in range(0,ensembles): 
	print en_no
	nc_f = text_content[en_no][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	year = nc_fid.variables['year'][41:71]
	TX = nc_fid.variables['TX'][41:71,:,:,:];
	TN = nc_fid.variables['TN'][41:71,:,:,:];
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
	for lat_index in range(192):
		for lon_index in range(288):
			HW_events = np.array([]);
			CS_events = np.array([]);
			if (~np.isnan(ocean_mask_CESM[lat_index,lon_index])):
				for iyear in range(len(year)):
				 	####HEAT WAVE
					HW_index = np.ones(365);
					TX_cache = TX_HW[iyear,:,lat_index,lon_index]; TN_cache = TN_HW[iyear,:,lat_index,lon_index];
					TXTN = (TX_cache+TN_cache)/2
					# TXTN = (TX[iyear,:,lat_index,lon_index]+TN[iyear,:,lat_index,lon_index])/2
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];HW_index[tag] = 0;				
					ranges = zeros_lookups(HW_index);
					for ino in range(np.shape(ranges)[0]):
						if (ranges[ino,1]-ranges[ino,0] < 3):
							HW_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN)- np.sum(HW_index)==0):
						HW_events=np.append(HW_events,np.nan);
					else:
						ranges = zeros_lookups(HW_index);
						for ino in range(np.shape(ranges)[0]):
							HW_events=np.append(HW_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
					####COLD SPELL
					CS_index = np.ones(365);
					TX_cache = TX_CS[iyear,:,lat_index,lon_index]; TN_cache = TN_CS[iyear,:,lat_index,lon_index];
					TXTN = (TX_cache+TN_cache)/2
					# TXTN = (TX[iyear,:,lat_index,lon_index]+TN[iyear,:,lat_index,lon_index])/2
					tag = [item for item in range(len(TXTN)) if (TX_cache[item]>0 and TN_cache[item]>0)];CS_index[tag] = 0;				
					ranges = zeros_lookups(CS_index);
					for ino in range(np.shape(ranges)[0]):
						if (ranges[ino,1]-ranges[ino,0] < 3):
							CS_index[ranges[ino,0]:ranges[ino,1]] = 1
					if (len(TXTN)- np.sum(CS_index)==0):
						CS_events=np.append(CS_events,np.nan)
					else:
						ranges = zeros_lookups(CS_index);
						for ino in range(np.shape(ranges)[0]):
							CS_events=np.append(CS_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
				HWIM[en_no,lat_index,lon_index] = stats.nanmean(HW_events);
				HWIS[en_no,lat_index,lon_index] = stats.nanstd(HW_events);
				CSIM[en_no,lat_index,lon_index] = stats.nanmean(CS_events);
				CSIS[en_no,lat_index,lon_index] = stats.nanstd(CS_events);
	nc_fid.close()
	

threshold_data ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/TmTn_percentile_calender_enMean.nc'
nc_fid = nc4.Dataset(threshold_data,mode='r')
lat = nc_fid.variables['lat']
lon = nc_fid.variables['lon']

file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/HWI_CSI_1961_1990_0595_CalDayThr.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4')
f.createDimension('en', ensembles)
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

ens = f.createVariable('en',np.float32, ('en'))
lats = f.createVariable('lat',np.float32, ('lat'))
lons = f.createVariable('lon',np.float32, ('lon'))
HWIMs = f.createVariable('HWIM',np.float32,('en','lat','lon'))
HWISs = f.createVariable('HWIS',np.float32,('en','lat','lon'))
CSIMs = f.createVariable('CSIM',np.float32,('en','lat','lon'))
CSISs = f.createVariable('CSIS',np.float32,('en','lat','lon'))

ens[:] = range(ensembles); 
lats[:] = lat; 
lons[:] = lon;
HWIMs[:] = HWIM;
HWISs[:] = HWIS;
CSIMs[:] = CSIM; 
CSISs[:] = CSIS; 

f.description = 'Heat wave index'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
HWIMs.standard_name = 'heat wave index mean'
HWISs.standard_name = 'heat wave index std'
CSIMs.standard_name = 'cold spell index mean'
CSISs.standard_name = 'cold spell index std'

nc_fid.close()
f.close()
