# -*- coding: utf-8 -*-
'''
This scipt is to compute the baseline heat wave intensity (HWI) AND COLD SPELL INDEX for the NCEP data
The baseline is difined as the the mean over 1961-1990
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
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio
## data readin 
ocean_mask_NCEP = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_NCEP_94_192.mat')['landocean']
ocean_mask_NCEP=np.flipud(ocean_mask_NCEP)
ocean_mask_NCEP[ocean_mask_NCEP==0]=np.nan;ocean_mask_NCEP[0:13,:]=np.nan
# plt.imshow(ocean_mask_NCEP);cbar = plt.colorbar();plt.show()


def zeros_lookups(data):
    # Create an array that is 1 where data is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(data, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

# TmTn 95percentiles for 1961-1990
threshold_data ='/scratch/local/s1667168/NCEP_TmTn_percentile_365.nc'
nc_fid = nc4.Dataset(threshold_data,mode='r')
lat = nc_fid.variables['lat']
lon = nc_fid.variables['lon']
TX95P = np.multiply(nc_fid.variables['TX95P'],ocean_mask_NCEP)-273.15
TN95P = np.multiply(nc_fid.variables['TN95P'],ocean_mask_NCEP)-273.15
TX5P = np.multiply(nc_fid.variables['TX5P'],ocean_mask_NCEP)-273.15
TN5P = np.multiply(nc_fid.variables['TN5P'],ocean_mask_NCEP)-273.15
# plt.imshow(TX90P,origin='lower');cbar = plt.colorbar();plt.show()
nc_fid.close();

HWIM= np.empty((94,192));HWIM[:]=np.nan;
HWIS= np.empty((94,192));HWIS[:]=np.nan;
CSIM= np.empty((94,192));CSIM[:]=np.nan;
CSIS= np.empty((94,192));CSIS[:]=np.nan;

year = range(1961,1991)
nc_f='/scratch/local/s1667168/NCEP_TXTN.2m.gauss.1948_2017.nc'
nc_fid = nc4.Dataset(nc_f,mode='r')
TX = nc_fid.variables['TX'][13:43,:,:,:]-273.15;  TX_HW = TX-TX95P; TX_CS = TX5P-TX; 
TN = nc_fid.variables['TN'][13:43,:,:,:]-273.15;  TN_HW = TN-TN95P; TN_CS = TN5P-TN;
nc_fid.close()

for lat_index in range(94):
	for lon_index in range(192):
		if (~np.isnan(TX95P[lat_index,lon_index])):
			HW_events = np.array([]);
			CS_events = np.array([]);
			for iyear in range(len(year)):
				####HEAT WAVE
				HW_index = np.ones(365);
				TX_cache = TX_HW[iyear,:,lat_index,lon_index]; TN_cache = TN_HW[iyear,:,lat_index,lon_index];
				# TXTN = (TX_cache+TN_cache)/2
				TXTN = (TX[iyear,:,lat_index,lon_index]+TN[iyear,:,lat_index,lon_index])/2
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
				# TXTN = (TX_cache+TN_cache)/2
				TXTN = (TX[iyear,:,lat_index,lon_index]+TN[iyear,:,lat_index,lon_index])/2
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
			HWIM[lat_index,lon_index] = stats.nanmean(HW_events);
			HWIS[lat_index,lon_index] = stats.nanstd(HW_events);
			CSIM[lat_index,lon_index] = stats.nanmean(CS_events);
			CSIS[lat_index,lon_index] = stats.nanstd(CS_events);	

threshold_data ='/scratch/local/s1667168/NCEP_TXTN.2m.gauss.1948_2017.nc'
nc_fid = nc4.Dataset(threshold_data,mode='r')
lat = nc_fid.variables['lat']
lon = nc_fid.variables['lon']

file_name = '/scratch/local/s1667168/NCEP_HWI_CSI_1961_1990_0595_abs.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4')
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
lats = f.createVariable('lat',np.float32, ('lat'))
lons = f.createVariable('lon',np.float32, ('lon'))
HWIMs = f.createVariable('HWIM',np.float32,('lat','lon'))
HWISs = f.createVariable('HWIS',np.float32,('lat','lon'))	
CSIMs = f.createVariable('CSIM',np.float32,('lat','lon'))
CSISs = f.createVariable('CSIS',np.float32,('lat','lon'))	

lats[:] = lat; 
lons[:] = lon;
HWIMs[:] = HWIM; 
HWISs[:] = HWIS; 
CSIMs[:] = CSIM; 
CSISs[:] = CSIS; 

f.description = 'heat wave and cold spell statistics from NCEP'
HWIMs.standard_name = 'heat wave index mean'
HWISs.standard_name = 'heat wave index std'
CSIMs.standard_name = 'cold spell index mean'
CSISs.standard_name = 'cold spell index std'

nc_fid.close()
f.close()
