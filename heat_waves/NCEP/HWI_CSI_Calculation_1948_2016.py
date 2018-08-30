# -*- coding: utf-8 -*-
'''
Calculate the HWI and CSI index statistics based on the baseline climatology (1961-1990)
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
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio
## data readin 
ocean_mask_NCEP = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_NCEP_94_192.mat')['landocean']
ocean_mask_NCEP=np.flipud(ocean_mask_NCEP)
ocean_mask_NCEP[ocean_mask_NCEP==0]=np.nan;ocean_mask_NCEP[0:13,:]=np.nan

# TmTn 95percentiles for 1961-1990
threshold_data ='/scratch/local/s1667168/NCEP_TmTn_percentile_365.nc'
nc_fid = nc4.Dataset(threshold_data,mode='r')
lat = nc_fid.variables['lat']
lon = nc_fid.variables['lon']
TX95P = np.multiply(nc_fid.variables['TX95P'],ocean_mask_NCEP)-273.15
TN95P = np.multiply(nc_fid.variables['TN95P'],ocean_mask_NCEP)-273.15
TX5P = np.multiply(nc_fid.variables['TX5P'],ocean_mask_NCEP)-273.15
TN5P = np.multiply(nc_fid.variables['TN5P'],ocean_mask_NCEP)-273.15
nc_fid.close();

# baseline data
baseline_data ='/scratch/local/s1667168/NCEP_HWI_CSI_1961_1990_0595_abs.nc'
nc_fid = nc4.Dataset(baseline_data,mode='r')
HWIM= nc_fid.variables['HWIM'][:]
HWIS= nc_fid.variables['HWIS'][:]
CSIM= nc_fid.variables['CSIM'][:]
CSIS= nc_fid.variables['CSIS'][:]
nc_fid.close()

# functions
def zeros_lookups(data):
    # Create an array that is 1 where data is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(data, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where abs diff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

# nc files for caculation
scenario = 'NCEP'
scenario_dic ={'NCEP':[1,0,69]}
ensemble = scenario_dic[scenario][0]; layer_s = scenario_dic[scenario][1]; layer_e = scenario_dic[scenario][2]



nc_f = '/scratch/local/s1667168/NCEP_TXTN.2m.gauss.1948_2017.nc'
nc_fid = nc4.Dataset(nc_f,mode='r')
year = nc_fid.variables['year'][layer_s:layer_e]
TX = nc_fid.variables['TX'][layer_s:layer_e,:,:,:]-273.15;  TX_HW = TX-TX95P; TX_CS = TX5P-TX; 
TN = nc_fid.variables['TN'][layer_s:layer_e,:,:,:]-273.15;  TN_HW = TN-TN95P; TN_CS = TN5P-TN;

HWI = np.zeros((len(year),94,192,21));HWI[:]=np.nan
CSI = np.zeros((len(year),94,192,21));CSI[:]=np.nan

for lat_index in range(94):
	for lon_index in range(192):
		if (~np.isnan(TX95P[lat_index,lon_index])):
			for iyear in range(len(year)):
				##########HEAT WAVES
				HW_events = np.array([]);RF_events = np.array([]);HW_intensity = np.array([]);event_no =0;HW_index = np.ones(365);
				TX_cache = TX_HW[iyear,:,lat_index,lon_index]; TN_cache = TN_HW[iyear,:,lat_index,lon_index];
				TXTN = (TX_cache+TN_cache)/2
				# TXTN = (TX[iyear,:,lat_index,lon_index]+TN[iyear,:,lat_index,lon_index])/2
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
					HWI_M = stats.nanmean(HW_scaled);HWI_X = np.nanmax(HW_scaled)
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
				# TXTN = (TX[iyear,:,lat_index,lon_index]+TN[iyear,:,lat_index,lon_index])/2
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
file_name = '/scratch/local/s1667168/HWI_CSI_'+scenario+'_abs.nc'
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

f.description = 'Heat wave and Cold spell statistics from NCEP'
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