# -*- coding: utf-8 -*-
"""
Created on Wed 07 March 2018
@author @Alcide.Zhao at Crew, KB, UoE
This is to pre-process the daily Temperature and precipitation into extreme induces!

Targets: Six WDGAR emission scenarios, 40yrs equilibrium of
Input: Daily maximum and minimum T at reference height (2m above the surface); 
	   Daily convective and large scale precipitation
output: T and P extremes: TXX,TNN;
						  rx5day, cdd, cwd, cddt, cwdt, r10, r20, sdii
		Mean and Std of T and P 
calculation periods: Northern Hemisphere summer (JJA) and calendar year
caveats: Percentile based upon induces are not taken into account for the reason 
		 that the baseline (1961-1990) transient simulations were not performed
		 
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import matplotlib.pyplot as plt


#########################
## Precp       Extremes # 
#########################
def precip_extreme_indeces(precp_year_data):
	size_data = np.shape(precp_year_data);
	# r10 r20  sdii
	ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
	ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;
	def rnnmm():
		#######################
		# variables initiation
		#######################
		rx5day = np.empty((size_data[1],size_data[2]));rx5day[:] = np.nan
		for row in range(size_data[1]):
			for colum in range(size_data[2]):
				if (~np.isnan(ocean_mask_CESM[row,colum])):
					prec_cache = precp_year_data[:,row,colum]	
					# rx5day
					rx5day[row, colum] = 0.
					for iday in range(2,len(prec_cache)-2):
						rx5day_cache = 0.
						for span in  range(iday-2 , iday+3):  # note here the definition of 5 consecive days 
							if (np.isnan(prec_cache[span])):
								continue
							else:
								rx5day_cache = rx5day_cache+prec_cache[span]
								if (rx5day_cache >= rx5day[row, colum]):
									rx5day[row, colum] = rx5day_cache
		rx5day[rx5day < 0] = np.nan
		return rx5day
	rx5day = rnnmm()
	mean_prep = np.nanmean(precp_year_data, axis = 0)
	std_prep = np.nanstd(precp_year_data, axis = 0)
	rx1day = np.nanmax(precp_year_data, axis = 0)
	return rx5day, mean_prep, std_prep,rx1day


def P_extremes(lon,lat,time,PRECP,interval,year_series):
	# year_series= np.unique([ iyear/10000 for iyear in map(int,time)]);
	rx5day = np.empty((len(year_series),len(lat),len(lon)));	rx5day[:] = np.nan
	mean_prep = np.empty((len(year_series),len(lat),len(lon)));	 mean_prep[:] = np.nan
	std_prep = np.empty((len(year_series),len(lat),len(lon)));	std_prep[:] = np.nan
	rx1day = np.empty((len(year_series),len(lat),len(lon)));	rx1day[:] = np.nan
	if interval == 'JJA':
		for iyear in year_series:	
			layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+601][0]  #June01
			layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+831][0]  #August 31
			PRECP_cache = PRECP[layer_b:layer_e+1,:,:]
			rx5day[iyear-year_series[0],:,:], mean_prep[iyear-year_series[0],:,:], std_prep[iyear-year_series[0],:,:],\
			rx1day[iyear-year_series[0],:,:]= precip_extreme_indeces(PRECP_cache)
		
	else:
		for iyear in year_series:
			## Because the model outputs are not exactly 365 days organized, few days records may be missed for the first or the last year
			if (iyear == year_series[0] and time[0]>= year_series[0] *10000+101):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+101][0]  #June01
			if (iyear == year_series[-1] and time[-1]<= year_series[-1] *10000+1231):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+1231][0]  #August 31
			PRECP_cache = PRECP[layer_b:layer_e+1,:,:]
			rx5day[iyear-year_series[0],:,:], mean_prep[iyear-year_series[0],:,:], std_prep[iyear-year_series[0],:,:],\
			rx1day[iyear-year_series[0],:,:]= precip_extreme_indeces(PRECP_cache)
	return 	rx5day,mean_prep,std_prep,rx1day


def data_readin(FREQ,scenario):
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		start_year =2000
		start =start_year*365
		ith=0	
		for iday in days:
			month_days =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
			calendar_days = np.array([0,31,59,90,120,151,181,212,243,273,304,334,365])
			total_days = int(iday) + start; 
			year = total_days//365; 
			remainder =  total_days%365; 
			if remainder ==0: year=year;month=1;day=1
			else:
				month = 1+[layer for layer in range(len(calendar_days)) if calendar_days[layer]<= remainder and calendar_days[layer+1]>remainder][0]
				day = remainder - calendar_days[month-1]+1
			date_int[ith] = year*10000+month*100+day; #print date_int[ith]
			ith=ith+1
		return date_int.astype(int)

	input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/CamOnly/'
	PC_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.day.PRECC.nc'  #convective precipitation
	PL_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.day.PRECL.nc'  #large-scale precipitation
	ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
	ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;

	##precipitation data
	nc_fid_pc = nc4.Dataset(PC_path,mode='r')
	lat = nc_fid_pc.variables['lat'][:]
	lon = nc_fid_pc.variables['lon'][:]
	days = nc_fid_pc.variables['time'][:]; time = day2datetime(scenario,days);#print time
	PC = np.multiply(nc_fid_pc.variables['PRECC'][:],ocean_mask_CESM)
	nc_fid_pc.close()

	nc_fid_pl = nc4.Dataset(PL_path,mode='r')	
	PL = np.multiply(nc_fid_pl.variables['PRECL'][:],ocean_mask_CESM)
	nc_fid_pl.close()

	PRECP = (PL+PC)*24*60*60*1000  # from m/s to mm/day ; precipitation =convective+large-scale
	PRECP[PRECP<=0] =np.nan 
	return lon,lat,time,PRECP

######
#main#
######
def main_body(scenario):
	FREQ = "day"  #"mon"
	lon,lat,time,PRECP = data_readin (FREQ,scenario)
	for interval in ['annual']: #'JJA',"annual"

		year_series = range(2000,2030)
		rx5day,PRECPM,PRECPSTD,RX1DAY =P_extremes(lon,lat,time,PRECP,interval,year_series)
		#output into netcdf
		output_path = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/CamOnly/analysis/'
		file_name_out = output_path+scenario+'.Precip.extremes.'+interval+'.nc'
		f = nc4.Dataset(file_name_out,'w', format='NETCDF4')
		f.createDimension('year', 30)
		f.createDimension('lat', len(lat))
		f.createDimension('lon', len(lon))

		years = f.createVariable('year',np.float32, ('year'))
		lats = f.createVariable('lat',np.float32, ('lat'))
		lons = f.createVariable('lon',np.float32, ('lon'))

		RX5DAYs= f.createVariable('RX5DAY',np.float32,('year','lat','lon'))
		PRECPMs= f.createVariable('PRECPM',np.float32,('year','lat','lon'))
		PRECPSTDs= f.createVariable('PRECPSTD',np.float32,('year','lat','lon'))
		RX1DAYs= f.createVariable('RX1DAY',np.float32,('year','lat','lon'))

		years[:] = range(30)
		RX5DAYs[:] = rx5day
		PRECPMs[:]=PRECPM
		PRECPSTDs[:] = PRECPSTD
		RX1DAYs[:] = RX1DAY
		lons[:] = lon
		lats[:] = lat

		RX5DAYs.long_name =  'maximum amount of precipitation over 5 days'
		PRECPMs.long_name =  'seasonal mean of daily precp'
		PRECPSTDs.long_name =  'seasonal std of daily precp'
		RX1DAYs.long_name =  'maximum amount of precipitation over 1 day'

		f.description = 'P extreme indices from the CESM-EDGAR experiments'
		f.institution = 'Alcide Zhao at the university of Edinburgh'
		f.close()
		
scenario = "F_EdgTech"  # F_1970 F_Edg70GO F_Edg70Oz F_Edg70T10SOZ F_EdgEne F_EdgRef F_EdgRef70A F_EdgTech
main_body(scenario)








