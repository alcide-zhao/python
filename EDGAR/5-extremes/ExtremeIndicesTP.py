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
## Temperature Extremes # 
#########################

## mean and std of TS and TXXX,TNN
def T_stats(ts_cache,tx_cache,tn_cache):
	txx = np.nanmax(tx_cache, axis=0)
	tnn = np.nanmin(tn_cache, axis=0)
	mean = np.nanmean(ts_cache, axis=0)	
	std = np.nanstd(ts_cache, axis=0)
	return txx,tnn,mean,std

def T_extremes(lon,lat,time,TS,TX,TN,interval,year_series):
	# year_series= np.unique([ iyear/10000 for iyear in map(int,time)]);
	txx=np.empty((len(year_series),len(lat),len(lon)));txx[:]=np.nan
	tnn=np.empty((len(year_series),len(lat),len(lon)));tnn[:]=np.nan
	mean=np.empty((len(year_series),len(lat),len(lon)));mean[:]=np.nan
	std=np.empty((len(year_series),len(lat),len(lon)));std[:]=np.nan
	if interval == 'JJA':
		for iyear in year_series:
			layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+601][0]  #June01
			layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+831][0]  #August 31
			ts_cache = TS[layer_b:layer_e+1,:,:]
			tx_cache = TN[layer_b:layer_e+1,:,:]
			tn_cache = TX[layer_b:layer_e+1,:,:]
			txx[iyear-year_series[0],:,:],tnn[iyear-year_series[0],:,:],mean[iyear-year_series[0],:,:],std[iyear-year_series[0],:,:] =\
				T_stats(ts_cache,tx_cache,tn_cache)
	else:
		for iyear in year_series:
			## Because the model outputs are not exactly 365 days organized, few days records may be missed for the first or the last year
			if (iyear == year_series[0] and time[0] >= year_series[0]*10000+101):
				layer_b = 0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+101][0]  #June01
			if (iyear == year_series[-1] and time[-1]<= year_series[-1] *10000+1231):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+1231][0]  #August 31
			# layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+101][0]
			# layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+1231][0]
			# print iyear,layer_b,layer_e
			ts_cache = TS[layer_b:layer_e+1,:,:]
			tx_cache = TN[layer_b:layer_e+1,:,:]
			tn_cache = TX[layer_b:layer_e+1,:,:]
			# print np.shape(tn_cache)
			txx[iyear-year_series[0],:,:],tnn[iyear-year_series[0],:,:],mean[iyear-year_series[0],:,:],std[iyear-year_series[0],:,:] =\
				T_stats(ts_cache,tx_cache,tn_cache)
	return 	txx,tnn,mean,std

#########################
## Precp       Extremes # 
#########################
def precip_extreme_indeces(precp_year_data):
	size_data = np.shape(precp_year_data);
	# r10 r20  sdii
	def rnnmm(r10t=10, r20t=20):
		#######################
		# variables initiation
		#######################
		r10 = np.empty((size_data[1],size_data[2]));r10[:] = np.nan
		r20 = np.empty((size_data[1],size_data[2]));r20[:] = np.nan	
		sdii = np.empty((size_data[1],size_data[2]));sdii[:] = np.nan
		rx5day = np.empty((size_data[1],size_data[2]));rx5day[:] = np.nan
		r1 = np.empty((size_data[1],size_data[2]));r1[:] = np.nan
		
		for row in range(size_data[1]):
			for colum in range(size_data[2]):
				prec_cache = precp_year_data[:,row,colum]	
				r1[row,colum]=len([value for value in prec_cache if value >= 1])
				r10[row,colum]=len([value for value in prec_cache if value >= r10t])
				r20[row,colum]=len([value for value in prec_cache if value >= r20t])
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
		r1 = r1.astype('float')
		r1[ r1 <= 0] = np.nan
		sdii = np.divide(total_precip, r1)
		return r10, r20, sdii, rx5day
							
	def cdd_cwd():
		#######################
		# variables initiation
		#######################
		cdd = np.empty((size_data[1],size_data[2]));cdd[:] = np.nan
		cwd = np.empty((size_data[1],size_data[2]));cwd[:] = np.nan
		for row in range(size_data[1]):
			for colum in range(size_data[2]):
				prec_cache = precp_year_data[:,row,colum]
				cdd[row, colum] = 0.
				cwd[row, colum] = 0.
				cdd_cache = 0.
				cwd_cache = 0.
				for iday in range(len(prec_cache)):
					if (np.isnan(prec_cache[iday])):
						cdd[row, colum] = cdd_cache
						cwd[row, colum] = cwd_cache
						cdd_cache = 0.
						cwd_cache = 0.
					elif (iday == 0 and prec_cache[iday] < 1):
						cdd_cache = cdd_cache +1.0
						cdd[row, colum] = cdd_cache
					elif (iday == 0 and prec_cache[iday] >= 1):
						cwd_cache = cwd_cache +1.0
						cwd[row, colum] = cwd_cache
					else:
						if ( prec_cache[iday] < 1 and prec_cache[iday-1] < 1):
							cdd_cache=cdd_cache+1
							if (cdd_cache >= cdd[row, colum]):
								cdd[row, colum] = cdd_cache
						elif ( prec_cache[iday] >= 1 and prec_cache[iday-1] >= 1):
							cwd_cache=cwd_cache+1
							if (cwd_cache >= cwd[row, colum]):
								cwd[row, colum]=cwd_cache
						elif ( prec_cache[iday] >= 1 and (prec_cache[iday-1] < 1 or np.isnan(prec_cache[iday-1]))):
							cdd_cache = 0
							cwd_cache = 1
							if (cwd_cache >= cwd[row, colum]):
								cwd[row, colum]=cwd_cache
						elif ( prec_cache[iday] < 1 and (prec_cache[iday-1] >= 1 or np.isnan(prec_cache[iday-1]))):
							cdd_cache = 1
							cwd_cache = 0
							if (cdd_cache >= cdd[row, colum]):
								cdd[row, colum]=cdd_cache
		return cdd,cwd
	total_precip = np.nansum(precp_year_data, axis = 0)
	mean_prep = np.nanmean(precp_year_data, axis = 0)
	std_prep = np.nanstd(precp_year_data, axis = 0)
	rx1day = np.nanmax(precp_year_data, axis = 0)
	r95 = np.nanpercentile(precp_year_data, 95,axis = 0)
	r99 = np.nanpercentile(precp_year_data, 99,axis = 0)
	r10, r20, sdii, rx5day = rnnmm(r10t=10, r20t=20)
	cdd,cwd = cdd_cwd()
	# plt.imshow(cdd);plt.show()
	return r10, r20,sdii, rx5day, cdd, cwd, mean_prep, std_prep,rx1day,r95,r99


def P_extremes(lon,lat,time,PRECP,interval,year_series):
	# year_series= np.unique([ iyear/10000 for iyear in map(int,time)]);
	r10 = np.empty((len(year_series),len(lat),len(lon)));	r10[:] = np.nan
	r20 = np.empty((len(year_series),len(lat),len(lon)));	r20[:] = np.nan
	sdii = np.empty((len(year_series),len(lat),len(lon)));	sdii[:] = np.nan
	rx5day = np.empty((len(year_series),len(lat),len(lon)));	rx5day[:] = np.nan
	cdd = np.empty((len(year_series),len(lat),len(lon)));	cdd[:] = np.nan
	cwd = np.empty((len(year_series),len(lat),len(lon)));	cwd[:] = np.nan
	mean_prep = np.empty((len(year_series),len(lat),len(lon)));	 mean_prep[:] = np.nan
	std_prep = np.empty((len(year_series),len(lat),len(lon)));	std_prep[:] = np.nan
	rx1day = np.empty((len(year_series),len(lat),len(lon)));	rx1day[:] = np.nan
	r95 = np.empty((len(year_series),len(lat),len(lon)));	 r95[:] = np.nan
	r99 = np.empty((len(year_series),len(lat),len(lon)));	r99[:] = np.nan
	if interval == 'JJA':
		for iyear in year_series:	
			layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+601][0]  #June01
			layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+831][0]  #August 31
			PRECP_cache = PRECP[layer_b:layer_e+1,:,:]
			r10[iyear-year_series[0],:,:], r20[iyear-year_series[0],:,:], sdii[iyear-year_series[0],:,:], rx5day[iyear-year_series[0],:,:], \
			cdd[iyear-year_series[0],:,:], cwd[iyear-year_series[0],:,:], mean_prep[iyear-year_series[0],:,:], std_prep[iyear-year_series[0],:,:],\
			rx1day[iyear-year_series[0],:,:], r95[iyear-year_series[0],:,:], r99[iyear-year_series[0],:,:]= precip_extreme_indeces(PRECP_cache)
		
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
			r10[iyear-year_series[0],:,:], r20[iyear-year_series[0],:,:], sdii[iyear-year_series[0],:,:], rx5day[iyear-year_series[0],:,:], \
			cdd[iyear-year_series[0],:,:], cwd[iyear-year_series[0],:,:], mean_prep[iyear-year_series[0],:,:], std_prep[iyear-year_series[0],:,:],\
			rx1day[iyear-year_series[0],:,:], r95[iyear-year_series[0],:,:], r99[iyear-year_series[0],:,:]= precip_extreme_indeces(PRECP_cache)
	return 	r10,r20,sdii,rx5day,cdd,cwd,mean_prep,std_prep,rx1day,r95,r99


def data_readin(exp,scenario):
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		if scenario.startswith('F_'):start_year =2000 
		else:start_year =2010
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
	
	if exp == 'SlowResp':
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
		TS_path_t = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TS.nc'
		TX_path_t = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMX.nc'
		TN_path_t = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMN.nc'
		PC_path_t = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECC.nc'  #convective precipitation
		PL_path_t = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECL.nc'  #large-scale precipitation
		
		##Temperature
		nc_fid_ts = nc4.Dataset(TS_path_t,mode='r')
		lat = nc_fid_ts.variables['lat'][:]
		lon = nc_fid_ts.variables['lon'][:]
		days = nc_fid_ts.variables['time'][:]; time = day2datetime(scenario,days);#print time
		TS_t = nc_fid_ts.variables['TS'][0:10936,:,:]-273.15
		nc_fid_ts.close()
		
		nc_fid_tn = nc4.Dataset(TN_path_t,mode='r')
		TN_t = nc_fid_tn.variables['TREFHTMN'][0:10936,:,:]-273.15
		nc_fid_tn.close()

		nc_fid_tx = nc4.Dataset(TX_path_t,mode='r')	
		TX_t = nc_fid_tx.variables['TREFHTMX'][0:10936,:,:]-273.15
		nc_fid_tx.close()

		##precipitation data
		nc_fid_pc = nc4.Dataset(PC_path_t,mode='r')
		PC_t = nc_fid_pc.variables['PRECC'][0:10936,:,:]
		nc_fid_pc.close()

		nc_fid_pl = nc4.Dataset(PL_path_t,mode='r')	
		PL_t = nc_fid_pl.variables['PRECL'][0:10936,:,:]
		nc_fid_pl.close()

		PRECP_t = (PL_t+PC_t)*24*60*60*1000  # from m/s to mm/day ; precipitation =convective+large-scale
		PRECP_t[PRECP_t<=0] =np.nan
		
		if scenario == 'T1970RCP': scenario ='F_1970'## different naming conventuon here
		else:
			scenario ='F_'+scenario
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/CamOnly/'
		TS_path_f = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TS.nc'
		TX_path_f = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMX.nc'
		TN_path_f = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMN.nc'
		PC_path_f = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECC.nc'  #convective precipitation
		PL_path_f = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECL.nc'  #large-scale precipitation
		##Temperature
		nc_fid_ts = nc4.Dataset(TS_path_f,mode='r')
		lat = nc_fid_ts.variables['lat'][:]
		lon = nc_fid_ts.variables['lon'][:]
		days = nc_fid_ts.variables['time'][:];
		TS_f = nc_fid_ts.variables['TS'][:]-273.15
		nc_fid_ts.close()
		
		nc_fid_tn = nc4.Dataset(TN_path_f,mode='r')
		TN_f = nc_fid_tn.variables['TREFHTMN'][:]-273.15
		nc_fid_tn.close()

		nc_fid_tx = nc4.Dataset(TX_path_f,mode='r')	
		TX_f = nc_fid_tx.variables['TREFHTMX'][:]-273.15
		nc_fid_tx.close()

		##precipitation data
		nc_fid_pc = nc4.Dataset(PC_path_f,mode='r')
		PC_f = nc_fid_pc.variables['PRECC'][:]
		nc_fid_pc.close()

		nc_fid_pl = nc4.Dataset(PL_path_f,mode='r')	
		PL_f = nc_fid_pl.variables['PRECL'][:]
		nc_fid_pl.close()

		PRECP_f = (PL_f+PC_f)*24*60*60*1000  # from m/s to mm/day ; precipitation =convective+large-scale
		PRECP_f[PRECP_f<=0] =np.nan
		TS = TS_t-TS_f
		TN = TN_t-TN_f
		TX = TX_t-TX_f
		PRECP = PRECP_t-PRECP_f	
	else:	
		if exp=='CamOnly':
			input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/CamOnly/'
			TS_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TS.nc'
			TX_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMX.nc'
			TN_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMN.nc'
			PC_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECC.nc'  #convective precipitation
			PL_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECL.nc'  #large-scale precipitation
		elif exp =='FullCp':
			input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
			TS_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TS.nc'
			TX_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMX.nc'
			TN_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.TREFHTMN.nc'
			PC_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECC.nc'  #convective precipitation
			PL_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.PRECL.nc'  #large-scale precipitation
			##Temperature
		nc_fid_ts = nc4.Dataset(TS_path,mode='r')
		lat = nc_fid_ts.variables['lat'][:]
		lon = nc_fid_ts.variables['lon'][:]
		days = nc_fid_ts.variables['time'][:]; time = day2datetime(scenario,days);#print time
		TS = nc_fid_ts.variables['TS'][:]-273.15
		nc_fid_ts.close()
		
		nc_fid_tn = nc4.Dataset(TN_path,mode='r')
		TN = nc_fid_tn.variables['TREFHTMN'][:]-273.15
		nc_fid_tn.close()

		nc_fid_tx = nc4.Dataset(TX_path,mode='r')	
		TX = nc_fid_tx.variables['TREFHTMX'][:]-273.15
		nc_fid_tx.close()

		##precipitation data
		nc_fid_pc = nc4.Dataset(PC_path,mode='r')
		PC = nc_fid_pc.variables['PRECC'][:]
		nc_fid_pc.close()

		nc_fid_pl = nc4.Dataset(PL_path,mode='r')	
		PL = nc_fid_pl.variables['PRECL'][:]
		nc_fid_pl.close()

		PRECP = (PL+PC)*24*60*60*1000  # from m/s to mm/day ; precipitation =convective+large-scale
		PRECP[PRECP<=0] =np.nan
		
	return lon,lat,time,TS,TN,TX,PRECP

######
#main#
######
def main_body(exp,scenario):
	lon,lat,time,TS,TN,TX,PRECP = data_readin(exp,scenario)
	for interval in ['annual']: #'JJA',"annual"
		if scenario.startswith('F_'):
			year_series = range(2000,2030)
		elif scenario=='T1970RCP':
			year_series = range(2020,2050)
		elif scenario=='EdgEne':
			year_series = range(2200,2230)
		elif scenario=='Edg70GO':
			year_series = range(2070,2100)
		else:
			year_series = range(2130,2160)
		txx,tnn,TSM,TSSTD = T_extremes(lon,lat,time,TS,TX,TN,interval,year_series)
		r10,r20,sdii,rx5day,cdd,cwd,PRECPM,PRECPSTD,RX1DAY,R95,R99 =P_extremes(lon,lat,time,PRECP,interval,year_series)
		#output into netcdf
		output_path = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/analysis/'+interval+'/'
		file_name_out = output_path+exp+'_'+scenario+'.TempPrecp.extremes.'+interval+'.nc'
		f = nc4.Dataset(file_name_out,'w', format='NETCDF4')
		f.createDimension('year', 30)
		f.createDimension('lat', len(lat))
		f.createDimension('lon', len(lon))

		years = f.createVariable('year',np.float32, ('year'))
		lats = f.createVariable('lat',np.float32, ('lat'))
		lons = f.createVariable('lon',np.float32, ('lon'))

		TXXs= f.createVariable('TXX',np.float32,('year','lat','lon'))
		TNNs= f.createVariable('TNN',np.float32,('year','lat','lon'))
		TSMs= f.createVariable('TSM',np.float32,('year','lat','lon'))
		TSSTDs= f.createVariable('TSSTD',np.float32,('year','lat','lon'))
		R10s= f.createVariable('R10',np.float32,('year','lat','lon'))
		R20s= f.createVariable('R20',np.float32,('year','lat','lon'))
		SDIIs= f.createVariable('SDII',np.float32,('year','lat','lon'))
		RX5DAYs= f.createVariable('RX5DAY',np.float32,('year','lat','lon'))
		CDDs= f.createVariable('CDD',np.float32,('year','lat','lon'))
		CWDs= f.createVariable('CWD',np.float32,('year','lat','lon'))
		PRECPMs= f.createVariable('PRECPM',np.float32,('year','lat','lon'))
		PRECPSTDs= f.createVariable('PRECPSTD',np.float32,('year','lat','lon'))
		RX1DAYs= f.createVariable('RX1DAY',np.float32,('year','lat','lon'))
		R95s= f.createVariable('R95',np.float32,('year','lat','lon'))
		R99s= f.createVariable('R99',np.float32,('year','lat','lon'))

		years[:] = range(30)
		TXXs[:] = txx
		TNNs[:] = tnn
		TSMs[:] = TSM
		TSSTDs[:] = TSSTD
		R10s[:] = r10
		R20s[:] = r20
		SDIIs[:] = sdii
		RX5DAYs[:] = rx5day
		CDDs[:] = cdd
		CWDs[:] = cwd
		PRECPMs[:]=PRECPM
		PRECPSTDs[:] = PRECPSTD
		RX1DAYs[:] = RX1DAY
		R95s[:]= R95
		R99s[:] = R99
		lons[:] = lon
		lats[:] = lat

		TXXs.long_name = 'Maximum of daily TX'
		TNNs.long_name =  'minimum of daily TN'
		TSMs.long_name =  'Seasonal mean of TS'
		TSSTDs.long_name =  'Seasonal std of TS'
		R10s.long_name =  'Count of days with precipitation ≥ 10 mm'
		R20s.long_name =  'Count of days with precipitation ≥ 20 mm'
		SDIIs.long_name =  'precipitation mean over wet (daily precp >1mm) days '
		RX5DAYs.long_name =  'maximum amount of precipitation over 5 days'
		CDDs.long_name=  'Maximum number of consecutive dry days '
		CWDs.long_name =  'Maximum number of consecutive wet days'
		PRECPMs.long_name =  'seasonal mean of daily precp'
		PRECPSTDs.long_name =  'seasonal std of daily precp'
		RX1DAYs.long_name =  'maximum amount of precipitation over 1 day'
		R95s.long_name =  '95th percentile of dialy precipitaiton'
		R99s.long_name =  '99th percentile of dialy precipitaiton'

		f.description = 'T and P extreme indices from the CESM-EDGAR experiments'
		f.institution = 'Alcide Zhao at the university of Edinburgh'
		f.close()

exp = 'SlowResp'	#SlowResp FullCp
scenario = "EdgEne" # "T1970RCP" "EdgRef" "Edg70GO" "Edg70Oz" "EdgEne" "EdgTech" "Edg70T10SOZ"
# exp = 'CamOnly'	#CamOnly
# scenario = "F_1970" # F_1970 F_Edg70GO F_Edg70Oz F_Edg70T10SOZ F_EdgEne F_EdgRef F_EdgTech

main_body(exp,scenario)








