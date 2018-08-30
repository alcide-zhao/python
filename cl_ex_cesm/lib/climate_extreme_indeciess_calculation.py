# -*- coding: utf-8 -*-
"""
Created on Friday 18 2016
This is to calculate the 27 ETCCDI indices.
The input should be daily precipitation TMAX, TMIN. 
NOTE: PRCP units = millimeters and Temperature units= degrees Celsius
"""

import os
import site
import math
import numpy as np
from scipy import stats

#leapyear judgement
# def leapyear(iyear):
	# leapyear = False
	# if (iyear % 400 == 0):
		# leapyear = True
	# else:
		# if (iyear % 100 == 0):
			# leapyear = False
		# else:
			# if (iyear % 4 == 0):
				# leapyear = True
			# else:
				# leapyear = False
				
"""				
# This is the data quality control section which examine the input data and 
# replae all missing and unreasonble values into NaNa the prepared data can 
# then be used for indecies calculation using corresponding functios.
# 
# a) daily precipitation amounts less than zero
# b) daily maximum temperature less than daily minimum temperature
# c) outliers in daily maximum and minimum temperature, defining as the mean
     plus or minus n times standard deviation of the value for the day

def quality_control(source_data)

"""



"""
this is to calculate the precipittion etreme indexies using rhe cliped one year precipitation data
input: the one year precipitation data (precp_year_data), user defined threshold of the  precip
output: 
	1) rx1day    : maximum 1-day precip amount
	2) rx5day    : maximum 5-day precip amount
	3) sdii      : annual precip amount devided by the number of wet days (> 1mm)
	4) r10       : (number of heavy precip days)annual count of days woth precip >= 10mm
	5) r20       : (number of very heavy precip days)annual count of days woth precip >= 20mm
	6) rnm       : annual count of days woth precip >= nmm
	7) cdd       : (consective dry days) maximum consecive days with rr<1mm
	8) cwd       : (consective wet days) maximum consecive days with rr>=1mm
	9) r95p      : (very wet days) annual total precp >95th percetile referring to 1961-1990
	10) r99p     : (extremely wet days) annual total precp >99th percetile referring to 1961-1990
	11) precptot : (annual total wet-day prrecip) annual total precp in wet days (rr>1mm)
"""

def precip_extreme_indeces(precp_year_data,r95p_threshold,r99p_threshold,rnt):
	total_precip = np.nansum(precp_year_data, axis = 0)  # note here the first dimenssion rpresents time series
	mean_prep = stats.nanmean(precp_year_data, axis = 0)
	std_prep = stats.nanstd(precp_year_data, axis = 0)
	rx1day =np.nanmax(precp_year_data,axis = 0)
	size_data = np.shape(precp_year_data)
	# r10 r20 rnm sdii
	def rnnmm(rnt, r1t=1, r10t=10, r20t=20):
		#######################
		# variables initiation
		#######################
		r10 = np.empty((size_data[1],size_data[2]))
		r10[:] = np.nan
		r20 = np.empty((size_data[1],size_data[2]))
		r20[:] = np.nan
		rnm = np.empty((size_data[1],size_data[2]))
		rnm[:] = np.nan
		sdii = np.empty((size_data[1],size_data[2]))
		sdii[:] = np.nan
		precptot = np.empty((size_data[1],size_data[2]))
		precptot[:] = np.nan
		rx5day = np.empty((size_data[1],size_data[2]))
		rx5day[:] = np.nan
		r95p = np.empty((size_data[1],size_data[2]))
		r95p[:] = np.nan
		r99p = np.empty((size_data[1],size_data[2]))
		r99p[:] = np.nan
		r1 = np.empty((size_data[1],size_data[2]))
		r1[:] = np.nan		

		for row in range(size_data[1]):
			for colum in range(size_data[2]):
				# r95p_threshold=np.percentile(prec_cache,95)
				prec_cache = precp_year_data[:,row,colum]
				r95p[row,colum] = np.nansum([value for value in prec_cache if value > r95p_threshold[row,colum] ])
				r99p[row,colum] = np.nansum([value for value in prec_cache if value > r99p_threshold[row,colum] ])	
				# r1 r10 r20 rnm and prectot
				r1[row,colum]=len([value for value in prec_cache if value >= r1t])
				precptot[row,colum] = np.nansum([value for value in prec_cache if value >= r1t])
				r10[row,colum]=len([value for value in prec_cache if value >= r10t])
				r20[row,colum]=len([value for value in prec_cache if value >= r20t])
				rnm[row,colum]=len([value for value in prec_cache if value >= rnt])
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
		rnm[rnm < 0] = np.nan
		rx5day[rx5day < 0] = np.nan
		r1 = r1.astype('float')
		r1[ r1 <= 0] = np.nan
		sdii = np.divide(total_precip, r1)
		return r10, r20, rnm, sdii, precptot, rx5day, rx1day, r95p, r99p
	
							
	def cdd_cwd():
		#######################
		# variables initiation
		#######################
		cdd = np.empty((size_data[1],size_data[2]))
		cdd[:] = np.nan
		cwd = np.empty((size_data[1],size_data[2]))
		cwd[:] = np.nan
		
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
	r10, r20, rnm, sdii, precptot, rx5day, rx1day, r95p, r99p= rnnmm(rnt, r1t=1, r10t=10, r20t=20)
	cdd,cwd = cdd_cwd()
	
	return r10, r20, rnm, sdii, precptot, rx5day, rx1day, r95p, r99p, cdd, cwd, total_precip, mean_prep, std_prep



"""
this is to calculate the temperature etreme indecies using rhe cliped one year precipitation data
input: the min and max daily temperature of the interested year
output: 
	1) Fd0      : (frost days) annual count when TN < 0
	2) su25     : (summer days) annual count when TX < =25
	3) id0   	: (ice days) annual count when TX < 0
	4) tr20		: (tropical nights ) annual count when TN > 20
	5) gsl 		: (growing season length ) Annual (1st Jan to 31st Dec in NH, 1st July to 30th June 
				  in SH) count between first span of at least 6 days with TG>5 ◦C and first span after 
				  July 1 (January 1 in SH) of 6 days with TG<5 ◦C
	6) txx 		: (max tmax) annual maximum of daily maximum
	7) tnx		: (max tmin) annual minimum of daily maximum
	8) txn		: (min tmax) annual maximum of daily minimum
	9) tn10p    : (cold night) percentage of days when TN < 10th percentile
	10)tx10p 	: (cool days) percentage of days when TX < 10th percentile
	11)tn90p    : (warm night) percentage of days when TN > 90th percentile
	12)tx90p 	: (warm days) percentage of days when TX > 90th percentile
	13)dtr 		: (diurnal temperature range) annual mean differance between TX and TN  
	14)tnn		: (min tmin) annual minimum of daily minimum
"""
def temperature_extreme_indeces(daily_minimum_year,daily_maximum_year,tn10p_min_threshold,tn90p_min_threshold,tx10p_max_threshold,tx90p_max_threshold):
	size_data = np.shape(daily_minimum_year)
	txx = np.nanmax(daily_maximum_year, axis=0)
	txn = np.nanmin(daily_maximum_year, axis=0)
	tnx = np.nanmax(daily_minimum_year, axis=0)
	tnn = np.nanmin(daily_minimum_year, axis=0)
	dtr = stats.nanmean(daily_maximum_year - daily_minimum_year, axis = 0)
	mean_tn = stats.nanmean(daily_minimum_year, axis = 0)
	mean_tx = stats.nanmean(daily_maximum_year, axis = 0)
	std_tx = stats.nanstd(daily_maximum_year, axis = 0)
	std_tn = stats.nanstd(daily_minimum_year, axis = 0)
	
	fd0 = np.empty((size_data[1],size_data[2]))
	fd0[:] = np.nan
	id0 = np.empty((size_data[1],size_data[2]))
	id0[:] = np.nan
	su25 = np.empty((size_data[1],size_data[2]))
	su25[:] = np.nan
	tr20 = np.empty((size_data[1],size_data[2]))
	tr20[:] = np.nan
	tn10p = np.empty((size_data[1],size_data[2]))
	tn10p[:] = np.nan
	tx90p = np.empty((size_data[1],size_data[2]))
	tx90p[:] = np.nan
	tn90p = np.empty((size_data[1],size_data[2]))
	tn90p[:] = np.nan	
	tx10p = np.empty((size_data[1],size_data[2]))
	tx10p[:] = np.nan	
	for row in range(size_data[1]):
		for colum in range(size_data[2]):
			min_temp_cache = daily_minimum_year[:,row,colum]
			max_temp_cache = daily_maximum_year[:,row,colum]
			fd0[row,colum] = len([value for value in min_temp_cache if value< 0])
			su25[row,colum] = len([value for value in max_temp_cache if value> 25])
			id0[row,colum] = len([value for value in max_temp_cache if value< 0])
			tr20[row,colum] = len ([value for value in min_temp_cache if value> 20])
			tn10p[row,colum] = 100*len([value for value in min_temp_cache if value< tn10p_min_threshold[row,colum]])/ float(len(min_temp_cache))
			tn90p[row,colum] = 100*len([value for value in min_temp_cache if value> tn90p_min_threshold[row,colum]])/ float(len(min_temp_cache))
			tx10p[row,colum] = 100*len([value for value in max_temp_cache if value< tx10p_max_threshold[row,colum]])/ float(len(max_temp_cache))
			tx90p[row,colum] = 100*len([value for value in max_temp_cache if value> tx90p_max_threshold[row,colum]])/ float(len(max_temp_cache))
	return mean_tx,mean_tn,std_tx,std_tn,txx, txn, tnx, tnn, dtr, fd0, su25, id0, tr20, tn10p, tn90p, tx10p, tx90p