"""
This script is to calculate the monthly mean of global T and the percentiles
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d
from scipy.stats import mannwhitneyu as man_test
from scipy.stats import ttest_ind as student_test

def get_global_monthly_T(FREQ='mon',index ='TS'):
	def box_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,mask):
		"""
		fill the range outside the box with 0
		"""
		lon = np.array(lon)
		lat = np.array(lat)
		colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
		colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
		row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
		row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
		if (colum_s> colum_e):
			cache = colum_e; colum_e = colum_s; colum_s = cache;
		if (row_s> row_e):
			cache = row_e; row_e = row_s; row_s = cache;
		mask[:,0:colum_s] =0; mask[:,colum_e:-1] =0
		# plt.imshow(mask,origin='lower');plt.show()
		mask[0:row_s,:] =0; mask[row_e:-1,:] =0
		# plt.imshow(mask,origin='lower');plt.show()
		return mask
			
	def AreaWeight(lon1,lon2,lat1,lat2):
		'''
		calculate the earth radius in m2
		'''
		radius = 6371000;
		area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
		(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
		# print np.nansum(np.nansum(area,axis=1),axis=0)
		return area

		
	def data_readin(variable,FREQ):	
		def day2datetime(scenario,days):
			"""
			# convert days from a reference into int datetime 
			# do not take leap years into account
			"""
			date_int = np.empty((len(days)));date_int[:]=np.nan
			if scenario =='T1970C': start_year =1970
			else: start_year =2010
			start =(start_year*365)
			ith=0	
			for iday in days:
				month_days =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
				calendar_days = np.array([0,31,59,90,120,151,181,212,243,273,304,334,365])
				total_days = int(iday) + start; 
				year = total_days//365; 
				remainder =  total_days%365
				if remainder ==0: year=year-1;month=12;day=31
				else: 
					month = 1+[layer for layer in range(len(calendar_days)) if calendar_days[layer]< remainder and calendar_days[layer+1]>=remainder][0]
					day = int(remainder - calendar_days[month-1])
					if day == 0: day = month_days[month-1]
				date_int[ith] = year*10000+month*100+day
				ith=ith+1
			return date_int.astype(int)
		def annual_and_monthly(scenario,time,data):
			# plt.imshow(data[3,:,:]);plt.show()		
			annual_month=np.empty((40,13,np.shape(data)[1],np.shape(data)[2]));annual_month[:]=np.nan
			# plt.imshow(annual_month[0,3,:,:]);plt.show()
			if scenario=='T1970RCP':
				year_series = range(2020,2050)
			elif scenario=='EdgEne':
				year_series = range(2200,2230)
			elif scenario=='Edg70GO':
				year_series = range(2070,2100)
			else:
				year_series = range(2130,2160)
			for iyear in year_series:
				## annual mean
				if (iyear == year_series[0] and time[0]//100 > year_series[0] *100+1):
					layer_b=0
				else:
					layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  #January
				if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
					layer_e=-2
				else:
					layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+12][0]   #december
				data_cache = data[layer_b:layer_e+1,:,:]
				annual_month[iyear-year_series[0],0,:,:] = np.nanmean(data_cache,axis=0)
				## monthly mean
				for imonth in range(1,13):
					layer = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+imonth][0]  #August 31
					annual_month[iyear-year_series[0],imonth,:,:] = data[layer,:,:]   #:len(layer)
			return annual_month	

		def data_netcdf(scenario,FREQ,variable):
			input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
			var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
			# print var_path
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			data = nc_fid.variables[variable][:]#-273.15
			# missing_value = nc_fid.variables[variable].missing_value; #print missing_value
			# fill_value = nc_fid.variables[variable]._FillValue; #print fill_value
			# data[data==missing_value] = np.nan;data[data==fill_value] = np.nan;
			nc_fid.close()
			if variable=='TS':
				data = data -273.15
			elif variable=='PRECL' or variable=='PRECC':
				data = data*24*60*60*1000  # m/s to mm/day
			var40map = annual_and_monthly(scenario,time,data)
			return lon,lat,var40map
		
		lon,lat,T1970 = data_netcdf('T1970RCP',FREQ,variable)
		lon,lat,Edg70GO = data_netcdf('Edg70GO',FREQ,variable)
		_,_,EdgRef = data_netcdf('EdgRef',FREQ,variable)
		_,_,EdgEne = data_netcdf('EdgEne',FREQ,variable)
		_,_,EdgTech = data_netcdf('EdgTech',FREQ,variable)
		return lon,lat,T1970,Edg70GO,EdgRef,EdgEne,EdgTech

		
	def print_domain_mean(variable,FREQ):
		"""
		This function block is to produce a weighted_mask for specific regions (either administrative or box)
		and then produce the spatial mean (weighted)
		"""

		def mask_weight(region_key,lon,lat):
			"""
			Read in the country mask
			interpolate it to the required resolution grids with lon_interp,lat_interp 
			
			"""
			lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
			lons,lats = np.meshgrid (lon,lat)
			area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
			if region_key == 'All':
				mask=area
				mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
			else:
				##OCEAN_MASKS FOR COUNTRIES
				ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
				lon_mask = ocean_mask['lon'][0,:];
				lat_mask = ocean_mask['lat'][0,:];
				box_region_dic={'Land':[0,360,-90,90],'ASIA':[65,145,5,45],'EUS':[265,280,30,50],'EA':[100,145,20,50],'SA':[65,100,5,30],'SESA':[295,315,-40,-25]}
				if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'Globe'):
					mask= ocean_mask[region_key][:]
				elif (region_key == 'Land' or region_key == 'ASIA' or region_key == 'EA' or region_key == 'SA' or region_key == 'SESA' or region_key == 'EUS'):
					mask= ocean_mask['Globe'][:]
					box = box_region_dic[region_key]
					mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
				else:
					print "no such region"
				# interpolate from 360*720 to 192*288
				mask[np.isnan(mask)]=0;	mask[mask>0]=1;
				f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
				mask = f(lon, lat);
				mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
				# weight each grid cell by its area weight against the total area
				mask=np.multiply(mask,area);  
				mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
				# print np.nansum(np.nansum(mask_weighted,axis=1),axis=0)
			return mask,mask_weighted	

		
		def dif_std_for_region(var1,var2,mask):
			"""
			for a regional dif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
			"""
			
			var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=3),axis=2)
			var2_domain_mean = np.nansum(np.nansum(np.multiply(mask,var2),axis=3),axis=2)
			dif =  np.nanmean(var1_domain_mean - var2_domain_mean,axis=0)
			P25 =  np.nanpercentile(var1_domain_mean - var2_domain_mean,25,axis=0)
			P75 =  np.nanpercentile(var1_domain_mean - var2_domain_mean,75,axis=0)
			return np.stack((dif,P25,P75),axis=1)

		lon,lat,T1970,Edg70GO,EdgRef,EdgEne,EdgTech = data_readin(variable,FREQ);
		region_dics ={'All':np.empty((4,13,3))}  #,'USA','SESA','EA','India'
		for region_key in region_dics:
			_,mask = mask_weight(region_key,lon,lat); #print mask
			region_dics[region_key][0,:,:] = dif_std_for_region(EdgRef,T1970,mask)  #2010-1970
			region_dics[region_key][1,:,:] = dif_std_for_region(Edg70GO,T1970,mask)  #BEoA
			region_dics[region_key][2,:,:] = dif_std_for_region(EdgRef,EdgEne,mask)  #ene
			region_dics[region_key][3,:,:] = dif_std_for_region(EdgRef,EdgTech,mask) ##tech
		return region_dics

  
	TS_dics = print_domain_mean(index,FREQ);
	return TS_dics
