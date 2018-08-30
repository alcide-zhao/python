"""
Thius is to test the relationship between Net radiation at TOA (FSNT-FLNT at TOM in this case)
namely, to calculate the radiative forcing and climate sensitivity
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import math


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


def data_readin(variable,FREQ):	
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		if scenario =='EdgRef70AP': start_year =1990
		elif scenario =='EdgTechP': start_year =2180
		else: start_year =2000
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
	def mon_mean2annual_mean(scenario,time,data):
		
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario=='EdgRef70AP':
			year_series = range(2020,2050)
		elif scenario=='EdgTechP':
			year_series = range(2210,2240)	
		else:
			year_series = range(2030,2060)
		# print scenario,iyear,layer_b,layer_e
		annual_mean=np.empty((len(year_series),192,288));annual_mean[:]=np.nan
		for iyear in year_series:
			# print scenario,iyear,time//100
			annual = np.empty((192,288));annual[:] = 0.0			
			if (iyear == year_series[0] and time[0]//100 >= year_series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  
			if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+12][0]
			
			data_cache = data[layer_b:layer_e+1,:,:]
			annual_mean[iyear-year_series[0],:,:] = stats.nanmean(data_cache,axis=0)
		return annual_mean

	def data_netcdf(scenario,FREQ,variable):
		def AreaWeight(lon1,lon2,lat1,lat2):
			'''
			calculate the earth radius in m2
			'''
			radius = 6371000;
			area = (np.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
			(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
			# print np.nansum(np.nansum(area,axis=1),axis=0)
			return area
		
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/full_chem/'
		var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]
		nc_fid.close()
		
		lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
		lons,lats = np.meshgrid (lon,lat)
		area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
		# if scenario == 'T1970': print time
		area = area/np.nansum(np.nansum(area,axis=1),axis=0)
		var = mon_mean2annual_mean(scenario,time,data)
		var_TS = np.nansum(np.nansum(np.multiply(var,area),axis=2),axis=1)
		
		return var_TS
	
	

	Edg1970 = data_netcdf('Edg1970',FREQ,variable)
	EdgRefP = data_netcdf('EdgRefP',FREQ,variable)
	EdgEneP = data_netcdf('EdgEneP',FREQ,variable)
	EdgTechP = data_netcdf('EdgTechP',FREQ,variable)
	EdgRef70AP = data_netcdf('EdgRef70AP',FREQ,variable)
	
	return Edg1970,EdgRefP,EdgEneP ,EdgTechP,EdgRef70AP
	
	
Edg1970_TS,EdgRefP_TS,EdgEneP_TS ,EdgTechP_TS,EdgRef70AP_TS = data_readin(variable='TS',FREQ='mon');   #,

print 'Total', np.nanmean(EdgRefP_TS[10:30]-Edg1970_TS[10:30],axis=0)
print 'AAs', np.nanmean(EdgRefP_TS[10:30]-EdgRef70AP_TS[10:30],axis=0)
print 'Ene', np.nanmean(EdgRefP_TS[10:30]-EdgEneP_TS[10:30],axis=0)
print 'Tech', np.nanmean(EdgRefP_TS[10:30]-EdgTechP_TS[10:30],axis=0)




# T1970_LW,Edg70GO_LW,Edg70Oz_LW,EdgRef_LW,EdgEne_LW,EdgTech_LW = data_readin(variable='FLNT',FREQ='mon');  #
# T1970_SW,Edg70GO_SW,Edg70Oz_SW,EdgRef_SW,EdgEne_SW,EdgTech_SW = data_readin(variable='FSNT',FREQ='mon');  #

fig = plt.figure(facecolor='White',figsize=[10,6]);plot_setup();pad= 5;
ax1 = plt.subplot(1,1,1);
ax1.annotate('(a) SAT (K)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
ax1.plot(range(0,len(Edg1970_TS)),Edg1970_TS,'-',color="g",lw=3,label ='1970')				
ax1.plot(range(0,len(EdgRefP_TS)),EdgRefP_TS,'-',color="B",lw=3,label ='2010')	
ax1.plot(range(0,len(EdgEneP_TS)),EdgEneP_TS,'-',color="R",lw=3,label ='STAG_ENE')	
ax1.plot(range(0,len(EdgTechP_TS)),EdgTechP_TS,'-',color="M",lw=3,label ='STAG_TECH')
ax1.plot(range(0,len(EdgRef70AP_TS)),EdgRef70AP_TS,'-',color="k",lw=3,label ='STAG_AAs')
legend = ax1.legend(shadow=False,ncol=2,loc ='lower right')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)
ax1.set_ylim([287,290]);
# ax1.set_ylim([-1.0,2.0]);ax1.set_xlim([-0.3,1.5])
# ax1.axvline(x=0,c='k')
# ax1.axhline(y=0,color='k',linewidth = 1);

plt.subplots_adjust(left=0.10, bottom=0.03, right=0.95, top=0.95, wspace=0.20, hspace=0.15); 
plt.savefig('SAT_evelution.png', format='png', dpi=1000)
