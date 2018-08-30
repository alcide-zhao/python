"""
This is to process the monthly fields into zonal mean and plot 
	currently, SAT ans precipitation
	the model internal variability (25-t5 th percentile) and mean are both plotted
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
from scipy import integrate
# from climlab import constants as const

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

def inferred_heat_transport( energy_in, lat_deg):
	'''Returns the inferred heat transport (in PW) by integrating the net energy imbalance from pole to pole.'''
	a = 6.373*10**6
	# print np.shape(energy_in)
	# plt.plot(lat_deg,energy_in)
	# plt.plot(np.flip(lat_deg),np.flip(energy_in))
	# print np.flip(energy_in,axis=0)
	imbal_ncep = np.average(energy_in, weights=np.cos(np.deg2rad(lat_deg)),axis=1); 
	for imember in range(np.shape(energy_in)[0]):
		energy_in[imember,:] = energy_in[imember,:] - imbal_ncep[imember]
	lat_rad = np.deg2rad(lat_deg)
	inter_S2N = ( 1E-15 * 2 * np.math.pi *a**2 * 
		integrate.cumtrapz(np.cos(lat_rad)*energy_in,x=lat_rad, initial=0.))
	# inter_N2S = ( 1E-15 * 2 * np.math.pi *a**2 * 
		# integrate.cumtrapz(np.cos(np.flip(lat_rad,axis=0))*np.flip(energy_in,axis=1),x=lat_rad, initial=0.))	
	return inter_S2N
			
def data_read_global_mean(scenario,variable):	
	"""
	calculate the earth radius in m2
	"""		
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
	def mask_weight(lon,lat):
		"""
		Read in the country mask
		interpolate it to the required resolution grids with lon_interp,lat_interp 
		
		"""
		def AreaWeight(lon1,lon2,lat1,lat2):
			'''
			calculate the earth radius in m2
			'''
			radius = 6371000;
			area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
			(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
			# print np.nansum(np.nansum(area,axis=1),axis=0)
			return area
		lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
		lons,lats = np.meshgrid (lon,lat)
		area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
		mask=area
		mask_weighted = np.divide(area,np.nansum(np.nansum(area,axis=1),axis=0))	
		return mask_weighted
	def mon_mean2annual_mean(scenario,time,data):
		annual_mean=np.empty((40,192,288));annual_mean[:]=np.nan
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario=='T1970RCP':
			year_series = range(2020,2050)
		elif scenario=='EdgEne':
			year_series = range(2200,2230)
		elif scenario=='Edg70GO':
			year_series = range(2070,2100)
		else:
			year_series = range(2130,2160)
		for iyear in year_series:
			
			if (iyear == year_series[0] and time[0]//100 >= year_series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  #June01
			if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+12][0]  #August 31
			data_cache = data[layer_b:layer_e+1,:,:]
			annual_mean[iyear-year_series[0],:,:] = stats.nanmean(data_cache,axis=0)
		return annual_mean

	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]
		if variable =='TS': data=data-273.15
		nc_fid.close()
		var40map = mon_mean2annual_mean(scenario,time,data)
		mask =mask_weight(lon,lat)
		global_mean = np.nansum(np.nansum(np.multiply(var40map,mask),axis=2),axis=1)
		return global_mean
	
	data = data_netcdf(scenario,variable)
	return data
	
def get_flux_calculate_Sensitivity(scenario1,scenario2,):
	ET = {'T1970RCP':0.45,'Edg70GO':0.38,'EdgRef':1.15,'EdgEne':1.44,'EdgTech':0.16}
	RF = {'T1970RCP':0.7,'Edg70GO':0.60,'EdgRef':1.25,'EdgEne':1.52,'EdgTech':1.12}
	# def get_mean_spread(data):
		# mean = round(np.nanmean(data,axis=0),2)
		# p25 = round(np.nanpercentile(data,25,axis=0),2)
		# p75 = round(np.nanpercentile(data,75,axis=0),2)
		# return mean
	
	impose = RF[scenario1]-RF[scenario2];
	normT = ET[scenario1]-ET[scenario2]; 
	variable = 'TS';normT =   np.nanmean(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))

	variable = 'FSNT'; FSNT =  np.nanmean(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'FSNSC';FSNSC=  np.nanmean(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'FLNT'; FLNT =  np.nanmean(-(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable)))
	variable = 'FLNSC';FLNSC=  np.nanmean( -(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable)))
	residual = np.nanmean(FSNT+FLNT);	
	
	FSNT =FSNT-residual; FLNT =FLNT+residual
	lw = round(((FLNT - impose))/normT,2)     # LW feedback
	lwclr = round(((FLNSC-impose))/normT,2)    # LW clearsky feedback
	lwclf = round((FLNT-FLNSC)/normT,2)      # lw cloud-forcing feedback
	sw = round((FSNT)/normT,2)                  # SW feedback
	swclr = round((FSNSC)/normT,2)             # SW  clearsky feedback
	swclf = round((FSNT-FSNSC)/normT ,2)     # LW cloud-forcing feedback
	return lw,lwclr,lwclf,sw,swclr,swclf,round(lw+sw,2),-round(impose/normT,2), round(residual/normT,2)

def get_flux_response(scenario1,scenario2):
	def get_mean_spread(data):
		mean = round(np.nanmean(data,axis=0),2)
		p25 = round(np.nanpercentile(data,25,axis=0),2)
		p75 = round(np.nanpercentile(data,75,axis=0),2)
		return np.stack((mean,p25,p75),axis=0)
	# variable = 'FSNT'; FSNT =  get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	# variable = 'FSNTC'; FSNTC = get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	# variable = 'FLNT'; FLNT = get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	# variable = 'FLNTC';FLNTC= get_mean_spread (data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'FSNS';FSNS=  get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'FSNSC';FSNSC=  get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'FLNS'; FLNS = get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'FLNSC';FLNSC=  get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'SWCF';SWCF=  get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	variable = 'LWCF';LWCF= get_mean_spread(data_read_global_mean(scenario1,variable)-data_read_global_mean(scenario2,variable))
	return np.stack((FSNS,FSNSC,FLNS,FLNSC,SWCF,FSNS-FSNSC,LWCF,FLNS-FLNSC),axis=1)
	
	
def bar_err(data,step,color):
	# print "------"+var+"--------"
	def autolabel(X_pos,values,height_lift):
		## Attach a text label above each bar displaying its height
		height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
		for i in range(len(height)):
			ax.text(X_pos[i],y_pos[i],'%4.2f' % height[i], ha='center', va='bottom',size=15)
	
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 2);
	ax.set_xlim([-0.5,15.5]);ax.set_xticks(np.arange(0.75,16.5,2));ax.set_xticklabels(('FSNS','FSNSC','FLNS','FLNSC','SWCF','FSNS-FSNSC','LWCF','FLNS-FLNSC'));
	ax.set_ylim([-1,1.0]);ax.set_yticks(np.arange(-1,1.01,0.5))
	ax.axvspan(1.5, 3.5, alpha=0.01, color='r',edgecolor='none');
	ax.axvspan(5.5, 7.5, alpha=0.01, color='r',edgecolor='none')
	ax.axvspan(9.5, 11.5, alpha=0.01, color='r',edgecolor='none')
	ax.axvspan(13.5, 15.5, alpha=0.01, color='r',edgecolor='none')
	

	X_pos=np.arange(0,15.5,2);dis=0.5
	y_value = data[0,:];yerr = np.stack((y_value-data[1,:],data[2,:]- y_value),axis=0);
	ax.bar(X_pos+step*dis, y_value, yerr=yerr, align='center',color=color, ecolor=color,capsize=0,alpha=0.7,width=0.5,lw=0)
	# autolabel(X_pos+step*dis,y_value,1.20)


fig = plt.figure(facecolor='White',figsize=[8,5]);
ax = plt.subplot(1,1,1);	
# scenario1='EdgRef';scenario2 = 'T1970RCP';RF=get_flux_response(scenario1,scenario2);
# bar_err(RF,0,'k')
scenario1='Edg70GO';scenario2 = 'T1970RCP';RF=get_flux_response(scenario1,scenario2);
bar_err(RF,0,'b')
scenario1='EdgRef';scenario2 = 'EdgEne';RF=get_flux_response(scenario1,scenario2);
bar_err(RF,1,'r')
scenario1='EdgRef';scenario2 = 'EdgTech';RF=get_flux_response(scenario1,scenario2);
bar_err(RF,2,'g')

plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.95, wspace=0.04, hspace=0.15); 
plt.savefig('energy_balance.png', format='png', dpi=1000)











