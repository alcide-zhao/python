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
from matplotlib.patches import Rectangle


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
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
		if scenario=='T1970C' or scenario=='T1970': start_year =1970
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
		
	def mon_mean2annual_mean(scenario,time,data):
		
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario=='T1970C' or scenario=='T1970':
			year_series = range(1970,2010)
		elif scenario=='T1970RCP':
			year_series = range(2010,2050)
		elif scenario=='EdgEne':
			year_series = range(2010,2230)
		elif scenario=='Edg70GO':
			year_series = range(2010,2101)
		elif scenario=='Edg70T10SOZ':
			year_series = range(2010,2160)
		elif scenario=='EdgRef70A':
			year_series = range(2010,2160)
		else:
			year_series = range(2010,2160)
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
		
		
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
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
	
	
	T1970 = data_netcdf('T1970',FREQ,variable);
	T1970C = data_netcdf('T1970C',FREQ,variable);
	T1970R = data_netcdf('T1970RCP',FREQ,variable);
	T1970 = np.concatenate([T1970,T1970C,T1970R],axis=0)   # two experiments in sequence
	
	Edg70GO = data_netcdf('Edg70GO',FREQ,variable);
	EdgRef = data_netcdf('EdgRef',FREQ,variable)
	Edg70Oz = data_netcdf('Edg70Oz',FREQ,variable)
	EdgEne = data_netcdf('EdgEne',FREQ,variable)
	EdgTech = data_netcdf('EdgTech',FREQ,variable)
	Edg70T10SOZ = data_netcdf('Edg70T10SOZ',FREQ,variable)
	EdgRef70A = data_netcdf('EdgRef70A',FREQ,variable)
	
	return T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne ,EdgTech,Edg70T10SOZ,EdgRef70A
	
	

def pp_regress(TS,SW,LW,windows=1):
	"""
	# This is to process SAT into net chanegs relative to the start year
	# Net radiation at TOM  NRF = (FSTN - FLTN)
	# Then regress NRF against dSAT
	# y =a*x+bn
	"""
	dNR_s =smooth(SW-LW,5)
	dTS_s = smooth(TS,5)- np.min(TS);    #  
	a,b,a_range,b_range,r2 = LinearRegression_Stats(dTS_s,dNR_s)
	return round(a,2),round(b,2),a_range,b_range,round(-1.0*b/a,2),TS,SW -LW

def TC_regress(TS):
	"""
	# This is to FIT SAT  change onto time to calculate the time constant
	"""
	# TS = TS[:-40]
	min = np.min(TS); max = np.max(TS)
	TS = (TS-min)/(max-min)
	def f(x,a,b):
		return (np.exp(-x/b))
	x=range(0,len(TS));
	popt, pcov = curve_fit(f, x, TS)
	a = popt[0]
	b = popt[1]
	return b,popt



T1970_TS,Edg70GO_TS,Edg70Oz_TS,EdgRef_TS,EdgEne_TS,EdgTech_TS,Edg70T10SOZ_TS,EdgRef70A_TS = data_readin(variable='TS',FREQ='mon');   #,

##################  Gregory plot ###################  
T1970_LW,Edg70GO_LW,Edg70Oz_LW,EdgRef_LW,EdgEne_LW,EdgTech_LW,Edg70T10SOZ_LW,EdgRef70A_LW = data_readin(variable='FLNT',FREQ='mon');  #
T1970_SW,Edg70GO_SW,Edg70Oz_SW,EdgRef_SW,EdgEne_SW,EdgTech_SW,Edg70T10SOZ_SW,EdgRef70A_SW = data_readin(variable='FSNT',FREQ='mon');  #	
	
T1970_s, T1970_i, T1970_r, T1970_p, T1970_std, T1970_TS, T1970_NR = pp_regress(T1970_TS,T1970_SW,T1970_LW)
print "E1970", T1970_s,T1970_i, T1970_r, T1970_p, T1970_std
EdgRef_s, EdgRef_i, EdgRef_r, EdgRef_p, EdgRef_std, EdgRef_TS, EdgRef_NR = pp_regress(EdgRef_TS,EdgRef_SW,EdgRef_LW)
print "EdgRef", EdgRef_s,EdgRef_i, EdgRef_r, EdgRef_p, EdgRef_std
Edg70GO_s, Edg70GO_i, Edg70GO_r, Edg70GO_p, Edg70GO_std, Edg70GO_TS, Edg70GO_NR = pp_regress(Edg70GO_TS,Edg70GO_SW,Edg70GO_LW)
print "Edg70GO", Edg70GO_s,Edg70GO_i, Edg70GO_r, Edg70GO_p, Edg70GO_std
Edg70Oz_s, Edg70Oz_i, Edg70Oz_r, Edg70Oz_p, Edg70Oz_std, Edg70Oz_TS, Edg70Oz_NR = pp_regress(Edg70Oz_TS,Edg70Oz_SW,Edg70Oz_LW)
print "Edg70Oz", Edg70Oz_s,Edg70Oz_i, Edg70Oz_r, Edg70Oz_p, Edg70Oz_std
Edg70GO_s, Edg70GO_i, Edg70GO_r, Edg70GO_p, Edg70GO_std, Edg70GO_TS, Edg70GO_NR = pp_regress(Edg70GO_TS,Edg70GO_SW,Edg70GO_LW)
print "Edg70GO", Edg70GO_s,Edg70GO_i, Edg70GO_r, Edg70GO_p, Edg70GO_std
Edg70Oz_s, Edg70Oz_i, Edg70Oz_r, Edg70Oz_p, Edg70Oz_std, Edg70Oz_TS, Edg70Oz_NR = pp_regress(Edg70Oz_TS,Edg70Oz_SW,Edg70Oz_LW)
print "Edg70Oz", Edg70Oz_s,Edg70Oz_i, Edg70Oz_r, Edg70Oz_p, Edg70Oz_std
Edg70T10SOZ_s, Edg70T10SOZ_i, Edg70T10SOZ_r, Edg70T10SOZ_p, Edg70T10SOZ_std, Edg70T10SOZ_TS, Edg70T10SOZ_NR = pp_regress(Edg70T10SOZ_TS,Edg70T10SOZ_SW,Edg70T10SOZ_LW)
print "Edg70T10SOZ", Edg70T10SOZ_s,Edg70T10SOZ_i, Edg70T10SOZ_r, Edg70T10SOZ_p, Edg70T10SOZ_std
EdgRef70A_s, EdgRef70A_i, EdgRef70A_r, EdgRef70A_p, EdgRef70A_std, EdgRef70A_TS, EdgRef70A_NR = pp_regress(EdgRef70A_TS,EdgRef70A_SW,EdgRef70A_LW)
print "EdgRef70A", EdgRef70A_s,EdgRef70A_i, EdgRef70A_r, EdgRef70A_p, EdgRef70A_std

fig = plt.figure(facecolor='White',figsize=[6,8]);plot_setup();pad= 5;
ax1 = plt.subplot(2,1,1);
ax1.annotate('(a) SAT (K)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
ax1.plot(range(1,len(T1970_TS)+1),smooth(T1970_TS,10),'-',color="k",lw=2,label ='B70')				
ax1.plot(range(1,len(Edg70GO_TS)+1),smooth(Edg70GO_TS,10),'-',color="g",lw=2,label ='SGO')	
ax1.plot(range(1,len(Edg70Oz_TS)+1),smooth(Edg70Oz_TS,10),'-',color="m",lw=2,label ='SOZ')	
ax1.plot(range(1,len(EdgRef_TS)+1),smooth(EdgRef_TS,10),'-',color="b",lw=2,label ='B10')	
ax1.plot(range(1,len(Edg70T10SOZ_TS)+1),smooth(Edg70T10SOZ_TS,10),'-',color="r",lw=2,label ='STO')	
# ax1.plot(range(1,len(EdgRef70A_TS)+1),smooth(EdgRef70A_TS,10),'-',color="r",lw=2,label ='STAG_AAs')

ax1.add_patch(Rectangle((120, 288.35), 30, 0.55,edgecolor='k', fill=None,alpha=1))
ax1.add_patch(Rectangle((90, 287.30), 30, 0.35, edgecolor='k', fill=None,alpha=1))
ax1.add_patch(Rectangle((60, 287.40), 30, 0.4, edgecolor='k', fill=None,alpha=1))


legend = ax1.legend(shadow=False,ncol=3,loc ='lower right')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)


ax1.plot(range(1,len(T1970_TS)+1),T1970_TS,'-',color="k",lw=0.5,alpha=1,label ='E1970')				
ax1.plot(range(1,len(Edg70GO_TS)+1),Edg70GO_TS,'-',color="g",lw=0.5,alpha=1,label ='STAG_G&O')	
ax1.plot(range(1,len(Edg70Oz_TS)+1),Edg70Oz_TS,'-',color="m",lw=0.5,alpha=1,label ='STAG_O3')	
ax1.plot(range(1,len(EdgRef_TS)+1),EdgRef_TS,'-',color="b",lw=0.5,alpha=1,label ='E2010')	
ax1.plot(range(1,len(Edg70T10SOZ_TS)+1),Edg70T10SOZ_TS,'-',color="r",lw=0.5,alpha=1,label ='STAG_TrO3')	
# ax1.plot(range(1,len(EdgRef70A_TS)+1),EdgRef70A_TS,'-',color="r",lw=1,alpha=1,label ='STAG_AAs')	

	
ax1.set_xlim([0,151]); ax1.set_xticks(np.arange(0,151,30)); #ax1.set_ylim([-0.3,1.0]);			
ax1.set_ylim([287,289]);			

print 'E1970',np.nanmean(T1970_TS[-30:-1]),np.nanmean(T1970_TS[-51:-11])
print 'STAG_G&O',np.nanmean(Edg70GO_TS[-30:-1]),np.nanmean(Edg70GO_TS[-51:-11])
print 'STAG_O3',np.nanmean(Edg70Oz_TS[-30:-1]),np.nanmean(Edg70Oz_TS[-51:-11])
print 'Edg70T10SOZ',np.nanmean(Edg70T10SOZ_TS[-30:-1]),np.nanmean(Edg70T10SOZ_TS[-51:-11])
print 'EdgRef70A',np.nanmean(EdgRef70A_TS[-30:-1]),np.nanmean(EdgRef70A_TS[-51:-11])
print 'EdgRef_TS',np.nanmean(EdgRef_TS[-30:-1]),np.nanmean(EdgRef_TS[-51:-11])


print 'E1970',np.nanmean(T1970_TS[-40:-1]),np.nanmean(T1970_TS[-51:-11])
ax1 = plt.subplot(2,1,2);
ax1.annotate('(b) Radiation imbalance '+r'($\mathrm{\/W\/m^{-2}}$)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
ax1.plot(range(1,len(T1970_NR)+1),smooth(T1970_NR,10),'-',color="k",lw=2,label ='E1970')				
ax1.plot(range(1,len(Edg70GO_NR)+1),smooth(Edg70GO_NR,10),'-',color="g",lw=2,label ='STAG_G&O')	
ax1.plot(range(1,len(Edg70Oz_NR)+1),smooth(Edg70Oz_NR,10),'-',color="m",lw=2,label ='STAG_O3')	
ax1.plot(range(1,len(EdgRef_NR)+1),smooth(EdgRef_NR,10),'-',color="b",lw=2,label ='E2010')	
ax1.plot(range(1,len(Edg70T10SOZ_NR)+1),smooth(Edg70T10SOZ_NR,10),'-',color="r",lw=2,label ='STAG_TrO3')	
# ax1.plot(range(1,len(EdgRef70A_NR)+1),smooth(EdgRef70A_NR,10),'-',color="r",lw=2,label ='STAG_AAs')
ax1.plot(range(1,len(T1970_NR)+1),T1970_NR,'-',color="k",lw=0.2,alpha=1,label ='E1970')				
ax1.plot(range(1,len(Edg70GO_NR)+1),Edg70GO_NR,'-',color="g",lw=0.2,alpha=1,label ='STAG_G&O')	
ax1.plot(range(1,len(Edg70Oz_NR)+1),Edg70Oz_NR,'-',color="m",lw=0.2,alpha=1,label ='STAG_O3')	
ax1.plot(range(1,len(EdgRef_NR)+1),EdgRef_NR,'-',color="b",lw=0.2,alpha=1,label ='E2010')	
ax1.plot(range(1,len(Edg70T10SOZ_NR)+1),Edg70T10SOZ_NR,'-',color="r",lw=0.2,alpha=1,label ='STAG_TrO3')	
# ax1.plot(range(1,len(EdgRef70A_NR)+1),EdgRef70A_NR,'-',color="r",lw=1,alpha=1,label ='STAG_AAs')	
ax1.axhline(y=0,color='k',linewidth = 1);
ax1.set_ylim([-0.5,1.5]);ax1.set_xlim([0,151]);ax1.set_xticks(np.arange(0,151,30));

plt.subplots_adjust(left=0.10, bottom=0.03, right=0.95, top=0.95, wspace=0.20, hspace=0.15); 
plt.savefig('TS_EVELUTION.png', format='png', dpi=1000)


# ax1.plot(range(1,len(EdgEne_TS)+1),smooth(EdgEne_TS,10),'-',color="r",lw=2,label ='STAG_ENE')	
# ax1.plot(range(1,len(EdgTech_TS)+1),smooth(EdgTech_TS,10),'-',color="m",lw=2,label ='STAG_TECH')
# ax1.plot(range(1,len(EdgEne_TS)+1),EdgEne_TS,'-',color="r",lw=0.5,alpha=1,label ='STAG_ENE')	
# ax1.plot(range(1,len(EdgTech_TS)+1),EdgTech_TS,'-',color="m",lw=0.5,alpha=1,label ='STAG_TECH')	
# ax1.plot(range(1,len(EdgEne_NR)+1),smooth(EdgEne_NR,10),'-',color="r",lw=2,label ='STAG_ENE')	
# ax1.plot(range(1,len(EdgTech_NR)+1),smooth(EdgTech_NR,10),'-',color="m",lw=2,label ='STAG_TECH')
# ax1.plot(range(1,len(EdgEne_NR)+1),EdgEne_NR,'-',color="r",lw=0.2,alpha=1,label ='STAG_ENE')	
# ax1.plot(range(1,len(EdgTech_NR)+1),EdgTech_NR,'-',color="m",lw=0.2,alpha=1,label ='STAG_TECH')		

###################   TIME CONSTANT ###################  
# T1970_s, T1970_i = TC_regress(T1970_SW-T1970_LW)
# print "E1970", T1970_s,T1970_i
# EdgRef_s, EdgRef_i = TC_regress(EdgRef_SW-EdgRef_LW)
# print "EdgRef", EdgRef_s,EdgRef_i
# Edg70GO_s, Edg70GO_i = TC_regress(Edg70GO_SW-Edg70GO_LW)
# print "Edg70GO", Edg70GO_s,Edg70GO_i
# Edg70Oz_s, Edg70Oz_i= TC_regress(Edg70Oz_SW-Edg70Oz_LW)
# print "Edg70Oz", Edg70Oz_s,Edg70Oz_i
# EdgEne_s, EdgEne_i = TC_regress(EdgEne_SW-EdgEne_LW)
# print "EdgEne", EdgEne_s,EdgEne_i
# EdgTech_s, EdgTech_i = TC_regress(EdgTech_SW-EdgTech_LW)
# print "EdgTech", EdgTech_s,EdgTech_i



"""				
ax1 = plt.subplot(3,1,2);
ax1.annotate('(b) FSNT-SLNT (W/m2)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
ax1.plot(range(1,len(T1970_NR)+1),smooth(T1970_NR,10),'-',color="g",lw=2,label ='E1970')				
ax1.plot(range(1,len(Edg70GO_NR)+1),smooth(Edg70GO_NR,10),'-',color="k",lw=2,label ='STAG_G&O')	
ax1.plot(range(1,len(Edg70Oz_NR)+1),smooth(Edg70Oz_NR,10),'-',color="m",lw=2,label ='STAG_O3')	
ax1.plot(range(1,len(EdgRef_NR)+1),smooth(EdgRef_NR,10),'-',color="b",lw=2,label ='E2010')	
ax1.plot(range(1,len(EdgEne_NR)+1),smooth(EdgEne_NR,10),'-',color="r",lw=2,label ='STAG_ENE')	
ax1.plot(range(1,len(EdgTech_NR)+1),smooth(EdgTech_NR,10),'-',color="m",lw=2,label ='STAG_TECH')
ax1.plot(range(1,len(T1970_NR)+1),T1970_NR,'-',color="g",lw=0.2,alpha=1,label ='E1970')				
ax1.plot(range(1,len(Edg70GO_NR)+1),Edg70GO_NR,'-',color="k",lw=0.2,alpha=1,label ='STAG_G&O')	
ax1.plot(range(1,len(Edg70Oz_NR)+1),Edg70Oz_NR,'-',color="m",lw=0.2,alpha=1,label ='STAG_O3')	
ax1.plot(range(1,len(EdgRef_NR)+1),EdgRef_NR,'-',color="b",lw=0.2,alpha=1,label ='E2010')	
ax1.plot(range(1,len(EdgEne_NR)+1),EdgEne_NR,'-',color="r",lw=0.2,alpha=1,label ='STAG_ENE')	
ax1.plot(range(1,len(EdgTech_NR)+1),EdgTech_NR,'-',color="m",lw=0.2,alpha=1,label ='STAG_TECH')		
ax1.axhline(y=0,color='k',linewidth = 1);
ax1.set_ylim([-0.5,2.0]);ax1.set_xlim([0,221]);ax1.set_xticks(np.arange(0,221,20));
		
ax1 = plt.subplot(3,1,3);
ax1.annotate('(c) (FSNT-SLNT) vs. SAT',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
x =np.arange(0,2,0.1)			
ax1.plot(T1970_TS-  np.min(T1970_TS),T1970_NR,'.',color="g",ms=5)
# print T1970_s*x+T1970_i
ax1.plot(x,T1970_s*x+T1970_i,'-',color="g",lw=2,label ='E1970')

ax1.plot(Edg70GO_TS-np.min(Edg70GO_TS),Edg70GO_NR,'.',color="k",ms=5)
ax1.plot(x,Edg70GO_s*x+Edg70GO_i,'-',color="k",lw=2,label ='STAG_G&O')

ax1.plot(Edg70Oz_TS-np.min(Edg70Oz_TS),Edg70Oz_NR,'.',color="m",ms=5)
ax1.plot(x,Edg70Oz_s*x+Edg70Oz_i,'-',color="m",lw=2,label ='STAG_O3')

ax1.plot(EdgRef_TS-np.min(EdgRef_TS),EdgRef_NR,'.',color="b",ms=5)
ax1.plot(x,EdgRef_s*x+EdgRef_i,'-',color="b",lw=2,label ='E2010')

ax1.plot(EdgEne_TS-np.min(EdgEne_TS),EdgEne_NR,'.',color="r",ms=5)
ax1.plot(x,EdgEne_s*x+EdgEne_i,'-',color="r",lw=2,label ='STAG_ENE')

ax1.plot(EdgTech_TS-np.min(EdgTech_TS),EdgTech_NR,'.',color="m",ms=5)
ax1.plot(x,EdgTech_s*x+EdgTech_i,'-',color="m",lw=2,label ='STAG_TECH')
ax1.set_ylim([-1.0,2.0]);ax1.set_xlim([0,1.5])
# ax1.axvline(x=0,c='k')
ax1.axhline(y=0,color='k',linewidth = 1);

plt.subplots_adjust(left=0.10, bottom=0.03, right=0.95, top=0.95, wspace=0.20, hspace=0.15); 
plt.savefig('climate_sensitivity_NTOAR_vs_SAT.png', format='png', dpi=1000)
"""