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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

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
			
def data_read_zonal_mean(variable,mask):	
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
		
	def mon_mean2annual_mean(scenario,time,data):
		annual_mean=np.empty((30,192,288));annual_mean[:]=np.nan
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
		ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask'][:]
		if mask == "land":
			ocean_mask [ocean_mask ==0 ] = np.nan
		elif mask == "ocean":
			ocean_mask [ocean_mask ==1 ] = np.nan
			ocean_mask [ocean_mask ==0 ] = 1
		elif mask == "All":
			ocean_mask [ocean_mask ==0 ] = 1
		annual_mean = np.multiply(annual_mean,ocean_mask)
		return annual_mean

	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
		var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		if variable =="VQ" or variable == "VT":
			data = np.nanmean(nc_fid.variables[variable][:,23:30,:,:],axis=1)  # 850hpa
		else:
			data = nc_fid.variables[variable][:]#-273.15
		nc_fid.close()
		var40map = mon_mean2annual_mean(scenario,time,data)
		zonal_mean = np.nanmean(var40map,axis=2)
		return lat,zonal_mean
	
	lat,Edg70GO = data_netcdf('Edg70GO',variable)
	_,T1970 = data_netcdf('T1970RCP',variable)
	_,EdgRef = data_netcdf('EdgRef',variable)
	_,Edg70Oz = data_netcdf('Edg70Oz',variable)
	_,EdgEne = data_netcdf('EdgEne',variable)
	_,EdgTech = data_netcdf('EdgTech',variable)
	return lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	
def zonal_mean_unceer(variable1,variable2,mask):
	if variable2 =='none':
		if variable1 == "SFSC":   # surface shortwave cloud forcing
			index = 'FSNS'; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_read_zonal_mean(index,mask='All');
			index = 'FLNS'; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_read_zonal_mean(index,mask='All');
			T1970 = (T1970_S-T1970_l) 
			Edg70GO = (Edg70GO_S-Edg70GO_l)
			EdgRef =( EdgRef_S-EdgRef_l)
			EdgEne =( EdgEne_S-EdgEne_l)
			EdgTech =( EdgTech_S-EdgTech_l)
		elif variable1 == "SFLC":   # surface longwave cloud forcing
			index = 'FSNSC'; lat,T1970_SC,Edg70GO_SC,Edg70Oz_SC,EdgRef_SC,EdgEne_SC,EdgTech_SC = data_read_zonal_mean(index,mask='All');
			index = 'FLNSC'; lat,T1970_lC,Edg70GO_lC,Edg70Oz_lC,EdgRef_lC,EdgEne_lC,EdgTech_lC = data_read_zonal_mean(index,mask='All');
			T1970 =  -(T1970_SC-T1970_lC)
			Edg70GO = -(Edg70GO_SC-Edg70GO_lC)
			EdgRef =-(EdgRef_SC-EdgRef_lC)
			EdgEne = -(EdgEne_SC-EdgEne_lC)
			EdgTech =-(EdgTech_SC-EdgTech_lC)
		elif variable1 == "ATMABS":  # atmsopheric absorption = Ftoa -Tsrf
			index = 'FSNT'; lat,T1970_ST,Edg70GO_ST,Edg70Oz_ST,EdgRef_ST,EdgEne_ST,EdgTech_ST = data_read_zonal_mean(index,mask='All');
			index = 'FLNT'; lat,T1970_lT,Edg70GO_lT,Edg70Oz_lT,EdgRef_lT,EdgEne_lT,EdgTech_lT = data_read_zonal_mean(index,mask='All');			
			index = 'FSNS'; lat,T1970_SS,Edg70GO_SS,Edg70Oz_SS,EdgRef_SS,EdgEne_SS,EdgTech_SS = data_read_zonal_mean(index,mask='All');
			index = 'FLNS'; lat,T1970_lS,Edg70GO_lS,Edg70Oz_lS,EdgRef_lS,EdgEne_lS,EdgTech_lS = data_read_zonal_mean(index,mask='All');
			T1970 = (T1970_ST-T1970_lT) -(T1970_SS-T1970_lS)
			Edg70GO = (Edg70GO_ST-Edg70GO_lT)-(Edg70GO_SS-Edg70GO_lS)
			EdgRef =( EdgRef_ST-EdgRef_lT)-( EdgRef_SS-EdgRef_lS)
			EdgEne =( EdgEne_ST-EdgEne_lT)-( EdgEne_SS-EdgEne_lS)
			EdgTech =( EdgTech_ST-EdgTech_lT)-( EdgTech_SS-EdgTech_lS)
		else:
			index = variable1; lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_read_zonal_mean(index,mask='All');
	elif variable1 == "HT":  # heat transport
		lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L,EdgEne_L,EdgTech_L = data_read_zonal_mean('FLNT',mask='All');
		lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_read_zonal_mean("FSNT",mask='All');
		lat,T1970_LHFLX,Edg70GO_LHFLX,Edg70Oz_LHFLX,EdgRef_LHFLX,EdgEne_LHFLX,EdgTech_LHFLX = data_read_zonal_mean('LHFLX',mask='All');
		lat,T1970_SHFLX,Edg70GO_SHFLX,Edg70Oz_SHFLX,EdgRef_SHFLX,EdgEne_SHFLX,EdgTech_SHFLX = data_read_zonal_mean("SHFLX",mask='All');
		lat,T1970_FLNS,Edg70GO_FLNS,Edg70Oz_FLNS,EdgRef_FLNS,EdgEne_FLNS,EdgTech_FLNS = data_read_zonal_mean('FLNS',mask='All');
		lat,T1970_FSNS,Edg70GO_FSNS,Edg70Oz_FSNS,EdgRef_FSNS,EdgEne_FSNS,EdgTech_FSNS = data_read_zonal_mean("FSNS",mask='All');
		lat,T1970_PRECSC,Edg70GO_PRECSC,Edg70Oz_PRECSC,EdgRef_PRECSC,EdgEne_PRECSC,EdgTech_PRECSC = data_read_zonal_mean('PRECSC',mask='All');
		lat,T1970_PRECSL,Edg70GO_PRECSL,Edg70Oz_PRECSL,EdgRef_PRECSL,EdgEne_PRECSL,EdgTech_PRECSL = data_read_zonal_mean("PRECSL",mask='All');		
		lat,T1970_QFLX,Edg70GO_QFLX,Edg70Oz_QFLX,EdgRef_QFLX,EdgEne_QFLX,EdgTech_QFLX = data_read_zonal_mean('QFLX',mask='All');
		lat,T1970_PRECC,Edg70GO_PRECC,Edg70Oz_PRECC,EdgRef_PRECC,EdgEne_PRECC,EdgTech_PRECC = data_read_zonal_mean("PRECC",mask='All');
		lat,T1970_PRECL,Edg70GO_PRECL,Edg70Oz_PRECL,EdgRef_PRECL,EdgEne_PRECL,EdgTech_PRECL = data_read_zonal_mean("PRECL",mask='All');
		if mask=='All':     # total heatg transport
			T1970 =  inferred_heat_transport(T1970_S-T1970_L,lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L,lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L,lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L,lat)
			EdgEne=  inferred_heat_transport(EdgEne_S-EdgEne_L,lat)
			EdgTech= inferred_heat_transport(EdgTech_S-EdgTech_L,lat)
		elif mask=='ocean':   #ocean heat transpor
			T1970 =  inferred_heat_transport(-(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(-(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(-(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(-(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
			EdgEne=  inferred_heat_transport(-(EdgEne_FLNS-EdgEne_FSNS +EdgEne_LHFLX+EdgEne_SHFLX +(EdgEne_PRECSC+EdgEne_PRECSL)/2),lat)
			EdgTech= inferred_heat_transport(-(EdgTech_FLNS-EdgTech_FSNS +EdgTech_LHFLX+EdgTech_SHFLX +(EdgTech_PRECSC+EdgTech_PRECSL)/2),lat)			
		elif mask=='atm':   # Atmospheric heat transpor
			T1970 =  inferred_heat_transport(T1970_S-T1970_L+(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L+(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L+(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L+(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
			EdgEne=  inferred_heat_transport(EdgEne_S-EdgEne_L+(EdgEne_FLNS-EdgEne_FSNS +EdgEne_LHFLX+EdgEne_SHFLX +(EdgEne_PRECSC+EdgEne_PRECSL)/2),lat)	
			EdgTech= inferred_heat_transport(EdgTech_S-EdgTech_L+(EdgTech_FLNS-EdgTech_FSNS +EdgTech_LHFLX+EdgTech_SHFLX +(EdgTech_PRECSC+EdgTech_PRECSL)/2),lat)	
	elif variable1 == "CF_SFC":   # cloud-sky net radiative flux at surface
		index = 'FSNS'; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_read_zonal_mean(index,mask='All');
		index = 'FLNS'; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_read_zonal_mean(index,mask='All');
		index = 'FSNSC'; lat,T1970_SC,Edg70GO_SC,Edg70Oz_SC,EdgRef_SC,EdgEne_SC,EdgTech_SC = data_read_zonal_mean(index,mask='All');
		index = 'FLNSC'; lat,T1970_lC,Edg70GO_lC,Edg70Oz_lC,EdgRef_lC,EdgEne_lC,EdgTech_lC = data_read_zonal_mean(index,mask='All');
		T1970 = (T1970_S-T1970_l) -(T1970_SC-T1970_lC)
		Edg70GO = (Edg70GO_S-Edg70GO_l)-(Edg70GO_SC-Edg70GO_lC)
		EdgRef =( EdgRef_S-EdgRef_l)-(EdgRef_SC-EdgRef_lC)
		EdgEne =( EdgEne_S-EdgEne_l)-(EdgEne_SC-EdgEne_lC)
		EdgTech =( EdgTech_S-EdgTech_l)-(EdgTech_SC-EdgTech_lC)
	elif variable1 == "SWCF":    # FOR CLOUD FORCING, LWCF is defined positive dowanward, SWCF is defined positive dowanward
		index = variable1; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_read_zonal_mean(index,mask='All');
		index = variable2; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_read_zonal_mean(index,mask='All');
		T1970 = T1970_S+T1970_l;
		Edg70GO = Edg70GO_S+Edg70GO_l;
		EdgRef = EdgRef_S+EdgRef_l;
		EdgEne = EdgEne_S+EdgEne_l; 
		EdgTech = EdgTech_S+EdgTech_l;		
	else:
		index = variable1; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_read_zonal_mean(index,mask='All');
		index = variable2; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_read_zonal_mean(index,mask='All');
		T1970 = T1970_S-T1970_l;
		Edg70GO = Edg70GO_S-Edg70GO_l;
		EdgRef = EdgRef_S-EdgRef_l;
		EdgEne = EdgEne_S-EdgEne_l; 
		EdgTech = EdgTech_S-EdgTech_l;
	def decompose_mean_std(var1,var2):
		diff = var1 - var2
		# print "diff", np.shape(diff)
		mean = np.nanmean(diff,axis=0)
		P25 = mean - np.abs(mean-np.nanpercentile(diff,25,axis=0))/1.8
		P75 = mean + np.abs(np.nanpercentile(diff,75,axis=0)-mean)/1.8
		return mean,P25,P75
	
	mean_TAs,P25_TAs,P75_TAs = 	decompose_mean_std(EdgRef,T1970)
	mean_AAs,P25_AAs,P75_AAs = 	decompose_mean_std(Edg70GO,T1970)
	mean_Ene,P25_Ene,P75_Ene = 	decompose_mean_std(EdgRef,EdgEne)
	mean_Tech,P25_Tech,P75_Tech = decompose_mean_std(EdgRef,EdgTech)
	
	ref1970,_,_ = decompose_mean_std(T1970,0)
	
	TAs = np.concatenate((mean_TAs,P25_TAs,P75_TAs),axis=0)
	AAs = np.concatenate((mean_AAs,P25_AAs,P75_AAs),axis=0)
	Ene = np.concatenate((mean_Ene,P25_Ene,P75_Ene),axis=0)
	Tech = np.concatenate((mean_Tech,P25_Tech,P75_Tech),axis=0)
	return lat,AAs,Ene,Tech,ref1970

def plot_zonal_mean_uncertainty(ax,x,y1,y2,y3,shading):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
	ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
	ax.set_xticklabels(('90S','70S','50S','30S','10S','EQ','10N','30N','50N','70N','90N'));
	ax.set_ylim([-3,3.0]);ax.set_yticks(np.arange(-3,3.01,1))
	# ax.plot(x,y1[0:192],'-',color="k",linewidth=3,label = 'Historical (2010-1970)')
	ax.plot(x,y1[0:192],'-',color="b",linewidth=3,label = 'Best Estimation')
	ax.plot(x,y2[0:192],'-',color="r",linewidth=3,label = 'Energy Consumption')
	ax.plot(x,y3[0:192],'-',color="g",linewidth=3,label = 'Technology Advancements')
	
	# ax.fill_between(x,y1[192:384],y1[384:576],facecolor="k",alpha=0.15)
	ax.fill_between(x,y1[192:384],y1[384:576],facecolor="b",alpha=0.15)
	ax.fill_between(x,y2[192:384],y2[384:576],facecolor="r",alpha=0.15)
	ax.fill_between(x,y3[192:384],y3[384:576],facecolor="g",alpha=0.15)
	if shading == 'midlat':ax.axvspan(-30, 60, alpha=0.2, color='grey',edgecolor='none')
	elif shading == 'Arctic':ax.axvspan(60, 90, alpha=0.2, color='grey',edgecolor='none')
	ax2 = ax.twinx();
	# # ax2.plot(x,y4,'-',color="k",linewidth=3, label = "Historical (2010-1970")
	# # ax2.plot(x,ref,'-',color="m",linewidth=4,label = 'E1970')
	# # ax2.fill_between(x,ref[192:384],ref[384:576],facecolor="k",alpha=0.15)
	ax2.set_xlim([-90,90]);ax2.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
	ax2.set_ylim([-3,3]);ax2.set_yticks(np.arange(-3,3.01,1))
	return ax
	
	
def plot_zonal_mean_uncertainty2(ax,x,y1,y2,y3,ref):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
	ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
	ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	# ax.set_ylim([-4,4.0]);ax.set_yticks(np.arange(-4,4.01,2))
	# ax.plot(x,y1[0:192],'-',color="k",linewidth=3,label = 'Historical (2010-1970)')
	ax.plot(x,y1[0:192],'-',color="b",linewidth=3)
	ax.plot(x,y2[0:192],'-',color="r",linewidth=3)
	ax.plot(x,y3[0:192],'-',color="g",linewidth=3)
	
	# ax.fill_between(x,y1[192:384],y1[384:576],facecolor="k",alpha=0.15)
	ax.fill_between(x,y1[192:384],y1[384:576],facecolor="b",alpha=0.15)
	ax.fill_between(x,y2[192:384],y2[384:576],facecolor="r",alpha=0.15)
	ax.fill_between(x,y3[192:384],y3[384:576],facecolor="g",alpha=0.15)
	
	ax2 = ax.twinx();
	# ax2.plot(x,y4,'-',color="k",linewidth=3, label = "Historical (2010-1970")
	ax2.plot(x,ref,'-',color="k",linewidth=4,label = 'E1970')
	ax2.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
	return ax,ax2


######################################################
## radiative fluxes at TOA,surface and atmospheric absorption
######################################################
"""
fig = plt.figure(facecolor='White',figsize=[20,16]);pad= 5   #;plot_setup();
colormap='RdBu';colormap = reverse_colourmap(colormap);
	
index1 ='FSNT'; index2='FLNT'; lat, AAs, Ene,Tech,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(3,3,1);
ax.annotate('(a) TOA net radiative flux '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)
ax= plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,shading='none')

index1 ='ATMABS'; index2='none'; lat, AAs, Ene,Tech,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(3,3,2);
ax.annotate('(b) Atmospheric absorption '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,shading='none')
legend = ax.legend(shadow=False,ncol=1,bbox_to_anchor=(2.2, 0.7),fontsize=25)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

index1 ='FSNS'; index2='FLNS'; lat, AAs, Ene,Tech,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(3,3,4);
ax.annotate('(c) Surface net radiative flux '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,shading='none')

index1 ='FSNSC'; index2='FLNSC'; lat, AAs, Ene,Tech,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(3,3,5);
ax.annotate('(d) Surface net clear-sky '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,shading='Arctic')


index1 ='CF_SFC'; index2='CF_SFC'; lat, AAs, Ene,Tech,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(3,3,6);
ax.annotate('(e) cloud forcing '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,shading='midlat')


###########################
##heat transport
###########################

ax = plt.subplot(3,3,7);
index1 ='HT'; index2='HT'	;lat, AAs, Ene,Tech,ref = zonal_mean_unceer(index1,index2,mask = 'All');

ax.annotate('(f) Total heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)									
ax5,ax6= plot_zonal_mean_uncertainty2(ax,lat, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.15,0.15]);ax5.set_yticks(np.arange(-0.15,0.16,0.05));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

ax = plt.subplot(3,3,8);
index1 ='HT'; index2='HT'	;lat, AAs, Ene,Tech,ref = zonal_mean_unceer(index1,index2,mask = 'atm');

ax.annotate('(g) Atmospheric heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)									
ax5,ax6= plot_zonal_mean_uncertainty2(ax,lat, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.15,0.15]);ax5.set_yticks(np.arange(-0.15,0.16,0.05));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)


ax = plt.subplot(3,3,9);
index1 ='HT'; index2='HT'	;lat, AAs, Ene,Tech,ref = zonal_mean_unceer(index1,index2,mask = 'ocean');

ax.annotate('(h) Oceanic heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)									
ax5,ax6= plot_zonal_mean_uncertainty2(ax,lat, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.15,0.15]);ax5.set_yticks(np.arange(-0.15,0.16,0.05));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)


plt.subplots_adjust(left=0.05, bottom=0.04, right=0.96, top=0.96, wspace=0.20, hspace=0.20); 
plt.savefig('RF_HT.png', format='png', dpi=300)


"""
######################################################
## surface radiative fluxes (also LW and SW) and Heat transport
######################################################

fig = plt.figure(facecolor='White',figsize=[20,6]);pad= 5   #;plot_setup();
colormap='RdBu';colormap = reverse_colourmap(colormap);
	

ax = plt.subplot(1,3,1);
# index1 ='FSNS'; index2='none';lat, AAs_FSNS, Ene_FSNS,Tech_FSNS,_ = zonal_mean_unceer(index1,index2,mask = 'All');
# index1 ='FLNS';index2='none'; lat, AAs_FLNS, Ene_FLNS,Tech_FLNS,_ = zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='FSNS';index2='none'; lat, AAs_SFSC, Ene_SFSC,Tech_SFSC,_ = zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='FLNS';index2='none'; lat, AAs_SFLC, Ene_SFLC,Tech_SFLC,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
ax.set_ylim([-5,5.0]);ax.set_yticks(np.arange(-5,5.01,1))
ax.annotate('(a) Surfae all-sky '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.plot(lat,AAs_SFSC[0:192],'-',color="b",linewidth=3,label = 'BEoA (2010 - 1970) SW')
ax.plot(lat,Ene_SFSC[0:192],'-',color="r",linewidth=3,label = 'Energy Comsumption SW')
ax.plot(lat,Tech_SFSC[0:192],'-',color="g",linewidth=3,label = 'Technology advancements SW')
ax.plot(lat,-AAs_SFLC[0:192],'--',color="b",linewidth=3,label = 'LW')
ax.plot(lat,-Ene_SFLC[0:192],'--',color="r",linewidth=3,label = 'LW')
ax.plot(lat,-Tech_SFLC[0:192],'--',color="g",linewidth=3,label = 'LW')
ax.axvspan(60, 90, alpha=0.2, color='grey',edgecolor='none')
ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
x=np.arange(-90,90); y =x/90
ax.plot(x,y,'-',color="w",linewidth=0)
ax.set_xlim([-90,90]);ax.set_ylim([-5,5]);
line_s = mlines.Line2D(x,y,ls='-',color='k',lw=2)
line_d = mlines.Line2D(x,y,ls='--',color='k',lw=2)
lines = [line_s,line_d]
labels = ['Shortwave','Longwave']
legend = plt.legend(lines,labels,ncol=1,loc ='lower left',labelspacing=1.5,markerscale =10)
legend.get_frame().set_facecolor('none');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)


ax2 = ax.twinx();
ax2.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax2.set_ylim([-5,5]);ax2.set_yticks(np.arange(-5,5.01,1))

ax = plt.subplot(1,3,2);
# index1 ='FSNS'; index2='none';lat, AAs_FSNS, Ene_FSNS,Tech_FSNS,_ = zonal_mean_unceer(index1,index2,mask = 'All');
# index1 ='FLNS';index2='none'; lat, AAs_FLNS, Ene_FLNS,Tech_FLNS,_ = zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='FSNSC';index2='none'; lat, AAs_SFSC, Ene_SFSC,Tech_SFSC,_ = zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='FLNSC';index2='none'; lat, AAs_SFLC, Ene_SFLC,Tech_SFLC,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
ax.set_ylim([-5,5.0]);ax.set_yticks(np.arange(-5,5.01,1))
ax.annotate('(b) Surfae clear-sky '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.plot(lat,AAs_SFSC[0:192],'-',color="b",linewidth=3,label = 'BEoA (2010 - 1970) SW')
ax.plot(lat,Ene_SFSC[0:192],'-',color="r",linewidth=3,label = 'Energy Comsumption SW')
ax.plot(lat,Tech_SFSC[0:192],'-',color="g",linewidth=3,label = 'Technology advancements SW')
ax.plot(lat,-AAs_SFLC[0:192],'--',color="b",linewidth=3,label = 'LW')
ax.plot(lat,-Ene_SFLC[0:192],'--',color="r",linewidth=3,label = 'LW')
ax.plot(lat,-Tech_SFLC[0:192],'--',color="g",linewidth=3,label = 'LW')
ax.axvspan(60, 90, alpha=0.2, color='grey',edgecolor='none')
ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
x=np.arange(-90,90); y =x/90
ax.plot(x,y,'-',color="w",linewidth=0)
ax.set_xlim([-90,90]);ax.set_ylim([-5,5]);
line_b = mlines.Line2D(x,y,ls='-',color='b',lw=2)
line_r = mlines.Line2D(x,y,ls='-',color='r',lw=2)
line_g = mlines.Line2D(x,y,ls='-',color='g',lw=2)
lines = [line_b,line_r,line_g]
labels = ['Best Estimation','Energy Consumption','Technology Advancements']
legend = plt.legend(lines,labels,ncol=1,loc ='lower left',labelspacing=1.5,markerscale =10)
legend.get_frame().set_facecolor('none');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)

ax2 = ax.twinx();
ax2.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax2.set_ylim([-5,5]);ax2.set_yticks(np.arange(-5,5.01,1))

ax = plt.subplot(1,3,3);
index1 ='SFSC';index2='none'; lat, AAs_SFSC, Ene_SFSC,Tech_SFSC,_ = zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='SFLC';index2='none'; lat, AAs_SFLC, Ene_SFLC,Tech_SFLC,_ = zonal_mean_unceer(index1,index2,mask = 'All');

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
ax.set_ylim([-3,3.0]);ax.set_yticks(np.arange(-3,3.01,1))
ax.annotate('(c) Surface cloud forcing '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.plot(lat,AAs_SFSC[0:192],'-',color="b",linewidth=3,label = 'BEoA SW')
ax.plot(lat,Ene_SFSC[0:192],'-',color="r",linewidth=3,label = 'ENE SW')
ax.plot(lat,Tech_SFSC[0:192],'-',color="g",linewidth=3,label = 'TECH SW')
ax.plot(lat,AAs_SFLC[0:192],'--',color="b",linewidth=3,label = 'BEoA LW')
ax.plot(lat,Ene_SFLC[0:192],'--',color="r",linewidth=3,label = 'ENE LW')
ax.plot(lat,Tech_SFLC[0:192],'--',color="g",linewidth=3,label = 'TECH LW')
ax.axvspan(-30, 60, alpha=0.2, color='grey',edgecolor='none')
ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);

x=np.arange(-90,90); y =x/90
ax.plot(x,y,'-',color="w",linewidth=0)
ax.set_xlim([-90,90]);ax.set_ylim([-3,3]);

ax2 = ax.twinx();
ax2.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax2.set_ylim([-3,3]);ax2.set_yticks(np.arange(-3,3.01,1))

plt.subplots_adjust(left=0.04, bottom=0.1, right=0.96, top=0.90, wspace=0.20, hspace=0.20); 
plt.savefig('zonal_mean_SW_VS_LW.png', format='png', dpi=300)

