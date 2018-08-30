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
    os.path.pardir,os.path.pardir,os.path.pardir,
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
			# print scenario,iyear
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
	_,Edg70T10SOZ = data_netcdf('Edg70T10SOZ',variable)
	return lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ
	
def zonal_mean_unceer(variable1,variable2,mask):
	if variable2 =='none':
		if variable1 == "SFSC":   # surface shortwave cloud forcing
			index = 'FSNS'; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,Edg70T10SOZ_S = data_read_zonal_mean(index,mask='All');
			index = 'FLNS'; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,Edg70T10SOZ_l = data_read_zonal_mean(index,mask='All');
			T1970 = (T1970_S-T1970_l) 
			Edg70GO = (Edg70GO_S-Edg70GO_l)
			EdgRef =( EdgRef_S-EdgRef_l)
			Edg70Oz =( Edg70Oz_S-Edg70Oz_l)
			Edg70T10SOZ =( Edg70T10SOZ_S-Edg70T10SOZ_l)

		elif variable1 == "SFLC":   # surface longwave cloud forcing
			index = 'FSNSC'; lat,T1970_SC,Edg70GO_SC,Edg70Oz_SC,EdgRef_SC,Edg70T10SOZ_SC  = data_read_zonal_mean(index,mask='All');
			index = 'FLNSC'; lat,T1970_LC,Edg70GO_LC,Edg70Oz_LC,EdgRef_LC,Edg70T10SOZ_LC  = data_read_zonal_mean(index,mask='All');
			T1970 =  -(T1970_SC-T1970_lC)
			Edg70GO = -(Edg70GO_SC-Edg70GO_lC)
			EdgRef =-(EdgRef_SC-EdgRef_lC)
			Edg70T10SOZ = -(Edg70T10SOZ_SC-Edg70T10SOZ_lC)
			Edg70Oz =( Edg70Oz_SC-Edg70Oz_lC)
		else:
			index = variable1; lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ = data_read_zonal_mean(index,mask='All');
	elif variable1 == "HT":  # heat transport
		lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L,Edg70T10SOZ_L = data_read_zonal_mean('FLNT',mask='All');
		lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,Edg70T10SOZ_S = data_read_zonal_mean("FSNT",mask='All');
		lat,T1970_LHFLX,Edg70GO_LHFLX,Edg70Oz_LHFLX,EdgRef_LHFLX,Edg70T10SOZ_LHFLX = data_read_zonal_mean('LHFLX',mask='All');
		lat,T1970_SHFLX,Edg70GO_SHFLX,Edg70Oz_SHFLX,EdgRef_SHFLX,Edg70T10SOZ_SHFLX = data_read_zonal_mean("SHFLX",mask='All');
		lat,T1970_FLNS,Edg70GO_FLNS,Edg70Oz_FLNS,EdgRef_FLNS,Edg70T10SOZ_FLNS = data_read_zonal_mean('FLNS',mask='All');
		lat,T1970_FSNS,Edg70GO_FSNS,Edg70Oz_FSNS,EdgRef_FSNS,Edg70T10SOZ_FSNS = data_read_zonal_mean("FSNS",mask='All');
		lat,T1970_PRECSC,Edg70GO_PRECSC,Edg70Oz_PRECSC,EdgRef_PRECSC,Edg70T10SOZ_PRECSC = data_read_zonal_mean('PRECSC',mask='All');
		lat,T1970_PRECSL,Edg70GO_PRECSL,Edg70Oz_PRECSL,EdgRef_PRECSL,Edg70T10SOZ_PRECSL = data_read_zonal_mean("PRECSL",mask='All');		
		lat,T1970_QFLX,Edg70GO_QFLX,Edg70Oz_QFLX,EdgRef_QFLX,Edg70T10SOZ_QFLX = data_read_zonal_mean('QFLX',mask='All');
		lat,T1970_PRECC,Edg70GO_PRECC,Edg70Oz_PRECC,EdgRef_PRECC,Edg70T10SOZ_PRECC = data_read_zonal_mean("PRECC",mask='All');
		lat,T1970_PRECL,Edg70GO_PRECL,Edg70Oz_PRECL,EdgRef_PRECL,Edg70T10SOZ_PRECL = data_read_zonal_mean("PRECL",mask='All');
		if mask=='All':     # total heatg transport
			T1970 =  inferred_heat_transport(T1970_S-T1970_L,lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L,lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L,lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L,lat)
			Edg70T10SOZ=  inferred_heat_transport(Edg70T10SOZ_S-Edg70T10SOZ_L,lat)
		elif mask=='ocean':   #ocean heat transpor
			T1970 =  inferred_heat_transport(-(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(-(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(-(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(-(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
			Edg70T10SOZ=  inferred_heat_transport(-(Edg70T10SOZ_FLNS-Edg70T10SOZ_FSNS +Edg70T10SOZ_LHFLX+Edg70T10SOZ_SHFLX +(Edg70T10SOZ_PRECSC+Edg70T10SOZ_PRECSL)/2),lat)			
		elif mask=='atm':   # Atmospheric heat transpor
			T1970 =  inferred_heat_transport(T1970_S-T1970_L+(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L+(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L+(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L+(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
			Edg70T10SOZ=  inferred_heat_transport(Edg70T10SOZ_S-Edg70T10SOZ_L+(Edg70T10SOZ_FLNS-Edg70T10SOZ_FSNS +Edg70T10SOZ_LHFLX+Edg70T10SOZ_SHFLX +(Edg70T10SOZ_PRECSC+Edg70T10SOZ_PRECSL)/2),lat)	
	elif variable1 == "SFCd":   # cloud-sky net radiative flux at surface
		index = 'FSNS'; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,Edg70T10SOZ_S = data_read_zonal_mean(index,mask='All');
		index = 'FLNS'; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,Edg70T10SOZ_l = data_read_zonal_mean(index,mask='All');
		index = 'FSNSC'; lat,T1970_SC,Edg70GO_SC,Edg70Oz_SC,EdgRef_SC,Edg70T10SOZ_SC = data_read_zonal_mean(index,mask='All');
		index = 'FLNSC'; lat,T1970_lC,Edg70GO_lC,Edg70Oz_lC,EdgRef_lC,Edg70T10SOZ_lC = data_read_zonal_mean(index,mask='All');
		T1970 = (T1970_S-T1970_l) -(T1970_SC-T1970_lC)
		Edg70GO = (Edg70GO_S-Edg70GO_l)-(Edg70GO_SC-Edg70GO_lC)
		EdgRef =( EdgRef_S-EdgRef_l)-(EdgRef_SC-EdgRef_lC)
		Edg70T10SOZ =( Edg70T10SOZ_S-Edg70T10SOZ_l)-(Edg70T10SOZ_SC-Edg70T10SOZ_lC)
		Edg70Oz =( Edg70Oz_S-Edg70Oz_l)-(Edg70Oz_SC-Edg70Oz_lC)
	elif variable1 == "SWCF":    # FOR CLOUD FORCING, LWCF is defined positive dowanward, SWCF is defined positive dowanward
		index = variable1; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,Edg70T10SOZ_S = data_read_zonal_mean(index,mask='All');
		index = variable2; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,Edg70T10SOZ_l = data_read_zonal_mean(index,mask='All');
		T1970 = T1970_S+T1970_l;
		Edg70GO = Edg70GO_S+Edg70GO_l;
		EdgRef = EdgRef_S+EdgRef_l;
		Edg70T10SOZ = Edg70T10SOZ_S+Edg70T10SOZ_l; 
		Edg70Oz = Edg70Oz_S+Edg70Oz_l;		
	else:
		index = variable1; lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,Edg70T10SOZ_S = data_read_zonal_mean(index,mask='All');
		index = variable2; lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,Edg70T10SOZ_l = data_read_zonal_mean(index,mask='All');
		T1970 = T1970_S-T1970_l;
		Edg70GO = Edg70GO_S-Edg70GO_l;
		EdgRef = EdgRef_S-EdgRef_l;
		Edg70T10SOZ = Edg70T10SOZ_S-Edg70T10SOZ_l; 
		Edg70Oz = Edg70Oz_S-Edg70Oz_l;
	def decompose_mean_std(var1,var2):
		diff = var1 - var2
		# print "diff", np.shape(diff)
		mean = np.nanmean(diff,axis=0)
		P25 = mean - np.abs(mean-np.nanpercentile(diff,25,axis=0))/1.8
		P75 = mean + np.abs(np.nanpercentile(diff,75,axis=0)-mean)/1.8
		return mean,P25,P75
	
	mean_TAs,P25_TAs,P75_TAs = 	decompose_mean_std(EdgRef,T1970)
	mean_AAs,P25_AAs,P75_AAs = 	decompose_mean_std(Edg70GO,T1970)
	mean_GHG,P25_GHG,P75_GHG = 	decompose_mean_std(Edg70Oz,Edg70GO)
	mean_TrO3,P25_TrO3,P75_TrO3 = 	decompose_mean_std(EdgRef,Edg70T10SOZ)
	mean_StO3,P25_StO3,P75_StO3 = 	decompose_mean_std(Edg70T10SOZ,Edg70Oz)
	ref1970,_,_ = decompose_mean_std(T1970,0)
	
	TAs = np.concatenate((mean_TAs,P25_TAs,P75_TAs),axis=0)
	AAs = np.concatenate((mean_AAs,P25_AAs,P75_AAs),axis=0)
	GHGs = np.concatenate((mean_GHG,P25_GHG,P75_GHG),axis=0)
	TrO3 = np.concatenate((mean_TrO3,P25_TrO3,P75_TrO3),axis=0)
	StO3 = np.concatenate((mean_StO3,P25_StO3,P75_StO3),axis=0)
	return lat,TAs,GHGs,AAs,TrO3,StO3,ref1970

# def plot_zonal_mean_uncertainty(ax,x,y1,y2,y3,y4,shading):
	# ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	# ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	# ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
	# ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
	# ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
	# ax.set_xticklabels(('90S','670S','30S','EQ','30N','60N','90N'));
	# ax.set_ylim([-5,5.0]);ax.set_yticks(np.arange(-5,5.01,2.5))
	# # ax.plot(x,y1[0:192],'-',color="k",linewidth=3,label = 'Historical (2010-1970)')
	# ax.plot(x,y1[0:192],'-',color="g",linewidth=3,label = 'GHGs')
	# ax.plot(x,y2[0:192],'-',color="r",linewidth=3,label = 'AAs')
	# ax.plot(x,y3[0:192],'-',color="m",linewidth=3,label = 'Trop. O3')
	# ax.plot(x,y4[0:192],'-',color="b",linewidth=3,label = 'Str. O3')
	# # ax.fill_between(x,y1[192:384],y1[384:576],facecolor="k",alpha=0.15)
	# # ax.fill_between(x,y1[192:384],y1[384:576],facecolor="g",alpha=0.15)
	# # ax.fill_between(x,y2[192:384],y2[384:576],facecolor="r",alpha=0.15)
	# # ax.fill_between(x,y3[192:384],y3[384:576],facecolor="m",alpha=0.15)
	# # ax.fill_between(x,y4[192:384],y4[384:576],facecolor="b",alpha=0.15)
	# if shading == 'midlat':ax.axvspan(-30, 60, alpha=0.2, color='grey',edgecolor='none')
	# elif shading == 'Arctic':ax.axvspan(60, 90, alpha=0.2, color='grey',edgecolor='none')
	# return ax
	
	
def plot_zonal_mean_uncertainty(ax,x,y0,y1,y2,y3,y4,ref):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
	ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
	ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	
	ax.plot(x,y0[0:192],'-',color="k",linewidth=3,label = 'All')
	ax.plot(x,y1[0:192],'-',color="g",linewidth=3,label = 'GHGs')
	ax.plot(x,y2[0:192],'-',color="r",linewidth=3,label = 'AAs')
	ax.plot(x,y3[0:192],'-',color="m",linewidth=3,label = 'Trop. O3')
	ax.plot(x,y4[0:192],'-',color="b",linewidth=3,label = 'Strat. O3')
	ax2 = ax.twinx();
	ax2.plot(x,ref,'--',color="k",linewidth=3,label = 'E1970')
	ax2.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
	return ax,ax2


fig = plt.figure(facecolor='White',figsize=[12,5]);pad= 5   #;plot_setup();
colormap='RdBu';colormap = reverse_colourmap(colormap);
#########################
#### LW versus SW    ####
#########################
"""
ax = plt.subplot(1,2,1);
index1 ='FSNT';index2='none'; lat,_,GHGs_SFSC,AAs_SFSC,TrO3_SFSC,StO3_SFSC,_= zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='FLNT';index2='none'; lat,_,GHGs_SFLC,AAs_SFLC,TrO3_SFLC,StO3_SFLC,_ = zonal_mean_unceer(index1,index2,mask = 'All');
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
ax.set_ylim([-6,6.0]);ax.set_yticks(np.arange(-6,6.01,2))
ax.annotate('(a) TOA SW and LW all-sky '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.plot(lat,GHGs_SFSC[0:192],'-',color="g",linewidth=1)
ax.plot(lat,AAs_SFSC[0:192],'-',color="r",linewidth=1)
ax.plot(lat,TrO3_SFSC[0:192],'-',color="m",linewidth=1)
ax.plot(lat,StO3_SFSC[0:192],'-',color="b",linewidth=1)
ax.plot(lat,-GHGs_SFLC[0:192],'--',color="g",linewidth=1)
ax.plot(lat,-AAs_SFLC[0:192],'--',color="r",linewidth=1)
ax.plot(lat,-TrO3_SFLC[0:192],'--',color="m",linewidth=1)
ax.plot(lat,-StO3_SFLC[0:192],'--',color="b",linewidth=1)

x=np.arange(-90,90); y =x/90
ax.plot(x,y,'-',color="w",linewidth=0)
ax.set_xlim([-90,90]);ax.set_ylim([-8,8]);
line_s = mlines.Line2D(x,y,ls='-',color='k',lw=2)
line_d = mlines.Line2D(x,y,ls='--',color='k',lw=2)
lines = [line_s,line_d]
labels = ['Shortwave','Longwave']
legend = plt.legend(lines,labels,ncol=1,loc ='lower left',labelspacing=1.5,markerscale =10)
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)

ax = plt.subplot(1,2,2);
index1 ='FSNS';index2='none'; lat,_,GHGs_SFSC,AAs_SFSC,TrO3_SFSC,StO3_SFSC,_= zonal_mean_unceer(index1,index2,mask = 'All');
index1 ='FLNS';index2='none'; lat,_,GHGs_SFLC,AAs_SFLC,TrO3_SFLC,StO3_SFLC,_= zonal_mean_unceer(index1,index2,mask = 'All');
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
ax.set_ylim([-6,6.0]);ax.set_yticks(np.arange(-6,6.01,2))
ax.annotate('(b) Surface SW and LW all-sky '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.plot(lat,GHGs_SFSC[0:192],'-',color="g",linewidth=1)
ax.plot(lat,AAs_SFSC[0:192],'-',color="r",linewidth=1)
ax.plot(lat,TrO3_SFSC[0:192],'-',color="m",linewidth=1)
ax.plot(lat,StO3_SFSC[0:192],'-',color="b",linewidth=1)
ax.plot(lat,-GHGs_SFLC[0:192],'--',color="g",linewidth=1)
ax.plot(lat,-AAs_SFLC[0:192],'--',color="r",linewidth=1)
ax.plot(lat,-TrO3_SFLC[0:192],'--',color="m",linewidth=1)
ax.plot(lat,-StO3_SFLC[0:192],'--',color="b",linewidth=1)
plt.subplots_adjust(left=0.05, bottom=0.08, right=0.95, top=0.90, wspace=0.20, hspace=0.30); 
plt.savefig('zonal_mean_SW_LW.png', format='png', dpi=300)





fig = plt.figure(facecolor='White',figsize=[12,8]);pad= 5   #;plot_setup();
colormap='RdBu';colormap = reverse_colourmap(colormap);
index1 ='FSNT'; index2='FLNT'; lat,TAs,GHGs,AAs,TrO3,StO3,ref1970= zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(2,2,1);
ax.annotate('(a) TOA net (SW-LW) all-sky radiation '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
ax5,ax6= plot_zonal_mean_uncertainty(ax,lat, TAs,GHGs,AAs,TrO3,StO3,ref1970)

index1 ='FSNS'; index2='FLNS'; lat,TAs,GHGs,AAs,TrO3,StO3,ref1970= zonal_mean_unceer(index1,index2,mask = 'All');
ax = plt.subplot(2,2,2);
ax.annotate('(b) Surface net (SW-LW) all-sky radiation '+r'$(W\/m^{-2})$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
ax5,ax6= plot_zonal_mean_uncertainty(ax,lat, TAs,GHGs,AAs,TrO3,StO3,ref1970)
# legend = ax.legend(shadow=False,ncol=1,loc ='lower left',fontsize=15)	 
# legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)
"""	
fig = plt.figure(facecolor='White',figsize=[7,10]);pad= 5   #;plot_setup();
colormap='RdBu';colormap = reverse_colourmap(colormap);

###########################
##heat transport
###########################

ax = plt.subplot(3,1,1);
index1 ='HT'; index2='HT'	;lat,TAs,GHGs,AAs,TrO3,StO3,ref1970= zonal_mean_unceer(index1,index2,mask = 'All');

ax.annotate('(a) Total heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)									
ax5,ax6= plot_zonal_mean_uncertainty(ax,lat, TAs,GHGs,AAs,TrO3,StO3,ref1970)
ax5.set_ylim([-0.16,0.16]);ax5.set_yticks(np.arange(-0.16,0.17,0.08));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

ax = plt.subplot(3,1,2);
index1 ='HT'; index2='HT'	;lat,TAs,GHGs,AAs,TrO3,StO3,ref1970= zonal_mean_unceer(index1,index2,mask = 'atm');

ax.annotate('(b) Atmospheric heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)									
ax5,ax6= plot_zonal_mean_uncertainty(ax,lat, TAs,GHGs,AAs,TrO3,StO3,ref1970)

ax5.set_ylim([-0.16,0.16]);ax5.set_yticks(np.arange(-0.16,0.17,0.08));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

ax = plt.subplot(3,1,3);
index1 ='HT'; index2='HT'	;lat,TAs,GHGs,AAs,TrO3,StO3,ref1970= zonal_mean_unceer(index1,index2,mask = 'ocean');

ax.annotate('(c) Oceanic heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)									
ax5,ax6= plot_zonal_mean_uncertainty(ax,lat, TAs,GHGs,AAs,TrO3,StO3,ref1970)
legend = ax.legend(shadow=False,ncol=3,loc ='lower left',fontsize=15)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

ax5.set_ylim([-0.16,0.16]);ax5.set_yticks(np.arange(-0.16,0.17,0.08));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)


plt.subplots_adjust(left=0.10, bottom=0.08, right=0.93, top=0.94, wspace=0.20, hspace=0.30); 
plt.savefig('Meridioanl_HT.png', format='png', dpi=300)


