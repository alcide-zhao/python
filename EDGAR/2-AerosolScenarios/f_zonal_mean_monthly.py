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
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
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
	
def zonal_mean_unceer(variable,mask):
	if variable == "precipitation":	
		lat,T1970_C,Edg70GO_C,Edg70Oz_C,EdgRef_C,EdgEne_C,EdgTech_C = data_read_zonal_mean('PRECL',mask);
		lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L,EdgEne_L,EdgTech_L = data_read_zonal_mean('PRECC',mask);
		T1970=(T1970_C+T1970_L)*24*60*60*1000  # m/s to mm/day
		Edg70GO=(Edg70GO_C+Edg70GO_L)*24*60*60*1000 
		Edg70Oz=(Edg70Oz_C+Edg70Oz_L)*24*60*60*1000 
		EdgRef=(EdgRef_C+EdgRef_L)*24*60*60*1000 
		EdgEne=(EdgEne_C+EdgEne_L)*24*60*60*1000 
		EdgTech=(EdgTech_C+EdgTech_L)*24*60*60*1000 
	elif variable == "TS":	
		# print variable,mask
		lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_read_zonal_mean(variable,mask); #,EdgEne,EdgTech
		T1970 = T1970-273.15
		Edg70GO = Edg70GO-273.15
		Edg70Oz = Edg70Oz-273.15
		EdgRef = EdgRef-273.15
		EdgEne = EdgEne-273.15
		EdgTech = EdgTech-273.15
	elif variable == "HT":  # heat transport
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
			EdgEne= inferred_heat_transport(EdgEne_S-EdgEne_L,lat)
			EdgTech=  inferred_heat_transport(EdgTech_S-EdgTech_L,lat)			
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
		elif mask=='latent':   # tmospheric latent heat transport,
			Lhvap = 2.5E6; rho_w = 1000.   # Latent heat of vaporization (J / kg)
			T1970 =  inferred_heat_transport((T1970_QFLX-rho_w*(T1970_PRECC+T1970_PRECL))*Lhvap,lat)  
			Edg70GO= inferred_heat_transport((Edg70GO_QFLX-rho_w*(Edg70GO_PRECC+Edg70GO_PRECL))*Lhvap,lat)
			Edg70Oz= inferred_heat_transport((Edg70Oz_QFLX-rho_w*(Edg70Oz_PRECC+Edg70Oz_PRECL))*Lhvap,lat)
			EdgRef=  inferred_heat_transport((EdgRef_QFLX-rho_w*(EdgRef_PRECC+EdgRef_PRECL))*Lhvap,lat)
			EdgEne=  inferred_heat_transport((EdgEne_QFLX-rho_w*(EdgEne_PRECC+EdgEne_PRECL))*Lhvap,lat)
			EdgTech= inferred_heat_transport((EdgTech_QFLX-rho_w*(EdgTech_PRECC+EdgTech_PRECL))*Lhvap,lat)			
		elif mask=='sensible':   # Atmospheric heat transpor
			Lhvap = 2.5E6 ;rho_w = 1000.  # Latent heat of vaporization (J / kg)
			T1970 =  inferred_heat_transport(T1970_S-T1970_L+(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2)-(T1970_QFLX-rho_w*(T1970_PRECC+T1970_PRECL))*Lhvap,lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L+(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2)-(Edg70GO_QFLX-rho_w*(Edg70GO_PRECC+Edg70GO_PRECL))*Lhvap,lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L+(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2)-(Edg70Oz_QFLX-rho_w*(Edg70Oz_PRECC+Edg70Oz_PRECL))*Lhvap,lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L+(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2)-(EdgRef_QFLX-rho_w*(EdgRef_PRECC+EdgRef_PRECL))*Lhvap,lat)			
			EdgEne=  inferred_heat_transport(EdgEne_S-EdgEne_L+(EdgEne_FLNS-EdgEne_FSNS +EdgEne_LHFLX+EdgEne_SHFLX +(EdgEne_PRECSC+EdgEne_PRECSL)/2)-(EdgEne_QFLX-rho_w*(EdgEne_PRECC+EdgEne_PRECL))*Lhvap,lat)	
			EdgTech= inferred_heat_transport(EdgTech_S-EdgTech_L+(EdgTech_FLNS-EdgTech_FSNS +EdgTech_LHFLX+EdgTech_SHFLX +(EdgTech_PRECSC+EdgTech_PRECSL)/2)-(EdgTech_QFLX-rho_w*(EdgTech_PRECC+EdgTech_PRECL))*Lhvap,lat)	
	else:
		lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_read_zonal_mean(variable,mask); #,EdgEne,EdgTech
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
	return lat, TAs,AAs, Ene,Tech,ref1970

def plot_zonal_mean_uncertainty(ax,x,y1,y2,y3,y4,ref):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
	ax.set_xticklabels(('90S','70S','50S','30S','10S','EQ','10N','30N','50N','70S','90N'));
	# ax.set_ylim([-4,4.0]);ax.set_yticks(np.arange(-4,4.01,2))
	ax.plot(x,y1[0:192],'-',color="k",linewidth=3,label = 'Historical (2010-1970)')
	ax.plot(x,y2[0:192],'-',color="b",linewidth=3,label = 'BEoA (2010-1970)')
	ax.plot(x,y3[0:192],'-',color="r",linewidth=3,label = 'Energy Consumption')
	ax.plot(x,y4[0:192],'-',color="g",linewidth=3,label = 'Technology Advancements')
	
	ax.fill_between(x,y1[192:384],y1[384:576],facecolor="k",alpha=0.15)
	ax.fill_between(x,y2[192:384],y2[384:576],facecolor="b",alpha=0.15)
	ax.fill_between(x,y3[192:384],y3[384:576],facecolor="r",alpha=0.15)
	ax.fill_between(x,y4[192:384],y4[384:576],facecolor="g",alpha=0.15)
	
	ax2 = ax.twinx();
	# ax2.plot(x,y4,'-',color="k",linewidth=3, label = "Historical (2010-1970")
	ax2.plot(x,ref,'-',color="m",linewidth=4,label = 'E1970')
	# ax2.fill_between(x,ref[192:384],ref[384:576],facecolor="k",alpha=0.15)
	ax2.set_xlim([-90,90]);ax2.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
	return ax,ax2

fig = plt.figure(facecolor='White',figsize=[15,12]);pad= 5   #;plot_setup();
colormap='RdBu';colormap = reverse_colourmap(colormap);

"""	
index ='TS'	 
lat, TAs, AAs, Ene,Tech,ref= zonal_mean_unceer(index,mask = 'All');

ax = plt.subplot(3,3,1);
ax.annotate(r'(a) $\mathrm{\mathsf{\Delta} SAT\/(land+ocean)\/(K)}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)			
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)

legend = ax2.legend(shadow=False,ncol=1,loc ='upper right')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

legend = ax1.legend(shadow=False,ncol=1,loc ='upper left')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

ax1.set_ylim([-2,2]);ax.set_yticks(np.arange(-2,2.1,1));
ax2.set_ylim([0,8]);ax2.set_yticks(np.arange(0,8.1,2))
align_yaxis(ax1,0,ax2,4)

lat, TAs, AAs, Ene,Tech,ref= zonal_mean_unceer(index,mask = 'land');

ax = plt.subplot(3,3,4);
ax.annotate(r'(b) $\mathrm{\mathsf{\Delta}SAT\/(land)\/(K)}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
					
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax1.set_ylim([-2,2]);ax.set_yticks(np.arange(-2,2.1,1));
ax2.set_ylim([0,8]);ax2.set_yticks(np.arange(0,8.1,2))
align_yaxis(ax1,0,ax2,4)


lat, TAs, AAs, Ene,Tech,ref= zonal_mean_unceer(index,mask = 'ocean');

ax = plt.subplot(3,3,7);
ax.annotate(r'(c) $\mathrm{\mathsf{\Delta}SAT\/(Ocean)\/(K)}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
	
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax1.set_ylim([-2,2]);ax.set_yticks(np.arange(-2,2.1,1));
ax2.set_ylim([0,8]);ax2.set_yticks(np.arange(0,8.1,2))
align_yaxis(ax1,0,ax2,4)


index ='precipitation';lat, TAs, AAs, Ene,Tech,ref =zonal_mean_unceer(index,mask = 'All');
ax = plt.subplot(3,3,2);
ax.annotate(r'(d) $\mathrm{\mathsf{\Delta}P\/(Land+Ocean)\/(mm\/day^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
				
ax3,ax4 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax3.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
ax4.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
align_yaxis(ax3,0,ax4,0)


index ='precipitation';lat, TAs, AAs, Ene,Tech,ref = zonal_mean_unceer(index,mask = 'land');
ax = plt.subplot(3,3,5);
ax.annotate(r'(e) $\mathrm{\mathsf{\Delta}P\/(Land)\/(mm\/day^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
				
ax3,ax4 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax3.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
ax4.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
align_yaxis(ax3,0,ax4,0)

index ='precipitation';lat, TAs, AAs, Ene,Tech,ref= zonal_mean_unceer(index,mask = 'ocean');
ax = plt.subplot(3,3,8);
ax.annotate(r'(f) $\mathrm{\mathsf{\Delta}P\/(Ocean)\/(mm\/day^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
				
ax3,ax4 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax3.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
ax4.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
align_yaxis(ax3,0,ax4,0)
"""
index ='HT'	;lat, TAs, AAs, Ene,Tech,ref = zonal_mean_unceer(index,mask = 'All');
ax = plt.subplot(3,2,1);
ax.annotate('(a) Total heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)									
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.3,0.30]);ax5.set_yticks(np.arange(-0.30,0.31,0.1));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

ax = plt.subplot(3,2,2);ax.axis('off')

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

x=np.arange(1,366); y =10*x/366
ax.plot(x,y,'-',color="w",linewidth=0)
ax.set_xlim([1,366]);ax.set_ylim([0,10]);

line_k = mlines.Line2D(x,y,ls='-',color='k',lw=2)
line_b = mlines.Line2D(x,y,ls='-',color='b',lw=2)
line_r = mlines.Line2D(x,y,ls='-',color='r',lw=2)
line_g = mlines.Line2D(x,y,ls='-',color='g',lw=2)
line_m = mlines.Line2D(x,y,ls='-',color='m',lw=2)

lines = [line_m,line_k,line_b,line_r,line_g]
labels = ['E1970','Historical(2010 -1970)','BEoA (2010 - 1970)','Energy consumption','Technology advancements']
legend = plt.legend(lines,labels,ncol=1,loc='best',labelspacing=1.5,markerscale =10)
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)


index ='HT'	;lat, TAs, AAs, Ene,Tech,ref = zonal_mean_unceer(index,mask = 'atm');
ax = plt.subplot(3,2,3);
ax.annotate('(b) Atmospheric heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)								
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.3,0.30]);ax5.set_yticks(np.arange(-0.30,0.31,0.1));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

index ='HT'	;lat, TAs, AAs, Ene,Tech,ref = zonal_mean_unceer(index,mask = 'ocean');
ax = plt.subplot(3,2,4);
ax.annotate('(c) Ocean heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)								
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.3,0.30]);ax5.set_yticks(np.arange(-0.30,0.31,0.1));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

index ='HT'	;lat, TAs, AAs, Ene,Tech,ref = zonal_mean_unceer(index,mask = 'sensible');
ax = plt.subplot(3,2,5);
ax.annotate('(d) Atmospheric dry static heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)								
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.3,0.30]);ax5.set_yticks(np.arange(-0.30,0.31,0.1));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)

index ='HT'	;lat, TAs, AAs, Ene,Tech,ref = zonal_mean_unceer(index,mask = 'latent');
ax = plt.subplot(3,2,6);
ax.annotate('(e) Atmospheric latent heat transport (PW)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)								
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax5.set_ylim([-0.3,0.30]);ax5.set_yticks(np.arange(-0.30,0.31,0.1));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,2));
align_yaxis(ax5,0,ax6,0)



plt.subplots_adjust(left=0.05, bottom=0.04, right=0.95, top=0.96, wspace=0.20, hspace=0.25); 
plt.savefig('zonal_mean_HT.png', format='png', dpi=1000)
