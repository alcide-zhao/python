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
			# # # annual = np.empty((192,288));annual[:] = 0.0
			# # # ### if ther is pnly 11 months, just average over the 11 months
			# # # if (iyear == year_series[0] and np.shape(data_cache)[0] <=11):
				# # # for i in range(12-np.shape(data_cache)[0],12):
					# # # annual = annual+data_cache[i-(12-np.shape(data_cache)[0]),:,:]*calendar_day[i]
				# # # annual_mean[iyear-year_series[0],:,:] = annual/np.sum(calendar_day[12-np.shape(data_cache)[0]:12])
			# # # elif (iyear == year_series[-1] and np.shape(data_cache)[0] <=11):
				# # # for i in range(0,np.shape(data_cache)[0]):
					# # # annual = annual+data_cache[i,:,:]*calendar_day[i]
				# # # annual_mean[iyear-year_series[0],:,:] = annual/np.sum(calendar_day[0:np.shape(data_cache)[0]])
			# # # else:
				# # # for i in range(0,12):
					# # # annual = annual+data_cache[i,:,:]*calendar_day[i]
				# # # annual_mean[iyear-year_series[0],:,:] = annual/365
		# mean_map = stats.nanmean(annual_mean,axis=0)
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
	# _,EdgEne = data_netcdf('EdgEne',variable)
	# _,EdgTech = data_netcdf('EdgTech',variable)
	return lat,T1970,Edg70GO,Edg70Oz,EdgRef  #,EdgEne,EdgTech
	
def zonal_mean_unceer(variable,mask):
	if variable == "precipitation":	
		lat,T1970_C,Edg70GO_C,Edg70Oz_C,EdgRef_C = data_read_zonal_mean('PRECL',mask);
		lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L = data_read_zonal_mean('PRECC',mask);
		T1970=(T1970_C+T1970_L)*24*60*60*1000  # m/s to mm/day
		Edg70GO=(Edg70GO_C+Edg70GO_L)*24*60*60*1000 
		Edg70Oz=(Edg70Oz_C+Edg70Oz_L)*24*60*60*1000 
		EdgRef=(EdgRef_C+EdgRef_L)*24*60*60*1000 
		# EdgEne=(EdgEne_C+EdgEne_L)
		# EdgTech=(EdgTech_C+EdgTech_L)
	elif variable == "TS":	
		print variable,mask
		lat,T1970,Edg70GO,Edg70Oz,EdgRef = data_read_zonal_mean(variable,mask); #,EdgEne,EdgTech
		T1970 = T1970-273.15
		Edg70GO = Edg70GO-273.15
		Edg70Oz = Edg70Oz-273.15
		EdgRef = EdgRef-273.15
	elif variable == "HT":  # heat transport
		lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L = data_read_zonal_mean('FLNT',mask='All');
		lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S = data_read_zonal_mean("FSNT",mask='All');
		lat,T1970_LHFLX,Edg70GO_LHFLX,Edg70Oz_LHFLX,EdgRef_LHFLX = data_read_zonal_mean('LHFLX',mask='All');
		lat,T1970_SHFLX,Edg70GO_SHFLX,Edg70Oz_SHFLX,EdgRef_SHFLX = data_read_zonal_mean("SHFLX",mask='All');
		lat,T1970_FLNS,Edg70GO_FLNS,Edg70Oz_FLNS,EdgRef_FLNS = data_read_zonal_mean('FLNS',mask='All');
		lat,T1970_FSNS,Edg70GO_FSNS,Edg70Oz_FSNS,EdgRef_FSNS = data_read_zonal_mean("FSNS",mask='All');
		lat,T1970_PRECSC,Edg70GO_PRECSC,Edg70Oz_PRECSC,EdgRef_PRECSC = data_read_zonal_mean('PRECSC',mask='All');
		lat,T1970_PRECSL,Edg70GO_PRECSL,Edg70Oz_PRECSL,EdgRef_PRECSL = data_read_zonal_mean("PRECSL",mask='All');		
		lat,T1970_QFLX,Edg70GO_QFLX,Edg70Oz_QFLX,EdgRef_QFLX = data_read_zonal_mean('QFLX',mask='All');
		lat,T1970_PRECC,Edg70GO_PRECC,Edg70Oz_PRECC,EdgRef_PRECC = data_read_zonal_mean("PRECC",mask='All');
		lat,T1970_PRECL,Edg70GO_PRECL,Edg70Oz_PRECL,EdgRef_PRECL = data_read_zonal_mean("PRECL",mask='All');
		if mask=='All':
			T1970 =  inferred_heat_transport(T1970_S-T1970_L,lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L,lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L,lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L,lat)
		elif mask=='ocean':
			T1970 =  inferred_heat_transport(-(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(-(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(-(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(-(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
		elif mask=='atm':
			T1970 =  inferred_heat_transport(T1970_S-T1970_L+(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L+(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L+(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L+(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
	else:
		lat,T1970,Edg70GO,Edg70Oz,EdgRef = data_read_zonal_mean(variable,mask); #,EdgEne,EdgTech
	
	def decompose_mean_std(var1,var2):
		diff = var1 - var2
		# print "diff", np.shape(diff)
		mean = np.nanmean(diff,axis=0)
		P25 = np.nanpercentile(diff,25,axis=0)
		P75 = np.nanpercentile(diff,75,axis=0)
		return mean,P25,P75
	
	mean_Total,P25_Total,P75_Total = decompose_mean_std(EdgRef,T1970)
	mean_GHG,P25_GHG,P75_GHG = 	decompose_mean_std(Edg70Oz,Edg70GO)
	mean_AAs,P25_AAs,P75_AAs = 	decompose_mean_std(Edg70GO,T1970)
	mean_O3,P25_O3,P75_O3 = decompose_mean_std(EdgRef,Edg70Oz)
	ref1970,_,_ = decompose_mean_std(T1970,0)
	# mean_GHG,P25_GHG,P75_GHG = 	decompose_mean_std(T1970,0)
	# mean_AAs,P25_AAs,P75_AAs = 	decompose_mean_std(Edg70GO,0)
	# mean_O3,P25_O3,P75_O3 = decompose_mean_std(Edg70Oz,0)
	# decompose_mean_std(EdgRef,EdgEne)
	# decompose_mean_std(EdgRef,EdgTech)
	Total = np.concatenate((mean_Total,P25_Total,P75_Total),axis=0)
	GHG = np.concatenate((mean_GHG,P25_GHG,P75_GHG),axis=0)
	AAs = np.concatenate((mean_AAs,P25_AAs,P75_AAs),axis=0)
	O3 = np.concatenate((mean_O3,P25_O3,P75_O3),axis=0)

	return lat,Total, GHG, AAs, O3,ref1970
	
def ECMWF_SAT_P(var,mask):	
	def read_and_interpolate(file_name,var,mask):
		var_path ='/exports/csce/datastore/geos/users/s1667168/obs/SAT/'+file_name
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['latitude'][:]
		lon = nc_fid.variables['longitude'][:]
		if var == "t2m":
			data = np.nanmean(nc_fid.variables['t2m'][:],axis=1)-273.15
		elif var == "tp":
			data = np.nanmean(nc_fid.variables['tp'][:],axis=1)*1000/30
		nc_fid.close()
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		annual_weighted_mean = np.empty(np.shape(data));annual_weighted_mean[:]=0
		
		for imonth in range(len(calendar_day)):
			annual_weighted_mean[imonth,:,:] = data[imonth,:,:]*calendar_day[imonth]
		annual_mean =np.nansum(annual_weighted_mean,axis=0)/365
		
		
		## interpolate the 14408720ocean_mask into the ECMF 1*1degreee grids
		ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/world.oceanmask.1440x720_0_360_-90_90.mat')['landoceanmask'][:]
		ocean_mask[np.isnan(ocean_mask)]=0;	ocean_mask[ocean_mask>0]=1;
		f = interp2d(np.arange(0,360,0.25), np.arange(-90,90,0.25), ocean_mask,kind='linear'); 
		ocean_mask = f(lon, lat);
		ocean_mask[ocean_mask >= 1] = 1;ocean_mask[ocean_mask < 1] = np.nan;
		if mask == "land":
			ocean_mask [ocean_mask ==0 ] = np.nan
		elif mask == "ocean":
			ocean_mask [ocean_mask ==1 ] = np.nan
			ocean_mask [ocean_mask ==0 ] = 1
		elif mask == "All":
			ocean_mask [ocean_mask ==0 ] = 1
		annual_mean = np.nanmean(np.multiply(annual_mean,ocean_mask),axis=1)
		return lat,annual_mean
	
	lat,M1970= read_and_interpolate('ecmwf_era20c_ts_p_1970.nc',var,mask)
	lat,M2010 =read_and_interpolate('ecmwf_era20c_ts_p_2010.nc',var,mask)
	diff = M2010-M1970
	return lat, diff 
	
def lENS_SAT_P(var,mask):
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask'][:]
	if mask == "land":
		ocean_mask [ocean_mask ==0 ] = np.nan
	elif mask == "ocean":
		ocean_mask [ocean_mask ==1 ] = np.nan
		ocean_mask [ocean_mask ==0 ] = 1
	elif mask == "All":
		ocean_mask [ocean_mask ==0 ] = 1
	
	def read_and_interpolate(file_name,var,ocean_mask):
		var_path ='/exports/csce/datastore/geos/users/s1667168/obs/SAT/'+file_name
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		if file_name == 'ensumble_mean_TS_192002_200512.nc':
			if var == "TS":data = nc_fid.variables['TS'][599:611,:,:]-273.15
			elif var == "precipitation":data = nc_fid.variables['P'][599:611,:,:]*24*60*60/1000
		elif file_name == 'ensumble_mean_TS_200602_210101.nc':
			if var == "TS":data = nc_fid.variables['rcp85'][47:59,:,:]-273.15
			elif var == "precipitation":data = nc_fid.variables['P'][47:59,:,:]*24*60*60/1000
		nc_fid.close()
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		annual_weighted_mean = np.empty(np.shape(data));annual_weighted_mean[:]=0
		
		for imonth in range(len(calendar_day)):
			annual_weighted_mean[imonth,:,:] = data[imonth,:,:]*calendar_day[imonth]
		annual_mean =np.nansum(annual_weighted_mean,axis=0)/365
		annual_mean = np.nanmean(np.multiply(annual_mean,ocean_mask),axis=1)
		return lat,annual_mean
	
	lat,LENS1970= read_and_interpolate('ensumble_mean_TS_192002_200512.nc',var,ocean_mask)
	lat,LENS2010 =read_and_interpolate('ensumble_mean_TS_200602_210101.nc',var,ocean_mask)
	diff = LENS2010 - LENS1970
	if var == 'TS': print 'LENS1970',np.average(LENS1970, weights=np.cos(np.deg2rad(lat)),axis=0)
	if var == 'TS': print 'LENS2010',np.average(LENS2010, weights=np.cos(np.deg2rad(lat)),axis=0)
	return lat,diff	

def plot_zonal_mean_uncertainty(ax,x,y1,y2,y3,y4,ref):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks(np.arange(-90,91,30));
	ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	# ax.set_ylim([-4,4.0]);ax.set_yticks(np.arange(-4,4.01,2))
	
	ax.plot(x,y1[0:192],'-',color="b",linewidth=3,label = 'Total')
	ax.plot(x,y2[0:192],'-',color="g",linewidth=3,label = 'GHGs')
	ax.plot(x,y3[0:192],'-',color="r",linewidth=3,label = 'AAs')
	ax.plot(x,y4[0:192],'-',color="m",linewidth=3,label = 'O3')
	ax.fill_between(x,y1[192:384],y1[384:576],facecolor="b",alpha=0.15)
	ax.fill_between(x,y2[192:384],y2[384:576],facecolor="g",alpha=0.15)
	ax.fill_between(x,y3[192:384],y3[384:576],facecolor="r",alpha=0.15)
	ax.fill_between(x,y4[192:384],y4[384:576],facecolor="m",alpha=0.15)
	
	ax2 = ax.twinx();
	ax2.plot(x,ref,'-',color="k",linewidth=3, label = "REF1970")
	ax2.set_xlim([-90,90]);ax2.set_xticks(np.arange(-90,91,30));
	return ax,ax2
# lev = np.array([3.64,7.59,14.36,24.61,39.27,54.60,72.01,87.82,103.31,121.55,142.89,168.23,197.91,232.83,273.91,322.24,379.10,445.99,524.68,609.77,691.39,763.40,820.86,859.53,887.02,912.64,936.20,957.49,976.33,992.56])
# lev_gradient = np.gradient(lev);  # print lev_gradient;
# for ilev in range(0,30):
	# UQ_TMP = UQI3150[ilev,:,:];VQ_TMP= VQI3150[ilev,:,:];
	# dvdlat = numerical_dif_2D(VQ_TMP,lonm, latm,ax=0); dudlon = numerical_dif_2D(UQ_TMP,lonm, latm,ax=1); 
	# QC_3150[ilev-0,:,:]=(dvdlat+dudlon)*lev_gradient[ilev]*factor_c2mmd*p_g
	# ##2081_2100
	# UQ_TMP = UQI8100[ilev,:,:];VQ_TMP = VQI8100[ilev,:,:];
	# dvdlat = numerical_dif_2D(VQ_TMP,lonm, latm,ax=0); dudlon = numerical_dif_2D(UQ_TMP,lonm, latm,ax=1); 
	# QC_8100[ilev-0,:,:]=(dvdlat+dudlon)*lev_gradient[ilev]*factor_c2mmd*p_g

	
fig = plt.figure(facecolor='White',figsize=[15,12]);plot_setup();pad= 5
colormap='RdBu';colormap = reverse_colourmap(colormap);

	
index ='TS'	 
lat,Total, GHG, AAs, O3,ref= zonal_mean_unceer(index,mask = 'All');
lat_ERA,TS_ERA = ECMWF_SAT_P(var = 't2m',mask = 'All')
lat_LENS,TS_LENS= lENS_SAT_P(var ='TS',mask = 'All')


print 'CESM: ', np.average(Total[0:192], weights=np.cos(np.deg2rad(lat)),axis=0)
print 'LENS: ', np.average(TS_LENS, weights=np.cos(np.deg2rad(lat_LENS)),axis=0)
print 'CESM1970',np.average(ref, weights=np.cos(np.deg2rad(lat_LENS)),axis=0) 


ax = plt.subplot(3,3,1);
ax.annotate(r'(a) $\mathrm{\mathsf{\Delta} SAT\/(land+ocean)\/(K)}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

				
# ax.plot(lat_ERA,TS_ERA,"k",linewidth=3,label='ERA-20CM')
ax.plot(lat_LENS,TS_LENS,"y",linewidth=3,label='CESM-LENS')				
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)

ax1.set_ylim([-2,8]);ax.set_yticks(np.arange(-2,8.1,2));
ax2.set_ylim([-50,30]);ax2.set_yticks(np.arange(-50,30.1,20))
align_yaxis(ax1,3,ax2,-10)

legend = ax2.legend(shadow=False,ncol=1,loc ='upper right')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

legend = ax1.legend(shadow=False,ncol=1,loc ='upper left')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)


lat,Total, GHG, AAs, O3,ref= zonal_mean_unceer(index,mask = 'land');
lat_ERA,TS_ERA = ECMWF_SAT_P(var = 't2m',mask = 'land')
lat_LENS,TS_LENS= lENS_SAT_P('TS',mask = 'land')

ax = plt.subplot(3,3,4);
ax.annotate(r'(b) $\mathrm{\mathsf{\Delta}SAT\/(land)\/(K)}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

# ax.plot(lat_ERA,TS_ERA,"k",linewidth=3,label='ERA-20CM')
ax.plot(lat_LENS,TS_LENS,"y",linewidth=3,label='CESM-LENS')						
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax1.set_ylim([-2,8]);ax.set_yticks(np.arange(-2,8.1,2));
ax2.set_ylim([-50,30]);ax2.set_yticks(np.arange(-50,30.1,20))
align_yaxis(ax1,3,ax2,-10)


lat,Total, GHG, AAs, O3,ref= zonal_mean_unceer(index,mask = 'ocean');
lat_ERA,TS_ERA = ECMWF_SAT_P(var = 't2m',mask = 'ocean')
lat_LENS,TS_LENS= lENS_SAT_P('TS',mask = 'ocean')


ax = plt.subplot(3,3,7);
ax.annotate(r'(c) $\mathrm{\mathsf{\Delta}SAT\/(Ocean)\/(K)}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

# ax.plot(lat_ERA,TS_ERA,"k",linewidth=3,label='ERA-20CM')
ax.plot(lat_LENS,TS_LENS,"y",linewidth=3,label='CESM-LENS')	
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax1.set_ylim([-2,8]);ax.set_yticks(np.arange(-2,8.1,2));
ax2.set_ylim([-50,30]);ax2.set_yticks(np.arange(-50,30.1,20))
align_yaxis(ax1,3,ax2,-10)


index ='precipitation';lat,Total, GHG, AAs, O3,ref= zonal_mean_unceer(index,mask = 'All');
ax = plt.subplot(3,3,2);
ax.annotate(r'(d) $\mathrm{\mathsf{\Delta}P\/(Land+Ocean)\/(mm\/day^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
ax3,ax4 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax3.set_ylim([-1.0,1.0]);ax.set_yticks(np.arange(-1,1.1,0.5));
ax4.set_ylim([0,8]);ax4.set_yticks(np.arange(0,8.1,2))
align_yaxis(ax3,0,ax4,4)

index ='precipitation';lat,Total, GHG, AAs, O3,ref= zonal_mean_unceer(index,mask = 'land');
ax = plt.subplot(3,3,5);
ax.annotate(r'(e) $\mathrm{\mathsf{\Delta}P\/(Land)\/(mm\/day^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
ax3,ax4 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax3.set_ylim([-1.0,1.0]);ax.set_yticks(np.arange(-1,1.1,0.5));
ax4.set_ylim([0,8]);ax4.set_yticks(np.arange(0,8.1,2))
align_yaxis(ax3,0,ax4,4)

index ='precipitation';lat,Total, GHG, AAs, O3,ref= zonal_mean_unceer(index,mask = 'ocean');
ax = plt.subplot(3,3,8);
ax.annotate(r'(f) $\mathrm{\mathsf{\Delta}P\/(Ocean)\/(mm\/day^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
ax3,ax4 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax3.set_ylim([-1.0,1.0]);ax.set_yticks(np.arange(-1,1.1,0.5));
ax4.set_ylim([0,8]);ax4.set_yticks(np.arange(0,8.1,2))
align_yaxis(ax3,0,ax4,4)



index ='HT'	;lat,Total, GHG, AAs, O3,ref = zonal_mean_unceer(index,mask = 'All');
ax = plt.subplot(3,3,3);
ax.annotate(r'(g) $\mathrm{\mathsf{\Delta}\/NHT\/(Ocean+Atmosphere)\/(PW)}$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)									
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax5.set_ylim([-0.4,0.4]);ax.set_yticks(np.arange(-0.4,0.41,0.2));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,3))
align_yaxis(ax5,0,ax6,0)

index ='HT'	;lat,Total, GHG, AAs, O3,ref = zonal_mean_unceer(index,mask = 'atm');
ax = plt.subplot(3,3,6);
ax.annotate(r'(h) $\mathrm{\mathsf{\Delta}\/NHT\/(Atmosphere)\/(PW)}$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)								
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax5.set_ylim([-0.4,0.4]);ax.set_yticks(np.arange(-0.4,0.41,0.2));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,3))
align_yaxis(ax5,0,ax6,0)


index ='HT'	;lat,Total, GHG, AAs, O3,ref = zonal_mean_unceer(index,mask = 'ocean');
ax = plt.subplot(3,3,9);
ax.annotate(r'(i) $\mathrm{\mathsf{\Delta}\/NHT\/(Ocean)\/(PW)}$',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)								
ax5,ax6 = plot_zonal_mean_uncertainty(ax,lat,Total, GHG, AAs, O3,ref)
ax5.set_ylim([-0.4,0.4]);ax.set_yticks(np.arange(-0.4,0.41,0.2));
ax6.set_ylim([-6,6]);ax6.set_yticks(np.arange(-6,6.1,3))
align_yaxis(ax5,0,ax6,0)


"""	
index ='VQ'	;lat,Total, GHG, AAs, O3,ref = zonal_mean_unceer(index);
ax = plt.subplot(6,2,4);
ax.annotate(r'(d) $\mathrm{\mathsf{\Delta}\/Northward\/Water\/Transport\/(m\/s^{-1} g\/kg^{-1})}$',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
ax7,ax8 = plot_zonal_mean_uncertainty(ax,lat,Total*1000, GHG*1000, AAs*1000, O3*1000,ref*1000)
ax7.set_ylim([-4,4]);ax.set_yticks(np.arange(-4,4.1,2));			
ax8.set_ylim([-20,20]);ax8.set_yticks(np.arange(-20,20.1,10))
align_yaxis(ax7,0,ax8,0)
"""			
plt.subplots_adjust(left=0.05, bottom=0.04, right=0.95, top=0.96, wspace=0.2, hspace=0.25); 
plt.savefig('./GHG_AEROSOL_OZONE/zonal_mean_TS_P_HT.png', format='png', dpi=1000)
