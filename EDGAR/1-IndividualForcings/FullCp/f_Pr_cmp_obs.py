"""
This is to compare the model SAT and P at 1970 and 2010 equilibrium with observations
	1. annual mean map difference
	2. zonal mean difference, especially, the tropics and the aArctic pole in SAT
	3. compare the model 1970 first and last 320 yrs of run
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import matplotlib.pyplot as plt
from datetime import date, timedelta
import math
from scipy.interpolate import interp2d  as interp2d


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *
def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area
	
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
	mask[:,:colum_s] =0; mask[:,colum_e:] =0
	mask[:row_s,:] =0; mask[row_e:,:] =0
	return mask

def mask_weight(region_key,lon,lat,reverse=False):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid(lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)

	##OCEAN_MASKS FOR COUNTRIES
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	box_region_dic={'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
	if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'GloLand'):
		mask= ocean_mask[region_key][:]
	elif  region_key in box_region_dic:
		mask= ocean_mask['All'][:]
		box = box_region_dic[region_key]
		mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
	else:
		print "error region name"
	
	# interpolate from 360*720 to 192*288
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon, lat);
	mask[mask >= 1] = 1;mask[mask < 1] = 0;
	mask[mask==0] = np.nan
	return mask


def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	lons[lons>180]-=360; 
	# calculate the origin of the map
	lon_0 = lons.mean(); 
	lat_0 = lats.mean(); 
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# s = map.pcolor(xi, yi, data)
	
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines(); #map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	# pm = np.ma.masked_not_equal(p_value, 1)
	# map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh

def CESM_readin(variable):	
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		if scenario =='T1970C': start_year =1970
		else: start_year =2010
		start =start_year*365
		ith=0	
		for iday in days:
			month_days =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
			calendar_days = np.array([31,59,90,120,151,181,212,243,273,304,334,365])
			total_days = int(iday) + start; 
			year = total_days//365; 
			remainder =  total_days%365
			if remainder ==0: year=year-1;month=12;day=31
			else : 
				month = 1+[layer for layer in range(len(calendar_days)) if calendar_days[layer]<= remainder and calendar_days[layer+1]>remainder][0]
				day = int(remainder - calendar_days[month-1])
				if day == 0: day = month_days[month-1]
			date_int[ith] = year*10000+month*100+day
			ith=ith+1
		return date_int.astype(int)
	def mon_mean2annual_mean(scenario,time,data):
		
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario=='T1970C' or scenario=='T1970':
			year_series = range(1970,1990)
		elif scenario=='T1970RCP':
			year_series = range(2020,2050)
		else:
			year_series = range(2130,2160)
		
		annual_mean=np.empty((len(year_series),192,288));annual_mean[:]=np.nan
		for iyear in year_series:
			annual = np.empty((192,288));annual[:] = 0.0			
			if (iyear == year_series[0] and time[0]//100 >= year_series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  
			if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+12][0]
			# print scenario,iyear,layer_b,layer_e
			data_cache = data[layer_b:layer_e+1,:,:]
			annual_mean[iyear-year_series[0],:,:] = stats.nanmean(data_cache,axis=0)
		mean_map = stats.nanmean(annual_mean,axis=0)	
		return mean_map

	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
		var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]*24*60*60*1000
		nc_fid.close()
		var = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var
	

	scenario = 'T1970RCP';  
	lon,lat,PRECC = data_netcdf(scenario,'PRECC'); lon,lat,PRECL =data_netcdf(scenario,'PRECL')
	T1970E20 = PRECC + PRECL
	scenario = 'EdgRef' ; 
	lon,lat,PRECC = data_netcdf(scenario,'PRECC'); lon,lat,PRECL =data_netcdf(scenario,'PRECL')
	T2010E20 = PRECC + PRECL
	return lon,lat,T1970E20,T2010E20

def lENS_PR():
	var_path ='/exports/csce/datastore/geos/users/s1667168/obs/Pr/'
	file_name = var_path+'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h0.PRECC.185001-200512.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	PRECC = nc_fid.variables['PRECC'][1320:,:,:]
	nc_fid.close()
	file_name = var_path+'b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h0.PRECL.185001-200512.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	PRECL = nc_fid.variables['PRECL'][1320:,:,:]
	nc_fid.close()
	PRECIP_6005 = (PRECC+PRECL)*1000*24*60*60

	file_name = var_path+'b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h0.PRECC.200601-208012.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	PRECC = nc_fid.variables['PRECC'][:166,:,:]
	nc_fid.close()
	file_name = var_path+'b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h0.PRECL.200601-208012.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	PRECL = nc_fid.variables['PRECL'][:166,:,:]
	nc_fid.close()
	PRECIP_0620 = (PRECC+PRECL)*1000*24*60*60
	PRECIP = np.concatenate((PRECIP_6005,PRECIP_0620),axis=0)
	M1970 = np.nanmean(PRECIP[36:216,:,:],axis=0)
	M2010 = np.nanmean(PRECIP[516:696,:,:],axis=0)
	return lon,lat,M1970,M2010
	

def GPCC_pr(lon_i,lat_i):
	def day2datetime(days):
		date_int = np.empty((len(days)));date_int[:]=np.nan;
		start = date(1800,1,1)
		ith = 0
		for day in days:
			delta = timedelta(int(day))
			offset = start + delta
			date_int[ith] = offset.year*10000+offset.month*100+offset.day
			ith=ith+1
		return date_int
	input_path ='/exports/csce/datastore/geos/users/s1667168/obs/Pr/'
	var_path = input_path+'GPCCprecip.mon.total.1x1.v7.nc'
	nc_fid = nc4.Dataset(var_path,mode='r')
	lat =nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	time1 = nc_fid.variables['time'][:];
	data1 = nc_fid.variables['precip'][:]
	
	input_path ='/exports/csce/datastore/geos/users/s1667168/obs/Pr/'
	var_path = input_path+'GPCCprecip.monitor.mon.total.1x1.v4.nc'
	nc_fid = nc4.Dataset(var_path,mode='r')
	lat =nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	time2 = nc_fid.variables['time'][:];
	data2 = nc_fid.variables['precip'][84:-2,:,:]
	data = np.concatenate((data1,data2),axis=0)
	time = np.concatenate((time1,time2),axis=0)
	data[data<0]= np.nan; data[data>3000]= np.nan; 
	data=data/30
	nc_fid.close()


	time=day2datetime(time);print time
	layer_b = [layer for layer in range(len(time)) if time[layer]//100 == 1963*100+1][0]
	layer_e = [layer for layer in range(len(time)) if time[layer]//100 == 1977*100+12][0]
	TS60_80 = data[layer_b:layer_e+1,:,:]
	layer_b = [layer for layer in range(len(time)) if time[layer]//100 == 2003*100+1][0]
	layer_e = [layer for layer in range(len(time)) if time[layer]//100 == 2017*100+12][0]
	TS00_17 = data[layer_b:layer_e+1,:,:]	
	def month_weighted_mean(years,data,lon,lat,lon_i,lat_i):
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		annual_weighted_mean = np.empty((len(years),np.shape(data)[1],np.shape(data)[2]));annual_weighted_mean[:]=0
		# monthly weigth
		for iyear in range(len(years)):
			for imonth in range(len(calendar_day)):
				annual_weighted_mean[iyear,:,:] = annual_weighted_mean[iyear,:,:]+ data[iyear*12+imonth,:,:]*calendar_day[imonth]
		annual_weighted_mean =annual_weighted_mean/365
		## annual mean 
		annual_mean = np.nanmean(annual_weighted_mean,axis = 0)
		f = interp2d(lon, lat, annual_mean,kind='linear'); 
		# annual_mean =f(lon_i, lat_i);
		return annual_mean
	M6080 = month_weighted_mean(range(1963,1978),TS60_80,lon,lat,lon_i,lat_i)
	M0017 = month_weighted_mean(range(2003,2018),TS00_17,lon,lat,lon_i,lat_i)
	return lon,lat,M6080,M0017	
	


	
lon,lat,T1970E20,T2010E20 = CESM_readin(variable='TS');
lon,lat1,LENS1970, LENS2010 = lENS_PR()
lon_N,lat_N,GPCC6080,GPCC0017 = GPCC_pr(lon,lat)
mask = mask_weight('GloLand',lon,lat,reverse=False)



fig = plt.figure(facecolor='White',figsize=[16.5,9.5]);pad= 5
colormap='BrBG';  #colormap = reverse_colourmap(colormap);
colorbar_min = -0.5 ;colorbar_max =0.5;
ax = plt.subplot(2,2,1);

ax.annotate('(a) CESM ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax, T2010E20 - T1970E20,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,2);

ax.annotate('(b) LENS ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,LENS2010 - LENS1970 ,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(2,2,3);
ax.annotate('(c) GPCC',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
colormesh1=spatial_figure(ax,GPCC0017-GPCC6080 ,lon_N,lat_N,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.18, 0.04, 0.64, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate('mm/day',xy=(1.05,-1.1), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)
# ax = plt.subplot(3,2,3);

# ax.annotate('(c) GISS',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=20)	
# spatial_figure(ax,GISS2010-GISS1970 ,lon_G,lat_G,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,4);
ax.annotate('(d) Zonal mean',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
plt.plot(lat,np.nanmean(T2010E20-T1970E20,axis=1),'k',lw=2,label ='CESM')

plt.plot(lat,np.nanmean(np.multiply(T2010E20-T1970E20,mask),axis=1),'k',lw=2,ls='--',label ='CESM land')
plt.plot(lat,np.nanmean(LENS2010-LENS1970,axis=1),'b',lw=2,label ='LENS')
plt.plot(lat,np.nanmean(np.multiply(LENS2010-LENS1970,mask),axis=1),'b',lw=2,ls='--',label ='LENS land')

GPCC = np.nanmean(GPCC0017-GPCC6080,axis=1); GPCC[GPCC == 0]=np.nan
plt.plot(lat_N,GPCC,'r',lw=2,label ='GPCC')
ax.set_xlim([-90,90]);ax.set_xticks((-90,-60,-30,0,30,60,90));
ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));

legend = ax.legend(shadow=False,ncol=2,loc ='upper left')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)


# plt.show()
plt.subplots_adjust(left=0.03, bottom=0.10, right=0.98, top=0.95, wspace=0.1, hspace=0.15); 
plt.savefig('Pr_CESM_GISS_NCEp_LENS_comp.png', format='png', dpi=1000)


