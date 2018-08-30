"""
This is to plot the PRECT+l  from the EDGAR six experiments
data inputs are forty years of monthly PRECT+l
	first step is to process the monthly PRECT+l into monthly mean
	second step is to process annual mean of monthly mean
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu as man_test
from scipy.stats import ttest_ind as student_test
from scipy.interpolate import interp2d  as interp2d

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


def spatial_figure(axs,data,p_value,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
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
	cmap = discrete_cmap(10,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh

		
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
			### if ther is pnly 11 months, just average over the 11 months
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
		# # # # mean_map = stats.nanmean(annual_mean,axis=0)
		return annual_mean

	def data_netcdf(scenario,FREQ,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]*24*60*60*1000
		nc_fid.close()
		var40map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var40map
	
	lon,lat,Edg70GO = data_netcdf('Edg70GO',FREQ,variable)
	_,_,T1970 = data_netcdf('T1970RCP',FREQ,variable)
	_,_,EdgRef = data_netcdf('EdgRef',FREQ,variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',FREQ,variable)
	_,_,EdgEne = data_netcdf('EdgEne',FREQ,variable)
	_,_,EdgTech = data_netcdf('EdgTech',FREQ,variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	
def print_domain_mean(variable1,variable2,FREQ):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
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
				print "error region name"
			# interpolate from 360*720 to 192*288
			mask[np.isnan(mask)]=0;	mask[mask>0]=1;
			f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
			mask = f(lon, lat);
			mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
			# weight each grid cell by its area weight against the total area
			mask=np.multiply(mask,area);  
			mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
			# print np.nansum(np.nansum(mask_weighted,axis=1),axis=0)
		return mask_weighted
		
	####calculate the difference and their significance
	def diff_sig(vairable1,variable2,ref):
		def mannwhitneyu_test(vairable1,variable2):
			p_threshold=0.10
			size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
			p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
			# from scipy.stats import mannwhitneyu as test
			for x in range(size[0]):
				for y in range(size[1]):
					cache1 = vairable1[:,x,y]
					cache2 = variable2[:,x,y]
					if (cache2==cache1).all(): p_value[x,y] = np.nan
					else: _,p_value[x,y] = man_test(cache1,cache2);
			p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
			return p_value
		##=======================================================
		##============= absolute or percent =====================
		##=======================================================
		# dif = np.divide(np.nanmean(vairable1[0:40,:,:],axis=0)-np.nanmean(variable2[0:40,:,:],axis=0),np.nanmean(ref[0:40,:,:],axis=0))*100  #np.nanmean(ref[0:30,:,:],axis=0)
		dif = np.nanmean(vairable1[0:40,:,:],axis=0)-np.nanmean(variable2[0:40,:,:],axis=0)
		sig = mannwhitneyu_test(vairable1[0:40,:,:], variable2[0:40,:,:])
		# sig_dif = np.multiply(dif,sig)
		return dif,sig
	
	lon,lat,T1970_C,Edg70GO_C,Edg70Oz_C,EdgRef_C,EdgEne_C,EdgTech_C = data_readin(variable1,FREQ);
	lon,lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L,EdgEne_L,EdgTech_L = data_readin(variable2,FREQ);
	T1970=(T1970_C+T1970_L)
	Edg70GO=(Edg70GO_C+Edg70GO_L)
	Edg70Oz=(Edg70Oz_C+Edg70Oz_L)
	EdgRef=(EdgRef_C+EdgRef_L)
	EdgEne=(EdgEne_C+EdgEne_L)
	EdgTech=(EdgTech_C+EdgTech_L)
	
	GHG,GHG_s = diff_sig(Edg70Oz,Edg70GO,T1970)
	OZONE,OZONE_s = diff_sig(EdgRef, Edg70Oz,T1970)
	AEROSOL,AEROSOL_s = diff_sig(Edg70GO,T1970,T1970)
	Total2010,Total2010_s = diff_sig(EdgRef,T1970,T1970)
	ENE,ENE_s= diff_sig(EdgRef,EdgEne,T1970)
	TECH,TECH_s= diff_sig(EdgRef,EdgTech,T1970)
	
	# def dif_std_for_region(var1,var2,region_key,ref):
		# """
		# for a rgional sif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		# """
		# mask = mask_weight(region_key,lon,lat);
		# var1_domain_mean = np.nanmean(np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1),axis=0);
		# var2_domain_mean =np.nanmean(np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1),axis=0);
		# ref_domain_mean = np.nanmean(np.nansum(np.nansum(np.multiply(mask,ref),axis=2),axis=1),axis=0);
		# dif =  var1_domain_mean - var2_domain_mean
		# dif_a = round(dif,4)
		# dif_p = round(100*dif/ref_domain_mean,4)
		# return dif_a,dif_p
	def dif_std_for_region(var1,var2,region_key,ref):
		"""
		for a rgional sif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
		var2_domain_mean =np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
		dif_a =  round(np.nanmean(var1_domain_mean - var2_domain_mean,axis=0),4)
		ref_domain_mean = np.nanmean(np.nansum(np.nansum(np.multiply(mask,ref),axis=2),axis=1),axis=0);
		_,p_value_a = student_test(var1_domain_mean , var2_domain_mean)
		p_value_a = round(p_value_a,4)
		var1_domain_mean_perc = 100*np.divide(var1_domain_mean,ref_domain_mean)
		var2_domain_mean_perc = 100*np.divide(var2_domain_mean,ref_domain_mean)
		dif_p =  round(np.nanmean(var1_domain_mean_perc - var2_domain_mean_perc,axis=0),4)
		_,p_value_p = student_test(var1_domain_mean_perc , var2_domain_mean_perc)		
		p_value_p = round(p_value_p,4)
		return dif_a,p_value_a,dif_p,p_value_p
	region_keys =['All','Globe','ASIA','Europe','USA']  #,'USA','SESA','EA','India'
	for region_key in region_keys:
		# print region_key
		mask = mask_weight(region_key,lon,lat); #print mask
		print "=========="+region_key+"=========="
		print 'Total2010 ',dif_std_for_region(EdgRef,T1970,region_key,T1970)
		print 'GHG       ',dif_std_for_region(Edg70Oz,Edg70GO,region_key,T1970)
		print 'AEROSOL   ',dif_std_for_region(Edg70GO,T1970,region_key,T1970)
		print 'OZONE     ',dif_std_for_region(EdgRef,Edg70Oz,region_key,T1970)
		print 'ENE      ',dif_std_for_region(EdgRef,EdgEne,region_key,T1970)
		print 'TECH     ',dif_std_for_region(EdgRef,EdgTech,region_key,T1970)
	return lon,lat,GHG,OZONE,AEROSOL,Total2010,ENE,TECH,GHG_s,OZONE_s,AEROSOL_s,Total2010_s,ENE_s,TECH_s


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
		annual_mean=np.empty((40,192,288));annual_mean[:]=np.nan
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario=='T1970RCP':
			year_series = range(2010,2050)
		elif scenario=='EdgEne':
			year_series = range(2190,2230)
		elif scenario=='Edg70GO':
			year_series = range(2060,2100)
		else:
			year_series = range(2140,2180)
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
		if mask=='All':
			T1970 =  inferred_heat_transport(T1970_S-T1970_L,lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L,lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L,lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L,lat)
			EdgEne= inferred_heat_transport(EdgEne_S-EdgEne_L,lat)
			EdgTech=  inferred_heat_transport(EdgTech_S-EdgTech_L,lat)			
		elif mask=='ocean':
			T1970 =  inferred_heat_transport(-(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(-(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(-(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(-(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
			EdgEne=  inferred_heat_transport(-(EdgEne_FLNS-EdgEne_FSNS +EdgEne_LHFLX+EdgEne_SHFLX +(EdgEne_PRECSC+EdgEne_PRECSL)/2),lat)
			EdgTech= inferred_heat_transport(-(EdgTech_FLNS-EdgTech_FSNS +EdgTech_LHFLX+EdgTech_SHFLX +(EdgTech_PRECSC+EdgTech_PRECSL)/2),lat)			
		elif mask=='atm':
			T1970 =  inferred_heat_transport(T1970_S-T1970_L+(T1970_FLNS-T1970_FSNS +T1970_LHFLX+T1970_SHFLX +(T1970_PRECSC+T1970_PRECSL)/2),lat)
			Edg70GO= inferred_heat_transport(Edg70GO_S-Edg70GO_L+(Edg70GO_FLNS-Edg70GO_FSNS +Edg70GO_LHFLX+Edg70GO_SHFLX +(Edg70GO_PRECSC+Edg70GO_PRECSL)/2),lat)
			Edg70Oz= inferred_heat_transport(Edg70Oz_S-Edg70Oz_L+(Edg70Oz_FLNS-Edg70Oz_FSNS +Edg70Oz_LHFLX+Edg70Oz_SHFLX +(Edg70Oz_PRECSC+Edg70Oz_PRECSL)/2),lat)
			EdgRef=  inferred_heat_transport(EdgRef_S-EdgRef_L+(EdgRef_FLNS-EdgRef_FSNS +EdgRef_LHFLX+EdgRef_SHFLX +(EdgRef_PRECSC+EdgRef_PRECSL)/2),lat)			
			EdgEne=  inferred_heat_transport(EdgEne_S-EdgEne_L+(EdgEne_FLNS-EdgEne_FSNS +EdgEne_LHFLX+EdgEne_SHFLX +(EdgEne_PRECSC+EdgEne_PRECSL)/2),lat)	
			EdgTech= inferred_heat_transport(EdgTech_S-EdgTech_L+(EdgTech_FLNS-EdgTech_FSNS +EdgTech_LHFLX+EdgTech_SHFLX +(EdgTech_PRECSC+EdgTech_PRECSL)/2),lat)	
	
	else:
		lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_read_zonal_mean(variable,mask); #,EdgEne,EdgTech
	def decompose_mean_std(var1,var2):
		diff = var1 - var2
		# print "diff", np.shape(diff)
		mean = np.nanmean(diff,axis=0)
		P25 = mean - np.abs(mean-np.nanpercentile(diff,25,axis=0))/1.5
		P75 = mean + np.abs(np.nanpercentile(diff,75,axis=0)-mean)
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
	ax.set_xlim([-90,90]);ax.set_xticks(np.arange(-90,91,30));
	ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	# ax.set_ylim([-4,4.0]);ax.set_yticks(np.arange(-4,4.01,2))
	# ax.plot(x,y1[0:192],'-',color="b",linewidth=3,label = 'Historical (2010-1970)')
	ax.plot(x,y2[0:192],'-',color="b",linewidth=3,label = 'BEoA (2010-1970)')
	ax.plot(x,y3[0:192],'-',color="r",linewidth=3,label = 'Energy Consumption')
	ax.plot(x,y4[0:192],'-',color="g",linewidth=3,label = 'Technology Advancements')
	
	# ax.fill_between(x,y1[192:384],y1[384:576],facecolor="b",alpha=0.15)
	ax.fill_between(x,y2[192:384],y2[384:576],facecolor="b",alpha=0.15)
	ax.fill_between(x,y3[192:384],y3[384:576],facecolor="r",alpha=0.15)
	ax.fill_between(x,y4[192:384],y4[384:576],facecolor="g",alpha=0.15)
	
	ax2 = ax.twinx();
	# ax2.plot(x,y4,'-',color="k",linewidth=3, label = "Historical (2010-1970")
	ax2.plot(x,y1[0:192],'-',color="k",linewidth=3,label = 'Historical (2010-1970)')
	ax2.fill_between(x,y1[192:384],y1[384:576],facecolor="k",alpha=0.15)
	ax2.set_xlim([-90,90]);ax2.set_xticks(np.arange(-90,91,30));
	return ax,ax2
	
lon,lat,GHG,OZONE,AEROSOL,Total2010,ENE,TECH,GHG_s,OZONE_s,AEROSOL_s,Total2010_s,ENE_s,TECH_s= print_domain_mean(variable1='PRECC',variable2='PRECL',FREQ='mon');


colormap='BrBG';colorbar_min=-0.4;colorbar_max=0.4;
# colormap='BrBG';colorbar_min=-15;colorbar_max=15;

	
fig2 = plt.figure(facecolor='White',figsize=[20.5,11]);pad= 5; #plot_setup()
ax = plt.subplot(2,2,1);
ax.annotate('(a) BEoA (2010 - 1970)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(2,2,2);
ax.annotate('(b) Energy Consumption',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot(2,2,3);
ax.annotate('(c) Technology Advancements',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
colormesh1=spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

index='precipitation'
lat, TAs, AAs, Ene,Tech,ref= zonal_mean_unceer(index,mask = 'All');

print np.corrcoef(Ene[96:192],Tech[96:192])
ax = plt.subplot(2,2,4);
ax.annotate('(d) zonal mean',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax1,ax2 = plot_zonal_mean_uncertainty(ax,lat, TAs, AAs, Ene,Tech,ref)
ax1.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
ax2.set_ylim([-0.6,0.6]);ax.set_yticks(np.arange(-0.6,0.61,0.2));
align_yaxis(ax1,0,ax2,0)

legend = ax2.legend(shadow=False,ncol=1,loc ='upper right',fontsize=20)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

legend = ax1.legend(shadow=False,ncol=1,loc ='upper left',fontsize=20)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

cbar_ax = fig2.add_axes([0.94, 0.07, 0.02, 0.83])
char = fig2.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))


cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/mm\/day^{-1}}$',xy=(2.5,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='bottom',rotation='vertical',fontsize=15)
# cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}}$ %',xy=(1.1,-0.9), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='bottom',rotation='horizontal',fontsize=15)	
				

# plt.subplots_adjust(left=0.02, bottom=0.05, right=0.90, top=0.92, wspace=0.08, hspace=0.15); 
# plt.savefig(index+'_3scen_plus_zonalmean.png', format='png', dpi=1000)
