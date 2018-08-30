"""
This is to plot the AODVIS at a wavelength of 550nm from the EDGAR six experiments
data inputs are forty years of monthly AODVIS
	first step is to process the monthly AODVIS into monthl mean
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
from scipy.interpolate import interp2d  as interp2d

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


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
	map.drawcoastlines(); map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(10,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_under('w') #cmap.set_over('r'); 
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True) #
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
		mean_map = np.nanmean(annual_mean,axis=0)
		return mean_map

	def data_netcdf(scenario,FREQ,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]
		ms_vale = nc_fid.variables[variable].missing_value; data[data == ms_vale] = np.nan
		FillValue = nc_fid.variables[variable]._FillValue; data[data == FillValue] = np.nan
		data[data>1] =np.nan;data[data<0] =np.nan;
		nc_fid.close()
		var_mean_map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var_mean_map
	
	lon,lat,Edg70GO = data_netcdf('Edg70GO',FREQ,variable)
	_,_,T1970 = data_netcdf('T1970C',FREQ,variable)
	_,_,EdgRef = data_netcdf('EdgRef',FREQ,variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',FREQ,variable)
	_,_,EdgEne = data_netcdf('EdgEne',FREQ,variable)
	_,_,EdgTech = data_netcdf('EdgTech',FREQ,variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	

	
def print_domain_mean(variable,FREQ):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weieghted)
	"""
	def AreaWeight(lon1,lon2,lat1,lat2):
		'''
		calculate the earth radius in m2
		'''
		radius = 6371000;
		area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
		(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
		# plt.imshow(area);plt.show()
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
		##OCEAN_MASKS FOR COUNTRIES
		ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
		lon_mask = ocean_mask['lon'][0,:];
		lat_mask = ocean_mask['lat'][0,:];
		box_region_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,60],'EA':[100,145,20,50],'SA':[65,100,5,30]}
		if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'Globe'):
			mask= ocean_mask[region_key][:]
		elif (region_key == 'ASIA' or region_key == 'EA' or region_key == 'SA'):
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
		# weight each grid cell by its area
		mask=np.multiply(mask,area);  
		mask=np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		return mask	
	
	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_readin(variable,FREQ);
	region_keys =['Globe','China','EA','Europe','India','USA']
	for region_key in region_keys:
		# print region_key
		mask = mask_weight(region_key,lon,lat)
		print "======"+region_key+"========"
		print 'T1970     ',round(np.nansum(np.nansum(np.multiply(mask,T1970),axis=1),axis=0),4)
		print 'Edg70GO   ',round(np.nansum(np.nansum(np.multiply(mask,Edg70GO),axis=1),axis=0),4)
		print 'Edg70Oz   ',round(np.nansum(np.nansum(np.multiply(mask,Edg70Oz),axis=1),axis=0),4)
		print 'EdgRef    ',round(np.nansum(np.nansum(np.multiply(mask,EdgRef),axis=1),axis=0),4)
		print 'EdgEne    ',round(np.nansum(np.nansum(np.multiply(mask,EdgEne),axis=1),axis=0),4)
		print 'EdgTech   ',round(np.nansum(np.nansum(np.multiply(mask,EdgTech),axis=1),axis=0),4)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	
lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = print_domain_mean(variable='AODVIS',FREQ='mon');

fig = plt.figure(facecolor='White',figsize=[12,7]);plot_setup();pad= 5;
colormap='RdYlBu';colormap = reverse_colourmap(colormap);


ax = plt.subplot(2,2,1);colorbar_min=0;colorbar_max=0.5;
ax.annotate('(a) REF1970',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
colormesh0 = spatial_figure(ax,T1970,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.06, 0.05, 0.40, 0.015])
char = fig.colorbar(colormesh0,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))

colorbar_min=-0.1;colorbar_max=0.1;

ax = plt.subplot(2,2,2);
ax.annotate('(b) REF2010 - REF1970',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,EdgRef-T1970,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

# ax = plt.subplot(3,2,3);
# ax.annotate('(c) STAG_GO2010',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)	
# spatial_figure(ax,Edg70GO-T1970,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
# ax = plt.subplot(3,2,4);
# ax.annotate('(d) STAG_ENE2010',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)	
# spatial_figure(ax,EdgEne-T1970,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot(2,2,3);
ax.annotate('(c) REF2010-STAG_ENE',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,EdgRef -EdgEne,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot(2,2,4);
ax.annotate('(d) REF2010-STAG_TECH',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
colormesh1=spatial_figure(ax,EdgRef-EdgTech,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


cbar_ax = fig.add_axes([0.55, 0.05, 0.40, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.02, bottom=0.08, right=0.98, top=0.95, wspace=0.04, hspace=0.15); 
plt.savefig('AODVIS550.png', format='png', dpi=1000)