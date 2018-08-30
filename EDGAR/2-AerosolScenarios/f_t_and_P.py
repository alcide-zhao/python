"""
This is to plot  TS and Pr together with its zonal mean from the EDGAR six experiments
data inputs are forty years of monthly TS
	first step is to process the monthly TS into monthly mean
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
from scipy.stats import mannwhitneyu as man_test
from scipy.stats import ttest_ind as student_test

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
	# plt.imshow(p_value);plt.show()
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
		return annual_mean

	def data_netcdf(scenario,FREQ,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]#-273.15
		if variable == 'TS':data= data-273.15
		elif variable == 'PRECL' or variable == 'PRECC': data= data*24*60*60*1000
		# elif variable == 'PSL':  data= data*100
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
	
def print_domain_mean(variable,FREQ):
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
	
	def mannwhitneyu_test(vairable1,variable2):
		p_threshold=0.05
		size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
		p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
		# from scipy.stats import ttest_ind as test
		for x in range(size[0]):
			for y in range(size[1]):
				cache1 = vairable1[:,x,y]
				cache2 = variable2[:,x,y]
				if (cache2==cache1).all(): p_value[x,y] = np.nan
				else: _,p_value[x,y] = man_test(cache1,cache2);
		p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
		return p_value
	####calculate the difference and their significance
	def diff_sig(vairable1,variable2):
		dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
		sig = mannwhitneyu_test(vairable1, variable2)
		# sig_dif = np.multiply(dif,sig)
		return dif,sig
		
	if variable == "Pr":	  # cloud forcing
		lon,lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_readin('PRECC',FREQ);
		lon,lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L,EdgEne_L,EdgTech_L = data_readin('PRECL',FREQ);
		T1970=(T1970_S+T1970_L)
		Edg70GO=(Edg70GO_S+Edg70GO_L) 
		Edg70Oz=(Edg70Oz_S+Edg70Oz_L)
		EdgRef=(EdgRef_S+EdgRef_L)
		EdgEne=(EdgEne_S+EdgEne_L)
		EdgTech=(EdgTech_S+EdgTech_L)
	else:
		lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_readin(variable,FREQ);
		
	# GHG,GHG_s = diff_sig(Edg70Oz,Edg70GO)
	# OZONE,OZONE_s = diff_sig(EdgRef, Edg70Oz)
	AEROSOL,AEROSOL_s = diff_sig(Edg70GO,T1970)
	# Total2010,Total2010_s = diff_sig(EdgRef,T1970)
	ENE,ENE_s= diff_sig(EdgRef,EdgEne)
	TECH,TECH_s= diff_sig(EdgRef,EdgTech)
	
	def dif_std_for_region(var1,var2,region_key):
		"""
		for a rgional sif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
		var2_domain_mean =np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
		dif =  round(np.nanmean(var1_domain_mean - var2_domain_mean,axis=0),4)
		_,p_value = student_test(var1_domain_mean , var2_domain_mean)
		p_value = round(p_value,4)
		return dif,p_value
	# region_keys =['All','Globe','ASIA','Europe','USA']  #,'USA','SESA','EA','India'
	# for region_key in region_keys:
		# # print region_key
		# mask = mask_weight(region_key,lon,lat); #print mask
		# print "=========="+region_key+"=========="
		# # print 'Total2010 ',dif_std_for_region(EdgRef,T1970,region_key)
		# # print 'GHG       ',dif_std_for_region(Edg70Oz,Edg70GO,region_key)
		# print 'AEROSOL   ',dif_std_for_region(Edg70GO,T1970,region_key)
		# print 'OZONE     ',dif_std_for_region(EdgRef,Edg70Oz,region_key)
		# print 'ENE      ',dif_std_for_region(EdgRef,EdgEne,region_key)
		# print 'TECH     ',dif_std_for_region(EdgRef,EdgTech,region_key)
	return lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s


def plot_zonal_mean_uncertainty(ax,x,y,xlim):
	ax.tick_params(axis='both', which='major', direction = 'in',left=True,right=True,bottom=True,top=True, pad=5)
	ax.axvline(x=0,color='k',linewidth = 1);ax.axhline(y=0,color='k',linewidth = 1);
	ax.set_ylim([-90,90]);ax.set_yticks(np.arange(-90,91,30));ax.set_yticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	ax.set_xlim([-xlim,xlim]);ax.set_xticks(np.arange(-xlim,xlim+.0000001,xlim/2.0))
	ax.plot(x,y,'-',color="k",linewidth=5)
	return ax

fig = plt.figure(facecolor='White',figsize=[27.5,16.5]);pad= 5;	

index ='TS'; colorbar_min=-1;colorbar_max=1;colormap='RdBu_r';
lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s = print_domain_mean(index,FREQ='mon');
################# TS ######################
ax = plt.subplot2grid((3, 8), (0, 0), colspan=3)
ax.annotate('Surface temperature '+r'($\mathrm{\mathsf{\Delta}\/K}$)',xy=(0.5,1.05), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)
ax.annotate('BEoA (2010 - 1970)',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)
ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot2grid((3, 8), (0, 3), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(AEROSOL,axis=1),192),lat,xlim=2)

ax = plt.subplot2grid((3, 8), (1, 0), colspan=3)
ax.annotate('Energy Consumption',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)
ax.annotate('(b)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot2grid((3, 8), (1, 3), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(ENE,axis=1),192),lat,xlim=2)

ax = plt.subplot2grid((3, 8), (2, 0), colspan=3)
ax.annotate('Technology Advancements',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)
ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
colormesh1 = spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot2grid((3, 8), (2, 3), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(TECH,axis=1),192),lat,xlim=2)

cbar_ax = fig.add_axes([0.06, 0.05, 0.42, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',
	extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

#################precipitation######################
index ='Pr'; colorbar_min=-0.4;colorbar_max=0.4;colormap='BrBG';
lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s = print_domain_mean(index,FREQ='mon');

ax = plt.subplot2grid((3, 8), (0, 4), colspan=3)
ax.annotate('Precipitation ' +r'($\mathrm{\mathsf{\Delta}\/mm\/day^{-1}}$)' ,xy=(0.5,1.05), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)
ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot2grid((3, 8), (0, 7), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(AEROSOL,axis=1),192),lat,xlim=0.2)

ax = plt.subplot2grid((3, 8), (1, 4), colspan=3)
ax.annotate('(e)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot2grid((3, 8), (1, 7), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(ENE,axis=1),192),lat,xlim=0.2)

ax = plt.subplot2grid((3, 8), (2, 4), colspan=3)
ax.annotate('(f)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
colormesh2 = spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot2grid((3, 8), (2, 7), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(TECH,axis=1),192),lat,xlim=0.2)
cbar_ax = fig.add_axes([0.54, 0.05, 0.42, 0.02])
char = fig.colorbar(colormesh2,orientation='horizontal',
	extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))


plt.subplots_adjust(left=0.02, bottom=0.10, right=0.99, top=0.95, wspace=0.08, hspace=0.15); 
plt.savefig('TP_ZONAL_MEAN.png', format='png', dpi=300)				