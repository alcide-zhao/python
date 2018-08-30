"""
This is to plot the TS  from the EDGAR six experiments
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
    os.path.pardir,os.path.pardir,os.path.pardir,
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
		
def data_readin(variable):	
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		start_year =2000
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

		year_series = range(2005,2029)
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
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/CamOnly/'
		var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]
		if variable == "TS": data=data-273.15
		elif variable == 'CLDTOT':data=data*100
		var40map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var40map
	
	if variable == 'precip':
		lon,lat,PRECC = data_netcdf('F_F_Edg70GO','PRECC');_,_,PRECL = data_netcdf('F_Edg70GO','PRECL');F_Edg70GO=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('F_1970','PRECC');_,_,PRECL = data_netcdf('F_1970','PRECL');T1970=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('F_EdgRef','PRECC');_,_,PRECL = data_netcdf('F_EdgRef','PRECL');F_EdgRef=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('F_Edg70Oz','PRECC');_,_,PRECL = data_netcdf('F_Edg70Oz','PRECL');F_Edg70Oz=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('F_Edg70T10SOZ','PRECC');_,_,PRECL = data_netcdf('F_Edg70T10SOZ','PRECL');F_Edg70T10SOZ=(PRECC+PRECL)*24*60*60*1000
	elif variable == 'ERFT':  # SW- - LW at TOA
		lon,lat,FSNT = data_netcdf('F_Edg70GO','FSNT');_,_,FLNT = data_netcdf('F_Edg70GO','FLNT');F_Edg70GO=(FSNT-FLNT)
		_,_,FSNT = data_netcdf('F_1970','FSNT');_,_,FLNT = data_netcdf('F_1970','FLNT');T1970=(FSNT-FLNT)
		_,_,FSNT = data_netcdf('F_EdgRef','FSNT');_,_,FLNT = data_netcdf('F_EdgRef','FLNT');F_EdgRef=(FSNT-FLNT)
		_,_,FSNT = data_netcdf('F_Edg70Oz','FSNT');_,_,FLNT = data_netcdf('F_Edg70Oz','FLNT');F_Edg70Oz=(FSNT-FLNT)
		_,_,FSNT = data_netcdf('F_Edg70T10SOZ','FSNT');_,_,FLNT = data_netcdf('F_Edg70T10SOZ','FLNT');F_Edg70T10SOZ=(FSNT-FLNT)
	else:
		lon,lat,F_Edg70GO = data_netcdf('F_Edg70GO',variable)
		_,_,T1970 = data_netcdf('F_1970',variable)
		_,_,F_EdgRef = data_netcdf('F_EdgRef',variable)
		_,_,F_Edg70Oz = data_netcdf('F_Edg70Oz',variable)
		_,_,F_Edg70T10SOZ = data_netcdf('F_Edg70T10SOZ',variable)
	return lon,lat,T1970,F_Edg70GO,F_Edg70Oz,F_EdgRef,F_Edg70T10SOZ
	
def spa_pat_reg_mean(variable):
	"""
	option: SpaPat or RegMean
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
		mask[:,:colum_s] =0; mask[:,colum_e:] =0
		mask[:row_s,:] =0; mask[row_e:,:] =0
		return mask

	def mask_weight(region_key,lon,lat):
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
		# plt.imshow(mask,origin='lower');plt.show()
		mask[mask >= 1] = 1;mask[mask < 1] = 0;
		# weight each grid cell by its area weight against the total area
		mask[mask==0] = np.nan
		mask=np.multiply(mask,area); 
		mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		# print np.nansum(np.nansum(mask_weighted,axis=1),axis=0)
		return mask_weighted
	
	def dif_std_for_region(var1,var2,mask):
		"""
		for a regional dif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		
		dif =  np.nansum(np.nansum(np.multiply(mask,np.nanmean(var1,axis=0) - np.nanmean(var2,axis=0)),axis=1),axis=0)
		var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
		var2_domain_mean = np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
		std = np.std(var1_domain_mean - var2_domain_mean);print std
		p25 =  np.abs(dif - np.nanpercentile(var1_domain_mean - var2_domain_mean,25,axis=0))/1.25
		p75 =  np.abs(np.nanpercentile(var1_domain_mean - var2_domain_mean,75,axis=0) - dif)/1.25
		# print dif, p25,p75
		return dif,p25,p75
		

	lon,lat,F_1970,F_Edg70GO,F_Edg70Oz,F_EdgRef,F_Edg70T10SOZ = data_readin(variable);
	region_dics ={'All':np.empty((5,3)),'TROPICS':np.empty((5,3)),'ASIA':np.empty((5,3)),'EUROPE':np.empty((5,3)),'US':np.empty((5,3)),'ARCTIC':np.empty((5,3))} 
	for region_key in region_dics:
		print region_key
		mask = mask_weight(region_key,lon,lat); 
		region_dics[region_key][0,:] = dif_std_for_region(F_EdgRef,F_1970,mask)
		region_dics[region_key][1,:] = dif_std_for_region(F_Edg70Oz,F_Edg70GO,mask)
		region_dics[region_key][2,:] = dif_std_for_region(F_Edg70GO,F_1970,mask)
		region_dics[region_key][3,:] = dif_std_for_region(F_EdgRef, F_Edg70T10SOZ,mask)
		region_dics[region_key][4,:] = dif_std_for_region(F_Edg70T10SOZ, F_Edg70Oz,mask)
	return region_dics
		
	
def plot_burden_bar(data):
	# print data
	# print "------"+var+"--------"
	def autolabel(X_pos,values,height_lift):
		## Attach a text label above each bar displaying its height
		height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
		ax.text(X_pos,y_pos,'%4.2f' % height, ha='center', va='bottom',size=10)
	ax.axhline(y=0,color='k',linewidth = 2);ax.set_xticklabels([])
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
	X_pos=0;dis=0.5
	rects1=ax.bar(X_pos, data[0,0], yerr=data[0,1:3].reshape((2,1)), align='center',color='k', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos,data[0,0],1.1)

	rects2=ax.bar(X_pos+dis, data[1,0], yerr=data[1,1:3].reshape((2,1)),  align='center',color='g', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis,data[1,0],1.1)

	rects3=ax.bar(X_pos+dis*2, data[2,0], yerr=data[2,1:3].reshape((2,1)),  align='center',color='r', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*2,data[2,0],1.1)
	
	rects4=ax.bar(X_pos+dis*3, data[3,0], yerr=data[3,1:3].reshape((2,1)),  align='center',color='m', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*3,data[3,0],1.1)
	
	rects5=ax.bar(X_pos+dis*4, data[4,0], yerr=data[4,1:3].reshape((2,1)),  align='center',color='b', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*4,data[4,0],1.1)
	return rects1,rects2,rects3,rects4,rects5
	
	

index ='ERFT';region_dics = spa_pat_reg_mean(index);	
fig = plt.figure(facecolor='White',figsize=[10,6]);pad= 5
ax = plt.subplot(2,3,1);
ax.annotate('(a) Globe',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
plot_burden_bar(region_dics['All'][:])
ax.set_ylim([-0.5,2.5]);

ax = plt.subplot(2,3,2);
ax.annotate('(b) Tropics (28S - 28N)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
plot_burden_bar(region_dics['TROPICS'][:])
ax.set_ylim([-0.5,2.5]);

ax = plt.subplot(2,3,3);
ax.annotate('(c) Arctic (60N - 90N)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
plot_burden_bar(region_dics['ARCTIC'][:])
ax.set_ylim([-0.5,2.5]);

ax = plt.subplot(2,3,4);
ax.annotate('(d) Asia',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
plot_burden_bar(region_dics['ASIA'][:])
ax.set_ylim([-2,6]);

ax = plt.subplot(2,3,5);
ax.annotate('(e) Europe',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
plot_burden_bar(region_dics['EUROPE'][:])
ax.set_ylim([-2,6]);

ax = plt.subplot(2,3,6);
ax.annotate('(f) USA',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
plot_burden_bar(region_dics['US'][:])
ax.set_ylim([-2,6]);

ax = fig.add_axes([0.12, 0.01, 0.80, 0.10])

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
k_patch = mpatches.Patch(color='k', alpha=0.7,label='All')
g_patch = mpatches.Patch(color='g', alpha=0.7,label='GHGs')
r_patch = mpatches.Patch(color='r', alpha=0.7,label='AAs')
m_patch = mpatches.Patch(color='m', alpha=0.7,label='Trop. O3')
b_patch = mpatches.Patch(color='b', alpha=0.7,label='Strat. O3')
lines = [k_patch,g_patch,r_patch,m_patch,b_patch]
labels = ['All','GHGs','AAs','Trop. O3','Strat. O3']
legend = plt.legend(lines,labels,ncol=5,loc='lower left',labelspacing=0.2)
legend.get_frame().set_facecolor('none');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)
ax.axis('off')

plt.subplots_adjust(left=0.08, bottom=0.12, right=0.95, top=0.95, wspace=0.25, hspace=0.25); 
plt.savefig('ERF_Regional_mean.png', format='png', dpi=1000)	
# plt.show()