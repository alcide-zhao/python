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
	# plt.imshow(mask,origin='lower');plt.show()
	mask[mask >= 1] = 1;mask[mask < 1] = 0;
	# weight each grid cell by its area weight against the total area
	if reverse:
		mask=1-mask
	mask[mask==0] = np.nan
	mask=np.multiply(mask,area); 
	mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
	return mask_weighted

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

	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
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
		lon,lat,PRECC = data_netcdf('Edg70GO','PRECC');_,_,PRECL = data_netcdf('Edg70GO','PRECL');Edg70GO=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('T1970RCP','PRECC');_,_,PRECL = data_netcdf('T1970RCP','PRECL');T1970=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('EdgRef','PRECC');_,_,PRECL = data_netcdf('EdgRef','PRECL');EdgRef=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('Edg70Oz','PRECC');_,_,PRECL = data_netcdf('Edg70Oz','PRECL');Edg70Oz=(PRECC+PRECL)*24*60*60*1000
		_,_,PRECC = data_netcdf('Edg70T10SOZ','PRECC');_,_,PRECL = data_netcdf('Edg70T10SOZ','PRECL');Edg70T10SOZ=(PRECC+PRECL)*24*60*60*1000
	else:
		lon,lat,Edg70GO = data_netcdf('Edg70GO',variable)
		_,_,T1970 = data_netcdf('T1970RCP',variable)
		_,_,EdgRef = data_netcdf('EdgRef',variable)
		_,_,Edg70Oz = data_netcdf('Edg70Oz',variable)
		_,_,Edg70T10SOZ = data_netcdf('Edg70T10SOZ',variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ
	
def print_domain_mean(variable):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""
	
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
		
	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ = data_readin(variable);
	Total,Total_s = diff_sig(EdgRef,T1970)
	GHG,GHG_s = diff_sig(Edg70Oz,Edg70GO)
	AEROSOL,AEROSOL_s = diff_sig(Edg70GO,T1970)
	TrO3,TrO3_s = diff_sig(EdgRef, Edg70T10SOZ)
	StO3,StO3_s = diff_sig(Edg70T10SOZ, Edg70Oz)
	REF = np.nanmean(T1970,axis=0)
	def dif_std_for_region(var1,var2,mask):
		"""
		for a rgional sif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
		var2_domain_mean =np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
		dif =  round(np.nanmean(var1_domain_mean - var2_domain_mean,axis=0),4)
		_,p_value = student_test(var1_domain_mean , var2_domain_mean)
		if p_value<0.05:sig=1
		else: sig=0
		return dif,sig
	
	# region_keys =['All','GloLand','ASIA','EUROPE','US']  #,'USA','SESA','EA','India'
	# for region_key in region_keys:
		# # print region_key
		# mask = mask_weight(region_key,lon,lat); #print mask
		# print "=========="+region_key+"=========="
		# print 'Total2010 ',dif_std_for_region(EdgRef,T1970,mask)
		# print 'GHG       ',dif_std_for_region(Edg70Oz,Edg70GO,mask)
		# print 'AEROSOL   ',dif_std_for_region(Edg70GO,T1970,mask)
		# print 'OZONE     ',dif_std_for_region(EdgRef,Edg70Oz,mask)
		# print 'Trop. O3  ',dif_std_for_region(EdgRef, Edg70T10SOZ,mask)
		# print 'Strat.O3  ',dif_std_for_region(Edg70T10SOZ, Edg70Oz,mask)	
	return lon,lat,REF,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s

	
def plot_zonal_mean_uncertainty(ax,y1,y2,x,ylim):
	ax.tick_params(axis='both', which='major', direction = 'in',left=True,right=True,bottom=True,top=True, pad=5)
	ax.axvline(x=0,color='k',linewidth = 1);ax.axhline(y=0,color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks(np.arange(-90,91,30));ax.set_xticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	ax.set_ylim([-ylim,ylim]);ax.set_yticks(np.arange(-xlim,xlim+.0000001,xlim/3.0))
	ax.plot(x,np.nanmean(y1,axis=1),'-',color="k",linewidth=2,label = 'Precip')
	ax.plot(x,np.nanmean(y2,axis=1),'-',color="r",linewidth=2,label = 'Scaled precip')
	return ax
	
index ='TS';lon,lat,REF_TS,Total_TS,GHG_TS,AEROSOL_TS,TrO3_TS,StO3_TS,Total_s_TS,GHG_s_TS,AEROSOL_s_TS,TrO3_s_TS,StO3_s_TS = print_domain_mean(index);		
index ='precip';lon,lat,REF_Pr,Total_Pr,GHG_Pr,AEROSOL_Pr,TrO3_Pr,StO3_Pr,Total_s_Pr,GHG_s_Pr,AEROSOL_s_Pr,TrO3_s_Pr,StO3_s_Pr = print_domain_mean(index);				

def global_mean(var):
	mask = mask_weight('All',lon,lat,reverse=False)
	mask = mask/np.nansum(np.nansum(mask,axis=1),axis=0)
	spatial_mean = round(np.nansum(np.nansum(np.multiply(var,mask),axis=1),axis=0),4)
	return str(spatial_mean)

fig = plt.figure(facecolor='White',figsize=[9,11]);pad= 5;colormap='BrBG';


####################ERF###########################
colorbar_min=-0.4;colorbar_max=0.4;
ax = plt.subplot(4,2,1);
ax.annotate('GHGs',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(a)                                         '+global_mean(GHG_Pr),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('Precipitation',xy=(0.5,1.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	

				
colormesh0 = spatial_figure(ax,GHG_Pr,GHG_s_Pr,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,3);
ax.annotate('AAs',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(b)                                        '+global_mean(AEROSOL_Pr),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,AEROSOL_Pr,AEROSOL_s_Pr,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,5);
ax.annotate('Trop. O3',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(c)                                         '+global_mean(TrO3_Pr),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TrO3_Pr,TrO3_s_Pr,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,7);
ax.annotate('Strat. O3',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(d)                                        '+global_mean(StO3_Pr),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,StO3_Pr,StO3_s_Pr,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.06, 0.06, 0.44, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/mm\/day^{-1}}$',xy=(0.5,-3.2), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)


####################TS###########################
from matplotlib.colors import LinearSegmentedColormap
colors = [(151.0/255, 151.0/255, 151.0/255), (240/255.0,128/255.0,128/255.0), (205/255.0,92/255.0,92/255.0), (220/255.0,20/255.0,60/255.0), (0.8, 0, 0)]  # R -> G -> B
colormap = LinearSegmentedColormap.from_list('my_list', colors, N=10)	
colorbar_min=0;colorbar_max=10;
def spatial_figure(axs,data,p_value,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
	lons[lons>180]-=360; 
	lon_0 = lons.mean(); 
	lat_0 = lats.mean(); 
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	map.drawcoastlines(); #map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	colormap.set_bad([1,1,1],alpha = 1.0); colormap.set_over((0.4, 0, 0)); colormap.set_under('white');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=colormap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh
	
def scale_P(dPr,dT,GHG_Pr,GHG_TS):
	def global_mean_Num(var):
		mask = mask_weight('All',lon,lat,reverse=False)
		spatial_mean = np.nansum(np.nansum(np.multiply(var,mask),axis=1),axis=0)
		return spatial_mean
	dPrP = np.divide(dPr,REF_Pr)*100;GHG_PrP = np.divide(GHG_Pr,REF_Pr)*100
	ratio = np.divide(np.divide(dPrP,dT),np.divide(GHG_PrP,GHG_TS))
	dT_mean = global_mean_Num(dT);GHG_TS_mean = global_mean_Num(GHG_TS)
	dPrPm = np.divide(global_mean_Num(dPr),global_mean_Num(REF_Pr))*100;
	GHG_PrPm = np.divide(global_mean_Num(GHG_Pr),global_mean_Num(REF_Pr))*100;
	Globalmean = str(round(np.divide(np.divide(dPrPm,dT_mean),np.divide(GHG_PrPm,GHG_TS_mean)),1))
	# Globalmean =str(round(global_mean_Num(ratio),1))
	
	return ratio,Globalmean

AEROSOL_PrS,AEROSOL_PrG= scale_P(AEROSOL_Pr,AEROSOL_TS,GHG_Pr,GHG_TS);
TrO3_PrS,TrO3_PrG= scale_P(TrO3_Pr,TrO3_TS,GHG_Pr,GHG_TS);
StO3_PrS,StO3_PrG= scale_P(StO3_Pr,StO3_TS,GHG_Pr,GHG_TS);

sig = GHG_s_TS;sig[:] =np.nan

ax = plt.subplot(4,2,2);
ax.annotate('Ratio of sensitivity (dP/dT)\n relative to GHGs',xy=(0.5,1.0), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	
ax.axis('off')			
ax = plt.subplot(4,2,4);
ax.annotate('(e)                                                '+AEROSOL_PrG,xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,AEROSOL_PrS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,6);
ax.annotate('(f)                                                 '+TrO3_PrG,xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TrO3_PrS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,8);
ax.annotate('(g)                                               '+StO3_PrG,xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,StO3_PrS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.54, 0.06, 0.44, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate('#',xy=(0.5,-3.2), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.92, wspace=0.1, hspace=0.20); 
plt.savefig('Pr_scaled_sensitivity.png', format='png', dpi=1000)
