"""
This is to plot the T and P extreme indices calculated from the EDGAR experiments
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
    os.path.pardir,os.path.pardir,
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
	mask[:,0:colum_s] =0; mask[:,colum_e:-1] =0
	# plt.imshow(mask,origin='lower');plt.show()
	mask[0:row_s,:] =0; mask[row_e:-1,:] =0
	# plt.imshow(mask,origin='lower');plt.show()
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
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	axs.set_ylim([-60,90])
	return colormesh
		
def data_readin_momthly(variable):	
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
	ref = T1970
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ,ref
	
def data_readin_extreme(variable,scale,exp):	
	def data_netcdf(scenario,scale,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/analysis/'+scale+'/'
		var_path = input_path+scenario+'.TempPrecp.extremes.'+scale+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		data = nc_fid.variables[variable][0:30,:,:]
		nc_fid.close()
		return lon,lat,data
	if exp =='FullCp':
		lon,lat,Edg70GO = data_netcdf('Edg70GO',scale,variable)
		_,_,T1970 = data_netcdf('T1970',scale,variable)
		_,_,EdgRef = data_netcdf('EdgRef',scale,variable)
		_,_,Edg70Oz = data_netcdf('Edg70Oz',scale,variable)
		_,_,Edg70T10SOZ = data_netcdf('Edg70T10SOZ',scale,variable)
	elif exp =='Fast':
		lon,lat,Edg70GO = data_netcdf('F_Edg70GO',scale,variable)
		_,_,T1970 = data_netcdf('F_T1970',scale,variable)
		_,_,EdgRef = data_netcdf('F_EdgRef',scale,variable)
		_,_,Edg70Oz = data_netcdf('F_Edg70Oz',scale,variable)
		_,_,Edg70T10SOZ = data_netcdf('F_Edg70T10SOZ',scale,variable)
	elif exp =='Slow':
		lon,lat,Edg70GO = data_netcdf('SlowResp_Edg70GO',scale,variable)
		_,_,T1970 = data_netcdf('SlowResp_T1970',scale,variable)
		_,_,EdgRef = data_netcdf('SlowResp_EdgRef',scale,variable)
		_,_,Edg70Oz = data_netcdf('SlowResp_Edg70Oz',scale,variable)
		_,_,Edg70T10SOZ = data_netcdf('SlowResp_Edg70T10SOZ',scale,variable)
	elif exp == 'FCminusFast':
		lon,lat,T_Edg70GO = data_netcdf('Edg70GO',scale,variable);lon,lat,F_Edg70GO = data_netcdf('F_Edg70GO',scale,variable); Edg70GO=T_Edg70GO-F_Edg70GO
		_,_,T_T1970 = data_netcdf('T1970',scale,variable);_,_,F_T1970 = data_netcdf('F_T1970',scale,variable);T1970=T_T1970-F_T1970
		_,_,T_EdgRef = data_netcdf('EdgRef',scale,variable);_,_,F_EdgRef = data_netcdf('F_EdgRef',scale,variable);EdgRef=T_EdgRef-F_EdgRef
		_,_,T_Edg70Oz = data_netcdf('Edg70Oz',scale,variable);_,_,F_Edg70Oz = data_netcdf('F_Edg70Oz',scale,variable);Edg70Oz=T_Edg70Oz-F_Edg70Oz
		_,_,T_Edg70T10SOZ = data_netcdf('Edg70T10SOZ',scale,variable);_,_,F_Edg70T10SOZ = data_netcdf('F_Edg70T10SOZ',scale,variable);Edg70T10SOZ=T_Edg70T10SOZ-F_Edg70T10SOZ
	_,_,ref = data_netcdf('T1970',scale,variable);
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ,ref
	
def individual_forcing_maps(exp,scale,var,percent):

	####calculate the difference and their significance
		
	def diff_sig(vairable1,variable2,ref,percent):
		def mannwhitneyu_test(vairable1,variable2):
			p_threshold=0.10
			size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
			p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
			from scipy.stats import mannwhitneyu as test
			for x in range(size[0]):
				for y in range(size[1]):
					cache1 = vairable1[:,x,y]
					cache2 = variable2[:,x,y]
					if np.array_equal(cache1,cache2): p_value[x,y] = np.nan;
					else: _,p_value[x,y] = test(cache1,cache2);
			p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
			return p_value
		if percent:
			dif = np.divide(np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0),np.nanmean(ref,axis=0))*100
		else:
			dif =np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)	
		sig = mannwhitneyu_test(vairable1, variable2)
		return dif,sig

	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ,ref = data_readin_extreme(var,scale,exp);
	lon,lat,T1970_T,Edg70GO_T,Edg70Oz_T,EdgRef_T,Edg70T10SOZ_T,ref_T = data_readin_extreme(var,scale,'FullCp');
	Total,Total_s = diff_sig(EdgRef,T1970,T1970_T,percent)
	GHG,GHG_s = diff_sig(Edg70Oz,Edg70GO,Edg70GO_T,percent)
	AEROSOL,AEROSOL_s = diff_sig(Edg70GO,T1970,T1970_T,percent)
	TrO3,TrO3_s = diff_sig(EdgRef, Edg70T10SOZ,Edg70T10SOZ_T,percent)
	StO3,StO3_s = diff_sig(Edg70T10SOZ, Edg70Oz,Edg70Oz_T,percent)
	return lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s,np.nanmean(ref,axis=0)



def plot_SpaPat(var,scale,percent):
	if var in ['TXX','TNN','TSM','TSSTD']:
		colormap='RdBu';colormap = reverse_colourmap(colormap);
		colorbar_min=-2;colorbar_max=2;
	elif var in ['RX5DAY','R10','CWD','CDD','SDII','RX1DAY','R95','R99','PRECPM','R20']:
		colormap='BrBG';
		if var == 'CDD': 
			colormap = reverse_colourmap(colormap);
		if var in ['RX5DAY','SDII','RX1DAY','R10','R20','PRECPM']:
			colorbar_min=-25;colorbar_max=25;
		elif var in ['R95','R99']:
			colorbar_min=-3;colorbar_max=3;
		else:
			colorbar_min=-5;colorbar_max=5;
	else:
		print "make sure of a right extreme variable"

	exp='FullCp';   lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s,ref= individual_forcing_maps(exp,scale,var,percent);
		
	fig = plt.figure(facecolor='White',figsize=[15,8]);pad= 5;
	ax = plt.subplot(3,3,1);
	ax.annotate('Total',xy=(0.5,1.1), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='center', va='baseline',rotation='horizontal',fontsize=20)
	ax.annotate('GHGs',xy=(-0.1,0.5), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='center',rotation='vertical',fontsize=15)	
	ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)			
	colormesh0 = spatial_figure(ax,GHG,GHG_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	ax = plt.subplot(3,3,4);
	ax.annotate('AAs',xy=(-0.1,0.5), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='center',rotation='vertical',fontsize=15)	
	ax.annotate('(b)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)					
	spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	ax = plt.subplot(3,3,7);
	ax.annotate('Trop. O3',xy=(-0.1,0.5), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='center',rotation='vertical',fontsize=15)	
	ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	colormesh0=spatial_figure(ax,TrO3,TrO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	# ax = plt.subplot(3,3,10);
	# ax.annotate('Strat. O3',xy=(-0.1,0.5), xytext=(0, pad),
					# xycoords='axes fraction', textcoords='offset points',
					# ha='left', va='center',rotation='vertical',fontsize=15)	
	# ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
					# xycoords='axes fraction', textcoords='offset points',
					# ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	# colormesh0=spatial_figure(ax,StO3,StO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	exp='Fast';  lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s,ref= individual_forcing_maps(exp,scale,var,percent);
		
	ax = plt.subplot(3,3,2);
	ax.annotate('Fast',xy=(0.5,1.1), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='center', va='baseline',rotation='horizontal',fontsize=20)	
	ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)			
	colormesh1 = spatial_figure(ax,GHG,GHG_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	ax = plt.subplot(3,3,5)
	ax.annotate('(e)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)					
	colormesh1 = spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	ax = plt.subplot(3,3,8);
	ax.annotate('(f)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	colormesh0=spatial_figure(ax,TrO3,TrO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	# ax = plt.subplot(3,3,11);
	# ax.annotate('(h)',xy=(0.02,1.01), xytext=(0, pad),
					# xycoords='axes fraction', textcoords='offset points',
					# ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	# colormesh0=spatial_figure(ax,StO3,StO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
	
	
	
	exp='FCminusFast';   lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s,ref= individual_forcing_maps(exp,scale,var,percent);
		
	ax = plt.subplot(3,3,3);
	ax.annotate('Total - Fast',xy=(0.5,1.1), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='center', va='baseline',rotation='horizontal',fontsize=20)	
	ax.annotate('(g)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	colormesh2 = spatial_figure(ax,GHG,GHG_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	ax = plt.subplot(3,3,6)
	ax.annotate('(h)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)					
	colormesh2 = spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	ax = plt.subplot(3,3,9);
	ax.annotate('(i)',xy=(0.02,1.01), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	colormesh2=spatial_figure(ax,TrO3,TrO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

	# ax = plt.subplot(3,3,12);
	# ax.annotate('(l)',xy=(0.02,1.01), xytext=(0, pad),
					# xycoords='axes fraction', textcoords='offset points',
					# ha='left', va='baseline',rotation='horizontal',fontsize=15)		
	# colormesh0=spatial_figure(ax,StO3,StO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
	
	cbar_ax = fig.add_axes([0.18, 0.04, 0.64, 0.015])
	char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

	plt.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.95, wspace=0.04, hspace=0.15); 
	plt.savefig(var+'_Fast_slow.png', format='png', dpi=1000)
	

#########################
#### Spatial Maps    ####
#########################
var = 'PRECPM'; scale='annual';
if var in ['TNN','TSM','TSSTD','CDD','R95','R99','CWD','TXX']:
	percent =False
elif var in ['RX5DAY','SDII','R10','RX1DAY','PRECPM','R20']: 
	percent =True
plot_SpaPat(var,scale,percent)     #TXX TNN PRECPM RX1DAY RX5DAY R95


#################################
####  regional mean ratios   ####
#################################
"""	
# # var_dic = ['TXX','TNN','PRECPM','RX5DAY','RX1DAY'];region_key='All';
# # for var in var_dic:
	# # print "------------"+var+'-------------------'
	# # PorT ='TS';print P_vs_T(var,PorT,region_key)
	# # if var =='RX5DAY': 
		# # PorT ='PRECIP';
		# # print "------------Against P-----------------"
def individual_forcing_RegMean(exp,scale,var,percent,region_key):
	###
	# This module is to calculate the percentage/absoluta changes in a variable that has to be an area-weighted mean 
	###
	def isolate_individual_forcing(var1,var2,baseline,percent):
		mask = mask_weight(region_key,lon,lat);
		if percent:
			baseline = np.nanmean(np.nansum(np.nansum(np.multiply(mask,baseline),axis=2),axis=1),axis=0);#print baseline
			diff = 100*np.divide(np.nansum(np.nansum(np.multiply(mask,var1-var2),axis=2),axis=1),baseline)
		else:
			diff = np.nansum(np.nansum(np.multiply(mask,var1-var2),axis=2),axis=1)
		# mean = np.nanmean(diff);p25=np.nanpercentile(diff,25);p75=np.nanpercentile(diff,75)
		return diff
	
	# lon,lat,T1970_T,Edg70GO_T,Edg70Oz_T,EdgRef_T,Edg70T10SOZ_T,ref_T = data_readin_extreme(var,scale,'FullCp');	
	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ,ref = data_readin_extreme(var,scale,exp);
	
	# Total = isolate_individual_forcing(EdgRef,T1970,ref,percent)
	GHG = isolate_individual_forcing(Edg70Oz,Edg70GO,ref,percent)
	AAs = isolate_individual_forcing(Edg70GO,T1970,ref,percent)
	O3 = isolate_individual_forcing(EdgRef, Edg70Oz,ref,percent)
	# StrO3 = isolate_individual_forcing(Edg70T10SOZ, Edg70Oz,ref,percent)
	return GHG,AAs,O3

def plot_cross(Total,fast,slow):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
	ax.axvline(x=2.0,color='k',linewidth = 0.2);
	ax.axvline(x=4,color='k',linewidth = 0.2);
	ax.axhline(y=0,color='k',linewidth = 1);
	x_pos = np.arange(0.5,6.0,2);
	y_mean=np.array([Total[0,0], fast[0,0],slow[0,0]])
	p25=np.array([Total[0,1], fast[0,1],slow[0,1]]);p75=np.array([Total[0,2], fast[0,2],slow[0,2]])
	ax.errorbar(x_pos, y_mean, [p25, p75], color='g',fmt='gx', lw=2,capsize=0,markeredgecolor='g',ms=10,mew=2,label='GHGs')
	x_pos = np.arange(1,6.0,2);
	y_mean=np.array([Total[1,0], fast[1,0],slow[1,0]])
	p25=np.array([Total[1,1], fast[1,1],slow[2,1]]);p75=np.array([Total[1,2], fast[1,2],slow[1,2]])
	ax.errorbar(x_pos, y_mean, [p25, p75], color='r',fmt='rx', lw=2,capsize=0,markeredgecolor='r',ms=10,mew=2,label='AAs')
	x_pos = np.arange(1.5,6.0,2);
	y_mean=np.array([Total[2,0], fast[2,0],slow[2,0]])
	p25=np.array([Total[2,1], fast[2,1],slow[2,1]]);p75=np.array([Total[2,2], fast[2,2],slow[2,2]])
	ax.errorbar(x_pos, y_mean, [p25, p75], color='b',fmt='bx', lw=2,capsize=0,markeredgecolor='b',ms=10,mew=2,label='Ozone')		
	
	ax.set_xticks(np.arange(1,6.0,2));ax.set_xticklabels(['Total','Fast','Slow'])
	ax.set_xlim([0,6.0]);
	return ax	
	
def P_vs_T(var,PorT,region_key):
	def get_ratios(exp,var,PorT,region_key):
		if exp == 'FCminusFast':
			if PorT == 'PRECIP':  ## is normalize against P, then p should in percent
				GHG,AAs,O3 = individual_forcing_RegMean('FullCp','annual','PRECPM',True,region_key)
			elif PorT == 'TS':  ## if normmalize against T, T is only in degree celcius
				GHG,AAs,O3 = individual_forcing_RegMean('FullCp','annual','TSM',False,region_key)
			Xvar = np.stack([GHG,AAs,O3],axis=0);
		else:
			if PorT == 'PRECIP':  ## is normalize against P, then p should in percent
				GHG,AAs,O3 = individual_forcing_RegMean('FullCp','annual','PRECPM',True,region_key)
			elif PorT == 'TS':  ## if normmalize against T, T is only in degree celcius
				GHG,AAs,O3 = individual_forcing_RegMean('FullCp','annual','TSM',False,region_key)
			Xvar = np.stack([GHG,AAs,O3],axis=0);
		if var in ['TXX','TNN']:  # for temperature related indices, we are more interested in the absolute chanegs
			GHG,AAs,O3 = individual_forcing_RegMean(exp,'annual',var,False,region_key)
		else:   # for precipitation related indices, percent change is mroe interesting
			GHG,AAs,O3 = individual_forcing_RegMean(exp,'annual',var,True,region_key)
		Yvar = np.stack([GHG,AAs,O3],axis=0);
		ratio = np.around(np.divide(Yvar,Xvar),4)
		return ratio
	
	ratio_FuCp = get_ratios('FullCp',var,PorT,region_key);
	ratio_FmiF = get_ratios('FCminusFast',var,PorT,region_key)
	# ratio_Fast = get_ratios('Fast',var,PorT,region_key)
	ratio_Fast = ratio_FuCp - ratio_FmiF
	mean = np.nanmean(ratio_FuCp,axis=1);p25=abs(mean-np.nanpercentile(ratio_FuCp,25,axis=1))/1.25;p75=abs(-mean+np.nanpercentile(ratio_FuCp,75,axis=1))/1.25
	Total = np.stack([mean,p25,p75],axis=1)
	mean = np.nanmean(ratio_FmiF,axis=1);p25=abs(mean-np.nanpercentile(ratio_FmiF,25,axis=1))/1.25;p75=abs(-mean+np.nanpercentile(ratio_FmiF,75,axis=1))/1.25
	Slow = np.stack([mean,p25,p75],axis=1)
	mean = np.nanmean(ratio_Fast,axis=1);p25=abs(mean-np.nanpercentile(ratio_Fast,25,axis=1))/1.25;p75=abs(-mean+np.nanpercentile(ratio_Fast,75,axis=1))/1.25
	fast = np.stack([mean,p25,p75],axis=1)
	return Total,fast,Slow		# print P_vs_T(var,PorT,region_key)
	
fig = plt.figure(facecolor='White',figsize=[10,8]);pad= 5;
region_key='GloLand'
ax = plt.subplot(3,2,3);
ax.annotate('(a) TXX against SAT (K/K)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
Total,fast,slow = P_vs_T(var='TXX',PorT='TS',region_key=region_key);ax=plot_cross(Total,fast,slow)
legend = ax.legend(shadow=False,ncol=1,bbox_to_anchor=(0.7,2.2))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

ax = plt.subplot(3,2,5);
ax.annotate('(b) TNN against SAT (K/K)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
Total,fast,slow = P_vs_T(var='TNN',PorT='TS',region_key=region_key);plot_cross(Total,fast,slow)

ax = plt.subplot(3,2,2);
ax.annotate('(c) Pr mean against SAT (%/K)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
Total,fast,slow = P_vs_T(var='PRECPM',PorT='TS',region_key=region_key);plot_cross(Total,fast,slow)

ax = plt.subplot(3,2,4);
ax.annotate('(d) RX5DAY against SAT (%/K)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
Total,fast,slow = P_vs_T(var='RX5DAY',PorT='TS',region_key=region_key);plot_cross(Total,fast,slow)

ax = plt.subplot(3,2,6);
ax.annotate('(e) RX5DAY against Pr mean (%/%)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
Total,fast,slow = P_vs_T(var='RX5DAY',PorT='PRECIP',region_key=region_key);ax = plot_cross(Total,fast,slow)

plt.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.95, wspace=0.20, hspace=0.40); 
plt.savefig(region_key+'sensitivity_fast_vs_slow.png', format='png', dpi=500)

"""