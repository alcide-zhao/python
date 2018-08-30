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
	
def box_clip(lon_S,lon_e,lat_S,lat_e,lon,lat,mask):
	"""
	fill the range outside the box with 0
	"""
	lon = np.array(lon)
	lat = np.array(lat)
	colum_S = [index for index in range(len(lon)) if np.abs(lon-lon_S)[index] == np.min(np.abs(lon-lon_S))][0]
	colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
	row_S = [index for index in range(len(lat)) if np.abs(lat-lat_S)[index] == np.min(np.abs(lat-lat_S))][0]
	row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
	if (colum_S> colum_e):
		cache = colum_e; colum_e = colum_S; colum_S = cache;
	if (row_S> row_e):
		cache = row_e; row_e = row_S; row_S = cache;
	mask[:,:colum_S] =0; mask[:,colum_e:] =0
	mask[:row_S,:] =0; mask[row_e:,:] =0
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
	box_region_Sic={'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
	if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'GloLand'):
		mask= ocean_mask[region_key][:]
	elif  region_key in box_region_Sic:
		mask= ocean_mask['All'][:]
		box = box_region_Sic[region_key]
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
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_over('k'); cmap.set_under('darkblue');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh
		
def data_readin(variable,region_key,exp):	
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		if scenario.startswith('F_'):
			start_year =2000
		elif scenario =='T1970C': 
			start_year =1970
		else: 
			start_year =2010
		start =(start_year*365)
		ith=0	
		for iday in days:
			month_Says =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
			calendar_Says = np.array([0,31,59,90,120,151,181,212,243,273,304,334,365])
			total_Says = int(iday) + start; 
			year = total_Says//365; 
			remainder =  total_Says%365
			if remainder ==0: year=year-1;month=12;day=31
			else: 
				month = 1+[layer for layer in range(len(calendar_Says)) if calendar_Says[layer]< remainder and calendar_Says[layer+1]>=remainder][0]
				day = int(remainder - calendar_Says[month-1])
				if day == 0: day = month_Says[month-1]
			date_int[ith] = year*10000+month*100+day
			ith=ith+1
		return date_int.astype(int)

	def mon_mean2annual_mean(scenario,time,spatial_mean,opt='AnuMean'):
		GloAnu_mean=np.empty((30));GloAnu_mean[:]=np.nan
		if scenario.startswith('F_'):
			year_Series = range(2000,2030)
		elif scenario=='T1970RCP':
			year_Series = range(2020,2050)
		elif scenario=='EdgEne':
			year_Series = range(2200,2230)
		elif scenario=='Edg70GO':
			year_Series = range(2070,2100)
		else:
			year_Series = range(2130,2160)
		for iyear in year_Series:
			if (iyear == year_Series[0] and time[0]//100 >= year_Series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  #June01
			if (iyear == year_Series[-1] and time[-1]//100 <= year_Series[-1] *100+12):
				data_cache = spatial_mean[layer_b:]#print data_cache
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+12][0]  #August 31
				data_cache = spatial_mean[layer_b:layer_e+1]#print data_cache
			if opt == 'AnuMean':
				GloAnu_mean[iyear-year_Series[0]] = np.nanmean(data_cache,axis=0);
			elif opt == 'AnuSum':
				month_Says =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
				GloAnu_mean[iyear-year_Series[0]] = np.nansum(np.multiply(data_cache,month_Says),axis=0);
		return GloAnu_mean

	def data_netcdf(exp,scenario,variable,region_key):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'+exp+'/'
		if variable == 'precip':
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.PRECC.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			PRECC = (nc_fid.variables['PRECC'][:])*24*60*60*1000
			nc_fid.close()
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.PRECL.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			PRECL = (nc_fid.variables['PRECL'][:])*24*60*60*1000			
			nc_fid.close()
			data = PRECC+PRECL
		else:
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
			# print var_path
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			data = nc_fid.variables[variable][:]
			nc_fid.close()
		if variable == "TS": data=data-273.15
		elif variable == 'CLDTOT':data=data*100
		elif variable == 'TGCLDLWP':data=data*10**(3)
		mask = mask_weight(region_key,lon,lat,reverse=False);
		spatial_mean = np.nansum(np.nansum(np.multiply(data,mask),axis=2),axis=1)
		if variable == 'precip':
			GloAnu_mean = mon_mean2annual_mean(scenario,time,spatial_mean,opt='AnuSum')
		else: GloAnu_mean = mon_mean2annual_mean(scenario,time,spatial_mean,opt='AnuMean')
		return GloAnu_mean
		
	if exp=='CamOnly':
		Edg70GO = data_netcdf(exp,'F_Edg70GO',variable,region_key);
		Ref1970 = data_netcdf(exp,'F_1970',variable,region_key);
		EdgRef = data_netcdf(exp,'F_EdgRef',variable,region_key);
		Edg70Oz = data_netcdf(exp,'F_Edg70Oz',variable,region_key);
		Edg70T10SOZ = data_netcdf(exp,'F_Edg70T10SOZ',variable,region_key);
		EdgEne = data_netcdf(exp,'F_EdgEne',variable,region_key);
		EdgTech = data_netcdf(exp,'F_EdgTech',variable,region_key);
	elif exp=='FullCp':
		Edg70GO = data_netcdf(exp,'Edg70GO',variable,region_key)
		Ref1970 = data_netcdf(exp,'T1970RCP',variable,region_key)
		EdgRef = data_netcdf(exp,'EdgRef',variable,region_key)
		Edg70Oz = data_netcdf(exp,'Edg70Oz',variable,region_key)
		Edg70T10SOZ = data_netcdf(exp,'Edg70T10SOZ',variable,region_key)
		EdgEne = data_netcdf(exp,'EdgEne',variable,region_key);
		EdgTech = data_netcdf(exp,'EdgTech',variable,region_key);		
	return Ref1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ,EdgEne,EdgTech
	
def print_FaSl_GloAnuMean(variable,region_key,percent):
	def decompose_forcings(Ref970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ,EdgEne,EdgTech,percent):
		Ref970_T,Edg70GO_T,Edg70Oz_T,EdgRef_T,Edg70T10SOZ_T,EdgEne_T,EdgTech_T = data_readin(variable,region_key,exp='FullCp');
		def mean_Std(var1,var2,ref,percent):
			diff = var1 - var2
			mean = np.nanmean(diff,axis=0)
			std = round(np.nanstd(diff,axis=0)/2,2)
			if percent:
				mean =100.0*mean/np.nanmean(ref,axis=0)
			return round(mean,2),std
		Total,Total_S = mean_Std(EdgRef,Ref970,Ref970_T,percent)
		GHG,GHG_S = mean_Std(Edg70Oz,Edg70GO,Edg70GO_T,percent)
		AAs,AAs_S = mean_Std(Edg70GO,Ref970,Ref970_T,percent)
		TrO3,TrO3_S = mean_Std(EdgRef, Edg70T10SOZ,Edg70T10SOZ_T,percent)
		StO3,StO3_S = mean_Std(Edg70T10SOZ, Edg70Oz,Edg70Oz_T,percent)
		Ene,Ene_S = mean_Std(EdgRef, EdgEne,EdgEne_T,percent)
		Tech,Tech_S = mean_Std(EdgRef, EdgTech,EdgTech_T,percent)		
		return np.stack((Total,GHG,AAs,TrO3,StO3,Ene,Tech),axis=0)
		
	Ref970_T,Edg70GO_T,Edg70Oz_T,EdgRef_T,Edg70T10SOZ_T,EdgEne_T,EdgTech_T = data_readin(variable,region_key,exp='FullCp');    # fully coupled total response
	Ref970_F,Edg70GO_F,Edg70Oz_F,EdgRef_F,Edg70T10SOZ_F,EdgEne_F,EdgTech_F = data_readin(variable,region_key,exp='CamOnly');    # local response as the dynamics are fixed
	
	Ref970_S = Ref970_T - Ref970_F;
	Edg70Oz_S = Edg70Oz_T - Edg70Oz_F;
	Edg70GO_S = Edg70GO_T - Edg70GO_F;
	EdgRef_S = EdgRef_T - EdgRef_F;
	Edg70T10SOZ_S = Edg70T10SOZ_T - Edg70T10SOZ_F;
	EdgEne_S = EdgEne_T - EdgEne_F;
	EdgTech_S = EdgTech_T - EdgTech_F;
	
	print "========================="+variable+"============================"
	print "           TOT   GHG   AAs  TrO3  StO3   Ene   Tech"
	print 'Fast : ',decompose_forcings(Ref970_F,Edg70GO_F,Edg70Oz_F,EdgRef_F,Edg70T10SOZ_F,EdgEne_F,EdgTech_F,percent)
	print 'Slow : ',decompose_forcings(Ref970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,Edg70T10SOZ_S,EdgEne_S,EdgTech_S,percent)
	print 'Total: ',decompose_forcings(Ref970_T,Edg70GO_T,Edg70Oz_T,EdgRef_T,Edg70T10SOZ_T,EdgEne_T,EdgTech_T,percent)
# VAR=['FSNT','FLNT','FSNS','FLNS','SHFLX','LHFLX','TS','precip']  #'FSNT','FLNT','FSNS','FLNS','SHFLX','LHFLX','TS',
VAR=['precip']
###############
##   main    ##
###############

for var in VAR:
	if var == 'precip':
		# print 'abs'
		# print_FaSl_GloAnuMean(var,'All',percent=False)
		print 'percent'
		print_FaSl_GloAnuMean(var,'All',percent=True)
	else:
		print_FaSl_GloAnuMean(var,'All',percent=False)	
		if var == 'TS':
			print 'GloLand'
			print_FaSl_GloAnuMean(var,'GloLand',percent=False)
		
