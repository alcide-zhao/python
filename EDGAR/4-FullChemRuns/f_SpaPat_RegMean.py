"""
This is to plot the spatial map and regional mean of air pollution variables from the full-chemistry runs
Input are daily variables, while it depends on the species what statisticcs to plot: 30 years mean of annual mean daily mean/max/percentile
	For ozone, the hourly surface concentration is processed into MDA8 (maximum of 8-hour running mean), 
		the units kg/kg can be converted to ppm with a wright of 10**(9)
	for PM2.5, all surface aerosol species have been summed up already, and the units is kg/kg can be converted to ug/m3 by 
		(kg/kg)*1.0e9/(1/1.185) is the conversion formula where 1.185 is the density of air under 25C and 1000hPa  
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
np.set_printoptions(threshold=np.inf)

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
	from matplotlib.colors import Normalize
	class MidpointNormalize(Normalize):
		def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
			self.midpoint = midpoint
			Normalize.__init__(self, vmin, vmax, clip)

		def __call__(self, value, clip=None):
			# I'm ignoring masked values and all kinds of edge cases to make a
			# simple example...
			x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
			return np.ma.masked_array(np.interp(value, x, y))
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
	cmap = discrete_cmap(8,colormap)
	cmap.set_bad('w');cmap.set_over('maroon'); cmap.set_under('darkgreen');  #darkgreen
	# cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_over('k'); cmap.set_under('green');  #darkgreen
	norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max) #
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True) #norm=norm
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='/', alpha=0.,lw=0.9,latlon=True)
	axs.set_ylim([-62,90])
	return colormesh
		
		
def data_readin(variable,stats_option):	
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		if scenario.startswith('F_'):
			start_year =2000
		elif scenario =='EdgRef70AP': start_year =1990
		elif scenario =='EdgTechP': start_year =2180
		elif scenario =='EdgRefP' or scenario =='EdgEneP'  or scenario =='Edg1970': start_year =2000
		elif scenario =='T1970C': 
			start_year =1970
		else: 
			start_year =2010
		start =(start_year*365)+1
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

	def daily_stats(scenario,time,data,stats_option):
		# stats_option determines eithr calculate annual mean, amx or percentile (must be an integer)
		annual_stats=np.empty((30,192,288));annual_stats[:]=np.nan
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario.startswith('F_') :
			year_series = range(2000,2030)
		elif scenario=='EdgRef70AP':
			year_series = range(2020,2050)
		elif scenario=='EdgTechP':
			year_series = range(2210,2240)	
		elif scenario =='EdgRefP' or scenario =='EdgEneP' or scenario =='Edg1970':
			year_series = range(2030,2060)
		elif scenario=='T1970RCP':
			year_series = range(2020,2050)
		elif scenario=='EdgEne':
			year_series = range(2200,2230)
		elif scenario=='Edg70GO':
			year_series = range(2070,2100)
		else:
			year_series = range(2130,2160)
		for iyear in year_series:
			if (iyear == year_series[0] and time[0] >= year_series[0] *10000+101):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer] == iyear *10000+101][0]  #June01
			if (iyear == year_series[-1] and time[-1] <= year_series[-1] *10000+1231):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]  == iyear*10000+1231][0]  #August 31
			data_cache = data[layer_b:layer_e+1,:,:]
			if stats_option == 'mean':
				annual_stats[iyear-year_series[0],:,:] = stats.nanmean(data_cache,axis=0)
			elif stats_option == 'max':
				annual_stats[iyear-year_series[0],:,:] = np.nanmax(data_cache,axis=0)
			else:
				annual_stats[iyear-year_series[0],:,:] = np.nanpercentile(data_cache,stats_option,axis=0)
		return annual_stats

	def data_netcdf(scenario,variable,stats_option):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/full_chem/'
		var_path = input_path+scenario+'/day/atm/'+scenario+'.atm.day.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		if variable in ['PM10_SRF','PM2.5.mass0.20','PM2.5_15']:
			data = nc_fid.variables['pm25'][:]*1.185*10**(9)  ## kgkg --> ug/m3
		elif variable == 'O3_SRF.MDA8':
			data = nc_fid.variables['O3_SRF'][:];
		else:
			data = nc_fid.variables[variable][:]
		if variable in ['CO_SRF','O3_SRF.MDA8']:
			data=data*10**(9)  ## mol/mol --> ppb
		elif variable =='NO_SRF':
			data=data*30/29*1.185*10**(9)   ## mol/mol--> ug/m3
		elif variable =='NO2_SRF':
			data=data*46/29*1.185*10**(9)   ## mol/mol--> ug/m3			
		annual_stats = daily_stats(scenario,time,data,stats_option)
		return lon,lat,annual_stats
		

	lon,lat,EdgRef70A = data_netcdf('EdgRef70AP',variable,stats_option)
	_,_,Ref1970 = data_netcdf('Edg1970',variable,stats_option)
	_,_,EdgRef = data_netcdf('EdgRefP',variable,stats_option)
	_,_,EdgEne = data_netcdf('EdgEneP',variable,stats_option)
	_,_,EdgTech = data_netcdf('EdgTechP',variable,stats_option)
	return lon,lat,Ref1970,EdgRef,EdgRef70A,EdgEne,EdgTech

		
def spa_pat_reg_mean(variable,MapOrReg,stats_option):
	if variable =='NOX_SRF':
		lon,lat,Ref1970_NO,EdgRef_NO,EdgRef70A_NO,EdgEne_NO,EdgTech_NO = data_readin('NO_SRF',stats_option);
		lon,lat,Ref1970_NO2,EdgRef_NO2,EdgRef70A_NO2,EdgEne_NO2,EdgTech_NO2 = data_readin('NO2_SRF',stats_option);
		Ref1970 = Ref1970_NO+Ref1970_NO2
		EdgRef = EdgRef_NO+EdgRef_NO2
		EdgRef70A = EdgRef70A_NO+EdgRef70A_NO2
		EdgEne = EdgEne_NO+EdgEne_NO2
		EdgTech = EdgTech_NO+EdgTech_NO2
		
	else:
		lon,lat,Ref1970,EdgRef,EdgRef70A,EdgEne,EdgTech = data_readin(variable,stats_option);	
	if MapOrReg == 'SpaPat':
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
			mask = mask_weight('GloLand',lon,lat,reverse=False);mask[mask>=0]=1
			dif = np.multiply(np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0),mask)
			sig = np.multiply(mannwhitneyu_test(vairable1, variable2),mask)
			# sig_dif = np.multiply(dif,sig)
			return dif,sig
		Tot,Tot_s =diff_sig(EdgRef,Ref1970)
		AAs,AAs_s = diff_sig(EdgRef,EdgRef70A)
		ENE,ENE_s= diff_sig(EdgRef,EdgEne)
		TECH,TECH_s= diff_sig(EdgRef,EdgTech)
		return lon,lat,Tot,AAs,ENE,TECH,Tot_s,AAs_s,ENE_s,TECH_s,np.nanmean(EdgRef,axis=0)
	
	elif MapOrReg == 'RegMean':
		def dif_std_for_region(var1,var2,mask):
			"""
			for a regional dif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
			"""
			if np.array_equal(var1,var2):
				dif =  np.nansum(np.nansum(np.multiply(mask,np.nanmean(var1,axis=0)),axis=1),axis=0)
				p25 = 0;p75=0;
			else:
				dif =  np.nansum(np.nansum(np.multiply(mask,np.nanmean(var1,axis=0) - np.nanmean(var2,axis=0)),axis=1),axis=0)
				var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
				var2_domain_mean = np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
				std = np.std(var1_domain_mean - var2_domain_mean);#print std
				p25 =  np.abs(dif - np.nanpercentile(var1_domain_mean - var2_domain_mean,25,axis=0))/1.25
				p75 =  np.abs(np.nanpercentile(var1_domain_mean - var2_domain_mean,75,axis=0) - dif)/1.25
				# print dif, p25,p75
			return dif,p25,p75

		region_dics ={'All':np.empty((4,3)),'GloLand':np.empty((4,3)),'ASIA':np.empty((4,3)),'EUROPE':np.empty((4,3)),'US':np.empty((4,3)),'ARCTIC':np.empty((4,3))} 
		for region_key in region_dics:
			# print region_key
			mask = mask_weight(region_key,lon,lat);
			plt.imshow(mask,origin='lower');plt.show()
			region_dics[region_key][0,:] = dif_std_for_region(EdgRef,Ref1970,mask)
			region_dics[region_key][3,:] = dif_std_for_region(EdgRef,EdgRef70A,mask)
			region_dics[region_key][1,:] = dif_std_for_region(EdgRef,EdgEne,mask)
			region_dics[region_key][2,:] = dif_std_for_region(EdgRef,EdgTech,mask)
		return region_dics


def global_mean(var):
	mask = mask_weight('All',lon,lat,reverse=False)
	spatial_mean = round(np.nansum(np.nansum(np.multiply(var,mask),axis=1),axis=0),2)
	return str(spatial_mean)

def produce_figures(MapOrReg):
	if MapOrReg =='SpaPat':
		####################################################################
		###   Spatial paattern of air pollution variables: CO_SRF,NOX_SRF,O3_SRF,pm25
		####################################################################
		for var in ['O3_SRF.MDA8','CO_SRF','NOX_SRF','PM2.5.mass0.20']:
			fig = plt.figure(facecolor='White',figsize=[9,13.5]);pad= 5;colormap='RdYlGn_r';
			if var =='CO_SRF':
				colorbar_min=-100;colorbar_max=100;stats_option = 'mean'
			elif var =='NOX_SRF':
				colorbar_min=-6;colorbar_max=6;stats_option = 'mean'
			elif var =='PM2.5.mass0.20':
				colorbar_min=-15;colorbar_max=15;stats_option = 'mean'
			elif var =='O3_SRF.MDA8':
				colorbar_min=-50;colorbar_max=50;stats_option = 'max'
			else:	
				print 'confirm the variable name'
				

			lon,lat,Tot,AAs,ENE,TECH,Tot_s,AAs_s,ENE_s,TECH_s,EdgRef = spa_pat_reg_mean(var,MapOrReg='SpaPat',stats_option=stats_option);		
			ax = plt.subplot(3,1,1);	
			ax.annotate('(a) Net total (2010 - 1970)',xy=(0.02,1.01), xytext=(0, pad),
							xycoords='axes fraction', textcoords='offset points',
							ha='left', va='baseline',rotation='horizontal',fontsize=20)	
			colormesh0 = spatial_figure(ax,Tot,Tot_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

			# ax = plt.subplot(2,2,2);
			# ax.annotate('(b) Best Estimations',xy=(0.02,1.01), xytext=(0, pad),
							# xycoords='axes fraction', textcoords='offset points',
							# ha='left', va='baseline',rotation='horizontal',fontsize=20)			
			# colormesh0 = spatial_figure(ax,AAs,AAs_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

			ax = plt.subplot(3,1,2);
			ax.annotate('(b) Energy Consumption',xy=(0.02,1.01), xytext=(0, pad),
							xycoords='axes fraction', textcoords='offset points',
							ha='left', va='baseline',rotation='horizontal',fontsize=20)					
			spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

			ax = plt.subplot(3,1,3);	
			ax.annotate('(c) Technology Advancements',xy=(0.02,1.01), xytext=(0, pad),
							xycoords='axes fraction', textcoords='offset points',
							ha='left', va='baseline',rotation='horizontal',fontsize=20)		
			colormesh1=spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

			cbar_ax = fig.add_axes([0.15, 0.03, 0.75, 0.015])
			char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2))

			plt.subplots_adjust(left=0.04, bottom=0.06, right=0.98, top=0.95, wspace=0.04, hspace=0.08); 
			plt.savefig(var+'_'+stats_option+'_spatil_maps.png', format='png', dpi=300)
			
			print var+'.finished'
	
	elif MapOrReg =='SpaPat_abs':
		####################################################################
		###   Spatial paattern of air pollution variables: CO_SRF,NOX_SRF,O3_SRF,pm25 in 2010
		####################################################################
		fig = plt.figure(facecolor='White',figsize=[18,11]);pad= 5;colormap='RdYlGn_r';
		
		ax = plt.subplot(2,2,1);	colorbar_min=0;colorbar_max=80;
		lon,lat,_,_,_,_,_,Tot_s,_,_,EdgRef = spa_pat_reg_mean('PM2.5.mass0.20',MapOrReg='SpaPat',stats_option='mean');
		mask = mask_weight('GloLand',lon,lat,reverse=False);mask[mask>=0]=1	
		EdgRef=np.multiply(EdgRef,mask);Tot_s[:] = np.nan
		ax.annotate(r'$(a) \mathrm{\/PM2.5\/(ug\/m^{-3})}$',xy=(0.02,1.01), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='baseline',rotation='horizontal',fontsize=20)	
		colormesh1 = spatial_figure(ax,EdgRef,Tot_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
		plt.colorbar(colormesh1,fraction=0.08,pad=0.1, orientation='horizontal',extend='both',ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2))
		# plt.show()
		
		ax = plt.subplot(2,2,2);colorbar_min=10;colorbar_max=80;
		lon,lat,_,_,_,_,_,_,_,_,EdgRef = spa_pat_reg_mean('O3_SRF.MDA8',MapOrReg='SpaPat',stats_option='mean');	
		EdgRef=np.multiply(EdgRef,mask)	
		ax.annotate(r'$(b) \mathrm{O3\/MDA8\/(ppb)}$',xy=(0.02,1.01), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='baseline',rotation='horizontal',fontsize=20)					
		colormesh1 = spatial_figure(ax,EdgRef,Tot_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
		plt.colorbar(colormesh1,fraction=0.08,pad=0.1,orientation='horizontal',extend='both',ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2))
		
		ax = plt.subplot(2,2,3);	colorbar_min=0;colorbar_max=15;
		lon,lat,_,_,_,_,_,_,_,_,EdgRef = spa_pat_reg_mean('NOX_SRF',MapOrReg='SpaPat',stats_option='mean');	
		EdgRef=np.multiply(EdgRef,mask)	
		ax.annotate(r'$(c) \mathrm{\/NO_x\/(ug\/m^{-3})}$',xy=(0.02,1.01), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='baseline',rotation='horizontal',fontsize=20)		
		colormesh1=spatial_figure(ax,EdgRef,Tot_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
		plt.colorbar(colormesh1,fraction=0.08,pad=0.1,orientation='horizontal',extend='both',ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2))
	
		ax = plt.subplot(2,2,4);colorbar_min=50;colorbar_max=350;		
		lon,lat,_,_,_,_,_,_,_,_,EdgRef = spa_pat_reg_mean('CO_SRF',MapOrReg='SpaPat',stats_option='mean');	
		EdgRef=np.multiply(EdgRef,mask)			
		ax.annotate(r'$(d) \mathrm{\/CO\/(ppb)}$',xy=(0.02,1.01), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='baseline',rotation='horizontal',fontsize=20)		
		colormesh1=spatial_figure(ax,EdgRef,Tot_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)		
		plt.colorbar(colormesh1,fraction=0.08,pad=0.1,orientation='horizontal',extend='both',ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2))
		
		plt.subplots_adjust(left=0.02, bottom=0.05, right=0.98, top=0.95, wspace=0.04, hspace=0.20); 
		plt.savefig('abs2010_spatil_maps.png', format='png', dpi=300)
		
	
	elif MapOrReg=='RegMean':

		####################################################################
		###   Regional Mean Bar of ERF, TS and Pr
		####################################################################
		def plot_bar_errbar(data):
			# print data
			# print "------"+var+"--------"
			def autolabel(X_pos,values,height_lift):
				## Attach a text label above each bar displaying its height
				height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
				for i in range(len(height)):
					ax.text(X_pos[i],y_pos[i],'%4.2f' % height[i], ha='center', va='bottom',size=10)
			
			ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
			ax.axhline(y=0,color='k',linewidth = 2);ax.axvline(x=1.5,color='k',linewidth = 1);
			ax.axvline(x=3.5,color='k',linewidth = 1);ax.axvline(x=5.5,color='k',linewidth = 1);
			ax.axvline(x=7.5,color='k',linewidth = 1);ax.axvline(x=9.5,color='k',linewidth = 1);
			
			X_pos=np.arange(0,12.0,2);dis=0.5
			
			y_value = np.array([data['All'][0,0],data['GloLand'][0,0],data['ARCTIC'][0,0],data['ASIA'][0,0],data['EUROPE'][0,0],data['US'][0,0]])
			yerr = np.stack((data['All'][0,1:3],data['GloLand'][0,1:3],data['ARCTIC'][0,1:3],data['ASIA'][0,1:3],data['EUROPE'][0,1:3],data['US'][0,1:3]),axis=1)	
			rects1=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='b', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
			autolabel(X_pos,y_value,1.1)

			y_value = np.array([data['All'][1,0],data['GloLand'][1,0],data['ARCTIC'][1,0],data['ASIA'][1,0],data['EUROPE'][1,0],data['US'][1,0]])
			yerr = np.stack((data['All'][1,1:3],data['GloLand'][1,1:3],data['ARCTIC'][1,1:3],data['ASIA'][0,1:3],data['EUROPE'][1,1:3],data['US'][1,1:3]),axis=1)	
			rects2=ax.bar(X_pos+dis*1, y_value, yerr=yerr, align='center',color='r', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
			autolabel(X_pos+dis*1,y_value,1.1)
			
			y_value = np.array([data['All'][2,0],data['GloLand'][2,0],data['ARCTIC'][2,0],data['ASIA'][2,0],data['EUROPE'][2,0],data['US'][2,0]])
			yerr = np.stack((data['All'][2,1:3],data['GloLand'][2,1:3],data['ARCTIC'][2,1:3],data['ASIA'][0,1:3],data['EUROPE'][2,1:3],data['US'][2,1:3]),axis=1)	
			rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr, align='center',color='g', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
			autolabel(X_pos+dis*2,y_value,1.1)

			return rects1,rects2,rects3

		fig = plt.figure(facecolor='White',figsize=[10,12]);pad= 5;
		ax = plt.subplot(4,1,1);
		var ='PM2.5.mass0.20';stats_option='mean'
		ERFT_dics = spa_pat_reg_mean(var,MapOrReg='RegMean',stats_option=stats_option);	
		ax.annotate( r'$\mathrm{\/PM2.5\/(ug\/m^{-3})}$',xy=(-0.08,0.5), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='center',rotation='vertical',fontsize=20)
		rects1,rects2,rects3 = plot_bar_errbar(ERFT_dics)
		ax.set_xlim([-0.5,11.5]);ax.set_xticks([]);
		legend1=ax.legend((rects1[0], rects2[0],rects3[0]),\
		('Net total (2010 - 1970)','Energy Consumption', 'Technology Advancements'), shadow=False,ncol=3,loc='lower left',fontsize = 13)
		legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)


		ax = plt.subplot(4,1,2);
		var ='O3_SRF.MDA8';stats_option='max'
		TS_dics =spa_pat_reg_mean(var,MapOrReg='RegMean',stats_option=stats_option);	
		ax.annotate( r'$\mathrm{\/O3\/MDA8\/(ppb)}$',xy=(-0.08,0.5), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='center',rotation='vertical',fontsize=20)
		rects1,rects2,rects3 = plot_bar_errbar(TS_dics)
		ax.set_xlim([-0.5,11.5]);ax.set_xticks([]);

		ax = plt.subplot(4,1,3);
		var ='NOX_SRF';stats_option='mean'
		Pr_dics =spa_pat_reg_mean(var,MapOrReg='RegMean',stats_option=stats_option);	
		ax.annotate( r'$\mathrm{\/NO_X\/(ug\/m^{-3})}$',xy=(-0.08,0.5), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='center',rotation='vertical',fontsize=20)
		rects1,rects2,rects3 = plot_bar_errbar(Pr_dics)
		ax.set_xlim([-0.5,11.5]);ax.set_xticks([]);

		ax = plt.subplot(4,1,4);
		var ='CO_SRF';stats_option='mean'
		Pr_dics = spa_pat_reg_mean(var,MapOrReg='RegMean',stats_option=stats_option);	
		ax.annotate( r'$\mathrm{\/CO\/(ug\/m^{-3})}$',xy=(-0.08,0.5), xytext=(0, pad),
						xycoords='axes fraction', textcoords='offset points',
						ha='left', va='center',rotation='vertical',fontsize=20)
		rects1,rects2,rects3 = plot_bar_errbar(Pr_dics)
		ax.set_xlim([-0.5,11.5]);ax.set_xticks(np.arange(0.5,11.5,2));ax.set_xticklabels(('Global','Land','Arctic','Asia','Europe','USA'),fontsize=20);

		plt.subplots_adjust(left=0.08, bottom=0.05, right=0.98, top=0.98, wspace=0.04, hspace=0.10); 
		plt.savefig('RegMean_PM2.5_O3_NOX_CO.png', format='png', dpi=300)

		
produce_figures(MapOrReg='SpaPat')  #SpaPat RegMean SpaPat_abs