"""
This is to calculate premature mortality based on Silva's formula. for PM2.5 and asurface ozone
Baseline moratlity rate is obtained from Global premature mortality due to anthropogenic outdoor air pollution and
the contribution of past climate change. 

population and the age 25+ fraction are from the NASA and UN DESA, respectively
"""
import scipy
import numpy as np
import netCDF4 as nc4
import scipy.io as sio
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp
from scipy.interpolate import interp2d  as interp2d
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from scipy.stats import mannwhitneyu as man_test

# from scipy.interpolate import interp2d
import os;import site
lib_path = os.path.join(
	os.path.realpath(
		os.path.dirname(__file__)
	), 
	os.path.pardir, 
	os.path.pardir,
)
site.addsitedir(lib_path)

from lib import *
def get_pollutant_concentrations(pollutant):
	def AreaWeight(lon1,lon2,lat1,lat2):
		'''
		calculate the earth radius in m2
		'''
		radius = 6371000;
		area = (np.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
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
		### This is because the boundary might not lo/lat might probably be ignored due to that the value does ntot necessarily match well
		if colum_e >=len(lon)-1:submap[:,:colum_s] =0;
		elif colum_e ==0:submap[:,colum_e:] =0
		else:submap[:,:colum_s] =0; submap[:,colum_e:] =0
		
		if row_s == 0:submap[row_e:,:] =0
		elif row_e >=len(lat)-1:submap[:row_s,:] =0; 
		else:submap[:row_s,:] =0; submap[row_e:,:] =0
	
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
		mask[np.isnan(mask)]=0;	mask[mask>0]=1;
		f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon, lat);
		mask[mask >= 1] = 1;mask[mask < 1] = 0;
		if reverse:
			mask=1-mask
		mask[mask==0] = np.nan
		mask=np.multiply(mask,area); 
		mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		# plt.imshow(mask_weighted,origin='lower');plt.show()
		return mask_weighted
		
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
					annual_stats[iyear-year_series[0],:,:] = np.nanmean(data_cache,axis=0)
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
			if variable in ['PM2.5.mass0.20','PM2.5_15']:
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

			
	def spa_pat(variable,stats_option):
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
	lon,lat,Tot,AAs,ENE,TECH,Tot_s,AAs_s,ENE_s,TECH_s,ref = spa_pat(pollutant,'mean')
	return lon,lat,Tot,AAs,ENE,TECH
	
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
	cmap = base = discrete_cmap(8, colormap)
	cmap.set_bad('w');cmap.set_over([129/255.0,0,0]); cmap.set_under([25/255.0,25/255.0,112/255.0]);  #darkgreen
	# cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_over('k'); cmap.set_under('green');  #darkgreen
	norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max) #
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True) #
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='/', alpha=0.,lw=0.9,latlon=True)
	axs.set_ylim([-62,90])
	return colormesh

def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in Km2
	'''
	radius = 6371;
	area = (np.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area
	
def box_fill(lon_s,lon_e,lat_s,lat_e,lon,lat,map,fill_value):
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
	submap = map;

	if colum_e >=len(lon)-1:submap[:,:colum_s] =0;
	elif colum_e ==0:submap[:,colum_e:] =0
	else:submap[:,:colum_s] =0; submap[:,colum_e:] =0
	
	if row_s == 0:submap[row_e:,:] =0
	elif row_e >=len(lat)-1:submap[:row_s,:] =0; 
	else:submap[:row_s,:] =0; submap[row_e:,:] =0
		
	submap = submap*fill_value
	# plt.imshow(submap);plt.show()
	return submap
	
def get_mask(region_key,lon_CESM,lat_CESM):
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask=ocean_mask[region_key][:];
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon_CESM,lat_CESM);
	mask[mask >= 1] = 1;mask[mask < 1] = 0;
	# plt.imshow(mask,origin='lower');plt.show()
	return lon_CESM,lat_CESM,mask	
	
def get_base_mortality(lon_CESM,lat_CESM):
	rergion_dic={'AUS':[110,155,-45,-11],'NA':[190,300,12,90],'SA':[278,326,-56,12],'GRL':[300,350,60,85],'RUS':[40,180,50,90],\
	'EUS1':[0,40,30,75],'EUS2':[350,360,30,75],\
	'MDE1':[40,75,30,50],'MDE2':[0,65,20,30],'MDE3':[340,360,20,30],\
	'AFC1':[0,52,-35,20],'AFC2':[340,360,-35,20],\
	'EA1':[100,145,20,50],'EA2':[75,100,30,50],\
	'SEA':[100,155,-11,20],'IND':[65,100,-5,30]}

	Mortality_y0 ={'AUS':[4.000,0.525,0.797],'NA':[4.673,0.654,1.127],'SA':[4.143,0.235,1.091],'GRL':[6.239,0.736,1.064],'RUS':[12.252,0.489,0.785],\
	'EUS1':[6.239,0.736,1.064],'EUS2':[6.239,0.736,1.064],\
	'MDE1':[5.062,0.143,0.621],'MDE2':[5.062,0.143,0.621],'MDE3':[5.062,0.143,0.621],\
	'AFC1':[6.338,0.083,1.346],'AFC2':[6.338,0.083,1.346],\
	'EA1':[5.971,0.620,1.930],'EA2':[5.971,0.620,1.930],\
	'SEA':[5.222,0.334,1.347],'IND':[6.768,0.147,2.006]}
	
	##OCEAN_MASKS FOR COUNTRIES
	lon_CESM,lat_CESM,GloLand = get_mask('GloLand',lon_CESM,lat_CESM)

	Mortality_y0_map ={'CardioPul':np.empty((192,288)),'Lung':np.empty((192,288)),'Respiratory':np.empty((192,288))}
	loop = -1  # get access to the dic Mortality_y0
	for diseas in Mortality_y0_map:
		loop = loop+1
		lon_CESM,lat_CESM,Global = get_mask('All',lon_CESM,lat_CESM)
		Global = Global*0;
		for region_key in rergion_dic:
			lon_CESM,lat_CESM,GloBase = get_mask('All',lon_CESM,lat_CESM)
			box = rergion_dic[region_key]
			mask = box_fill(box[0],box[1],box[2],box[3],lon_CESM,lat_CESM,GloBase,Mortality_y0[region_key][loop]/1000)
			Global = Global+mask
		Mortality_y0_map[diseas]=np.multiply(Global,GloLand)
		# plt.imshow(Mortality_y0_map[diseas],origin='lower');plt.show()
	return 	Mortality_y0_map
	
def get_CESM_gc_area():
	path_CESM = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_PS_200602_210101.nc'
	CESM_fid = nc4.Dataset(path_CESM,mode='r')
	lat_CESM = CESM_fid.variables['lat'][:]
	lon_CESM = CESM_fid.variables['lon'][:]
	lons,lats = np.meshgrid (lon_CESM,lat_CESM)
	lat_res = lat_CESM[1] - lat_CESM[0];lon_res = lon_CESM[1] - lon_CESM[0];
	CESM_gc_area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	return 	lon_CESM,lat_CESM,CESM_gc_area	
	
def get_exposed_population(lon_CESM,lat_CESM,CESM_gc_area):	
	def regrid_pop2cesm(lon,data):
		def regulate_grids(data):
			data_rregulargrid = np.empty((720,1440));data_rregulargrid[:] = 0;
			data_rregulargrid[20:600,:]=data;
			## shift -180 - 180 to 0-360, 90 - -90 to -90 - 90
			shift = np.zeros([720,1440]); shift[:,0:720] = data_rregulargrid[:,720:1440];shift[:,720:1440] = data_rregulargrid[:,0:720]
			return np.flipud(shift)
		
		pop_regulated = regulate_grids(data)  #  regulate the original data to a regular global  range grids
		x = lon;y=np.arange(-89.875,90,0.25);
		xout,yout=np.meshgrid(lon,np.arange(-89.875,90,0.0625))
		original_interp = interp(pop_regulated, x, y, xout, yout, checkbounds=False, masked=False, order=1)   # interpolation into a regular gids to avoid boundary issues
		
		# area weighted population into persons per grids
		lons,lats = np.meshgrid (lon,np.arange(-89.875,90,0.0625))
		lat_res = 0.0625;lon_res = 0.25;
		area = getArea(lons,lons+lon_res,lats,lats+lat_res)
		Original_griddata = np.multiply(original_interp,area)
		
		## remap into CESM grids
		cesm_griddata = np.zeros((192,288));
		x_cesm=range(0,288);y_cesm=range(0,192)
		for row in range(0,192):
			for column in range(0,288):
				cesm_griddata[row,column] =np.nansum(np.nansum(Original_griddata[row*15:(row+1)*15,column*5:(column+1)*5],axis=1),axis=0)
		cesm_griddata[cesm_griddata<=0]=np.nan
		cesm_griddata = np.divide(cesm_griddata,CESM_gc_area)*1000   # from km2 to 1000 km2
		return cesm_griddata   # Population in each gridcell
		
	def get_pop_25plus_fraction(lon_CESM,lat_CESM):
		rergion_dic={'AUS':[110,155,-45,-11],'SA':[278,326,-56,10],'NA':[190,300,30,90],'CAM':[244,277,10,30],\
		'GRL':[300,350,60,85],'RUS':[40,180,50,90],\
		'EUS1':[0,40,48,75],'EUS2':[350,360,48,75],\
		'SEU1':[0,40,20,48],'SEU2':[350,360,20,48],\
		'MDE1':[40,75,30,50],'MDE2':[0,65,20,30],'MDE3':[340,360,20,30],\
		'AFC1':[0,52,-35,20],'AFC2':[340,360,-35,20],\
		'EA1':[100,145,20,50],'EA2':[75,100,30,50],\
		'SEA':[100,155,-11,20],'IND':[65,100,-5,30]}
        
		## The 25+ fraction is from the (World Population Prospects, 2010 Revision, http://esa.un.org/unpd/wpp/Excel-Data/population.htm
		frac_25plus ={'AUS':66.6,'NA':66.1,'SA':55.7,'CAM':49.8,'GRL':69.5,'RUS':70.9,\
		'EUS1':69.5,'EUS2':69.5,\
		'SEU1':74.4,'SEU2':74.4,\
		'MDE1':74.4,'MDE2':74.4,'MDE3':74.4,\
		'AFC1':48.8,'AFC2':48.8,\
		'EA1':66.2,'EA2':50.2,\
		'SEA':53.7,'IND':49.1}
		
		lon_CESM,lat_CESM,Global = get_mask('All',lon_CESM,lat_CESM)
		Global = Global*0
		for region_key in rergion_dic:
			lon_CESM,lat_CESM,GloBase = get_mask('All',lon_CESM,lat_CESM)
			box = rergion_dic[region_key]
			# print region_key
			mask = box_fill(box[0],box[1],box[2],box[3],lon_CESM,lat_CESM,GloBase,frac_25plus[region_key])
			Global = Global+mask
		lon_CESM,lat_CESM,GloLand = get_mask('GloLand',lon_CESM,lat_CESM)
		frac_25plus_map=np.multiply(Global,GloLand)/100
		return frac_25plus_map
	
	data_path = '/exports/csce/datastore/geos/users/s1667168/obs/gpw_v4_e_atotpopbt_dens_15_min.nc'
	nc_fid = nc4.Dataset(data_path,mode='r')
	lon = nc_fid.variables['longitude'][:]
	pop = nc_fid.variables['Population Density, v4.10 (2000, 2005, 2010, 2015, 2020): 15 arc-minutes'][2,:,:];  # year 2010
	pop[pop<0]=0
	cesm_gridpop = regrid_pop2cesm(lon,pop)
	fraction = get_pop_25plus_fraction(lon_CESM,lat_CESM);
	exposed_pop= np.multiply(cesm_gridpop,fraction);
	# plt.imshow(fraction,origin='lower');plt.show()
	# plt.imshow(exposed_pop,origin='lower');plt.show()
	return exposed_pop

		
def get_mortality(pollutant,stats):
	def health_impact_unction(disease,DelPollutant):
		Beta = np.log(RR_dic[disease])/10
		PopTimesy0 = np.multiply(Mortality_y0[disease][:],exposed_pop)
		Attributable_fraction = 1-np.exp(-Beta*DelPollutant)
		mortality = np.multiply(Attributable_fraction,PopTimesy0)
		mortality[mortality==0]=np.nan
		print 'Mortality_y0', Mortality_y0[disease][:,0]
		print 'exposed_pop', exposed_pop[:,0]
		print 'mortality', mortality[:,0]
		return mortality
	lon_CESM,lat_CESM,CESM_gc_area = get_CESM_gc_area()	
	Mortality_y0 =get_base_mortality(lon_CESM,lat_CESM)
	exposed_pop = get_exposed_population(lon_CESM,lat_CESM,CESM_gc_area)
	lon,lat,Tot,AAs,ENE,TECH = get_pollutant_concentrations(pollutant)
	if stats =='mean':
		RR_dic ={'CardioPul':1.128,'Lung':1.142,'Respiratory':1.040}
	elif stats =='min':
		RR_dic ={'CardioPul':1.077,'Lung':1.057,'Respiratory':1.013}
	elif stats =='max':
		RR_dic ={'CardioPul':1.182,'Lung':1.234,'Respiratory':1.067}
	
	if pollutant =='PM2.5.mass0.20':
		M_Tot = health_impact_unction('CardioPul',Tot)+health_impact_unction('Lung',Tot)
		M_AAs = health_impact_unction('CardioPul',AAs)+health_impact_unction('Lung',AAs)
		M_ENE = health_impact_unction('CardioPul',ENE)+health_impact_unction('Lung',ENE)
		M_TECH = health_impact_unction('CardioPul',TECH)+health_impact_unction('Lung',TECH)
	elif pollutant =='O3_SRF.MDA8':
		M_Tot = health_impact_unction('Respiratory',Tot)
		M_AAs = health_impact_unction('Respiratory',AAs)
		M_ENE = health_impact_unction('Respiratory',ENE)
		M_TECH = health_impact_unction('Respiratory',TECH)
	return lon,lat,M_Tot,M_AAs,M_ENE,M_TECH

##########
## MAIN  #
##########	

colors = [(51/255.0,161/255.0,201/255.0),(154/255.0,192/255.0,205/255.0),
	 (201/255.0,201/255.0,201/255.0), (255/255.0,236/255.0,139/255.0),
	 (255/255.0,255/255.0,0), (255/255.0,128/255.0,0),
	 (205/255.0,97/255.0,3/255.0),(255/255.0,0,0),
	 ]
from matplotlib.colors import LinearSegmentedColormap
colormap = LinearSegmentedColormap.from_list('my_list', colors, N=8)
fig = plt.figure(facecolor='White',figsize=[18,13.5]);pad= 5;
		
var = 'O3_SRF.MDA8';colorbar_min=-3;colorbar_max=9;stats='mean'
lon,lat,M_Tot,M_AAs,M_ENE,M_TECH = get_mortality(var,stats);

	
ax = plt.subplot(3,2,1);	
ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)
ax.annotate('Ozone',xy=(0.5,1.10), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='center',rotation='horizontal',fontsize=20)	
ax.annotate('Net total (2010 - 1970)',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)	
colormesh0 = spatial_figure(ax,M_Tot,M_Tot,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

# ax = plt.subplot(2,2,2);
# ax.annotate('(b) Best Estimations',xy=(0.02,1.01), xytext=(0, pad),
				# xycoords='axes fraction', textcoords='offset points',
				# ha='left', va='baseline',rotation='horizontal',fontsize=20)			
# colormesh0 = spatial_figure(ax,AAs,AAs_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,3);
ax.annotate('(b) ',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Energy Consumption',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)					
spatial_figure(ax,M_ENE,M_ENE,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,5);	
ax.annotate('(c) ',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)
ax.annotate('Technology Advancements',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)					
colormesh1=spatial_figure(ax,M_TECH,M_TECH,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.08, 0.03, 0.43, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)) #

var = 'PM2.5.mass0.20';colorbar_min=-100;colorbar_max=300;
lon,lat,M_Tot,M_AAs,M_ENE,M_TECH = get_mortality(var,stats);		
ax = plt.subplot(3,2,2);	
ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('PM2.5',xy=(0.5,1.10), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='center',rotation='horizontal',fontsize=20)		
colormesh0 = spatial_figure(ax,M_Tot,M_Tot,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot(3,2,4);
ax.annotate('(e) ',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)				
spatial_figure(ax,M_ENE,M_ENE,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
ax = plt.subplot(3,2,6);	
ax.annotate('(f) ',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)				
colormesh1=spatial_figure(ax,M_TECH,M_TECH,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.55, 0.03, 0.43, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)) #

plt.subplots_adjust(left=0.08, bottom=0.06, right=0.98, top=0.95, wspace=0.06, hspace=0.08); 
plt.savefig('pm_ozone_mortality_'+stats+'.png', format='png', dpi=300)


