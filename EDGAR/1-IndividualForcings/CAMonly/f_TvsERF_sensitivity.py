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
	return mask	

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
		
def data_readin(variable):	
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
		if scenario.startswith('F_'):
			year_series = range(2005,2029)
		elif scenario=='T1970RCP':
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

	def data_netcdf(exp,scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'+exp+'/'
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
	
	if variable == 'ERFT':  # ERF at TOA
		exp='CamOnly'
		lon,lat,FSNT = data_netcdf(exp,'F_Edg70GO','FSNT');_,_,FLNT = data_netcdf(exp,'F_Edg70GO','FLNT');Edg70GO=(FSNT-FLNT)
		_,_,FSNT = data_netcdf(exp,'F_1970','FSNT');_,_,FLNT = data_netcdf(exp,'F_1970','FLNT');Ref1970=(FSNT-FLNT)
		_,_,FSNT = data_netcdf(exp,'F_EdgRef','FSNT');_,_,FLNT = data_netcdf(exp,'F_EdgRef','FLNT');EdgRef=(FSNT-FLNT)
		_,_,FSNT = data_netcdf(exp,'F_Edg70Oz','FSNT');_,_,FLNT = data_netcdf(exp,'F_Edg70Oz','FLNT');Edg70Oz=(FSNT-FLNT)
		_,_,FSNT = data_netcdf(exp,'F_Edg70T10SOZ','FSNT');_,_,FLNT = data_netcdf(exp,'F_Edg70T10SOZ','FLNT');Edg70T10SOZ=(FSNT-FLNT)
	elif variable == 'TS': 
		exp='FullCp'   #  CamOnly
		lon,lat,Edg70GO = data_netcdf(exp,'Edg70GO',variable)
		_,_,Ref1970 = data_netcdf(exp,'T1970RCP',variable)
		_,_,EdgRef = data_netcdf(exp,'EdgRef',variable)
		_,_,Edg70Oz = data_netcdf(exp,'Edg70Oz',variable)
		_,_,Edg70T10SOZ = data_netcdf(exp,'Edg70T10SOZ',variable)
	return lon,lat,Ref1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ
	
def spa_pat_reg_mean(variable,option):
	"""
	option: SpaPat or RegMean
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""

	def mannwhitneyu_test(vairable1,variable2,p_threshold):
		# p_threshold=0.10
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
		sig = mannwhitneyu_test(vairable1, variable2,p_threshold=0.1)
		# sig_dif = np.multiply(dif,sig)
		return dif,sig
	
	def GradientTvsERF(y1,y2,x1,x2,mask,GloMask):
		####calculate the difference and their significance
		def mask_sig(vairable1,variable2,mask,p_threshold):
			dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
			sig = mannwhitneyu_test(vairable1, variable2,p_threshold)
			mask_sig = np.multiply(sig,mask)
			mask_sig_weighted = mask_sig/np.nansum(np.nansum(mask_sig,axis=1),axis=0)
			dif =  np.multiply(dif,mask_sig_weighted)
			return dif
		"""
		for a rgional sif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		global_dx = np.nansum(np.nansum(mask_sig(x1,x2,GloMask,p_threshold = 0.01),axis=1),axis=0)
		dy = np.nansum(np.nansum(mask_sig(y1,y2,mask,p_threshold = 0.2),axis=1),axis=0)
		dx = np.nansum(np.nansum(mask_sig(x1,x2,mask,p_threshold = 0.10),axis=1),axis=0)
		GloGra = np.divide(dy,global_dx)  # gradient of T versus global mean of ERF
		RegGra = np.divide(dy,dx)  # gradient of T versus regional mean of ERF
		return GloGra,RegGra
		
	if option == 'RegMean':	
		region_dics ={'All':np.empty((5,2)),'TROPICS':np.empty((5,2)),'ASIA':np.empty((5,2)),'EUROPE':np.empty((5,2)),'US':np.empty((5,2)),'ARCTIC':np.empty((5,2))} 
		lon,lat,Ref1970_TS,Edg70GO_TS,Edg70Oz_TS,EdgRef_TS,Edg70T10SOZ_TS = data_readin('TS');
		lon,lat,Ref1970_ERFT,Edg70GO_ERFT,Edg70Oz_ERFT,EdgRef_ERFT,Edg70T10SOZ_ERFT = data_readin('ERFT');
		
		GloMask = mask_weight('All',lon,lat);
		for region_key in region_dics:
			mask = mask_weight(region_key,lon,lat); 
			region_dics[region_key][0,:] = GradientTvsERF(EdgRef_TS,Ref1970_TS,EdgRef_ERFT,Ref1970_ERFT,mask,GloMask)
			region_dics[region_key][1,:] = GradientTvsERF(Edg70Oz_TS,Edg70GO_TS,Edg70Oz_ERFT,Edg70GO_ERFT,mask,GloMask)
			region_dics[region_key][2,:] = GradientTvsERF(Edg70GO_TS,Ref1970_TS,Edg70GO_ERFT,Ref1970_ERFT,mask,GloMask)
			region_dics[region_key][3,:] = GradientTvsERF(EdgRef_TS, Edg70T10SOZ_TS,EdgRef_ERFT, Edg70T10SOZ_ERFT,mask,GloMask)
			region_dics[region_key][4,:] = GradientTvsERF(Edg70T10SOZ_TS, Edg70Oz_TS,Edg70T10SOZ_ERFT, Edg70Oz_ERFT,mask,GloMask)
		return region_dics
	elif option == 'SpaPat':
		lon,lat,Ref970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ = data_readin(variable);
		Total,Total_s = diff_sig(EdgRef,Ref970)
		GHG,GHG_s = diff_sig(Edg70Oz,Edg70GO)
		AEROSOL,AEROSOL_s = diff_sig(Edg70GO,Ref970)
		TrO3,TrO3_s = diff_sig(EdgRef, Edg70T10SOZ)
		StO3,StO3_s = diff_sig(Edg70T10SOZ, Edg70Oz)
		return lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s


def bar_plots(data,NoneGlobe=True):
	def autolabel(X_pos,values,height_lift):
		## Attach a text label above each bar displaying its height
		height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
		ax.text(X_pos,y_pos,'%4.2f' % height, ha='center', va='bottom',size=10)
	
	def bar_block(X_pos_Star,y_value):
		color_list=['k','g','r','m','b']
		X_pos=X_pos_Star;dis=0.5
		for item in range(5):
			ax.bar(X_pos+dis*item, y_value[item], align='center',color=color_list[item], ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
			# autolabel(X_pos+dis*item,y_value[item],1.05)
	ax.set_ylim([-0.5,1.5]);ax.set_yticks(np.arange(-0.5,1.51,0.5))
	ax.axhline(y=0,color='k',linewidth = 2);ax.set_xticklabels([])
	ax.axhline(y=0.36,ls='--',color='k',linewidth = 2);
	ax.axhline(y=0.43,ls='--',color='b',linewidth = 2);
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
	### versus global mean ERF
	y_value = data[:,0];
	X_pos_Star=0.5;bar_block(X_pos_Star,y_value)
	### versus regional mean ERF
	if NoneGlobe:
		ax.set_xlim([0,6]);
		ax.axvspan(3, 6, alpha=0.2, color='y',edgecolor='none');
		y_value = data[:,1];
		X_pos_Star=3.5;bar_block(X_pos_Star,y_value)

"""
####################################################################
###   Spatial paattern of ERF, TS response & scaled diffeences   ###
####################################################################
index ='TS';lon,lat,Total_TS,GHG_TS,AEROSOL_TS,TrO3_TS,StO3_TS,Total_s_TS,GHG_s_TS,AEROSOL_s_TS,TrO3_s_TS,StO3_s_TS = spa_pat_reg_mean(index,option='SpaPat');		
index ='ERFT';lon,lat,Total_EF,GHG_EF,AEROSOL_EF,TrO3_EF,StO3_EF,Total_s_EF,GHG_s_EF,AEROSOL_s_EF,TrO3_s_EF,StO3_s_EF = spa_pat_reg_mean(index,option='SpaPat');				

def global_mean(var):
	mask = mask_weight('All',lon,lat,reverse=False)
	mask = mask/np.nansum(np.nansum(mask,axis=1),axis=0)
	spatial_mean = round(np.nansum(np.nansum(np.multiply(var,mask),axis=1),axis=0),2)
	return str(spatial_mean)

fig = plt.figure(facecolor='White',figsize=[14,11]);pad= 5;colormap='RdBu_r';

####################ERF###########################
colorbar_min=-4;colorbar_max=4;
ax = plt.subplot(4,3,1);
ax.annotate('GHGs',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(a)                                                 '+global_mean(GHG_EF),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('ERF',xy=(0.5,1.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	

				
colormesh0 = spatial_figure(ax,GHG_EF,GHG_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,4);
ax.annotate('AAs',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(b)                                             '+global_mean(AEROSOL_EF),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,AEROSOL_EF,AEROSOL_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,7);
ax.annotate('Trop. O3',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(c)                                              '+global_mean(TrO3_EF),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TrO3_EF,TrO3_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,10);
ax.annotate('Strat. O3',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(d)                                             '+global_mean(StO3_EF),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,StO3_EF,StO3_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.06, 0.06, 0.27, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/W\/m^{-2}}$',xy=(0.5,-3.2), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)


####################TS###########################
colorbar_min=-1.5;colorbar_max=1.5;

ax = plt.subplot(4,3,2);
ax.annotate('(e)                                              '+global_mean(GHG_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('T. response',xy=(0.5,1.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	
				
colormesh0 = spatial_figure(ax,GHG_TS,GHG_s_TS,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,5)
ax.annotate('(f)                                                '+global_mean(AEROSOL_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,AEROSOL_TS,AEROSOL_s_TS,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,8);
ax.annotate('(g)                                               '+global_mean(TrO3_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TrO3_TS,TrO3_s_TS,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,11);
ax.annotate('(h)                                              '+global_mean(StO3_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,StO3_TS,StO3_s_TS,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.38, 0.06, 0.27, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/K}$',xy=(0.5,-3.2), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)


####################Scaled diff###########################
colorbar_min=-2;colorbar_max=2;

def scale_T(ERF,dT):
	def global_mean_Num(var):
		mask = mask_weight('All',lon,lat,reverse=False)
		mask = mask/np.nansum(np.nansum(mask,axis=1),axis=0)
		spatial_mean = np.nansum(np.nansum(np.multiply(var,mask),axis=1),axis=0)
		return spatial_mean
	GHG_EF_mean = global_mean_Num(GHG_EF);GHG_TS_mean = global_mean_Num(GHG_TS)
	EF_MEAN = global_mean_Num(ERF);
	# return dT -(EF_MEAN/GHG_EF_mean)*GHG_TS_mean
	# return dT - np.multiply(GHG_TS,np.divide(ERF,GHG_EF))
	return dT - np.multiply(GHG_TS_mean,np.divide(ERF,GHG_EF_mean))
	# return -GHG_TS_mean+np.multiply(GHG_EF_mean,np.divide(ERF,dT))


AEROSOL_TS= scale_T(AEROSOL_EF,AEROSOL_TS);
TrO3_TS= scale_T(TrO3_EF,TrO3_TS);
StO3_TS= scale_T(StO3_EF,StO3_TS);

sig = AEROSOL_s_EF; sig[:]=np.nan

ax = plt.subplot(4,3,3);
ax.annotate('Scaled T difference',xy=(0.5,1.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	
ax.axis('off')
# ax.annotate('(i)                                                   '+global_mean(GHG_TS),xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=15)			
# spatial_figure(ax,AEROSOL_TS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,6);
ax.annotate('(i)                                                '+global_mean(AEROSOL_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,AEROSOL_TS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,9);
ax.annotate('(j)                                               '+global_mean(TrO3_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TrO3_TS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,12);
ax.annotate('(k)                                               '+global_mean(StO3_TS),xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
colormesh1=spatial_figure(ax,StO3_TS,sig,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


cbar_ax = fig.add_axes([0.70, 0.06, 0.27, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/K}$',xy=(0.5,-3.2), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)

plt.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.92, wspace=0.1, hspace=0.20); 
plt.savefig('ERF_TS_SCALED_DIF.png', format='png', dpi=1000)

		
"""		
############################
###   T versus ERF bar   ###
############################		
		
		
index ='TS';region_dics = spa_pat_reg_mean(index,option='RegMean');	
print region_dics['All'][:]
fig = plt.figure(facecolor='White',figsize=[10,6]);pad= 5


ax = plt.subplot(2,3,1);
ax.annotate('(a) Globe',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
bar_plots(region_dics['All'],NoneGlobe=False)

# ax.set_ylim([-2,2]);


ax = plt.subplot(2,3,2);
ax.annotate('(b) Tropics (28S - 28N)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
bar_plots(region_dics['TROPICS'])
# ax.set_ylim([-2,2]);

ax = plt.subplot(2,3,3);
ax.annotate('(c) Arctic (70N - 90N)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
bar_plots(region_dics['ARCTIC'])
# ax.set_ylim([-2,2]);

ax = plt.subplot(2,3,4);
ax.annotate('(d) Asia',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
bar_plots(region_dics['ASIA'])
# ax.set_ylim([-1,7]);

ax = plt.subplot(2,3,5);
ax.annotate('(e) Europe',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
bar_plots(region_dics['EUROPE'])
# ax.set_ylim([-1,7]);

ax = plt.subplot(2,3,6);
ax.annotate('(f) USA',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
bar_plots(region_dics['US'])
# ax.set_ylim([-1,7]);

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


plt.subplots_adjust(left=0.05, bottom=0.12, right=0.95, top=0.95, wspace=0.25, hspace=0.25); 
plt.savefig('TvsERFsensitivity.png', format='png', dpi=1000)	

