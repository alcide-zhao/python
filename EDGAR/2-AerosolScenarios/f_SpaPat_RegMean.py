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
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_over('k'); cmap.set_under('darkblue');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh
		
		
def data_readin(variable,exp):	
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
		elif variable == 'TGCLDLWP':data=data*10**(3)
		var40map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var40map
		
	if exp=='CamOnly':
		if variable == 'ERFT':
			lon,lat,var1 = data_netcdf(exp,'F_Edg70GO','FSNT');lon,lat,var2 = data_netcdf(exp,'F_Edg70GO','FLNT');Edg70GO =var1-var2
			_,_,var1 = data_netcdf(exp,'F_1970','FSNT');_,_,var2 = data_netcdf(exp,'F_1970','FLNT');Ref1970 =var1-var2
			_,_,var1 = data_netcdf(exp,'F_EdgRef','FSNT');_,_,var2 = data_netcdf(exp,'F_EdgRef','FLNT');EdgRef =var1-var2
			_,_,var1 = data_netcdf(exp,'F_EdgEne','FSNT');_,_,var2 = data_netcdf(exp,'F_EdgEne','FLNT');EdgEne =var1-var2
			_,_,var1 = data_netcdf(exp,'F_EdgTech','FSNT');	_,_,var2 = data_netcdf(exp,'F_EdgTech','FLNT');EdgTech =var1-var2
		else:
			lon,lat,Edg70GO = data_netcdf(exp,'F_Edg70GO',variable)
			_,_,Ref1970 = data_netcdf(exp,'F_1970',variable)
			_,_,EdgRef = data_netcdf(exp,'F_EdgRef',variable)
			_,_,EdgEne = data_netcdf(exp,'F_EdgEne',variable)
			_,_,EdgTech = data_netcdf(exp,'F_EdgTech',variable)

	elif exp=='FullCp':
		if variable == 'precip':
			lon,lat,var1 = data_netcdf(exp,'Edg70GO','PRECC');lon,lat,var2 = data_netcdf(exp,'Edg70GO','PRECL');Edg70GO =(var1+var2)*24*60*60*1000
			_,_,var1 = data_netcdf(exp,'T1970RCP','PRECC');_,_,var2 = data_netcdf(exp,'T1970RCP','PRECL');Ref1970 =(var1+var2)*24*60*60*1000
			_,_,var1 = data_netcdf(exp,'EdgRef','PRECC');_,_,var2 = data_netcdf(exp,'EdgRef','PRECL');EdgRef =(var1+var2)*24*60*60*1000
			_,_,var1 = data_netcdf(exp,'EdgEne','PRECC');_,_,var2 = data_netcdf(exp,'EdgEne','PRECL');EdgEne =(var1+var2)*24*60*60*1000
			_,_,var1 = data_netcdf(exp,'EdgTech','PRECC');	_,_,var2 = data_netcdf(exp,'EdgTech','PRECL');EdgTech =(var1+var2)*24*60*60*1000	
		else:
			lon,lat,Edg70GO = data_netcdf(exp,'Edg70GO',variable)
			_,_,Ref1970 = data_netcdf(exp,'T1970RCP',variable)
			_,_,EdgRef = data_netcdf(exp,'EdgRef',variable)
			_,_,EdgEne = data_netcdf(exp,'EdgEne',variable)
			_,_,EdgTech = data_netcdf(exp,'EdgTech',variable)
	return lon,lat,Ref1970,Edg70GO,EdgRef,EdgEne,EdgTech

		
def spa_pat_reg_mean(variable,exp,option):
	"""
	option: SpaPat or RegMean
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""
	
	lon,lat,Ref1970,Edg70GO,EdgRef,EdgEne,EdgTech = data_readin(variable,exp);	
	if option == 'SpaPat':
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
		AEROSOL,AEROSOL_s = diff_sig(Edg70GO,Ref1970)
		ENE,ENE_s= diff_sig(EdgRef,EdgEne)
		TECH,TECH_s= diff_sig(EdgRef,EdgTech)
		return lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s
	elif option == 'RegMean':
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

		region_dics ={'All':np.empty((4,3)),'TROPICS':np.empty((4,3)),'ASIA':np.empty((4,3)),'EUROPE':np.empty((4,3)),'US':np.empty((4,3)),'ARCTIC':np.empty((4,3))} 
		for region_key in region_dics:
			# print region_key
			mask = mask_weight(region_key,lon,lat); 
			region_dics[region_key][0,:] = dif_std_for_region(Edg70GO,Ref1970,mask)
			region_dics[region_key][1,:] = dif_std_for_region(EdgRef,EdgEne,mask)
			region_dics[region_key][2,:] = dif_std_for_region(EdgRef,EdgTech,mask)
			region_dics[region_key][3,:] = dif_std_for_region(Ref1970,Ref1970,mask)
		return region_dics
		

####################################################################
###   Spatial paattern of AOD, ERF, and TS
####################################################################


def global_mean(var):
	mask = mask_weight('All',lon,lat,reverse=False)
	spatial_mean = round(np.nansum(np.nansum(np.multiply(var,mask),axis=1),axis=0),2)
	return str(spatial_mean)

fig = plt.figure(facecolor='White',figsize=[9,8]);pad= 5;colormap='RdBu_r';

################### AOD ###########################
index ='AODVIS';lon,lat,AEROSOL_AOD,ENE_AOD,TECH_AOD,AEROSOL_s_AOD,ENE_s_AOD,TECH_s_AOD = spa_pat_reg_mean(index,exp='FullCp',option='SpaPat');		
colormap='RdYlBu_r';colorbar_min=-0.05;colorbar_max=0.05;

ax = plt.subplot(3,2,1);
ax.annotate('Best\nEstimation',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('AOD (#)',xy=(0.5,1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)	

				
colormesh0 = spatial_figure(ax,AEROSOL_AOD,AEROSOL_s_AOD,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,3);
ax.annotate('Energy\nConsumption',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(b)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,ENE_AOD,ENE_s_AOD,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,5);
ax.annotate('Technology\nAdvancements',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TECH_AOD,TECH_s_AOD,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.08, 0.04, 0.44, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))

				
				
################### ERF ###########################
index ='ERFT';lon,lat,AEROSOL_EF,ENE_EF,TECH_EF,AEROSOL_s_EF,ENE_s_EF,TECH_s_EF = spa_pat_reg_mean(index,exp='CamOnly',option='SpaPat');	

colorbar_min=-4;colorbar_max=4;colormap ='RdBu_r'

ax = plt.subplot(3,2,2);
ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('ERF ' +r'($\mathrm{W\/m^{-2}}$)',xy=(0.5,1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)	
				
colormesh0 = spatial_figure(ax,AEROSOL_EF,AEROSOL_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,4)
ax.annotate('(e)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)					
spatial_figure(ax,ENE_EF,ENE_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,6);
ax.annotate('(f)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)		
colormesh1=spatial_figure(ax,TECH_EF,TECH_s_EF,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.54, 0.04, 0.44, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))

		
plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.92, wspace=0.04, hspace=0.20); 
plt.savefig('SpaPat_AOD_ERF.png', format='png', dpi=1000)



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
	
	y_value = np.array([data['All'][0,0],data['TROPICS'][0,0],data['ARCTIC'][0,0],data['ASIA'][0,0],data['EUROPE'][0,0],data['US'][0,0]])
	yerr = np.stack((data['All'][0,1:3],data['TROPICS'][0,1:3],data['ARCTIC'][0,1:3],data['ASIA'][0,1:3],data['EUROPE'][0,1:3],data['US'][0,1:3]),axis=1)	
	rects1=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='b', ecolor='b',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos,y_value,1.3)

	y_value = np.array([data['All'][1,0],data['TROPICS'][1,0],data['ARCTIC'][1,0],data['ASIA'][1,0],data['EUROPE'][1,0],data['US'][1,0]])
	yerr = np.stack((data['All'][1,1:3],data['TROPICS'][1,1:3],data['ARCTIC'][1,1:3],data['ASIA'][0,1:3],data['EUROPE'][1,1:3],data['US'][1,1:3]),axis=1)	
	rects2=ax.bar(X_pos+dis*1, y_value, yerr=yerr, align='center',color='r', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*1,y_value,1.3)
	
	y_value = np.array([data['All'][2,0],data['TROPICS'][2,0],data['ARCTIC'][2,0],data['ASIA'][2,0],data['EUROPE'][2,0],data['US'][2,0]])
	yerr = np.stack((data['All'][2,1:3],data['TROPICS'][2,1:3],data['ARCTIC'][2,1:3],data['ASIA'][0,1:3],data['EUROPE'][2,1:3],data['US'][2,1:3]),axis=1)	
	rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr, align='center',color='g', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*2,y_value,1.3)

	return rects1,rects2,rects3

fig = plt.figure(facecolor='White',figsize=[10,7]);pad= 5;
ax = plt.subplot(3,1,1);
index ='ERFT';
ERFT_dics = spa_pat_reg_mean(index,exp='CamOnly',option='RegMean');	
ax.annotate( r'$\mathrm{\/ERF\/(W\/m^{-2})}$',xy=(-0.08,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=10)
rects1,rects2,rects3 = plot_bar_errbar(ERFT_dics)
ax.set_xlim([-0.5,11.5]);ax.set_xticks([]);
legend1=ax.legend((rects1[0], rects2[0],rects3[0]),\
('Best Estimation','Energy Consumption', 'Technology Advancements'), shadow=False,ncol=1,loc='upper left',fontsize = 10)
legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)


ax = plt.subplot(3,1,2);
index ='TS'
TS_dics = spa_pat_reg_mean(index,exp='FullCp',option='RegMean');
ax.annotate( r'$\mathrm{\/SAT\/(K)}$',xy=(-0.08,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=10)
rects1,rects2,rects3 = plot_bar_errbar(TS_dics)
ax.set_xlim([-0.5,11.5]);ax.set_xticks(np.arange(0.5,11.5,2));ax.set_xticklabels(('Global','Tropics','Arctic','Asia','Europe','USA'),fontsize=10);
ax.set_xlim([-0.5,11.5]);ax.set_xticks([]);

ax = plt.subplot(3,1,3);
index ='precip'
Pr_dics = spa_pat_reg_mean(index,exp='FullCp',option='RegMean');
ax.annotate( r'$\mathrm{\/Pr\/(mm\/day^{-1})}$',xy=(-0.08,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=10)
rects1,rects2,rects3 = plot_bar_errbar(Pr_dics)
ax.set_xlim([-0.5,11.5]);ax.set_xticks(np.arange(0.5,11.5,2));ax.set_xticklabels(('Global','Tropics','Arctic','Asia','Europe','USA'),fontsize=10);

plt.subplots_adjust(left=0.08, bottom=0.05, right=0.98, top=0.95, wspace=0.04, hspace=0.08); 
plt.savefig('RegMean_ERF_TS_Pr.png', format='png', dpi=300)

"""


####################################################################
###  Sensitivity: T vs ERF and Pr vs T
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
	
	y_value = np.array([data['All'][0,1],data['TROPICS'][0,1],data['ARCTIC'][0,1],data['ASIA'][0,1],data['EUROPE'][0,1],data['US'][0,1]])
	# yerr = np.stack((data['All'][0,1:3],data['TROPICS'][0,1:3],data['ARCTIC'][0,1:3],data['ASIA'][0,1:3],data['EUROPE'][0,1:3],data['US'][0,1:3]),axis=1)	
	rects1=ax.bar(X_pos, y_value, align='center',color='k', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos,y_value,1.3)

	y_value = np.array([data['All'][1,1],data['TROPICS'][1,1],data['ARCTIC'][1,1],data['ASIA'][1,1],data['EUROPE'][1,1],data['US'][1,1]])
	# yerr = np.stack((data['All'][1,1:3],data['TROPICS'][1,1:3],data['ARCTIC'][1,1:3],data['ASIA'][0,1:3],data['EUROPE'][1,1:3],data['US'][1,1:3]),axis=1)	
	rects2=ax.bar(X_pos+dis*1, y_value,  align='center',color='r', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*1,y_value,1.3)
	
	y_value = np.array([data['All'][2,1],data['TROPICS'][2,1],data['ARCTIC'][2,1],data['ASIA'][2,1],data['EUROPE'][2,1],data['US'][2,1]])
	# yerr = np.stack((data['All'][2,1:3],data['TROPICS'][2,1:3],data['ARCTIC'][2,1:3],data['ASIA'][0,1:3],data['EUROPE'][2,1:3],data['US'][2,1:3]),axis=1)	
	rects3=ax.bar(X_pos+dis*2, y_value,  align='center',color='g', ecolor='K',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*2,y_value,1.3)

	return rects1,rects2,rects3

	
index ='ERFT';ERFT_dics = spa_pat_reg_mean(index,exp='CamOnly',option='RegMean');	
index ='TS';TS_dics = spa_pat_reg_mean(index,exp='FullCp',option='RegMean');
index ='precip';Pr_dics = spa_pat_reg_mean(index,exp='FullCp',option='RegMean');
Effecacy_dics ={'All':np.empty((3,2)),'TROPICS':np.empty((3,2)),'ASIA':np.empty((3,2)),'EUROPE':np.empty((3,2)),'US':np.empty((3,2)),'ARCTIC':np.empty((3,2))}
Hydro_dics ={'All':np.empty((3,2)),'TROPICS':np.empty((3,2)),'ASIA':np.empty((3,2)),'EUROPE':np.empty((3,2)),'US':np.empty((3,2)),'ARCTIC':np.empty((3,2))} 

for region_key in Effecacy_dics:
	## first three rows, dT over global mean dERF, bottmom three rows regional dT over regional dERF 
	Effecacy_dics[region_key][0,0] = round(np.divide(TS_dics[region_key][0,0],ERFT_dics['All'][0,0]),2)
	Effecacy_dics[region_key][1,0] = round(np.divide(TS_dics[region_key][1,0],ERFT_dics['All'][1,0]),2)
	Effecacy_dics[region_key][2,0] = round(np.divide(TS_dics[region_key][2,0],ERFT_dics['All'][2,0]),2)
	Effecacy_dics[region_key][0,1] = round(np.divide(TS_dics[region_key][0,0],ERFT_dics[region_key][0,0]),2)
	Effecacy_dics[region_key][1,1] = round(np.divide(TS_dics[region_key][1,0],ERFT_dics[region_key][1,0]),2)
	Effecacy_dics[region_key][2,1] = round(np.divide(TS_dics[region_key][2,0],ERFT_dics[region_key][2,0]),2)
	## first three rows, dPr over global mean dT, bottmom three rows regional dPr over regional dT
	# # The precipitaiton is calculated into percent chaneg relative to 1970
	Hydro_dics[region_key][0,0] = round(np.divide(100*Pr_dics[region_key][0,0]/Pr_dics[region_key][3,0],TS_dics['All'][0,0]),2)
	Hydro_dics[region_key][1,0] = round(np.divide(100*Pr_dics[region_key][1,0]/Pr_dics[region_key][3,0],TS_dics['All'][1,0]),2)
	Hydro_dics[region_key][2,0] = round(np.divide(100*Pr_dics[region_key][2,0]/Pr_dics[region_key][3,0],TS_dics['All'][2,0]),2)
	Hydro_dics[region_key][0,1] = round(np.divide(100*Pr_dics[region_key][0,0]/Pr_dics[region_key][3,0],TS_dics[region_key][0,0]),2)
	Hydro_dics[region_key][1,1] = round(np.divide(100*Pr_dics[region_key][1,0]/Pr_dics[region_key][3,0],TS_dics[region_key][1,0]),2)
	Hydro_dics[region_key][2,1] = round(np.divide(100*Pr_dics[region_key][2,0]/Pr_dics[region_key][3,0],TS_dics[region_key][2,0]),2)	


fig = plt.figure(facecolor='White',figsize=[10,7.5]);pad= 5;
ax = plt.subplot(2,1,1);
ax.annotate( r'$\mathrm{\mathsf{\Delta}T}$'+'/' r'$\mathrm{\mathsf{\Delta}ERF}$',xy=(-0.08,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=10)
rects1,rects2,rects3 = plot_bar_errbar(Effecacy_dics)
ax.set_xlim([-0.5,11.5]);ax.set_xticks([]);
legend1=ax.legend((rects1[0], rects2[0],rects3[0]),\
('BEoA ( 2010 -1970)','Energy Consumption', 'Technology Advancements'), shadow=False,ncol=1,loc='upper left',fontsize = 10)
legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)


ax = plt.subplot(2,1,2);
ax.annotate( r'$\mathrm{\mathsf{\Delta}Pr}$'+'/' r'$\mathrm{\mathsf{\Delta}T}$',xy=(-0.08,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=10)
rects1,rects2,rects3 = plot_bar_errbar(Hydro_dics)
ax.set_xlim([-0.5,11.5]);ax.set_xticks(np.arange(0.5,11.5,2));ax.set_xticklabels(('Global','Tropics','Arctic','Asia','Europe','USA'),fontsize=10);

plt.subplots_adjust(left=0.08, bottom=0.05, right=0.98, top=0.95, wspace=0.04, hspace=0.08); 
plt.savefig('Effecacy_HydroSensitivity.png', format='png', dpi=300)
"""