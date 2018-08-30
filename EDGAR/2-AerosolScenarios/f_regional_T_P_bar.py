"""
This script is to produce the bar plot of SSAT and P over regional scales
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
	
	
def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area

	
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
		# missing_value = nc_fid.variables[variable].missing_value; #print missing_value
		# fill_value = nc_fid.variables[variable]._FillValue; #print fill_value
		# data[data==missing_value] = np.nan;data[data==fill_value] = np.nan;
		nc_fid.close()
		if variable=='TS':
			data = data -273.15
		elif variable=='PRECL' or variable=='PRECC':
			data = data*24*60*60*1000  # m/s to mm/day
		var40map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var40map
	
	lon,lat,T1970 = data_netcdf('T1970RCP',FREQ,variable)
	lon,lat,Edg70GO = data_netcdf('Edg70GO',FREQ,variable)
	_,_,EdgRef = data_netcdf('EdgRef',FREQ,variable)
	_,_,EdgEne = data_netcdf('EdgEne',FREQ,variable)
	_,_,EdgTech = data_netcdf('EdgTech',FREQ,variable)
	return lon,lat,T1970,Edg70GO,EdgRef,EdgEne,EdgTech

	
def print_domain_mean(variable,FREQ):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""

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
		print colum_s,colum_e,row_s,row_e
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
		box_region_dic={'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'EA':[100,145,20,50],'SA':[65,100,5,30],'SESA':[295,315,-40,-25],'ARCTIC':[0,360,70,90]}
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
		mask[mask >= 1] = 1;mask[mask < 1] = np.nan;
		# weight each grid cell by its area weight against the total area
		mask=np.multiply(mask,area);  
		mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		return mask_weighted	

	
	def dif_std_for_region(var1,var2,mask):
		"""
		for a regional dif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		
		var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
		var2_domain_mean =np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
		dif =  np.nanmean(var1_domain_mean - var2_domain_mean,axis=0)
		p25 =  np.abs(dif - np.nanpercentile(var1_domain_mean - var2_domain_mean,25,axis=0))/1.25
		p75 =  np.abs(np.nanpercentile(var1_domain_mean - var2_domain_mean,75,axis=0) - dif)/1.25
		return dif,p25,p75
	
	if variable == 'precipitation'	:
		lon,lat,T1970_C,Edg70GO_C,EdgRef_C,EdgEne_C,EdgTech_C = data_readin('PRECL',FREQ);
		lon,lat,T1970_L,Edg70GO_L,EdgRef_L,EdgEne_L,EdgTech_L = data_readin('PRECC',FREQ);
		T1970=(T1970_C+T1970_L)
		Edg70GO=(Edg70GO_C+Edg70GO_L)
		EdgRef=(EdgRef_C+EdgRef_L)
		EdgEne=(EdgEne_C+EdgEne_L)
		EdgTech=(EdgTech_C+EdgTech_L) 
	else:
		lon,lat,T1970,Edg70GO,EdgRef,EdgEne,EdgTech = data_readin(variable,FREQ);
	region_dics ={'All':np.empty((4,3)),'ASIA':np.empty((4,3)),'Europe':np.empty((4,3)),'USA':np.empty((4,3)),'ARCTIC':np.empty((4,3))}  #,'USA','SESA','EA','India'
	for region_key in region_dics:
		_,mask = mask_weight(region_key,lon,lat); #print mask
		region_dics[region_key][0,:] = dif_std_for_region(EdgRef,T1970,mask)  #BEoA
		region_dics[region_key][1,:] = dif_std_for_region(Edg70GO,T1970,mask)  #BEoA
		region_dics[region_key][2,:] = dif_std_for_region(EdgRef,EdgEne,mask)  #ene
		region_dics[region_key][3,:] = dif_std_for_region(EdgRef,EdgTech,mask) ##tech
	return region_dics


FREQ='mon'

def plot_bar(data):
	# print data
	# print "------"+var+"--------"
	def autolabel(X_pos,values,height_lift):
		## Attach a text label above each bar displaying its height
		height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
		for i in range(len(height)):
			ax.text(X_pos[i],y_pos[i],'%4.2f' % height[i], ha='center', va='bottom',size=10)
	
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 2);ax.axvline(x=2,color='k',linewidth = 1);
	ax.axvline(x=4.5,color='k',linewidth = 1);ax.axvline(x=7.0,color='k',linewidth = 1);
	ax.axvline(x=9.5,color='k',linewidth = 1);
	# ax.set_xlim([-0.5,9.5]);ax.set_xticks(np.arange(0.75,10,2.5));ax.set_xticklabels(('Global','Asia','Europe','USA'),fontsize=15);	
	# ax.set_ylim([-6,6.0]);ax.set_yticks(np.arange(-6,6.01,2))
	
	X_pos=np.arange(0,12.5,2.5);dis=0.5
	
	y_value = np.array([data['All'][0,0],data['ASIA'][0,0],data['Europe'][0,0],data['USA'][0,0],data['ARCTIC'][0,0]])
	yerr = np.stack((data['All'][0,1:3],data['ASIA'][0,1:3],data['Europe'][0,1:3],data['USA'][0,1:3],data['ARCTIC'][0,1:3]),axis=1)	
	print np.shape(yerr)
	rects1=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='K', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos,y_value,1.05)

	y_value = np.array([data['All'][1,0],data['ASIA'][1,0],data['Europe'][1,0],data['US'][1,0],data['ARCTIC'][1,0]])
	yerr = np.stack((data['All'][1,1:3],data['ASIA'][1,1:3],data['Europe'][1,1:3],data['USA'][1,1:3],data['ARCTIC'][1,1:3]),axis=1)	
	rects2=ax.bar(X_pos+dis*1, y_value, yerr=yerr, align='center',color='b', ecolor='b',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*1,y_value,1.05)
	
	y_value = np.array([data['All'][2,0],data['ASIA'][2,0],data['Europe'][2,0],data['USA'][2,0],data['ARCTIC'][2,0]])
	yerr = np.stack((data['All'][2,1:3],data['ASIA'][2,1:3],data['Europe'][2,1:3],data['USA'][2,1:3],data['ARCTIC'][2,1:3]),axis=1)	
	rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr, align='center',color='r', ecolor='r',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*2,y_value,1.05)

	y_value = np.array([data['All'][3,0],data['ASIA'][3,0],data['Europe'][3,0],data['USA'][3,0],data['ARCTIC'][3,0]])
	yerr = np.stack((data['All'][3,1:3],data['ASIA'][3,1:3],data['Europe'][3,1:3],data['USA'][3,1:3],data['ARCTIC'][3,1:3]),axis=1)	
	rects4=ax.bar(X_pos+dis*3, y_value, yerr=yerr, align='center',color='g', ecolor='g',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*3,y_value,1.05)
	
	return rects1,rects2,rects3,rects4
	
fig = plt.figure(facecolor='White',figsize=[10,5]);plot_setup();pad= 3;

ax = plt.subplot(2,1,1);
ax.annotate( r'$\mathrm{\mathsf{\Delta}\/SAT\/(K)}$',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)
index ='TS'	  
TS_dics = print_domain_mean(index,FREQ);
rects1,rects2,rects3,rects4 = plot_burden_bar(TS_dics)
ax.set_ylim([-4.5,4.5]);ax.set_xlim([-0.5,12]);ax.set_xticklabels(('','','',''));
# ax.set_yticklabels(np.arange(-4.5,4.51,1.5));
legend1=ax.legend((rects1[0], rects2[0],rects3[0], rects4[0]),\
('Historical (2010 - 1970)', 'BEoA ( 2010 -1970)','Energy Consumption', 'Technology Advancements'), shadow=False,ncol=1,loc='lower left')
legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)


ax = plt.subplot(2,1,2);
ax.annotate( r'$\mathrm{\mathsf{\Delta}\/P\/(mm\/day^{-1})}$',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='vertical',fontsize=15)

index ='precipitation'	  
P_dics = print_domain_mean(index,FREQ);
rects1,rects2,rects3,rects4 = plot_burden_bar(P_dics)
ax.set_xlim([-0.5,12]);ax.set_xticks(np.arange(0.75,12.5,2.5));ax.set_xticklabels(('Global','Asia','Europe','USA','Arctic'),fontsize=12);


plt.subplots_adjust(left=0.10, bottom=0.10, right=0.98, top=0.95, wspace=0.04, hspace=0.08); 
plt.savefig('Regional_response.png', format='png', dpi=300)