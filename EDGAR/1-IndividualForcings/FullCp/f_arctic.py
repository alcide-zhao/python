"""
This is to plot the TS  from the EDGAR six experiments
data inputs are forty years of monthly TS
	first step is to process the monthly TS into monthly mean
	second step is to process annual mean of monthly mean
"""
import site
import os
import numpy as np
np.set_printoptions(threshold='nan')
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


		
def data_readin(variable):	
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		# if scenario =='T1970C': start_year =1970
		start_year =0
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
		
		
	def annual_and_monthly(scenario,time,data):
		# plt.imshow(data[3,:,:]);plt.show()		
		annual_month=np.empty((30,13,np.shape(data)[1],np.shape(data)[2]));annual_month[:]=np.nan
		# plt.imshow(annual_month[0,3,:,:]);plt.show()
		if scenario=='T1970RCP':
			year_series = range(2020,2050)
		elif scenario=='EdgEne':
			year_series = range(2200,2230)
		elif scenario=='Edg70GO':
			year_series = range(2070,2100)
		else:
			year_series = range(2130,2160)
		for iyear in year_series:
			# print iyear*100+12,time//100 
			if (iyear == year_series[0] and time[0]//100 > year_series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  #January
			# print layer_b
			if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+12][0]   #december
			data_cache = data[layer_b:layer_e+1,:,:]
			annual_month[iyear-year_series[0],0,:,:] = np.nanmean(data_cache,axis=0)
			# plt.imshow(annual_month[iyear-year_series[0],0,:,:]);plt.show()
		for imonth in range(1,13):
			layer = [layer for layer in range(len(time)) if (time[layer]%10000)/100  == imonth]  #August 31
			annual_month[0:len(layer),imonth,:,:] = data[layer,:,:]   #:len(layer)
			# plt.imshow(np.nanmean(data[layer,:,:],axis=0));plt.show()
		return annual_month	
	
	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		if variable == "albice":
			var_path = input_path+scenario+'/mon/ice/'+scenario+'.ice.mon.aice.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['TLAT'][:]
			lon = nc_fid.variables['TLON'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			aice = nc_fid.variables['aice'][:]
			ms_value = nc_fid.variables['aice'].missing_value
			fi_value = nc_fid.variables['aice']._FillValue
			aice[aice == ms_value]=np.nan;aice[aice==fi_value]=np.nan;
			nc_fid.close()
			var_path = input_path+scenario+'/mon/ice/'+scenario+'.ice.mon.fswdn.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			fswdn = nc_fid.variables['fswdn'][:]
			ms_value = nc_fid.variables['fswdn'].missing_value
			fi_value = nc_fid.variables['fswdn']._FillValue
			fswdn[fswdn == ms_value]=np.nan;fswdn[fswdn==fi_value]=np.nan;
			nc_fid.close()
			var_path = input_path+scenario+'/mon/ice/'+scenario+'.ice.mon.fswup.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			fswup = nc_fid.variables['fswup'][:]
			ms_value = nc_fid.variables['fswup'].missing_value
			fi_value = nc_fid.variables['fswup']._FillValue
			fswup[fswup == ms_value]=np.nan;fswup[fswup==fi_value]=np.nan;
			nc_fid.close()
			data =np.multiply(aice,np.divide(fswup,fswdn))
		else:
			var_path = input_path+scenario+'/mon/ice/'+scenario+'.ice.mon.'+variable+'.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['TLAT'][:]
			lon = nc_fid.variables['TLON'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			data = nc_fid.variables[variable][:]
			ms_value = nc_fid.variables[variable].missing_value
			fi_value = nc_fid.variables[variable]._FillValue
			data[data == ms_value]=np.nan;data[data==fi_value]=np.nan;
			nc_fid.close()
		if variable == 'aice':
			data =data*0.01   # from percent to fraction
			data[data>=0.15]=1; #data[data<0.15]=0 #15% threshold tp account for sea ice content rather than area
		elif variable == 'Tsfc':
			data = np.ma.masked_where(data==-1.836000,data)
		elif variable == 'Tref':
			data=data-273.15
		annual_monthly = annual_and_monthly(scenario,time,data)
		return lon,lat,annual_monthly
	
	lon,lat,Edg70GO = data_netcdf('Edg70GO',variable)
	_,_,T1970 = data_netcdf('T1970RCP',variable)
	_,_,EdgRef = data_netcdf('EdgRef',variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',variable)
	_,_,Edg70T10SOZ = data_netcdf('Edg70T10SOZ',variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ
	
def annual_spatial_and_mon_vary(variable,pole):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""
	def spatial_diff_sig(vairable1,variable2):
		"""
		####calculate spatial the difference and their significance
		"""
		def mannwhitneyu_test(vairable1,variable2):
			p_threshold=0.15
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
		dif = np.nanmean(vairable1[:,0,:,:]-variable2[:,0,:,:],axis=0)
		sig = mannwhitneyu_test(vairable1[:,0,:,:], variable2[:,0,:,:])
		# sig_dif = np.multiply(dif,sig)
		
		return dif,sig
	
	def TimeSeries_diff_std(var1,var2):
		"""
		# calculate  the difference and their significance of time series
		# in this case the 12 months and the annual mean
		"""
		diff = var1 - var2
		# print "diff", np.shape(diff)
		mean = np.nanmean(diff,axis=0)
		P25 = mean - np.abs(mean-np.nanpercentile(diff,25,axis=0))/1.25
		P75 = mean + np.abs(np.nanpercentile(diff,75,axis=0)-mean)/1.25
		return mean,P25,P75
	def get_area_mask():
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+'cice_tarea.nc'
		nc_fid = nc4.Dataset(var_path,mode='r')
		area = nc_fid.variables['tarea'][:]*10**(-6)  # from m2 to km2
		nc_fid.close()
		#mask
		var_path = input_path+'cice_tmask.nc'
		nc_fid = nc4.Dataset(var_path,mode='r')
		mask = nc_fid.variables['tmask'][:] # from m2 to km2
		TLAT = nc_fid.variables['TLAT'][:]; 
		TLAT[abs(TLAT)<60] = np.nan; TLAT[abs(TLAT)>=60] = 1;
		mask[mask<1]=np.nan
		mask= np.multiply(TLAT,mask)
		nc_fid.close()
		mask_area = np.multiply(area,mask)
		return mask_area
	
	def monthly_ts(data,pole,area):
		def area_mean_or_sum(datain,pole,area,option):
			"""
			# calculate either area  weighted mean or sum of the valid points, 
			# only over the points that sea ice exists ot albedo value exists
			# The Arctic and Antarctic is different from each other
			"""
			if pole =="Arctic":
				data = datain[:,:,192:384,:]; area = area[192:384,:]
			elif pole =="Antarctic":
				data = datain[:,:,0:192,:]; area = area[0:192,:]
			masked = np.empty((40,13));masked[:]=np.nan
			for iyear in range(np.shape(data)[0]):
				for imonth in range(np.shape(data)[1]):
					cache = data[iyear,imonth,:,:];
					area[np.isnan(cache)]= np.nan
					# plt.imshow(area);plt.show()
					if option =='area_mean':  # weight the data with area
						area = area/np.nansum(np.nansum(area,axis=1),axis=0)
					masked [iyear,imonth] = np.nansum(np.nansum(np.multiply(cache,area),axis=1),axis=0)
			return masked
		if variable == 'alidr' or variable == 'alvdr' or variable == 'Tsfc' or variable == 'Tref' or variable == 'albice':
			result = area_mean_or_sum(data,pole,area,option='area_mean') 
		elif variable == 'hi' or variable == 'hs':# ice/snow volume from m to 1000 km in order to get 1000km3
			result = area_mean_or_sum(data*10**(-6),pole,area,option='area_sum')  
		elif variable == 'fs': # snow area
			result = area_mean_or_sum(data,pole,area,option='area_sum')  
		elif variable == 'aice': # ice area, it was in percent so need tobe divided by 100  and then 10**(-6) to convert into 10**(6) KM2
			result = area_mean_or_sum(data*10**(-6),pole,area,option='area_sum')  
		return result
	
	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ = data_readin(variable);
	area = get_area_mask()
	T1970_ts = monthly_ts(T1970,pole,area); Edg70GO_ts = monthly_ts(Edg70GO,pole,area);
	Edg70Oz_ts = monthly_ts(Edg70Oz,pole,area);EdgRef_ts = monthly_ts(EdgRef,pole,area);
	Edg70T10SOZ_ts = monthly_ts(Edg70T10SOZ,pole,area);
	
	##spatial dif and sig
	GHG,GHG_s = spatial_diff_sig(Edg70Oz,Edg70GO)
	AEROSOL,AEROSOL_s = spatial_diff_sig(Edg70GO,T1970)
	Total,Total_s = spatial_diff_sig(EdgRef,T1970)
	OZONE,OZONE_s = spatial_diff_sig(EdgRef, Edg70Oz)
	TrO3,TrO3_s = spatial_diff_sig(EdgRef, Edg70T10SOZ)
	StO3,StO3_s = spatial_diff_sig(Edg70T10SOZ, Edg70Oz)	
	
	
	ref = np.nanmean(T1970[:,0,:,:],axis=0)
	
	spatial_dif = np.stack((ref,Total,GHG,AEROSOL,TrO3,StO3),axis=0)
	spatial_sig = np.stack((Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s),axis=0)
	
	## monthly sif and sig
	Total2010,Total2010_25,Total2010_75 = TimeSeries_diff_std(EdgRef_ts,T1970_ts)
	GHG,GHG_25,GHG_75 = TimeSeries_diff_std(Edg70Oz_ts,Edg70GO_ts)
	AEROSOL,AEROSOL_25,AEROSOL_75 = TimeSeries_diff_std(Edg70GO_ts,T1970_ts)
	TrO3,TrO3_25,TrO3_75 = TimeSeries_diff_std(EdgRef_ts, Edg70T10SOZ_ts)
	StO3,StO3_25,StO3_75 = TimeSeries_diff_std(Edg70T10SOZ_ts, Edg70Oz_ts)
	ref1970,ref1970_25,ref1970_75 = TimeSeries_diff_std(T1970_ts,0)


	TimeSeries_dif = np.stack((ref1970,Total2010,GHG,AEROSOL,TrO3,StO3),axis=0)
	TimeSeries_P25 = np.stack((ref1970_25,Total2010_25,GHG_25,AEROSOL_25,TrO3_25,StO3_25),axis=0)
	TimeSeries_P75 = np.stack((ref1970_75,Total2010_75,GHG_75,AEROSOL_75,TrO3_75,StO3_75),axis=0)
	return lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75

	
def polar_figure(ax,data,p_value,lons,lats,colormap,colorbar_min,colorbar_max,pole='Arctic'): #c_bad,c_under,c_over,c_number=20,
	# setup north polar stereographic basemap.
	# The longitude lon_0 is at 6-o'clock, and the
	# latitude circle boundinglat is tangent to the edge
	# of the map at lon_0. Default value of lat_ts
	# (latitude of true scale) is pole.
	if pole=='Arctic': 
		projection = 'npstere';boundinglat = 70
	elif pole=='Antarctic': 
		projection = 'spstere';boundinglat = -70
	
	m = Basemap(projection=projection,boundinglat=boundinglat,lon_0=0,resolution='l',round=True)
	m.drawcoastlines()
	m.fillcontinents(color='grey',lake_color='aqua')
	# draw parallels and meridians.
	m.drawparallels(np.arange(-80.,81.,10.))
	m.drawmeridians(np.arange(-180.,181.,30.))
	# m.drawmapboundary(fill_color='aqua')
	# m.bluemarble() # The NASA 'Blue Marble' image
	
	# lon, lat = np.meshgrid(lons, lats);
	xi, yi = m(lon, lat)
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(16,colormap)
	# cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	colormesh = m.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)
	pm = np.ma.masked_not_equal(p_value, 1)
	# m.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	plt.colorbar(orientation='vertical',extend='both',shrink=0.8)
	return colormesh

def plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=12,color='k',ls = '--',linewidth = 1);
	ax.set_xlim([0,15.5]);ax.set_xticks(np.array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,14]));
	ax.set_xticklabels(('Jan','Feb','Mar','April','May','June','July','Aug','Sep','Oct','Nov','Dec','Annual'));
	ax.axvspan(12, 15.5, alpha=0.1, color='y',edgecolor='none');
	x=np.arange(0.5,12)
	##line
	ax.plot(x,TimeSeries_dif[1,1:13],'-',color="b",linewidth=3,label = 'E2010-E1970')
	ax.plot(x,TimeSeries_dif[2,1:13],'-',color="g",linewidth=3,label = 'GHGs')
	ax.plot(x,TimeSeries_dif[3,1:13],'-',color="r",linewidth=3,label = 'AAs')
	ax.plot(x,TimeSeries_dif[4,1:13],'-',color="m",linewidth=3,label = 'Trop. O3')
	ax.plot(x,TimeSeries_dif[5,1:13],'-',color="y",linewidth=3,label = 'Strat. O3')
	
	## range
	ax.fill_between(x,TimeSeries_P25[1,1:13],TimeSeries_P75[1,1:13],facecolor="b",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[2,1:13],TimeSeries_P75[2,1:13],facecolor="g",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[3,1:13],TimeSeries_P75[3,1:13],facecolor="r",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[4,1:13],TimeSeries_P75[4,1:13],facecolor="m",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[5,1:13],TimeSeries_P75[5,1:13],facecolor="y",alpha=0.15)
	## bar
	ax.errorbar(np.array([12.5]), TimeSeries_dif[1,0], yerr=TimeSeries_P75[1,0]-TimeSeries_P25[1,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='b')
	ax.errorbar(np.array([13]), TimeSeries_dif[2,0], yerr=TimeSeries_P75[2,0]-TimeSeries_P25[2,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='g')
	ax.errorbar(np.array([14]), TimeSeries_dif[3,0], yerr=TimeSeries_P75[3,0]-TimeSeries_P25[3,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='r')
	ax.errorbar(np.array([13.5]), TimeSeries_dif[4,0], yerr=TimeSeries_P75[4,0]-TimeSeries_P25[4,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='m')
	ax.errorbar(np.array([14.5]), TimeSeries_dif[5,0], yerr=TimeSeries_P75[5,0]-TimeSeries_P25[5,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='y')


	
	ax2 = ax.twinx();
	ax2.plot(x,TimeSeries_dif[0,1:13],'-',color="k",linewidth=3, label = "E1970")
	ax2.fill_between(x,TimeSeries_P25[0,1:13],TimeSeries_P75[0,1:13] ,facecolor="k",alpha=0.15)
	ax2.errorbar(np.array([15]),TimeSeries_dif[0,0], yerr=TimeSeries_P75[0,0]-TimeSeries_P25[0,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='k')	
	
	ax.set_xlim([0,15.5]);ax.set_xticks(np.array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,13.5]));
	return ax,ax2

	
fig = plt.figure(facecolor='White',figsize=[16,10]);plot_setup();pad= 5
colormap='RdBu_r';


pole = "Arctic"

index = 'albice'  # fs hs hi alvdr alidr aice
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,1);
ax.annotate('(a) Visible Direct Albedo (%)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal',fontsize=10)
			
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-25,5]);ax.set_yticks(np.arange(-25,5.1,5));
ax2.set_ylim([25,75]);ax2.set_yticks(np.arange(25,75.1,10))
align_yaxis(ax1,-10,ax2,50)

index = 'Tref' 
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,2);'Annual peak intensity\n($^\circ$C)'
ax.annotate('(b)  Reference Height (ocean only) Temperature ($^\circ$C)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)			
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-1,8]);ax.set_yticks(np.arange(-1,8.1,1));
ax2.set_ylim([-25,5]);ax2.set_yticks(np.arange(-25,5.1,5))
align_yaxis(ax1,3.5,ax2,-10)
legend = ax2.legend(shadow=False,ncol=1,loc ='upper right',bbox_to_anchor=(0.80, 0.98))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

legend = ax1.legend(shadow=False,ncol=1,loc ='upper left')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

index = 'aice' 
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,3);
ax.annotate(r'(c) Sea Ice Extent $\mathrm{(10^{6} km^{2})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-5,1]);ax.set_yticks(np.arange(-5,1.1,1));
ax2.set_ylim([7,14.0]);ax2.set_yticks(np.arange(7,14.1,1))
align_yaxis(ax1,-2,ax2,10.5)

index = 'hi'  
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,4);
ax.annotate(r'(d) Sea Ice Volume $\mathrm{(10^{3} km^{3})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-20,2]);ax.set_yticks(np.arange(-20,2.1,2));
ax2.set_ylim([22,38]);ax2.set_yticks(np.arange(22,38.1,2))
align_yaxis(ax1,-9,ax2,30)


plt.subplots_adjust(left=0.05, bottom=0.06, right=0.95, top=0.94, wspace=0.15, hspace=0.15); 
plt.savefig('./'+pole+'_albedo.png', format='png', dpi=1000)


"""
pole = "Antarctic"

index = 'alvdr'  # fs hs hi alvdr alidr aice
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,1);
# ax.annotate('Arctic',xy=(0.5,1.05), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='baseline',rotation='horizontal',fontsize=15)
ax.annotate('(a) Visible Direct Albedo (%)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal',fontsize=10)
			
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
# ax1.set_ylim([-20,5]);ax.set_yticks(np.arange(-20,5.1,5));
# ax2.set_ylim([0,70]);ax2.set_yticks(np.arange(0,70.1,10))
# align_yaxis(ax1,-7.5,ax2,35)

legend = ax2.legend(shadow=False,ncol=1,loc ='upper right',bbox_to_anchor=(0.80, 0.98))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

legend = ax1.legend(shadow=False,ncol=1,loc ='upper left')	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

index = 'Tref'  # fs hs hi alvdr alidr aice
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,2);
ax.annotate('(b) Reference Height (2 m) SAT (degree C)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)			
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
# ax1.set_ylim([-1,8]);ax.set_yticks(np.arange(-1,8.1,1));
# ax2.set_ylim([-25,5]);ax2.set_yticks(np.arange(-25,5.1,5))
# align_yaxis(ax1,3.5,ax2,-10)

index = 'aice'  # fs hs hi alvdr alidr aice
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,3);
ax.annotate(r'(c) Sea Ice Extent $\mathrm{(10^{6} km^{2})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
# ax1.set_ylim([-5,1]);ax.set_yticks(np.arange(-5,1.1,1));
# ax2.set_ylim([7,14.0]);ax2.set_yticks(np.arange(7,14.1,1))
# align_yaxis(ax1,-2,ax2,10.5)

index = 'hi'  # fs hs hi alvdr alidr aice
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,4);
ax.annotate(r'(d) Sea Ice Volume $\mathrm{(10^{3} km^{3})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
# ax1.set_ylim([-20,2]);ax.set_yticks(np.arange(-20,2.1,2));
# ax2.set_ylim([22,38]);ax2.set_yticks(np.arange(22,38.1,2))
# align_yaxis(ax1,-9,ax2,30)
"""



