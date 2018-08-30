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
    os.path.pardir,os.path.pardir,
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
	def annual_and_monthly(scenario,time,data):
		# plt.imshow(data[3,:,:]);plt.show()		
		annual_month=np.empty((40,13,np.shape(data)[1],np.shape(data)[2]));annual_month[:]=np.nan
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
				layer = [layer for layer in range(len(time)) if (time[layer]//100  == iyear*100+imonth)]
				annual_month[iyear-year_series[0],imonth,:,:] = data[layer,:,:]
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
		elif variable == 'iceextent':
			var_path = input_path+scenario+'/mon/ice/'+scenario+'.ice.mon.aice.nc'			
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['TLAT'][:]
			lon = nc_fid.variables['TLON'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			data = nc_fid.variables['aice'][:]
			ms_value = nc_fid.variables['aice'].missing_value
			fi_value = nc_fid.variables['aice']._FillValue
			data[data == ms_value]=np.nan;data[data==fi_value]=np.nan;
			nc_fid.close()
			data =data*0.01   # from percent to fraction
			data[data>=0.15]=1; #data[data<0.15]=0 #15% threshold tp account for sea ice content rather than area
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
			if variable == 'Tsfc':
				data = np.ma.masked_where(data==-1.836000,data)
			elif variable == 'Tref':
				data=data-273.15
		annual_monthly = annual_and_monthly(scenario,time,data)
		# print annual_monthlyk
		return lon,lat,annual_monthly
	
	lon,lat,Edg70GO = data_netcdf('Edg70GO',variable);
	_,_,T1970 = data_netcdf('T1970RCP',variable)
	_,_,EdgRef = data_netcdf('EdgRef',variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',variable)
	_,_,EdgEne = data_netcdf('EdgEne',variable)
	_,_,EdgTech = data_netcdf('EdgTech',variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	
def annual_spatial_and_mon_vary(variable,pole):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""
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
		TLAT[abs(TLAT)<70] = np.nan; TLAT[abs(TLAT)>=70] = 1;
		mask[mask<1]=np.nan
		mask= np.multiply(TLAT,mask)
		nc_fid.close()
		mask_area = np.multiply(area,mask)
		return mask_area


	def spatial_diff_sig(vairable1,variable2):
		"""
		####calculate spatial the difference and their significance
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
		dif = np.nanmean(vairable1[:,0,:,:]-variable2[:,0,:,:],axis=0)
		sig = mannwhitneyu_test(vairable1[:,0,:,:], variable2[:,0,:,:])
		return dif,sig
	
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
					masked[iyear,imonth] = np.nansum(np.nansum(np.multiply(cache,area),axis=1),axis=0)
			return masked
		if variable == 'alidr' or variable == 'alvdr' or variable == 'Tsfc' or variable == 'Tref' or variable == 'albice':
			result = area_mean_or_sum(data,pole,area,option='area_mean') 
		elif variable == 'hi' or variable == 'hs':# ice/snow volume from m to 1000 km in order to get 1000km3
			result = area_mean_or_sum(data*10**(-6),pole,area,option='area_sum')  
		elif variable == 'iceextent': # ice extent 10**(-6) to convert into 10**(6) KM2
			result = area_mean_or_sum(data*10**(-6),pole,area,option='area_sum')  
		return result

	def TimeSeries_diff_std(var1,var2):
		"""
		# calculate  the difference and their significance of time series
		# in this case the 12 months and the annual mean
		"""
		diff = var1 - var2
		# print "diff", np.shape(diff)
		mean = np.nanmean(diff,axis=0)
		P25 = mean - np.abs(mean-np.nanpercentile(diff,25,axis=0))/1.5
		P75 = mean + np.abs(np.nanpercentile(diff,75,axis=0)-mean)/1.5
		return mean,P25,P75	
	
	area = get_area_mask(); 
	lon,lat,T1970,Edg70GO,_,EdgRef,EdgEne,EdgTech = data_readin(variable);
	T1970_ts = monthly_ts(T1970,pole,area);EdgRef_ts = monthly_ts(EdgRef,pole,area);
	Edg70GO_ts = monthly_ts(Edg70GO,pole,area);EdgRef_ts = monthly_ts(EdgRef,pole,area);
	EdgTech_ts = monthly_ts(EdgTech,pole,area);EdgEne_ts = monthly_ts(EdgEne,pole,area);
	print T1970_ts
	##spatial dif and sig
	AEROSOL,AEROSOL_s = spatial_diff_sig(Edg70GO,T1970)
	Ene,Ene_s = spatial_diff_sig(EdgRef,EdgEne)
	Tech,Tech_s = spatial_diff_sig(EdgRef, EdgTech)	
	ref = np.nanmean(EdgRef[:,0,:,:],axis=0)
	spatial_dif = np.stack((ref,AEROSOL,Ene,Tech),axis=0)
	spatial_sig = np.stack((AEROSOL_s,Ene_s,Tech_s),axis=0)
	
	## monthly mean fields
	AEROSOL,AEROSOL_25,AEROSOL_75 = TimeSeries_diff_std(Edg70GO_ts,T1970_ts);
	Ene,Ene_25,Ene_75 = TimeSeries_diff_std(EdgRef_ts,EdgEne_ts); 
	Tech,Tech_25,Tech_75 = TimeSeries_diff_std(EdgRef_ts, EdgTech_ts)
	TA,TA_25,TA_75 = TimeSeries_diff_std(EdgRef_ts,T1970_ts)
	ref1970,ref1970_25,ref1970_75 = TimeSeries_diff_std(T1970_ts,0)
	
	TimeSeries_dif = np.stack((ref1970,TA,AEROSOL,Ene,Tech),axis=0)
	TimeSeries_P25 = np.stack((ref1970_25,TA_25,AEROSOL_25,Ene_25,Tech_25),axis=0)
	TimeSeries_P75 = np.stack((ref1970_75,TA_75,AEROSOL_75,Ene_75,Tech_75),axis=0)
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
	
	m = Basemap(projection=projection,boundinglat=boundinglat,lon_0=0,resolution='h',round=True)
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
	cmap = discrete_cmap(10,colormap)
	# cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	colormesh = m.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)
	pm = np.ma.masked_not_equal(p_value, 1)
	m.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='.', alpha=0.,lw=1,vmin=colorbar_min, vmax=colorbar_max)
	# plt.colorbar(orientation='vertical',extend='both',shrink=0.8,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
	return colormesh

def plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=12,color='k',ls = '--',linewidth = 1);
	ax.set_xlim([0,15]);ax.set_xticks(np.array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,13.5]));
	ax.set_xticklabels(('Jan','Feb','Mar','April','May','June','July','Aug','Sep','Oct','Nov','Dec','Annual'));
	ax.axvspan(12, 15, alpha=0.1, color='y',edgecolor='none');
	x=np.arange(0.5,12)
	##line
	# ax.plot(x,TimeSeries_dif[1,1:13],'-',color="k",linewidth=3,label = 'E1970')
	ax.plot(x,TimeSeries_dif[2,1:13],'-',color="b",linewidth=3,label = 'BEoA (2010-1970)')
	ax.plot(x,TimeSeries_dif[3,1:13],'-',color="r",linewidth=3,label = 'Energy Consumption')
	ax.plot(x,TimeSeries_dif[4,1:13],'-',color="g",linewidth=3,label = 'Technology Advancements')
	# ax.plot(x,TimeSeries_dif[4,1:13],'-',color="m",linewidth=3,label = 'O3')
	## range
	# ax.fill_between(x,TimeSeries_P25[1,1:13],TimeSeries_P75[1,1:13],facecolor="k",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[2,1:13],TimeSeries_P75[2,1:13],facecolor="b",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[3,1:13],TimeSeries_P75[3,1:13],facecolor="r",alpha=0.15)
	ax.fill_between(x,TimeSeries_P25[4,1:13],TimeSeries_P75[4,1:13],facecolor="g",alpha=0.15)
	# ax.fill_between(x,TimeSeries_P25[4,1:13],TimeSeries_P75[4,1:13],facecolor="m",alpha=0.15)
	## bar
	# ax.errorbar(np.array([12.5]), TimeSeries_dif[1,0], yerr=TimeSeries_P75[1,0]-TimeSeries_P25[1,0],fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='k')
	ax.errorbar(np.array([12.5]), TimeSeries_dif[2,0], yerr=(TimeSeries_P75[2,0]-TimeSeries_P25[2,0])/2,fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='b')
	ax.errorbar(np.array([13]), TimeSeries_dif[3,0], yerr=(TimeSeries_P75[3,0]-TimeSeries_P25[3,0])/2,fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='r')
	ax.errorbar(np.array([14]), TimeSeries_dif[4,0], yerr=(TimeSeries_P75[4,0]-TimeSeries_P25[4,0])/2,fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='g')	

	ax2 = ax.twinx();
	ax2.plot(x,TimeSeries_dif[1,1:13],'-',color="k",linewidth=4, label = "Historical (2010-1970)")
	ax2.fill_between(x,TimeSeries_P25[1,1:13],TimeSeries_P75[1,1:13] ,facecolor="k",alpha=0.15)
	ax2.errorbar(np.array([14.5]),TimeSeries_dif[1,0], yerr=(TimeSeries_P75[1,0]-TimeSeries_P25[1,0])/2,fmt='+w', lw=6,capsize=0,markeredgecolor='w',ms=7,c='k')	
	print np.nanmean(TimeSeries_dif[1,1:13]), TimeSeries_dif[1,0]
	ax2.set_xlim([0,15]);ax2.set_xticks(np.array([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,13.5]));
	return ax,ax2

pole = "Arctic"; season = 'sep';colormap='RdBu_r';
fig = plt.figure(facecolor='White',figsize=[18,12]);pad= 5  #plot_setup()

index = 'Tref' 
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(2,2,1);'Annual peak intensity\n($^\circ$C)'
ax.annotate('(a) Arctic mean (70-90N) surface air temperature (K)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)			
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-4,4]);ax.set_yticks(np.arange(-4,4.1,1));
ax2.set_ylim([1,10]);ax2.set_yticks(np.arange(1,10.1,1))
align_yaxis(ax1,0,ax2,5.5)


ax = plt.subplot(2,2,3);
index = 'iceextent' 
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);
ax.annotate(r'(b) Sea ice extent $\mathrm{(10^{6} km^{2})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
				
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-2,2]);ax.set_yticks(np.arange(-2,2.1,0.5));
ax2.set_ylim([-5,0.0]);ax2.set_yticks(np.arange(-5,0.1,0.5))
align_yaxis(ax1,0,ax2,-2.5)
legend = ax2.legend(shadow=False,ncol=1,loc ='lower right',fontsize=15)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

legend = ax1.legend(shadow=False,ncol=1,loc ='lower left',fontsize=15)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

ax = plt.subplot(2,2,4);
index = 'hi'  
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);
ax.annotate(r'(c) Sea ice volume $\mathrm{(10^{3} km^{3})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
				
ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
ax1.set_ylim([-6,6]);ax.set_yticks(np.arange(-6,6.1,1));
ax2.set_ylim([-19,-14.0]);ax2.set_yticks(np.arange(-19,-13.9,0.5))
align_yaxis(ax1,0,ax2,-16.5)
plt.subplots_adjust(left=0.03, bottom=0.08, right=0.96, top=0.95, wspace=0.15, hspace=0.20); 
plt.savefig('./'+pole+'_seaice_seasonal_variation.png', format='png', dpi=1000)

# index = 'Tref'
# lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);
# TimeSeries_P25 = TimeSeries_dif - np.abs(TimeSeries_dif-TimeSeries_P25)/2.0
# TimeSeries_P75 = TimeSeries_dif + np.abs(TimeSeries_dif-TimeSeries_P75)/1.5
# from get_global_monthly_T import get_global_monthly_T
# TS_dics =get_global_monthly_T(FREQ='mon',index ='TS')
# for i in range (1,5):   ##  normalize against global mean warming to see the arctic amplificaiton
	# norm = TS_dics['All'][i-1,:,0]; print np.shape(norm)
	# TimeSeries_dif[i,:] = np.divide(TimeSeries_dif[i,:],norm)
	# TimeSeries_P25[i,:] = np.divide(TimeSeries_P25[i,:],norm)
	# TimeSeries_P75[i,:] = np.divide(TimeSeries_P75[i,:],norm)

# ax = plt.subplot(2,2,2);'Annual peak intensity\n($^\circ$C)'
# ax.annotate('(b) Normalized Arctic mean surface air temperature',xy=(0.02,1.02), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=15)			
# ax1,ax2 = plot_annual_month_mean_uncertainty(ax,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75)
# ax1.axhline(y=1,color='r',linewidth = 2);

# print np.round(TimeSeries_dif,1)
# ax1.set_ylim([-10,50]);ax.set_yticks(np.arange(-10,50.10,10));
# ax2.set_ylim([-10,50]);ax2.set_yticks(np.arange(-10,50.10,10));
# align_yaxis(ax1,0,ax2,0)

"""

#####spatial patterns

fig = plt.figure(facecolor='White',figsize=[16,15]);pad= 5  #plot_setup();
colormap='RdBu_r';  #colormap = reverse_colourmap(colormap);


index ='Tref'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT

if season=='sep':
	colorbar_min=-3;colorbar_max=3;	
else:
	colorbar_min=-1.8;colorbar_max=1.8;
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(3,3,1);
ax.annotate('(a)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('BEoA (2010-1970)',xy=(0.5,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	
# ax.annotate(r'CDNC ($\mathrm{\mathsf{10^{10}\/m^{-2}}}$)',xy=(-0.1,0.5), xytext=(0, pad),
ax.annotate('2m Air temperature (K)',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)	
polar_figure(ax,spatial_dif[1,:,:],spatial_sig[0,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(3,3,2);
ax.annotate('(d)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Energy Consumption',xy=(0.5,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)		
polar_figure(ax,spatial_dif[2,:,:],spatial_sig[1,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(3,3,3);
ax.annotate('(g)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Technology advancements',xy=(0.5,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)		
colormesh1 = polar_figure(ax,spatial_dif[3,:,:],spatial_sig[2,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

cbar_ax = fig.add_axes([0.92, 0.67, 0.015, 0.28])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))


index ='aice'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT

if season=='sep':
	colorbar_min=-30;colorbar_max=30;
else:
	colorbar_min=-10;colorbar_max=10;		
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(3,3,4);
ax.annotate('(b)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Gridcell sea ice area (%)',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)	
polar_figure(ax,spatial_dif[1,:,:],spatial_sig[0,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(3,3,5);
ax.annotate('(e)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
polar_figure(ax,spatial_dif[2,:,:],spatial_sig[1,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(3,3,6);
ax.annotate('(h)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = polar_figure(ax,spatial_dif[3,:,:],spatial_sig[2,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

cbar_ax = fig.add_axes([0.92, 0.35, 0.015, 0.28])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

index ='hi'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT

if season=='sep':
	colorbar_min=-0.8;colorbar_max=0.8;
else:	
	colorbar_min=-0.6;colorbar_max=0.6;
lon,lat,spatial_dif,spatial_sig,TimeSeries_dif,TimeSeries_P25,TimeSeries_P75  = annual_spatial_and_mon_vary(index,pole);

ax = plt.subplot(3,3,7);
ax.annotate('(c)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Gridcell mean ice thickness (m)',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)
polar_figure(ax,spatial_dif[1,:,:],spatial_sig[0,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(3,3,8);
ax.annotate('(f)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
polar_figure(ax,spatial_dif[2,:,:],spatial_sig[1,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(3,3,9);
ax.annotate('(i)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = polar_figure(ax,spatial_dif[3,:,:],spatial_sig[2,:,:],lon,lat,colormap,colorbar_min,colorbar_max)

cbar_ax = fig.add_axes([0.92, 0.035, 0.015, 0.28])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.08, bottom=0.03, right=0.90, top=0.96, wspace=0.02, hspace=0.02); 
plt.savefig('./'+pole+'_seaice_spatial_sep.png', format='png', dpi=1000)
"""


