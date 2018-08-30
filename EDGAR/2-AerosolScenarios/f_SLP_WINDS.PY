"""
This is to plot the slp overlapped with WIND850 from the EDGAR six experiments
	interpolation is needed to make the plots legible
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

def spatial_scaler_vector(ax,lons,lats,colormap,colorbar_min,colorbar_max,scaler,vector_u,vector_v,qk_scale=None,qk_caption=None,qk_is = False):
	
	
	speed = np.sqrt(np.power(AA_U,2)+np.power(AA_V,2))
	vector_u =np.ma.masked_where(speed<0.02,vector_u)
	vector_v =np.ma.masked_where(speed<0.02,vector_v)
	lons[lons>180]-=360; 
	lon_0 = lons.mean();lat_0 = lats.mean()
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 20; lat_bin = 15
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=ax)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	cmap = discrete_cmap(10,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over([0,0,0]); cmap.set_under('b')
	masked_obj = np.ma.masked_where(np.isnan(scaler), scaler)
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)  #norm=norm,
	map.drawcoastlines(); 
	Q = map.quiver(xi,yi,vector_u,vector_v,scale=None,latlon=True)
	if qk_is:
		qk = plt.quiverkey(Q,1.3,-0.2,U=qk_scale, label=qk_caption, coordinates='axes',labelpos='E') #10 m/s ,coordinates='data'   1.12,0.1 1.16,0.05  0.80,-0.20
		qk.text.set_backgroundcolor('White');qk.text.set_fontsize(20);
		qk.text.set_verticalalignment('center')
	return colormesh

		
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
		if variable == 'U' or variable == 'V': 
			data= data[:,23,:,:]#*24*60*60/1000  # from m/s to km/hr
		nc_fid.close()
		var40map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var40map
	
	####calculate the difference
	def dif_interp(vairable1,variable2,lon,lat):
		dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
		lon_interp = np.arange(0,360,5);lat_interp= np.arange(-90,90.1,5);
		f = interp2d(lon,lat,dif,kind='linear')
		dif_interp = f(lon_interp, lat_interp)
		return lon_interp, lat_interp,dif_interp
		
	lon,lat,Edg70GO = data_netcdf('Edg70GO',FREQ,variable)
	_,_,T1970 = data_netcdf('T1970RCP',FREQ,variable)
	_,_,EdgRef = data_netcdf('EdgRef',FREQ,variable)
	_,_,EdgEne = data_netcdf('EdgEne',FREQ,variable)
	_,_,EdgTech = data_netcdf('EdgTech',FREQ,variable)

	lon_in,lat_in,AEROSOL = dif_interp(Edg70GO,T1970,lon,lat)
	lon_in,lat_in,ENE= dif_interp(EdgRef,EdgEne,lon,lat)
	lon_in,lat_in,TECH= dif_interp(EdgRef,EdgTech,lon,lat)
	
	return lon_in,lat_in,AEROSOL,ENE,TECH

def plot_zonal_mean_uncertainty(ax,x,y,xlim):
	ax.tick_params(axis='both', which='major', direction = 'in',left=True,right=True,bottom=True,top=True, pad=5)
	ax.axvline(x=0,color='k',linewidth = 1);ax.axhline(y=0,color='k',linewidth = 1);
	ax.set_ylim([-90,90]);ax.set_yticks(np.arange(-90,91,30));ax.set_yticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	ax.set_xlim([-xlim,xlim]);ax.set_xticks(np.arange(-xlim,xlim+.0000001,xlim/2.0))
	ax.plot(x,y,'-',color="k",linewidth=5)
	return ax
	
index ='PSL';lon,lat,AA_PSL,ENE_PSL,TECH_PSL = data_readin(index,FREQ='mon');
index ='U';lon,lat,AA_U,ENE_U,TECH_U = data_readin(index,FREQ='mon'); 
index ='V';lon,lat,AA_V,ENE_V,TECH_V = data_readin(index,FREQ='mon'); 


fig = plt.figure(facecolor='White',figsize=[13.5,16.5]);pad= 5;	
colorbar_min=-100;colorbar_max=100;colormap='RdBu_r';

ax = plt.subplot2grid((3, 4), (0, 0), colspan=3)
ax.annotate('Sea level pressure '+r'($\mathrm{\mathsf{\Delta}\/Pa}$)' + ' + Wind850 '+r'($\mathrm{\mathsf{\Delta}\/m\/s^{-1}}$)',xy=(0.4,1.05), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)
ax.annotate('BEoA (2010 - 1970)',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)
ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,AA_PSL,AA_U,AA_V,qk_scale=None,qk_caption=None,qk_is = False)
ax = plt.subplot2grid((3, 4), (0, 3), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(AA_PSL,axis=1),37),lat,xlim=100)

ax = plt.subplot2grid((3, 4), (1, 0), colspan=3)
ax.annotate('Energy Consumption',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)
ax.annotate('(b)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,ENE_PSL,AA_U,ENE_V,qk_scale=None,qk_caption=None,qk_is = False)
ax = plt.subplot2grid((3, 4), (1, 3), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(ENE_PSL,axis=1),37),lat,xlim=100)

ax = plt.subplot2grid((3, 4), (2, 0), colspan=3)
ax.annotate('Technology Advancements',xy=(-0.1,0.5), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='center',rotation='vertical',fontsize=20)
ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='left', va='baseline',rotation='horizontal',fontsize=20)	
qk_caption = r'$\mathrm{\mathsf{1\/m\/s^{-1}}}$'; qk_scale=1
colormesh1 = spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,TECH_PSL,TECH_U,TECH_V,qk_scale,qk_caption=qk_caption,qk_is = True)
ax = plt.subplot2grid((3, 4), (2, 3), colspan=1)
plot_zonal_mean_uncertainty(ax,np.reshape(np.nanmean(TECH_PSL,axis=1),37),lat,xlim=100)

cbar_ax = fig.add_axes([0.10, 0.03, 0.75, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate('Pa',xy=(1.10,-1.1), xytext=(0, pad),
					xycoords='axes fraction', textcoords='offset points',
					ha='center', va='bottom',rotation='horizontal',fontsize=15)
plt.subplots_adjust(left=0.04, bottom=0.10, right=0.97, top=0.95, wspace=0.09, hspace=0.15); 
# plt.show()
plt.savefig('SLP_Wind850.png', format='png', dpi=1000)	
		
				
