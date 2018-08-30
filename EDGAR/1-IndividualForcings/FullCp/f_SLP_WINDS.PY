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
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

def spatial_scaler_vector(ax,lons,lats,colormap,colorbar_min,colorbar_max,scaler,vector_u,vector_v,qk_scale=None,qk_caption=None,qk_is = False):
	speed = np.sqrt(np.power(vector_u,2)+np.power(vector_v,2))
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
		qk = plt.quiverkey(Q,1.3,0.5,U=qk_scale, label=qk_caption, coordinates='axes',labelpos='E') #10 m/s ,coordinates='data'   1.12,0.1 1.16,0.05  0.80,-0.20
		qk.text.set_backgroundcolor('White');qk.text.set_fontsize(20);
		qk.text.set_verticalalignment('center')
	return colormesh

		
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
		for iyear in year_series[1:-1]:
			if (iyear == year_series[0] and time[0]//100 >= year_series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+12][0]  
			if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100  == (1+iyear)*100+2][0]
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
		
	lon,lat,Edg70GO = data_netcdf('Edg70GO',variable)
	_,_,T1970 = data_netcdf('T1970RCP',variable)
	_,_,EdgRef = data_netcdf('EdgRef',variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',variable)
	_,_,Edg70T10SOZ = data_netcdf('Edg70T10SOZ',variable)

	lon_in,lat_in,AEROSOL = dif_interp(Edg70GO,T1970,lon,lat)
	_,_,Total = dif_interp(EdgRef,T1970,lon,lat)
	_,_,GHG = dif_interp(Edg70Oz,Edg70GO,lon,lat)
	_,_,TrO3 = dif_interp(EdgRef, Edg70T10SOZ,lon,lat)
	_,_,StO3 = dif_interp(Edg70T10SOZ, Edg70Oz,lon,lat)
	
	return lon_in,lat_in,Total,GHG,AEROSOL,TrO3,StO3

def plot_zonal_mean_uncertainty(ax,x,y,xlim):
	ax.tick_params(axis='both', which='major', direction = 'in',left=True,right=True,bottom=True,top=True, pad=5)
	ax.axvline(x=0,color='k',linewidth = 1);ax.axhline(y=0,color='k',linewidth = 1);
	ax.set_ylim([-90,90]);ax.set_yticks(np.arange(-90,91,30));ax.set_yticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	ax.set_xlim([-xlim,xlim]);ax.set_xticks(np.arange(-xlim,xlim+.0000001,xlim/2.0))
	ax.plot(x,y,'-',color="k",linewidth=5)
	return ax
	
index ='PSL';lon,lat,Total_PSL,GHG_PSL,AEROSOL_PSL,TrO3_PSL,StO3_PSL = data_readin(index);
index ='U';lon,lat,Total_U,GHG_U,AEROSOL_U,TrO3_U,StO3_U = data_readin(index); 
index ='V';lon,lat,Total_V,GHG_V,AEROSOL_V,TrO3_V,StO3_V = data_readin(index); 


fig = plt.figure(facecolor='White',figsize=[12.5,11.5]);pad= 5
colorbar_min=-300;colorbar_max=300;colormap='RdBu_r';

ax = plt.subplot(3,2,1);
ax.annotate('(a) All',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,Total_PSL,Total_U,Total_V,qk_scale=None,qk_caption=None,qk_is = False)


ax = plt.subplot(3,2,2);
ax.annotate('(b) GHGs',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,GHG_PSL,GHG_U,GHG_V,qk_scale=None,qk_caption=None,qk_is = False)


ax = plt.subplot(3,2,3);
ax.annotate('(c) AAs',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,AEROSOL_PSL,AEROSOL_U,AEROSOL_V,qk_scale=None,qk_caption=None,qk_is = False)



ax = plt.subplot(3,2,4);
ax.annotate('(d) Trop. O3',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,TrO3_PSL,TrO3_U,TrO3_V,qk_scale=None,qk_caption=None,qk_is = False)
# plt.show()

ax = plt.subplot(3,2,5);
ax.annotate('(e) Strat. O3',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
qk_caption = r'$\mathrm{\mathsf{1\/m\/s^{-1}}}$'; qk_scale=1
colormesh1 = spatial_scaler_vector(ax,lon,lat,colormap,colorbar_min,colorbar_max,StO3_PSL,StO3_U,StO3_V,qk_scale=qk_scale,qk_caption=qk_caption,qk_is = True)

cbar_ax = fig.add_axes([0.06, 0.04, 0.74, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',
	extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate('Pa',xy=(1.10,-1.1), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)
# ax = plt.subplot(3,2,6);
# ax.annotate('(f) Zonal mean OF SLP (Pa)',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=15)	
# plot_zonal_mean_uncertainty(ax,lat,Total_PSL,GHG_PSL,AEROSOL_PSL,TrO3_PSL,StO3_PSL,ylim_l=-100,ylim_u=100)
# legend = ax.legend(shadow=False,ncol=2,loc ='upper left')	 
# legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)


plt.subplots_adjust(left=0.04, bottom=0.10, right=0.97, top=0.95, wspace=0.09, hspace=0.15); 
plt.savefig('SLP_Wind850_DJF.png', format='png', dpi=1000)
# plt.show()	
		
				
