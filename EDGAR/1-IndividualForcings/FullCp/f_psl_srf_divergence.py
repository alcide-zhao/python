"""
Calculate the 850hpa wind divergence and plot it as shades with surface level pressure overlapped as contours
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d  as interp1d
# from scipy.stats import mannwhitneyu as man_test
# from scipy.stats import ttest_ind as student_test
# from scipy import integrate
# from climlab import constants as const

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


	
def data_read_annual_mean(variable,levs):	
	"""
	calculate the earth radius in m2
	"""	
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
		
		
	def mon_mean2annual_mean(scenario,time,data,levs):
		
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
			# print iyear,layer_b
			if levs:  # 4d DATA
				annual_mean=np.empty((40,30,192,288));annual_mean[:]=np.nan
				data_cache = data[layer_b:layer_e+1,:,:,:]
				annual_mean[iyear-year_series[0],:,:,:] = stats.nanmean(data_cache,axis=0)
			else:
				annual_mean=np.empty((40,192,288));annual_mean[:]=np.nan
				data_cache = data[layer_b:layer_e+1,:,:]
				annual_mean[iyear-year_series[0],:,:] = stats.nanmean(data_cache,axis=0)
		annual_mean = np.nanmean(annual_mean,axis=0)
		return annual_mean

	def data_netcdf(scenario,FREQ,variable,levs):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]# 850hpa
		VarAnMean = mon_mean2annual_mean(scenario,time,data,levs)
		if levs: lev = nc_fid.variables['lev'][:]
		else : lev =np.nan	
		nc_fid.close()
		return lev,lat,lon,VarAnMean
		
	FREQ = 'mon'
	lev,lat,lon,Edg70GO = data_netcdf('Edg70GO',FREQ,variable,levs)
	lev,lat,lon,T1970 = data_netcdf('T1970RCP',FREQ,variable,levs)
	lev,lat,lon,EdgRef = data_netcdf('EdgRef',FREQ,variable,levs)
	lev,lat,lon,Edg70Oz = data_netcdf('Edg70Oz',FREQ,variable,levs)
	#lev,lat,lon,EdgEne = data_netcdf('EdgEne',FREQ,variable,levs)
	#lev,lat,lon,EdgTech = data_netcdf('EdgTech',FREQ,variable,levs)
	return lev,lat,lon,T1970,Edg70GO,Edg70Oz,EdgRef  #,EdgEne,EdgTech

	
def hybrid2p_interp(PS,DataIn,midpoint=True,interface=False):
	"""
	This is to convert the CESM hybrid vertical coordinate to the pressure system
	datain must be an array of 3, 4, or 5 dimensions. Needs to contain a level dimension in hybrid coordinates. 
	The order of the dimensions is specific. The three rightmost dimensions must be 
	level x lat x lon [e.g. T(time,lev,lat,lon)]. 
	The order of the level dimension must be top-to-bottom
	"""
	
	PO = 1000   #surface reference pressure
	if midpoint:
		hya = np.array([ 0.00364346569404006, 0.00759481964632869, 0.0143566322512925,\
			0.0246122200042009, 0.0382682997733355, 0.0545954797416925,\
			0.0720124505460262, 0.0878212302923203, 0.103317126631737,\
			0.121547240763903, 0.142994038760662, 0.168225079774857,\
			0.178230673074722, 0.170324325561523, 0.161022908985615,\
			0.150080285966396, 0.137206859886646, 0.122061938047409,\
			0.104244712740183, 0.0849791541695595, 0.0665016956627369,\
			0.0501967892050743, 0.037188658490777, 0.028431948274374,\
			0.0222089774906635, 0.016407382208854, 0.0110745579004288,\
			0.00625495356507599, 0.00198940909467638, 0])
		hyb = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0196774136275053,\
			0.062504293397069, 0.112887907773256, 0.172161616384983,\
			0.241894043982029, 0.323930636048317, 0.420442461967468,\
			0.524799540638924, 0.624887734651566, 0.713207691907883,\
			0.783669710159302, 0.831102818250656, 0.864811271429062,\
			0.896237164735794, 0.92512384057045, 0.951230525970459,\
			0.974335998296738, 0.992556095123291])
	elif interface:
		hya = np.array([ 0.00225523952394724, 0.00503169186413288, 0.0101579474285245,\
			0.0185553170740604, 0.0306691229343414, 0.0458674766123295,\
			0.0633234828710556, 0.0807014182209969, 0.0949410423636436,\
			0.11169321089983, 0.131401270627975, 0.154586806893349,\
			0.181863352656364, 0.17459799349308, 0.166050657629967,\
			0.155995160341263, 0.14416541159153, 0.130248308181763,\
			0.113875567913055, 0.0946138575673103, 0.0753444507718086,\
			0.0576589405536652, 0.0427346378564835, 0.0316426791250706,\
			0.0252212174236774, 0.0191967375576496, 0.0136180268600583,\
			0.00853108894079924, 0.00397881818935275, 0, 0])
		hyb = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0393548272550106,\
			0.0856537595391273, 0.140122056007385, 0.204201176762581,\
			0.279586911201477, 0.368274360895157, 0.47261056303978,\
			0.576988518238068, 0.672786951065063, 0.753628432750702,\
			0.813710987567902, 0.848494648933411, 0.881127893924713,\
			0.911346435546875, 0.938901245594025, 0.963559806346893,\
			0.985112190246582, 1])
	Pnew = np.empty((len(hya),np.shape(PS)[0],np.shape(PS)[1])); Pnew[:] = np.nan
	
	for ilev in range(np.shape(hya)[0]):
		Pnew[ilev,:,:] = PO * hya[ilev]+ hyb[ilev] *PS/100
	
	p_interp=np.array([850])
	D850 = np.empty((len(lat),len(lon)));D850[:] = np.nan
	for ilat in range(len(lat)):
		for ilon in range(len(lon)):
			f = interp1d(Pnew[:,ilat,ilon], DataIn[:,ilat,ilon],bounds_error=False,fill_value = 0); 
			D850[ilat,ilon] = f(p_interp[:]); 
	return  D850     


def decompose_field(VAR,levs=False):
	if levs:
		variable='PS'; lev,lat,lon,PS_T1970,PS_Edg70GO,PS_Edg70Oz,PS_EdgRef = data_read_annual_mean(variable,levs=False)
		variable=VAR; lev,lat,lon,VAR_T1970,VAR_Edg70GO,VAR_Edg70Oz,VAR_EdgRef = data_read_annual_mean(variable,levs=True)
		T1970 = hybrid2p_interp(PS_T1970,VAR_T1970,midpoint=True,interface=False)
		Edg70GO = hybrid2p_interp(PS_Edg70GO,VAR_Edg70GO,midpoint=True,interface=False)
		Edg70Oz = hybrid2p_interp(PS_Edg70Oz,VAR_Edg70Oz,midpoint=True,interface=False)
		EdgRef = hybrid2p_interp(PS_EdgRef,VAR_EdgRef,midpoint=True,interface=False)
	else:
		variable='PSL'; lev,lat,lon,T1970,Edg70GO,Edg70Oz,EdgRef = data_read_annual_mean(variable,levs=False)
	total = EdgRef  - T1970
	GHG =   Edg70Oz - Edg70GO
	AAs =   Edg70GO - T1970
	O3 =    EdgRef  - Edg70Oz
	return lev,lat,lon,total,GHG,AAs,O3


def div_winds():
	def div_from_uv(U,V,lon,lat):
		def numerical_dif_2D(dependetnt,lon,lat,ax):
			if (ax == 1):
				lon_interval = 111.320*1000*np.cos(np.pi*lat/180)
				variable_gradient = np.gradient(lon)[1]*lon_interval
				dependetnt_gradient = np.gradient(dependetnt)[1]
			elif (ax == 0):
				variable_gradient = np.gradient(lat)[0]*110.574*1000
				dependetnt_gradient = np.gradient(dependetnt)[0]
			deriviative = np.divide(dependetnt_gradient,variable_gradient)
			return deriviative
		lonm, latm = np.meshgrid(lon, lat)
		dvdlat = numerical_dif_2D(V,lonm, latm,ax=0); 
		dudlon = numerical_dif_2D(U,lonm, latm,ax=1); 
		div = (dudlon+dvdlat)  # raise by a factor of 10exp(6)
		print np.nanmean(div,axis=(0,1))
		print np.nanmax(div,axis=(0,1))
		print np.nanmin(div,axis=(0,1))
		return div
	VAR = 'U';  lev,lat,lon,U850_total,U850_GHG,U850_AAs,U850_O3 = decompose_field(VAR,levs=True)
	VAR = 'V';  lev,lat,lon,V850_total,V850_GHG,V850_AAs,V850_O3 = decompose_field(VAR,levs=True)
	print np.nanmin(U850_total,axis=(0,1)),np.nanmax(U850_total,axis=(0,1))
	print np.nanmin(V850_total,axis=(0,1)),np.nanmax(V850_total,axis=(0,1))
	total = div_from_uv(U850_total,V850_total,lon,lat)
	GHG = div_from_uv(U850_GHG,V850_GHG,lon,lat)
	AAs = div_from_uv(U850_AAs,V850_AAs,lon,lat)
	O3 = div_from_uv(U850_O3,V850_O3,lon,lat)
	return 	lev,lat,lon,total,GHG,AAs,O3

	
	
def spatial_figure(axs,shade,contour,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): 
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
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
	
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	map.drawcoastlines(); #map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(shade), shade)
	cmap = discrete_cmap(16,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	CS1 = map.contourf(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)#
	clevs = np.arange(np.nanmin(contour[:]),np.nanmax(contour[:]),(np.nanmax(contour[:])-np.nanmin(contour[:]))/10)
	CS2 = map.contour(xi,yi,contour,clevs,linewidths=0.5,colors='k',animated=True)
	return CS1

VAR = 'PSL';  lev,lat,lon,PSL_total,PSL_GHG,PSL_AAs,PSL_O3 = decompose_field(VAR)
lev,lat,lon,div_total,div_GHG,div_AAs,div_O3 = div_winds()


fig = plt.figure(facecolor='White',figsize=[12,7]);plot_setup();pad= 5
colorbar_min=-10;colorbar_max=10;
colormap='RdBu';colormap = reverse_colourmap(colormap);


ax = plt.subplot(2,2,1);
ax.annotate('(a) 2010-1970',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,div_total,PSL_total,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,2);
ax.annotate('(b) GHGs',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,div_GHG,PSL_GHG,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,3);
ax.annotate('(c) AAs',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,div_AAs,PSL_AAs,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,4);
ax.annotate('(d) O3',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
colormesh1=spatial_figure(ax,div_O3,PSL_O3,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


cbar_ax = fig.add_axes([0.20, 0.03, 0.60, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.02, bottom=0.08, right=0.98, top=0.95, wspace=0.04, hspace=0.15); 
plt.savefig('./GHG_AEROSOL_OZONE/PSL_WINDdiv_Decompose.png', format='png', dpi=1000)










