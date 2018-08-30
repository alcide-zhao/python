"""
This is to process the monthly fields into zonal mean and plot 
	currently, SAT ans precipitation
	the model internal variability (25-t5 th percentile) and mean are both plotted
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
from scipy import integrate
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

# from climlab import constants as const

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

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
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh
			
def data_read(variable):	
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

	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
		var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		if variable =="VQ" or variable == "VT":
			data = np.nanmean(nc_fid.variables[variable][:,23:30,:,:],axis=1)  # 850hpa
		else:
			data = nc_fid.variables[variable][:]#-273.15
		nc_fid.close()
		var40map = mon_mean2annual_mean(scenario,time,data)
		return lon,lat,var40map
	
	lon,lat,Edg70GO = data_netcdf('Edg70GO',variable)
	_,_,T1970 = data_netcdf('T1970RCP',variable)
	_,_,EdgRef = data_netcdf('EdgRef',variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',variable)
	_,_,EdgEne = data_netcdf('EdgEne',variable)
	_,_,EdgTech = data_netcdf('EdgTech',variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	
def spa_dif_sig(variable1,variable2):
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

	def decompose_forcings(Ref970,Edg70GO,EdgEne,EdgRef,EdgTech):
		def diff_sig(vairable1,variable2):
			dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
			sig = mannwhitneyu_test(vairable1, variable2,p_threshold=0.1)
			# sig_dif = np.multiply(dif,sig)
			return dif,sig
		BEoA,BEoA_s = diff_sig(Edg70GO,Ref970)
		Ene,Ene_s = diff_sig(EdgRef,EdgEne)
		Tech,Tech_s = diff_sig(EdgRef,EdgTech)
		return np.stack((BEoA,Ene,Tech,BEoA_s,Ene_s,Tech_s),axis=0)
		
	if variable2 =='none':
		if variable1 == 'albedo':
			index = 'FSNS'; lon,lat,T1970_N,Edg70GO_N,Edg70Oz_N,EdgRef_N,EdgEne_N,EdgTech_N = data_read(index);
			index = 'FSDS'; lon,lat,T1970_D,Edg70GO_D,Edg70Oz_D,EdgRef_D,EdgEne_D,EdgTech_D = data_read(index);
			Ref970 = 100*np.divide(-(T1970_N-T1970_D),T1970_D);
			Edg70GO = 100*np.divide(-(Edg70GO_N-Edg70GO_D),Edg70GO_D);
			EdgRef = 100*np.divide(-(EdgRef_N-EdgRef_D),EdgRef_D);
			EdgEne = 100*np.divide(-(EdgEne_N-EdgEne_D),EdgEne_D); 
			EdgTech = 100*np.divide(-(EdgTech_N-EdgTech_D),EdgTech_D);			
		else:
			index = variable1;lon,lat,Ref970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_read(index);		
	else:
		index = variable1; lon,lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_read(index);
		index = variable2; lon,lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_read(index);
		Ref970 = T1970_S-T1970_l;
		Edg70GO = Edg70GO_S-Edg70GO_l;
		EdgRef = EdgRef_S-EdgRef_l;
		EdgEne = EdgEne_S-EdgEne_l; 
		EdgTech = EdgTech_S-EdgTech_l;		
	dics = decompose_forcings(Ref970,Edg70GO,EdgEne,EdgRef,EdgTech)
	return lon,lat,dics


fig = plt.figure(facecolor='White',figsize=[11,9]);pad= 5;
colorbar_min=-5;colorbar_max=5;colormap='RdBu_r'; 

#####################################################################
####  spatial patterns of Clear/scloudy sky radiFlux at TOA/sfc  ####
#####################################################################
"""
variable1 = 'FSNTOAC';variable2 = 'none';lon,lat,dics  = spa_dif_sig(variable1,variable2)
ax = plt.subplot(3,2,1);
ax.annotate('Best\nEstimation',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('TOA clear-sky',xy=(0.5,1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)			
colormesh0 = spatial_figure(ax,dics[0,:,:],dics[3,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,3);
ax.annotate('Energy\nConsumption',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)			
colormesh0 = spatial_figure(ax,dics[1,:,:],dics[4,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,5);
ax.annotate('Technology\nAdvancements',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=15)	
ax.annotate('(e)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
colormesh0 = spatial_figure(ax,dics[2,:,:],dics[5,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


variable1 = 'FSNSC';variable2 = 'none';lon,lat,dics  = spa_dif_sig(variable1,variable2)

ax = plt.subplot(3,2,2);
ax.annotate('(b)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.annotate('Surface clear-sky',xy=(0.5,1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)			
colormesh0 = spatial_figure(ax,dics[0,:,:],dics[3,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,4);
ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)			
colormesh1 = spatial_figure(ax,dics[1,:,:],dics[4,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,6);
ax.annotate('(f)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
colormesh1 = spatial_figure(ax,dics[2,:,:],dics[5,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


# variable1 = 'FSNS';variable2 = 'FSNSC';lon,lat,dics  = spa_dif_sig(variable1,variable2)

# ax = plt.subplot(3,3,3);
# ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=15)	
# ax.annotate('Cloud forcing',xy=(0.5,1.1), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='baseline',rotation='horizontal',fontsize=15)			
# colormesh0 = spatial_figure(ax,dics[0,:,:],dics[3,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

# ax = plt.subplot(3,3,6);
# ax.annotate('(f)',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=15)			
# colormesh0 = spatial_figure(ax,dics[1,:,:],dics[4,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

# ax = plt.subplot(3,3,9);
# ax.annotate('(i)',xy=(0.02,1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=15)				
# colormesh1 = spatial_figure(ax,dics[2,:,:],dics[5,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.10, 0.03, 0.75, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/W\/m^{-2}}$',xy=(1.10,-1.1), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)


plt.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.92, wspace=0.1, hspace=0.15); 
plt.savefig('clearsky_TOA_SFC.png', format='png', dpi=1000)


"""
#####################################################################
####                          surface albedo                     ####
#####################################################################
	
def plot_zonal_mean_uncertainty(ax,x1,y,xlim):
	ax.tick_params(axis='both', which='major', direction = 'in',left=True,right=True,bottom=True,top=True, pad=5)
	ax.axvline(x=0,color='k',linewidth = 1);ax.axhline(y=0,color='k',linewidth = 1);
	ax.set_ylim([-90,90]);ax.set_yticks(np.arange(-90,91,30));ax.set_yticklabels(('90S','60S','30S','EQ','30N','60N','90N'));
	ax.set_xlim([-xlim,xlim]);ax.set_xticks(np.arange(-xlim,xlim+.0000001,xlim/1.0))
	ax.plot(x1,y,'-',color="k",linewidth=2)
	ax.yaxis.tick_right()
	return ax
	
fig = plt.figure(facecolor='White',figsize=[7.5,9.00]);pad= 5;
colorbar_min=-4;colorbar_max=4;colormap='RdBu_r'; xlim=4

variable1 = 'albedo';variable2 = 'none';lon,lat,dics  = spa_dif_sig(variable1,variable2)
ax = plt.subplot2grid((3, 4), (0, 0), colspan=3)
ax.annotate('(a) Best Estimation',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)			
colormesh0 = spatial_figure(ax,dics[0,:,:],dics[3,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot2grid((3, 4), (0, 3), colspan=1)
x1 = np.reshape(np.nanmean(dics[0,:,:],axis=1),192)
plot_zonal_mean_uncertainty(ax,x1,lat,xlim)

ax = plt.subplot2grid((3, 4), (1, 0), colspan=3)
ax.annotate('(b) Energy Consumption',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)			
colormesh0 = spatial_figure(ax,dics[1,:,:],dics[4,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot2grid((3, 4), (1, 3), colspan=1)
x1 = np.reshape(np.nanmean(dics[1,:,:],axis=1),192)
plot_zonal_mean_uncertainty(ax,x1,lat,xlim)


ax = plt.subplot2grid((3, 4), (2, 0), colspan=3)
ax.annotate('(c) Technology Advancements',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
colormesh1 = spatial_figure(ax,dics[2,:,:],dics[5,:,:],lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot2grid((3, 4), (2, 3), colspan=1)
x1 = np.reshape(np.nanmean(dics[2,:,:],axis=1),192)
plot_zonal_mean_uncertainty(ax,x1,lat,xlim)


cbar_ax = fig.add_axes([0.10, 0.03, 0.75, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate('%',xy=(1.10,-1.1), xytext=(0, pad),
				xycoords='axes fraction', textcoords='offset points',
				ha='center', va='bottom',rotation='horizontal',fontsize=15)

plt.subplots_adjust(left=0.04, bottom=0.08, right=0.92, top=0.96, wspace=0.1, hspace=0.15); 
plt.savefig('surface_albedo.png', format='png', dpi=1000)

