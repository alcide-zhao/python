import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import site
import os
import numpy as np
from scipy import stats
import scipy
from scipy.interpolate import interp2d
from lib import discrete_cmap
import netCDF4 as nc4
from matplotlib import rcParams
import scipy.io as spio
import numpy.ma as ma


#######################################################
# 1. time series                                      #
#######################################################
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/time_series/ocean_masked/'
file_name = 'time_series_Precip_extremes_global_2006_2100.nc'
time_s = 2006
time_e = 2100

nc_fid = nc4.Dataset(input_path+file_name,mode='r')
sub_plot=1
plt.figure(1, facecolor='White',figsize=[20,16])
time_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
# ['rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot']:
for att_name in [ 'rx1day']:
	time_series_range = np.empty([3,time_e-time_s+1])
	time = nc_fid.variables['time'][:]
	att_value = nc_fid.variables[att_name][:]
	time_series_range[0,:] = np.nanmax(att_value,axis = 0)
	time_series_range[2,:] = np.nanmin(att_value,axis = 0)	
	time_series_range[1,:] = stats.nanmean(att_value,axis = 0)
	time_dic[att_name] = time_series_range
	ax1=plt.subplot(1,1,sub_plot)
	ax1.plot(time,time_series_range[1,:],'-',color="red",linewidth=2,label='RCP8.5_FixA')
	ax1.fill_between(time,time_series_range[0,:],time_series_range[2,:], where=time_series_range[0,:]>=time_series_range[1,:], color="red", alpha = 0.2) 
	sub_plot = sub_plot+1

plt.hold(True)
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/time_series/ocean_masked/'
file_name = 'time_series_Precip_extremes_global_2006_2080.nc'
time_s = 2006
time_e = 2080

nc_fid = nc4.Dataset(input_path+file_name,mode='r')
sub_plot=1
plt.figure(1, facecolor='White',figsize=[20,16])
time_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
for att_name in [ 'rx1day']:
	time_series_range = np.empty([3,time_e-time_s+1])
	time = nc_fid.variables['time'][:]
	att_value = nc_fid.variables[att_name][:]
	time_series_range[0,:] = np.nanmax(att_value,axis = 0)
	time_series_range[2,:] = np.nanmin(att_value,axis = 0)	
	time_series_range[1,:] = stats.nanmean(att_value,axis = 0)
	time_dic[att_name] = time_series_range
	ax1=plt.subplot(1,1,sub_plot)
	ax1.plot(time,time_series_range[1,:],'-',color="blue",linewidth=2,label='RCP8.5')
	ax1.fill_between(time,time_series_range[0,:],time_series_range[2,:], where=time_series_range[0,:]>=time_series_range[1,:], color="blue", alpha = 0.2) 
	sub_plot = sub_plot+1

plt.hold(True)
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/time_series/ocean_masked/'
file_name = 'time_series_Precip_extremes_global_2081_2100.nc'
time_s = 2081
time_e = 2100

nc_fid = nc4.Dataset(input_path+file_name,mode='r')
sub_plot=1
plt.figure(1, facecolor='White',figsize=[20,4])
unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}
time_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
# ['rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot']:
for att_name in [ 'rx1day']:
	time_series_range = np.empty([3,time_e-time_s+1])
	time = nc_fid.variables['time'][:]
	att_value = nc_fid.variables[att_name][:]
	time_series_range[0,:] = np.nanmax(att_value,axis = 0)
	time_series_range[2,:] = np.nanmin(att_value,axis = 0)	
	time_series_range[1,:] = stats.nanmean(att_value,axis = 0)
	time_dic[att_name] = time_series_range
	ax1=plt.subplot(1,1,sub_plot)
	ax1.plot(time,time_series_range[1,:],'-',color="blue",linewidth=2)
	ax1.fill_between(time,time_series_range[0,:],time_series_range[2,:], where=time_series_range[0,:]>=time_series_range[1,:], color="blue", alpha = 0.2) 
	ax1.set_title(att_name.title(),fontsize=30)
	ax1.set_xlabel('Year',fontsize=30)
	ax1.set_xlim([2000,2100])
	ax1.set_ylim([4,7.5])
	ax1.set_ylabel(unit_dic[att_name].title(),fontsize=30)
	ax1.set(aspect=10)
	sub_plot = sub_plot+1
	ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	ax1.tick_params(labelsize=25,axis='both', which='major', pad=30)
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(-1.11, 0.94),fontsize=20)	

plt.tight_layout()
# plt.show()
plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/figures/rx1day_time_series_global.tif')
# #######################################################
# # 1. spatial features                                 #
# #######################################################


def range_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,data):
	"""
	clip the data based on given range
	"""
	lon = np.array(lon)
	lat = np.array(lat)
	colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
	colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
	row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
	row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
	lon_clipped = lon[colum_s:colum_e+1]
	lat_clipped = lat[row_s:row_e+1]
	data_clipped = data[:,row_s:row_e+1,colum_s:colum_e+1]
	return lon_clipped, lat_clipped, data_clipped
	
	
def spatial_figure(lons,lats,plot_title,colorbar_unit,objective,colorbar_min,colorbar_max):
	# calculate the origin of the map
	lon_0 = lons.mean()
	lat_0 = lats.mean()

	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=0.,llcrnrlat=-90.,urcrnrlon=360.,urcrnrlat=90.)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, objective)
	# Add Grid Lines
	map.drawparallels(np.arange(-90., 90., 60.), labels=[1,0,0,0],linewidth=0.0, fontsize=10)
	map.drawmeridians(np.arange(0., 360., 60.), labels=[0,0,0,1],linewidth=0.0, fontsize=10)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines()
	# map.drawstates()
	# map.drawcountries()
	# Add Colorbar
	cmap = discrete_cmap(20,'hsv')
	cmap.set_bad('w',alpha = 1.0)
	cmap.set_over('r')
	cmap.set_under('k')
	# cmap.set_over((1., 0., 0.))
	# cmap.set_under((0., 0., 1.))
	masked_obj = np.ma.masked_where(np.isnan(objective), objective)
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)
	cbar = map.colorbar(colormesh, location='bottom', pad="20%",extend='both')
	cbar.set_label(colorbar_unit,fontsize=10)
	# Add Title
	plt.title(plot_title)
	
	
oceanmask=spio.loadmat('//home/s1667168/coding/python_scripts/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/'
att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}

figure_no = 1
for att_name in ['r10']:
	colorbar_unit = unit_dic[att_name]
	plt.figure(figure_no, facecolor='White',figsize=[12,8])
	
	ax1=plt.subplot(3,3,1)
	
	file_name = 'ensumble_mean_Precip_extremes_global_2006_2080.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lats = nc_fid.variables['lat'][:]
	lons = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[att_name][:]
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	att_0615_R = stats.nanmean(att_value[0:9,:,:],axis=0)
	att_0615_R = np.multiply(att_0615_R,oceanmask)
	figure_title = att_name+'_2006-2015'+'_RCP8.5'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_0615_R,None,None)
	nc_fid.close()
	
	ax1=plt.subplot(3,3,2)
	file_name = 'ensumble_mean_Precip_extremes_global_2081_2100.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lats = nc_fid.variables['lat'][:]
	lons = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[att_name][:]
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	att_9000_R = stats.nanmean(att_value[10:19,:,:],axis=0)
	att_9000_R = np.multiply(att_9000_R,oceanmask)
	figure_title = att_name+'_2091-2100'+'_RCP8.5'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_9000_R,None,None)
	
	
	ax1=plt.subplot(3,3,3)
	att_diff = att_9000_R - att_0615_R
	att_diff = np.multiply(att_diff,oceanmask)
	figure_title = att_name+'_DIFF'+'_RCP8.5'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_diff,-8,9)
	nc_fid.close()
	
	file_name = 'ensumble_mean_Precip_extremes_global_2006_2100_fixA.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lats = nc_fid.variables['lat'][:]
	lons = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[att_name][:]
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	att_0615 = stats.nanmean(att_value[0:9,:,:],axis=0)
	att_9000 = stats.nanmean(att_value[85:94,:,:],axis=0)
	att_diff = att_9000 - att_0615
	att_0615 = np.multiply(att_0615,oceanmask)
	att_9000 = np.multiply(att_9000,oceanmask)
	att_diff = np.multiply(att_diff,oceanmask)
	nc_fid.close()
	
	ax1=plt.subplot(3,3,4)
	figure_title = att_name+'_2006-2015_'+'RCP8.5_FixA'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_0615,None,None)
	
	
	ax1=plt.subplot(3,3,5)
	figure_title = att_name+'_2091-2100_'+'RCP8.5_FixA'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_9000,None,None)
	
	
	ax1=plt.subplot(3,3,6)
	figure_title = att_name+'_DIFF_'+'RCP8.5_FixA'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_diff,-8,9)
	
	att_diff = att_0615_R -att_0615
	ax1=plt.subplot(3,3,7)
	figure_title = att_name+'_2006-2015_'+'RCP8.5-RCP8.5_FixA'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_diff,None,None)
	
	att_diff = att_9000_R -att_9000
	ax1=plt.subplot(3,3,9)
	figure_title = att_name+'_2091-2100_'+'-CP8.5-RCP8.5_FixA'
	spatial_figure(lons,lats,figure_title,colorbar_unit,att_9000,None,None)
	
	
	
	
	plt.show()
	# plt.savefig(input_path+att_name+'.tif')
	figure_no+ figure_no+1
	
