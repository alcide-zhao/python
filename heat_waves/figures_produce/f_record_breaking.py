# -*- coding: utf-8 -*-
'''
Global mao of HW and CS with magnitude greater than 4(8)
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d
import math

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


scenario_dic ={'his':[21,0,45],'rcp85':[14,0,95],'fixa':[14,0,95]}

def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,p_value,p_value1,tb_lef=True,tb_bot=True ): #c_bad,c_under,c_over,c_number=20,
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
	# lon_bin = int((lon_e -lon_b)/5)
	# lat_bin = int((lat_e -lat_b)/5)
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, data)
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='.',color='k',zorder=5)
	sc2 = map.scatter(xi[p_value1==1], yi[p_value1==1],p_value1[p_value1==1],marker='.',color='k',zorder=5)
	# Add Grid Lines
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines();# map.drawstates(); map.drawcountries()
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_under('w')#cmap.set_over('r')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min-0.01, vmax=colorbar_max,latlon=True) #
	ax.set_ylim([-60,90])
	return colormesh

def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# plt.imshow(area);plt.show()
	return area

def mask_match(country_key,lon,lat):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask= ocean_mask[country_key][:]
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	# print np.unique(mask)
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
	mask = f(lon, lat);
	mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid (lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
	mask=np.multiply(mask,area);  #crop the interested region
	mask=np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
	return mask

def get_baseline():
	scenario = 'his'; #'fixa','his'
	ensemble = scenario_dic[scenario][0]; layer_s = 15;layer_e = scenario_dic[scenario][2]  #scenario_dic[scenario][1]; 
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'+scenario+'/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + '/file_list.txt', "r")
	text_content = text_file.readlines()
	HWI_max_his = np.empty((ensemble,layer_e-layer_s,192,288))
	for en_no in range(0,ensemble):
		nc_f = text_content[en_no][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		HWI_max_his[en_no,:,:,:] = nc_fid.variables['HWI'][layer_s:layer_e,:,:,2]
	HWI_max_his=np.reshape(HWI_max_his,[ensemble*(layer_e-layer_s),192,288])
	# mean=np.nanmean(HWI_max_his,axis=0)
	# std=np.nanstd(HWI_max_his,axis=0)
	# HisHWImax = mean +3*std
	
	scenario = 'rcp85'; #'fixa','his'
	ensemble = scenario_dic[scenario][0]; layer_s = 0;layer_e = 10  #scenario_dic[scenario][1]; 
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'+scenario+'/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + '/file_list.txt', "r")
	text_content = text_file.readlines()
	HWI_max_rcp = np.empty((ensemble,layer_e-layer_s,192,288))
	for en_no in range(0,ensemble):
		nc_f = text_content[en_no][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		HWI_max_rcp[en_no,:,:,:] = nc_fid.variables['HWI'][layer_s:layer_e,:,:,2]
	HWI_max_rcp=np.reshape(HWI_max_rcp,[ensemble*(layer_e-layer_s),192,288])
	hwi_cache= np.concatenate((HWI_max_his,HWI_max_rcp),axis=0)
	print np.shape(hwi_cache)
	mean=np.nanmean(hwi_cache,axis=0)
	std=np.nanstd(hwi_cache,axis=0)
	HisHWImax = mean +5*std
	HisHWImax = np.nanmax(hwi_cache,axis=0)
	return HisHWImax

def get_RecordBreak(scenario,year_s,year_e):
	ensemble = scenario_dic[scenario][0]; layer_s = year_s-2006; layer_e = year_e-2006
	HisHWImax = get_baseline()
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'+scenario+'/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + '/file_list.txt', "r")
	text_content = text_file.readlines()
	RecordBreakNo=np.zeros((192,288))
	for en_no in range(0,ensemble):
		nc_f = text_content[en_no][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		lon= nc_fid.variables['lon'][:];lat=nc_fid.variables['lat'][:];
		cache = nc_fid.variables['HWI'][layer_s:layer_e,:,:,2]-HisHWImax
		cache[cache>0]=1;cache[cache<=0]=np.nan;
		RecordBreakNo = RecordBreakNo+ np.nansum(cache,axis=0)
	RecordBreakRatio =100.0* RecordBreakNo/(year_e-year_s)/ensemble
	RecordBreakRatio[RecordBreakRatio==0]=np.nan
	return lon,lat,RecordBreakRatio

# lon,lat,rcp3150 = get_RecordBreak('rcp85',2031,2051)
# lon,lat,fix3150 = get_RecordBreak('fixa',2031,2051)
# AA3150=rcp3150-fix3150
lon,lat,rcp8100 = get_RecordBreak('rcp85',2051,2101)
lon,lat,fix8100 = get_RecordBreak('fixa',2051,2101)
AA8100=rcp8100-fix8100

country_list=['Globe','Australia','Europe','China','USA','India','Brazil','Southern African']
for country in country_list:
	mask = mask_match(country,lon,lat)
	print country
	print np.nansum(np.nansum(np.multiply(fix8100,mask),axis=1),axis=0),np.nansum(np.nansum(np.multiply(rcp8100,mask),axis=1),axis=0),\
		np.nansum(np.nansum(np.multiply(AA8100,mask),axis=1),axis=0)
	
fig = plt.figure(facecolor='White',figsize=[12,15]);plot_setup();pad= 5;
colormap='RdYlBu';colormap = reverse_colourmap(colormap)
# colormap='Reds';#colormap = reverse_colourmap(colormap)
ax = plt.subplot(3,1,1);colorbar_min=0;colorbar_max=100;
p_value = np.zeros((192,288)); p_value[rcp8100>95]=1;
p_value1 = np.zeros((192,288)); p_value1[rcp8100>80]=1;
ax.annotate('',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(a) RCP8.5',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,rcp8100,lon,lat,colormap,colorbar_min,colorbar_max,p_value,p_value1,tb_lef=False,tb_bot=False)


ax = plt.subplot(3,1,2);
p_value = np.zeros((192,288)); p_value[fix8100>80]=1;
ax.annotate('(b) RCP8.5_FixA',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
colormesh1 = spatial_figure(ax,fix8100,lon,lat,colormap,colorbar_min,colorbar_max,p_value,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.93, 0.4, 0.015, 0.52])
char = fig.colorbar(colormesh1,orientation='vertical',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.set_ylabel('%', rotation=90)
ax = plt.subplot(3,1,3);colorbar_min=0;colorbar_max=30;
p_value = np.zeros((192,288)); p_value[fix8100>100]=1;
ax.annotate('(c) RCP8.5 - RCP8.5_FixA',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
colormesh1 = spatial_figure(ax,AA8100,lon,lat,colormap,colorbar_min,colorbar_max,p_value,p_value,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.93, 0.06, 0.015, 0.25])
char = fig.colorbar(colormesh1,orientation='vertical',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.ax.set_ylabel('%', rotation=90)
					
plt.subplots_adjust(left=0.02, bottom=0.05, right=0.96, top=0.95, wspace=0.04, hspace=0.2); 
plt.savefig('Fig4_HW_RecordBreaking.png', format='png', dpi=1000)
