# -*- coding: utf-8 -*-
'''
This is to plot the boxplots of HW and CS characteristics over 7 different countries (regions)
His(1986-2005), rcp(2081-2100) and FiXA (2081-2100) are compared
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d
from scipy.stats import mannwhitneyu as man_test
from scipy.stats import ttest_ind as student_test

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

def spatial_figure(axs,data,p_value,lons,lats,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False): #c_bad,c_under,c_over,c_number=20,
	from matplotlib import colors as mcolors
	from matplotlib import cm as cm
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
	# s = map.pcolor(xi, yi, data)
	
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines(); #map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	N = len(levels)-1
	cmap = discrete_cmap(N,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_over('k'); cmap.set_under('m');
	norm = mcolors.BoundaryNorm(levels, ncolors=N, clip=False)
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,norm=norm,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='/', alpha=0.,lw=0.9,latlon=True)
	axs.set_ylim([-60, 90])
	return colormesh

def HWI_deff_sig_pp(variable):
	def date_readin(scenario,directory):
		file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp'+directory  #_InterOsi
		nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_'+scenario+'.nc',mode='r')
		lon = nc_fid.variables['lon'][:]
		lat = nc_fid.variables['lat'][:]
		year = nc_fid.variables['year'][:]
		data = nc_fid.variables[variable][1,-21:-1,:,:]   # 0-mean 1-median
		nc_fid.close()
		return lon,lat,data
	
	def mannwhitneyu_test(vairable1,variable2):
		p_threshold=0.01
		size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
		p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
		# from scipy.stats import ttest_ind as test
		for x in range(size[0]):
			for y in range(size[1]):
				cache1 = vairable1[:,x,y]
				cache2 = variable2[:,x,y]
				if np.array_equal(cache1,cache2):p_value[x,y] = np.nan
				else:_,p_value[x,y] = student_test(cache1,cache2);	 # 		 man_test
		p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
		return p_value
		
	####calculate the difference and their significance
	def diff_sig(vairable1,variable2):
		dif = np.nanmean(vairable1-variable2,axis=0)
		sig = mannwhitneyu_test(vairable1, variable2)
		# sig_dif = np.multiply(dif,sig)
		# print np.shape(dif),np.shape(sig)
		return dif,sig
		
	lon,lat,his = date_readin('his','/CalDayThr/')
	lon,lat,rcp = date_readin('rcp85','/CalDayThr/')
	lon,lat,fix = date_readin('fixa','/CalDayThr/')
	
	GHG_D,GHG_S = diff_sig(fix,his)
	AAS_D,AAS_S = diff_sig(rcp,fix)
	
	return  lon,lat,GHG_D,GHG_S,AAS_D,AAS_S


fig = plt.figure(facecolor='White',figsize=[15,12]);plot_setup();pad= 5;colormap='Spectral_r';

### INTENSITY
lon,lat,GHG_D,GHG_S,AAS_D,AAS_S = HWI_deff_sig_pp(variable='HWInX')
colorbar_min=0;colorbar_max=4;
levels= [0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.5,2.0,2.5,3.0,4]
ax = plt.subplot(4,2,1);
ax.annotate('Annual peak intensity\n($^\circ$C)',xy=(-0.08, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('GHG increase',xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')				
ax.annotate('(a)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
spatial_figure(ax,GHG_D,GHG_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,2);
ax.annotate('Anthropogenic aerosol reduciton',xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')				
ax.annotate('(e)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
colormesh = spatial_figure(ax,AAS_D,AAS_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.92, 0.74, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',extend='both',cax=cbar_ax,ticks=levels,spacing ='uniform')  #ticks=np.round(np.arange(0,1.2,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)

## DURATION
lon,lat,GHG_D,GHG_S,AAS_D,AAS_S = HWI_deff_sig_pp(variable='HWDuX')
colorbar_min=0;colorbar_max=50;
levels= [0,1,2,3,4,5,6,8,10,12,14,16,18,20,30,50]
ax = plt.subplot(4,2,3);
ax.annotate('Anuual maximum duration\n'+r'($\mathrm{day\/event^{-1}}$)',xy=(-0.08, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')				
ax.annotate('(b)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
spatial_figure(ax,GHG_D,GHG_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,4);				
ax.annotate('(f)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
colormesh = spatial_figure(ax,AAS_D,AAS_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.92, 0.505, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',extend='both',cax=cbar_ax,ticks=levels,spacing ='uniform')  #ticks=np.round(np.arange(0,1.2,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)


## FREQUENCY
lon,lat,GHG_D,GHG_S,AAS_D,AAS_S = HWI_deff_sig_pp(variable='HWI_NO')
colorbar_min=-5;colorbar_max=30;
levels= [-5,0,1,2,3,4,5,6,7,8,9,10,12,16,20,30]
ax = plt.subplot(4,2,5);
ax.annotate('Anuual frequency\n'+r'($\mathrm{events\/year^{-1}}$)',xy=(-0.08, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')			
ax.annotate('(c)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
spatial_figure(ax,GHG_D,GHG_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,6);				
ax.annotate('(g)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
colormesh = spatial_figure(ax,AAS_D,AAS_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.92, 0.27, 0.01, 0.21])

char = fig.colorbar(colormesh,orientation='vertical',extend='both',cax=cbar_ax,ticks=levels,spacing ='uniform')  #ticks=np.round(np.arange(0,1.2,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)


## FREQUENCY
lon,lat,GHG_D,GHG_S,AAS_D,AAS_S = HWI_deff_sig_pp(variable='HWDC')
colorbar_min=0;colorbar_max=300;
levels= [0,5,10,15,20,25,30,35,40,60,80,100,150,200,300]
ax = plt.subplot(4,2,7);
ax.annotate('Total hot days\n'+r'($\mathrm{days\/year^{-1}}$)',xy=(-0.08, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')			
ax.annotate('(d)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
spatial_figure(ax,GHG_D,GHG_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,2,8);				
ax.annotate('(h)',xy=(0.02,0.87), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
colormesh = spatial_figure(ax,AAS_D,AAS_S,lon,lat,colormap,colorbar_min,colorbar_max,levels,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.92, 0.035, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',extend='both',cax=cbar_ax,ticks=levels,spacing ='uniform')  #ticks=np.round(np.arange(0,1.2,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)


plt.subplots_adjust(left=0.08, bottom=0.03, right=0.90, top=0.95, wspace=0.05, hspace=0.08); 
plt.savefig('figureS1spatial_In_Du_Fr_Td.png', format='png', dpi=1000)

