# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 2017
This is to produce the time evelution of 10-years averaged spatial precipitation extremes 
Here using the rx5day
@author: Alcide.Zhao
"""
import netCDF4 as nc4
import numpy as np
import scipy.io as spio
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import os; import site
lib_path = os.path.join(
	os.path.realpath(
		os.path.dirname(__file__)
	), 
	os.path.pardir, 
	os.path.pardir,
    # 'lib'
)
site.addsitedir(lib_path)
from lib import *

oceanmask=spio.loadmat('//home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}

region = rergion_dic['ASIA'][:]
#############################
# climatology references
#############################
# precipitation indices climatology : reference time period is 1971-2000 
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
value = np.multiply(nc_fid.variables['total_precip'][66:85,:,:],1)/92
lon_CESM =  nc_fid.variables['lon'][:]
lat_CESM =  nc_fid.variables['lat'][:]
nc_fid.close()
value_ref = stats.nanmean(value,axis=0)
# print rx5day_ref
# #######################################################
# # 1. spatial features                                 #
# #######################################################


att_name='total_precip'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
# series if the starting points of the time slices

file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

time = nc_fid.variables['time'][:]
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
att_value = nc_fid.variables[att_name][:]/92-value_ref
att_value = np.multiply(att_value,1)
lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
nc_fid.close()
_,_,ocean_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)


file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
# time = nc_fid.variables['time'][:]
# lat = nc_fid.variables['lat'][:]
# lon = nc_fid.variables['lon'][:]
att_value = nc_fid.variables[att_name][:]/92-value_ref
att_value = np.multiply(att_value,1)
lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
nc_fid.close()

######################
# mean of time slice
######################

att_0615_R =stats.nanmean(att_clipped_r85[0:10,:,:],axis=0)
att_0615_F =stats.nanmean(att_clipped_r85F[0:10,:,:],axis=0)
att_0615_D = att_0615_R -att_0615_F
att_3150_R = stats.nanmean(att_clipped_r85[25:45,:,:],axis=0)
att_3150_F = stats.nanmean(att_clipped_r85F[25:45,:,:],axis=0)
att_3150_D = att_3150_R -att_3150_F
att_8100_R = stats.nanmean(att_clipped_r85[75:95,:,:],axis=0)
att_8100_F = stats.nanmean(att_clipped_r85F[75:95,:,:],axis=0)
att_8100_D = att_8100_R -att_8100_F


# now plotting
# fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(7.9, 4), facecolor='White');plot_setup();
fig = plt.figure(facecolor='White',figsize=(7.9, 4));plot_setup();	
		
# rcp85 and fixa
colormap ='BrBG'; p_value = np.zeros((np.shape(att_0615_R))); colorbar_min=-1.5;  colorbar_max = 1.5;
pad= 5 
				
ax = plt.subplot(2,3,1)
spatial_figure_norm(ax,att_3150_R,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=False )
ax.annotate('2031-2050',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('RCP8.5',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
				
ax = plt.subplot(2,3,4)
ax.annotate('2081-2100',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
spatial_figure_norm(ax,att_8100_R,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=True )
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax=plt.subplot(2,3,2)
spatial_figure_norm(ax,att_3150_F,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=False,tb_bot=False )
ax.annotate('RCP8.5_FixA',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
	
ax=plt.subplot(2,3,5)
colormesh1 =spatial_figure_norm(ax,att_8100_F,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=False,tb_bot=True )
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
#   Rcp and Fixa share the same colorbar while the diff plots use its own colorbar
cbar_ax = fig.add_axes([0.12, 0.10, 0.55, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal', cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$')

# rcp85 minus fixa
# from scipy.stats import ttest_ind as test	
p_threshold=0.05
def mannwhitneyu_test(data1,data2):
	size = np.array([np.shape(data2)[1],np.shape(data2)[2]]); 
	p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
	from scipy.stats import mannwhitneyu as test
	for x in range(size[0]):
		for y in range(size[1]):
			cache1 = data1[:,x,y]
			cache2 = data2[:,x,y]
			_,p_value[x,y] = test(cache1,cache2);
			# print p_value[x,y]
	p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
	# p_value=np.multiply(ocean_mask,p_value)
	return p_value


p_value_thres = 0.075;colorbar_min=-1.5;  colorbar_max = 1.5 #BrBG seismic

ax=plt.subplot(2,3,3)
p_value=mannwhitneyu_test(att_clipped_r85[25:45,:,:], att_clipped_r85F[25:45,:,:])
# p_value[p_value>=p_value_thres] =np.nan; p_value[p_value<p_value_thres] =1
spatial_figure_norm(ax,att_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=False,tb_bot=False )
ax.annotate('RCP8.5-RCP8.5_FixA',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')	
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(2,3,6)
p_value=mannwhitneyu_test(att_clipped_r85[75:95,:,:], att_clipped_r85F[75:95,:,:])
# p_value[p_value>=p_value_thres] =np.nan; p_value[p_value<p_value_thres] =1
colormesh2 =spatial_figure_norm(ax,att_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=False,tb_bot=True )
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

cbar_ax = fig.add_axes([0.70, 0.10, 0.27, 0.02])
char = fig.colorbar(colormesh2, orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$')


plt.subplots_adjust(left=0.1, bottom=0.16, right=0.98, top=0.95, wspace=0.05, hspace=0.05);
plt.savefig('Fig8.png', format='png', dpi=1000)
plt.show()