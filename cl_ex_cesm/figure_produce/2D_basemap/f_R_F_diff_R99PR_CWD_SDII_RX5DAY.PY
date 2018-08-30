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

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}

region = rergion_dic['ASIA'][:]
oceanmask=spio.loadmat('//home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
nc_fid.close()
_,_,ocean_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)



# #######################################################
# # 1. spatial features                                 #
# #######################################################

def mannwhitneyu_test(data1,data2):
	p_threshold=0.05
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
	p_value=np.multiply(ocean_mask,p_value)
	return p_value
	
def spatial_diff_sig(region, variable):
	from scipy.stats import ttest_ind as test	
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	#fixa
	file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	#RCP8.5
	file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	
	att_0615_R =stats.nanmean(att_clipped_r85[0:10,:,:],axis=0)
	att_0615_F =stats.nanmean(att_clipped_r85F[0:10,:,:],axis=0)
	att_0615_D = att_0615_R -att_0615_F
	att_3150_R = stats.nanmean(att_clipped_r85[25:45,:,:],axis=0)
	att_3150_F = stats.nanmean(att_clipped_r85F[25:45,:,:],axis=0)
	att_3150_D = att_3150_R -att_3150_F
	att_8100_R = stats.nanmean(att_clipped_r85[75:95,:,:],axis=0)
	att_8100_F = stats.nanmean(att_clipped_r85F[75:95,:,:],axis=0)
	att_8100_D = att_8100_R -att_8100_F
	att_0615_P=mannwhitneyu_test(att_clipped_r85[0:10,:,:], att_clipped_r85F[0:10,:,:])
	att_3150_P=mannwhitneyu_test(att_clipped_r85[25:45,:,:], att_clipped_r85F[25:45,:,:])
	att_8100_P=mannwhitneyu_test(att_clipped_r85[75:95,:,:], att_clipped_r85F[75:95,:,:])
	
	return lons,lats,att_0615_D,att_3150_D,att_8100_D,att_0615_P,att_3150_P,att_8100_P


# INDICES REPRODUCING
lons,lats,cwd_0615_D,cwd_3150_D,cwd_8100_D,cwd_0615_P,cwd_3150_P,cwd_8100_P = spatial_diff_sig(region, variable = 'cwd')
_,_,sdii_0615_D,sdii_3150_D,sdii_8100_D,sdii_0615_P,sdii_3150_P,sdii_8100_P = spatial_diff_sig(region, variable = 'r10')
_,_,rx5day_0615_D,rx5day_3150_D,rx5day_8100_D,rx5day_0615_P,rx5day_3150_P,rx5day_8100_P = spatial_diff_sig(region, variable = 'rx5day')
lons,lats,r95p_0615_D,r95p_3150_D,r95p_8100_D,r95p_0615_P,r95p_3150_P,r95p_8100_P = spatial_diff_sig(region, variable = 'r95p')

# now plotting
"""
def spatial_diff_r95pr(region):
	from scipy.stats import ttest_ind as test	
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	#fixa
	file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = np.divide(nc_fid.variables['r95p'][:],nc_fid.variables['total_precip'][:])*100
	att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	#RCP8.5
	file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	att_value = np.divide(nc_fid.variables['r95p'][:],nc_fid.variables['total_precip'][:])*100
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	
	att_0615_R =stats.nanmean(att_clipped_r85[0:10,:,:],axis=0)
	att_0615_F =stats.nanmean(att_clipped_r85F[0:10,:,:],axis=0)
	att_0615_D = att_0615_R -att_0615_F
	att_3150_R = stats.nanmean(att_clipped_r85[25:45,:,:],axis=0)
	att_3150_F = stats.nanmean(att_clipped_r85F[25:45,:,:],axis=0)
	att_3150_D = att_3150_R -att_3150_F
	att_8100_R = stats.nanmean(att_clipped_r85[75:95,:,:],axis=0)
	att_8100_F = stats.nanmean(att_clipped_r85F[75:95,:,:],axis=0)
	att_8100_D = att_8100_R -att_8100_F
	att_0615_P=mannwhitneyu_test(att_clipped_r85[0:10,:,:], att_clipped_r85F[0:10,:,:])
	att_3150_P=mannwhitneyu_test(att_clipped_r85[25:45,:,:], att_clipped_r85F[25:45,:,:])
	att_8100_P=mannwhitneyu_test(att_clipped_r85[75:95,:,:], att_clipped_r85F[75:95,:,:])
	
	return lons,lats,att_0615_D,att_3150_D,att_8100_D,att_0615_P,att_3150_P,att_8100_P

lons,lats,r95p_0615_D,r95p_3150_D,r95p_8100_D,r95p_0615_P,r95p_3150_P,r95p_8100_P = spatial_diff_r95pr(region)

"""

##################
####plotting
##################
# fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(6.5, 7.65), facecolor='White');plot_setup();
fig = plt.figure(facecolor='White',figsize=(6.5, 7.65));plot_setup();	pad= 5 
#rx5day
colormap ='BrBG'; colorbar_min=-25;  colorbar_max = 25;
# ax = plt.subplot(4,2,1)
# spatial_figure_norm(ax,rx5day_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,rx5day_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('2006-2015',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')
# ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)		

ax = plt.subplot(4,2,1)
spatial_figure_norm(ax,rx5day_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,rx5day_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('RX5DAY (mm)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('2031-2050',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax = plt.subplot(4,2,2)
ax.annotate('2081-2100',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
colormesh1 = spatial_figure_norm(ax,rx5day_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,rx5day_8100_P, tb_lef=False,tb_bot=False )
cbar_ax = fig.add_axes([0.92, 0.74, 0.01, 0.21])
char = fig.colorbar(colormesh1,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))


#SDII
colormap = 'BrBG'; colorbar_min=-2.5;  colorbar_max = 2.5;
# ax=plt.subplot(4,2,4)
# spatial_figure_norm(ax,sdii_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,sdii_0615_P, tb_lef=True,tb_bot=False )

# ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)				
ax=plt.subplot(4,2,3)
spatial_figure_norm(ax,sdii_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,sdii_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('R10 (days)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)	
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(4,2,4)
colormesh2=spatial_figure_norm(ax,sdii_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,sdii_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.92, 0.505, 0.01, 0.21])
char = fig.colorbar(colormesh2,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
# char.set_label()


# CWD
colormap ='BrBG'; colorbar_min=-5;  colorbar_max = 5;
# ax=plt.subplot(4,2,7)
# spatial_figure_norm(ax,cwd_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,cwd_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(4,2,5)
spatial_figure_norm(ax,cwd_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,cwd_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('CWD (days)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(4,2,6)
colormesh3 =spatial_figure_norm(ax,cwd_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,cwd_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
#   Rcp and Fixa share the same colorbar while the diff plots use its own colorbar
cbar_ax = fig.add_axes([0.92, 0.27, 0.01, 0.21])
char = fig.colorbar(colormesh3,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
# char.set_label('day')

#r95rp
colormap ='BrBG'; colorbar_min=-35;  colorbar_max = 35;
# ax = plt.subplot(4,2,10)
# spatial_figure_norm(ax,r95p_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,r95p_0615_P, tb_lef=True,tb_bot=True )
# ax.annotate('(j)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax = plt.subplot(4,2,7)
spatial_figure_norm(ax,r95p_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,r95p_3150_P, tb_lef=True,tb_bot=True )
ax.annotate('R95P (mm)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax = plt.subplot(4,2,8)
colormesh4 = spatial_figure_norm(ax,r95p_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,r95p_8100_P, tb_lef=False,tb_bot=True )
ax.annotate('(h)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.92, 0.04, 0.01, 0.21])
char = fig.colorbar(colormesh4,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))


plt.subplots_adjust(left=0.1, bottom=0.03, right=0.9, top=0.95, wspace=0.01, hspace=0.07);
plt.savefig('Fig9.png', format='png', dpi=1000)
plt.show()