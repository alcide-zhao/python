'''
This is to compare the HWI and CSI from NCEP and CESM ensemble mean ()
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as io
from scipy.interpolate import interp2d  as interp2d

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/ano/'

NCEP_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_NCEP.nc',mode='r')
NCEP_lon = NCEP_fid.variables['lon'][:]
NCEP_lat = NCEP_fid.variables['lat'][:]
NCEP_year = NCEP_fid.variables['year'][:]
NCEP_HWDC = stats.nanmean(NCEP_fid.variables['HWDC'][38:58,:,:],axis=0)
NCEP_HWDuM = stats.nanmean(NCEP_fid.variables['HWDuM'][38:58,:,:],axis=0)
NCEP_HWI_NO = stats.nanmean(NCEP_fid.variables['HWI_NO'][38:58,:,:],axis=0)
NCEP_HWInM = stats.nanmean(NCEP_fid.variables['HWInM'][38:58,:,:],axis=0)
NCEP_fid.close()


CESM_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_his.nc',mode='r')
CESM_lon = CESM_fid.variables['lon'][:]
CESM_lat = CESM_fid.variables['lat'][:]
CESM_year = CESM_fid.variables['year'][:]
CESM_HWDC = stats.nanmean(CESM_fid.variables['HWDC'][1,66:86,:,:],axis=0)
CESM_HWDuM = stats.nanmean(CESM_fid.variables['HWDuM'][1,66:86,:,:],axis=0)
CESM_HWI_NO = stats.nanmean( CESM_fid.variables['HWI_NO'][1,66:86,:,:],axis=0)
CESM_HWInM = stats.nanmean(CESM_fid.variables['HWInM'][1,66:86,:,:],axis=0)
CESM_fid.close()

# def relative_bias(CESM,NCEP):
	# LON = np.arange(0,360,2.5); LAT = np.arange(-90,90,2.5)
	# NCEP[np.isnan(NCEP)] =0
	# f1 = interp2d(NCEP_lon,NCEP_lat,NCEP)
	# NCEP_in = f1(LON, LAT);NCEP_in[NCEP_in==0]=np.nan;
	# CESM[np.isnan(CESM)] =0
	# f1 = interp2d(CESM_lon,CESM_lat,CESM)
	# CESM_in = f1(LON, LAT);CESM_in[CESM_in==0]=np.nan;
	# bias_a = stats.nanmean(stats.nanmean(CESM_in-NCEP_in,axis=1),axis=0)
	# bias_r = stats.nanmean(stats.nanmean(np.divide(CESM_in-NCEP_in,NCEP_in),axis=1),axis=0)*100
	# # CESM = stats.nanmean(stats.nanmean(CESM,axis=1),axis=0)
	# # NCEP = stats.nanmean(stats.nanmean(NCEP,axis=1),axis=0)
	# # bias_a = CESM-NCEP
	# # bias_r = (np.divide(CESM,NCEP)-1)*100
	# bias=np.array([bias_a,bias_r])
	# return bias

bias = relative_bias(CESM_HWInM,NCEP_HWInM); print bias
bias = relative_bias(CESM_HWDuM,NCEP_HWDuM); print bias
bias = relative_bias(CESM_HWI_NO,NCEP_HWI_NO); print bias
bias = relative_bias(CESM_HWDC,NCEP_HWDC); print bias

"""
def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,p_value,color_slice,tb_lef=False,tb_bot=False): #c_bad,c_under,c_over,c_number=20,
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
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='.',color='k',zorder=10)
	
	# Add Grid Lines
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines();# map.drawstates(); map.drawcountries()
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(color_slice,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); cmap.set_under('b')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min-0.01, vmax=colorbar_max,latlon=True) #
	return colormesh

### plot
fig = plt.figure(facecolor='White',figsize=[8.5,8]);plot_setup();
colormap ='RdBu';  colormap=reverse_colourmap(colormap);pad= 5;

p_value = np.zeros((np.shape(CESM_HWDC)))
# ax = plt.subplot(4,2,1);cb_min = 0;cb_max =4;
ax = plt.subplot(4,2,1);cb_min = -6;cb_max =0;
spatial_figure(ax,CESM_HWInM,CESM_lon,CESM_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('Intensity\n($^\circ$C)',xy=(-0.10, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('CESM1 ensemble mean',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('(a)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax = plt.subplot(4,2,3);cb_min = 3;cb_max =6
spatial_figure(ax,CESM_HWDuM,CESM_lon,CESM_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('Duration\n(days/event)',xy=(-0.10, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
				
ax = plt.subplot(4,2,5);cb_min = 0;cb_max =3
spatial_figure(ax,CESM_HWI_NO,CESM_lon,CESM_lat,colormap,cb_min,cb_max,p_value,color_slice=12)
ax.annotate('Frequency\n(events/year)',xy=(-0.10, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(e)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						

ax = plt.subplot(4,2,7);cb_min = 0;cb_max =16
spatial_figure(ax,CESM_HWDC,CESM_lon,CESM_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('Total days\n(days/year)',xy=(-0.10, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(g)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		

#### NCEP				
p_value = np.zeros((np.shape(NCEP_HWDC)))
# ax = plt.subplot(4,2,2);cb_min = 0;cb_max =4
ax = plt.subplot(4,2,2);cb_min = -6;cb_max =0;
colormesh=spatial_figure(ax,NCEP_HWInM,NCEP_lon,NCEP_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('NCEP/NCAR reanalysis',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('(b)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.93, 0.735, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.arange(-6,0.1,2))	
		
ax = plt.subplot(4,2,4);cb_min = 3;cb_max =6
colormesh=spatial_figure(ax,NCEP_HWDuM,NCEP_lon,NCEP_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('(d)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
cbar_ax = fig.add_axes([0.93, 0.50, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.arange(3,6.1,1))	
				
ax = plt.subplot(4,2,6);cb_min = 0;cb_max =3
colormesh=spatial_figure(ax,NCEP_HWI_NO,NCEP_lon,NCEP_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('(f)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
cbar_ax = fig.add_axes([0.93, 0.27, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.arange(0,3.1,1))	

ax = plt.subplot(4,2,8);cb_min = 0;cb_max =12
colormesh=spatial_figure(ax,NCEP_HWDC,NCEP_lon,NCEP_lat,colormap,cb_min,cb_max,p_value,color_slice=12 )
ax.annotate('(h)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.93, 0.035, 0.01, 0.21])
char = fig.colorbar(colormesh,orientation='vertical',cax=cbar_ax,extend='both',ticks=np.arange(0,12.1,4))					
				
plt.subplots_adjust(left=0.1, bottom=0.03, right=0.92, top=0.95, wspace=0.01, hspace=0.05); 
plt.savefig('CS_INDICES_CESM_NCEP_86_05.png', format='png', dpi=1000)
plt.show()

"""









