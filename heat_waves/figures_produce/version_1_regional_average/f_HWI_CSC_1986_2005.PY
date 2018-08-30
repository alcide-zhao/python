mport matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats


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

import scipy.io as sio
from matplotlib.pyplot import *
import matplotlib
import numpy as np

cdict = {'red': ((0., 0.85,0.85),
                 (1.0/6.0, 1,1),
                 (2.0/6.0, 1, 1),
                 (3.0/6.0, 0.80, 0.80),
				 (4.0/6.0,1.0,1.0),
				 (5.0/6.0,0.82, 0.82),
				 (6.0/6.0,0.80, 0.80)),
         'green': ((0., 0.85,0.85),
                   (1.0/6.0, 1,1),
                   (2.0/6.0, 0.41, 0.41),
                   (3.0/6.0, 0.52, 0.52),
				   (4.0/6.0, 0,0),
				   (5.0/6.0,0.71, 0.71),
				   (6.0/6.0, 0.13, 0.13)),
         'blue': ((0., 0.85,0.85),
				  (1.0/6.0, 0,0), 
                  (2.0/6.0, 0.71, 0.71),
		          (3.0/6.0, 0.25, 0.25),
				  (4.0/6.0, 0,0),
				  (5.0/6.0,0.55, 0.55),
				  (6.0/6.0, 0.13, 0.13))}
HWM_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,7)  # 7 slots in the colorlist


## baseline data
baseline_data ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/HWI_CSI_1961_1990_9010_anomoly.nc'
nc_fid = nc4.Dataset(baseline_data,mode='r')
RFWY= stats.nanmean(nc_fid.variables['RFWY'],axis=0);

nc_fid.close()

file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/HWI_CSI_EnMe_fixa.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
HWI = nc_fid.variables['HWI'][:]
HWI = nc_fid.variables['HWI'][:]


fig = plt.figure(facecolor='White',figsize=(8, 12));plot_setup();pad= 5 


#0-annual mean;            1- annual minimum;         2-annual maximum;            3-annual counts;
#4-no. of mild(<=-1);      5-no. of nornal(-1-1);     6-no. of moderate(1-2);      7-no. of severe(2-4);
#8-no. of V. extreme(4-8); 9-no. of V. extreme(8-16); 10-no. of S. extreme(16-32); 11-no. of U. extreme(>=32)
#12-RF mean    13-RF min   14-RF max 

# print lon
# data_r =np.nansum(HWI[75:95,:,:,3],axis=0)-np.nansum(HWI[75:95,:,:,4],axis=0)\
# -np.nansum(HWI[75:95,:,:,5],axis=0)-np.nansum(HWI[75:95,:,:,6],axis=0)\
# -np.nansum(HWI[75:95,:,:,7],axis=0)-np.nansum(HWI[75:95,:,:,8],axis=0)#-np.nansum(HWI[75:95,:,:,9],axis=0)
# HW magnitude 
ax = plt.subplot(3,1,1);colormap='jet';colorbar_min=-2;colorbar_max=5
data = np.nansum(HWI[0:15,:,:,4],axis=0); 
# data[(data<=1)] = -0.5;data[(data>1) & (data<=2)] = 0.5;
# data[(data>2) & (data<=4)] = 1.5;data[(data>4) & (data<=8)] = 2.5;
# data[(data>8) & (data<=16)] = 3.5;data[(data>16) & (data<=32)] = 4.5;
# data[(data>32) & (data<=64)] = 5.5;data[(data>64)] = 6.5;

p_value=np.empty((np.size(data)));p_value[:]=np.nan
colormesh1= spatial_figure(ax,data,lon,lat,colormap,colorbar_min,colorbar_max,p_value,c_number=17, c_bad=[1,1,1],c_under=[1,1,1],c_over=[0,0,0],tb_lef=False,tb_bot=False );  #0.65,0.16,0.16  [0.12,0.56,1]
ax.annotate('(a) Maximum of Annual Maximum HWMI (2060-2080)',xy=(0.02, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

cbar_ax = fig.add_axes([0.95, 0.67, 0.01, 0.28])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks= np.array([-2,0,1,2,4,8,16,32])); 
# char.ax.set_yticklabels(['1','2','4','8','16','32','64'])

## HW counts
ax = plt.subplot(3,1,2);colormap='jet';colorbar_min=0;colorbar_max=10
data =stats.nanmean(HWI[55:75,:,:,3],axis=0)#-np.nansum(HWI[55:75,:,:,4],axis=0)\
#-np.nansum(HWI[55:75,:,:,5],axis=0)-np.nansum(HWI[75:95,:,:,6],axis=0)\
#-np.nansum(HWI[75:95,:,:,7],axis=0)-np.nansum(HWI[75:95,:,:,8],axis=0)#-np.nansum(HWI[75:95,:,:,9],axis=0)
ax.annotate('(b) Annual mean no. of HW events during 2060-2080',xy=(0.02, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
p_value=np.empty((np.size(data)));p_value[:]=np.nan
colormesh1= spatial_figure(ax,data,lon,lat,colormap,colorbar_min,colorbar_max,p_value,c_number=10, c_bad=[1,1,1],c_under=[1,1,1],c_over=[0,0,0],tb_lef=False,tb_bot=False  );
cbar_ax = fig.add_axes([0.95, 0.36, 0.01, 0.28])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks= np.arange(0,10.1,2)); 

ax = plt.subplot(3,1,3);colormap='jet';colorbar_min=-1;colorbar_max=1
data =1-np.divide(np.nanmin(HWI[55:75,:,:,12],axis=0),RFWY)
ax.annotate('(c) The minimum of RF 2060-2080',xy=(0.02, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
p_value=np.empty((np.size(data)));p_value[:]=np.nan
colormesh1= spatial_figure(ax,data,lon,lat,colormap,colorbar_min,colorbar_max,p_value,c_number=10, c_bad=[1,1,1],c_under=[1,1,1],c_over=[0,0,0],tb_lef=False,tb_bot=False  );
cbar_ax = fig.add_axes([0.95, 0.05, 0.01, 0.28])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks= np.arange(-1.0,1.1,0.4)); 

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.1, hspace=0.1);
plt.savefig('Fig1.png', format='png', dpi=1000)
plt.show()


