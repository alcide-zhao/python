# -*- coding: utf-8 -*-
'''
Global mao of HWI and CSI gradients against per degree of warming for GHG and AAs
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as sio

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



file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/' #intereannual_osci/
#rcp85
nc_fid = nc4.Dataset(file_path+'HWI_CSI_1961_1990_0595_ano.nc',mode='r')
lon= nc_fid.variables['lon'][:];lat= nc_fid.variables['lat'][:];
HWIM = stats.nanmean(nc_fid.variables['HWIM'][:],axis=0)
CSIM = stats.nanmean(nc_fid.variables['CSIM'][:],axis=0)
HWIS = stats.nanmean(nc_fid.variables['HWIS'][:],axis=0)
CSIS = stats.nanmean(nc_fid.variables['CSIS'][:],axis=0)
nc_fid.close()

fig = plt.figure(facecolor='White',figsize=[12,7]);plot_setup();pad= 5;
colorbar_min=0;colorbar_max=20;
colormap='RdBu';colormap = reverse_colourmap(colormap)
ax = plt.subplot(2,2,1);p_value = np.zeros((np.shape(HWIM)));
ax.annotate('mean',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('Heat Waves',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(a)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,HWIM,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)

ax = plt.subplot(2,2,3);
ax.annotate('std',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
colormesh1 = spatial_figure(ax,HWIS,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.1, 0.05, 0.40, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))

colorbar_min=10;colorbar_max=30;
ax = plt.subplot(2,2,2);
ax.annotate('Cold Spells',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(b)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
spatial_figure(ax,CSIM,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,4);
ax.annotate('(d)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
colormesh1=spatial_figure(ax,CSIS,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.56, 0.05, 0.40, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))


plt.subplots_adjust(left=0.08, bottom=0.10, right=0.98, top=0.92, wspace=0.05, hspace=0.05); 
plt.savefig('HW_CS_MEAN_STD_ANO.png', format='png', dpi=1000)
