# -*- coding: utf-8 -*-
"""
This is to show the variation of Aerosol burdens and AOD variations 
in 2036-2045 and 2090-2100 with respct to current day (2006-2015) levels
@author: Alcide.Zhao
"""

import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import os
from scipy import stats
import site
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
# from scipy.interpolate import interp2d

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


rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]} #  'SA':[75,98,25,32]

region = rergion_dic['SA'][:]
import scipy.io as sio
file_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
# burdens
burden_file = file_path+'Aerosol_Burden_AOD_JJA_2006_2100.mat'
data = sio.loadmat(burden_file)
time_B = data['time'][0,:];lon = data['lon'][0,:];lat = data['lat'][0,:];
BC_RCP85 = data['BC_RCP85'];BC_FIXA = data['BC_FIXA'];
OC_RCP85 = data['OC_RCP85'];OC_FIXA = data['OC_FIXA'];
SO4_RCP85 = data['SO4_RCP85'];SO4_FIXA = data['SO4_FIXA'];
AOD_RCP85 = data['AOD_RCP85'];AOD_FIXA = data['AOD_FIXA'];
diff=BC_RCP85-BC_FIXA # np.divide(BC_RCP85-BC_FIXA,BC_FIXA)*100
lons,lats,BC = range_clip(region[0],region[1],region[2],region[3],lon,lat,diff)
diff=OC_RCP85-OC_FIXA # np.divide(OC_RCP85-OC_FIXA,OC_FIXA)*100
lons,lats,OC = range_clip(region[0],region[1],region[2],region[3],lon,lat,diff)
diff= SO4_RCP85-SO4_FIXA#np.divide(SO4_RCP85-SO4_FIXA,SO4_FIXA)*100
lons,lats,SO4 = range_clip(region[0],region[1],region[2],region[3],lon,lat,diff)
diff= AOD_RCP85-AOD_FIXA
lons,lats,AOD = range_clip(region[0],region[1],region[2],region[3],lon,lat,diff)


AOD_3645_D = stats.nanmean(AOD[75:95,:,:],axis=0)
print stats.nanmean(stats.nanmean(AOD_3645_D,axis=1),axis=0)


# BC_0615_D = stats.nanmean(BC[0:10,:,:],axis=0)*10**6
# OC_0615_D = stats.nanmean(OC[0:10,:,:],axis=0)*10**6
# SO4_0615_D = stats.nanmean(SO4[0:10,:,:],axis=0)*10**6
# AOD_0615_D = stats.nanmean(AOD[0:10,:,:],axis=0)
BC_3645_D = stats.nanmean(BC[25:45,:,:],axis=0)#*10**6
OC_3645_D = stats.nanmean(OC[25:45,:,:],axis=0)#*10**6
SO4_3645_D = stats.nanmean(SO4[25:45,:,:],axis=0)#*10**6
AOD_3645_D = stats.nanmean(AOD[25:45,:,:],axis=0)
BC_9100_D = stats.nanmean(BC[75:95,:,:],axis=0)#*10**6
OC_9100_D = stats.nanmean(OC[75:95,:,:],axis=0)#*10**6
SO4_9100_D = stats.nanmean(SO4[75:95,:,:],axis=0)#*10**6
AOD_9100_D = stats.nanmean(AOD[75:95,:,:],axis=0)
print stats.nanmean(stats.nanmean(BC_3645_D,axis=1),axis=0)
print stats.nanmean(stats.nanmean(OC_3645_D,axis=1),axis=0)
print stats.nanmean(stats.nanmean(SO4_3645_D,axis=1),axis=0)
# now plotting
# fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(8.5, 3.5), facecolor='White');plot_setup();
fig = plt.figure(facecolor='White',figsize=(8.5, 3.5));plot_setup();pad= 5 
p_value = np.zeros((np.shape(BC_3645_D))); colormap ='RdBu';colormap= reverse_colourmap(colormap)	 
# BC Burden
cb_min=-2.0; cb_max = 2.0;
cb_min=0; cb_max = 50;
# ax = plt.subplot(2,4,1)
# spatial_figure_norm(ax,BC_0615_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
# ax.annotate('2006-2015',xy=(0.5, 1.02), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',fontsize=10)

# ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax = plt.subplot(2,4,1)
spatial_figure_norm(ax,BC_3645_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
ax.annotate(r'BC ($\mathrm{\mathsf{mg\/m^{-2}}}$)',xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',fontsize=10)
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('2031-2050',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',fontsize=10,rotation = 90)
ax = plt.subplot(2,4,5)
colormesh1=spatial_figure_norm(ax,BC_9100_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=True)
ax.annotate('2081-2100',xy=(-0.25, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',fontsize=10,rotation = 90)
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

cbar_ax = fig.add_axes([0.09, 0.07, 0.19, 0.02])
char = fig.colorbar(colormesh1, cax=cbar_ax,extend='both',orientation='horizontal',ticks=np.round(np.arange(cb_min,cb_max+0.01,(cb_max-cb_min)/4),2))


# OC Burden		
cb_min=-8; cb_max = 8;	
cb_min=0; cb_max = 50;
# ax=plt.subplot(2,4,3)
# spatial_figure_norm(ax,OC_0615_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)
# ax.annotate('OC',xy=(-0.2, 0.5), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',fontsize=10,rotation = 90)
# ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax = plt.subplot(2,4,2)
spatial_figure_norm(ax,OC_3645_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('OC ($\mathrm{\mathsf{mg\/m^{-2}}}$)',xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',fontsize=10,)
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(2,4,6)
colormesh2 =spatial_figure_norm(ax,OC_9100_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.32, 0.07, 0.195, 0.02])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both',orientation='horizontal',ticks=np.round(np.arange(cb_min,cb_max+0.01,(cb_max-cb_min)/4),2))

# SO4
cb_min=-20; cb_max = 20;	
cb_min=0; cb_max = 50;
# ax=plt.subplot(2,4,7)
# spatial_figure_norm(ax,SO4_0615_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=False)

# ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax=plt.subplot(2,4,3)
spatial_figure_norm(ax,SO4_3645_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('SO4 ($\mathrm{\mathsf{mg\/m^{-2}}}$)',xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',fontsize=10)
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(2,4,7)
colormesh3 =spatial_figure_norm(ax,SO4_9100_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

cbar_ax = fig.add_axes([0.55, 0.07, 0.19, 0.02])
char = fig.colorbar(colormesh3, cax=cbar_ax,extend='both',orientation='horizontal',ticks=np.round(np.arange(cb_min,cb_max+0.01,(cb_max-cb_min)/4),2))
# AOD550
colormap= 'RdBu';colormap= reverse_colourmap(colormap)	
cb_min=-0.25; cb_max = 0.25;	
# ax=plt.subplot(2,4,10)
# spatial_figure_norm(ax,AOD_0615_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=True,tb_bot=True)

# ax.annotate('(j)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(2,4,4)
spatial_figure_norm(ax,AOD_3645_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=False)
ax.annotate('AOD',xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',fontsize=10)
ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(2,4,8)
colormesh4=spatial_figure_norm(ax,AOD_9100_D,lons,lats,colormap,cb_min,cb_max,p_value,tb_lef=False,tb_bot=True)		
ax.annotate('(h)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
cbar_ax = fig.add_axes([0.775, 0.07, 0.19, 0.02])
# char = fig.colorbar(colormesh4,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.25,0.4)*(cb_max-cb_min)+cb_min,2))
char = fig.colorbar(colormesh4, cax=cbar_ax,extend='both',orientation='horizontal',ticks=np.round(np.arange(cb_min,cb_max+0.01,(cb_max-cb_min)/4),2))


plt.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.93, wspace=0.1, hspace=0.00001);
plt.savefig('Fig5.png', format='png', dpi=1000)
plt.show()