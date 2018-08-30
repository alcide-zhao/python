"""
This is to plot the ozone concentration difference between 2010 and 1970 
the 1970 data is from ozone_1.9x2.5_L26_1850-2005_c090803.nc
the 2010 data is from ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202.nc
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

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *



def get_ozone():
	input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/ozone/'
	ozone1970 = input_path+'ozone_1.9x2.5_L26_1850-2005_c090803.nc'
	nc_fid = nc4.Dataset(ozone1970,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	lev = nc_fid.variables['lev'][:]
	o1970 = nc_fid.variables['O3'][156:168,:,:,:]
	nc_fid.close()
	ozone2010 = input_path+'ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202.nc'
	nc_fid = nc4.Dataset(ozone2010,mode='r')
	o2010 = nc_fid.variables['O3'][48:60,:,:,:]
	nc_fid.close()
	ozone_lat_alt_cross = stats.nanmean(stats.nanmean(o2010-o1970,axis=3),axis=0)
	ozone_lat_alt_cross = np.ma.masked_invalid(ozone_lat_alt_cross)
	ozone_lat_alt_cross = np.ma.masked_equal(ozone_lat_alt_cross,0)
	return lat,lev,ozone_lat_alt_cross
	

lat,lev,ozone_lat_alt_cross = get_ozone();
x,y= np.meshgrid(lat,lev)

print lev[:15]
fig = plt.figure(facecolor='White',figsize=[6,4]);plot_setup();pad= 5;cmap = 'RdBu_r';colorbar_min=-0.2;colorbar_max=0.2;
ax =plt.subplot2grid((3, 1), (0, 0), rowspan=1)
cmap = discrete_cmap(16,cmap)
colormesh1 = plt.pcolormesh(x[0:13,:],np.flipud(y[0:13,:]),np.flipud(ozone_lat_alt_cross[0:13,:])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([0,192.539935]);ax.set_xticks([]);
plt.gca().invert_yaxis()
cbar_ax = fig.add_axes([0.90, 0.67, 0.015, 0.27])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.25)*(colorbar_max-colorbar_min)+colorbar_min,2))
ax.annotate(r'Pressure (hPa)',xy=(-0.1,-0.35), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)
				
ax =plt.subplot2grid((3, 1), (1, 0), rowspan=2);colorbar_min=-0.02;colorbar_max=0.02;
cmap = discrete_cmap(16,cmap)
colormesh1 = plt.pcolormesh(x[13:-1,:],np.flipud(y[13:-1,:]),np.flipud(ozone_lat_alt_cross[13:-1,:])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([226.513265,1000]);ax.set_xticks(range(-90,91,30));ax.set_ylabel('');ax.set_xlabel('Lattitude (degree)');
plt.gca().invert_yaxis()
cbar_ax = fig.add_axes([0.90, 0.12, 0.015, 0.52])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.25)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.1, bottom=0.10, right=0.88, top=0.95, wspace=0.04, hspace=0.01); 
plt.savefig('Ozone.png', format='png', dpi=600)







