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

def get_TS():
	##surface temperature
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_TS_200602_210101.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	R85=nc_fid.variables['rcp85'];
	FiX=nc_fid.variables['rcp85_fixA'];
	
	rcp85_T=np.empty((95,192,288)); rcp85_T[:]=np.nan
	FixA_T=np.empty((95,192,288)); FixA_T[:]=np.nan
	for iyear in range(0,95):
		rcp85_T[iyear,:,:] = stats.nanmean(R85[iyear*12:(iyear+1)*12,:,:],axis=0)
		FixA_T[iyear,:,:] = stats.nanmean(FiX[iyear*12:(iyear+1)*12,:,:],axis=0)
	Aerosol_T = rcp85_T-FixA_T
	nc_fid.close();
	del rcp85_T
	return FixA_T,Aerosol_T

	
def HWI_CSI_gradient_pp():
	"""
	GET THE SPATIAL MAPS OF GRADIENT AGAIST PER DEGREE OF WARMING
	"""
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	#rcp85
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_rcp85.nc',mode='r')
	HWI_X_R = nc_fid.variables['HWI_X'][1,:,:,:]
	CSI_X_R = nc_fid.variables['CSI_X'][1,:,:,:]
	nc_fid.close()
	#fixa
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_fixa.nc',mode='r')
	HWI_X_F = nc_fid.variables['HWI_X'][1,:,:,:]
	CSI_X_F =  nc_fid.variables['CSI_X'][1,:,:,:]
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	nc_fid.close()
	## aerosol
	HWI_X_A= HWI_X_R-HWI_X_F; CSI_X_A= CSI_X_R-CSI_X_F;
	del HWI_X_R,CSI_X_R
	GHG_T,AA_T = get_TS();
	HW_GHG=np.empty((192,288));HW_GHG[:]=np.nan;
	CS_GHG=np.empty((192,288));CS_GHG[:]=np.nan;
	HW_AA=np.empty((192,288));HW_AA[:]=np.nan;
	CS_AA=np.empty((192,288));CS_AA[:]=np.nan;
	
	for row in range(len(lat)):
		for column in range(len(lon)):
			if (~np.isnan(HWI_X_A[0,row,column])):
				x=GHG_T[:,row,column]; 
				y = HWI_X_F[:,row,column]; mask = ~np.isnan(x) & ~np.isnan(y); HW_GHG[row,column], _, _, _, _ = stats.linregress(x[mask],y[mask]);
				y = CSI_X_F[:,row,column]; mask = ~np.isnan(x) & ~np.isnan(y); CS_GHG[row,column], _, _, _, _ = stats.linregress(x[mask],y[mask]);
				x=AA_T[:,row,column]; 
				y = HWI_X_A[:,row,column]; mask = ~np.isnan(x) & ~np.isnan(y); HW_AA[row,column], _, _, _, _ = stats.linregress(x[mask],y[mask]);
				y = CSI_X_A[:,row,column]; mask = ~np.isnan(x) & ~np.isnan(y); CS_AA[row,column], _, _, _, _ = stats.linregress(x[mask],y[mask]);
	return lon,lat,HW_GHG,CS_GHG,HW_AA,CS_AA


lon,lat,HW_GHG,CS_GHG,HW_AA,CS_AA =HWI_CSI_gradient_pp()


fig = plt.figure(facecolor='White',figsize=[12,7]);plot_setup();pad= 5;

colormap='RdBu';colormap = reverse_colourmap(colormap)
ax = plt.subplot(2,2,1);colorbar_min=0;colorbar_max=10;p_value = np.zeros((np.shape(HW_GHG)))
ax.annotate('GHG',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('Heat Waves',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(a)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,HW_GHG,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)

ax = plt.subplot(2,2,3);
ax.annotate('Aerosol',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
colormesh1 = spatial_figure(ax,HW_AA,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.1, 0.05, 0.40, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))


colorbar_min=-0.08;colorbar_max=0.05;
ax = plt.subplot(2,2,2);
ax.annotate('Cold Spells',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(b)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
spatial_figure(ax,CS_GHG,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)


ax = plt.subplot(2,2,4);
ax.annotate('(d)',xy=(0.02,0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
colormesh1=spatial_figure(ax,CS_AA,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.56, 0.05, 0.40, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))


plt.subplots_adjust(left=0.08, bottom=0.10, right=0.98, top=0.92, wspace=0.05, hspace=0.05); 
plt.savefig('gradient_bar_local.png', format='png', dpi=1000)
