# -*- coding: utf-8 -*-
'''
Global mao of HW and CS with magnitude greater than 4(8)
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

ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;ocean_mask_CESM[0:27,:]=np.nan

def total_no_slice(value):
	"""
	count the number of HW with magnitude exceeding 4 over 20yr time slice
	"""
	no_4160=np.multiply(np.nansum(value[25:45,:,:],axis=0),ocean_mask_CESM);no_4160[no_4160==0]=np.nan
	no_8100=np.multiply(np.nansum(value[75:95,:,:],axis=0),ocean_mask_CESM);no_8100[no_8100==0]=np.nan
	return no_4160,no_8100
	
def HWI_gradient_pp():

	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	#his
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_his.nc',mode='r')
	his_no = np.nansum((nc_fid.variables['HWCV'][1,:,:,:]+nc_fid.variables['HWCP'][1,:,:,:]+nc_fid.variables['HWCU'][1,:,:,:])[66:86,:,:],axis=0) #
	nc_fid.close()	
	#rcp85
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_rcp85.nc',mode='r')
	rcp_no = nc_fid.variables['HWCV'][1,:,:,:]+nc_fid.variables['HWCP'][1,:,:,:]+nc_fid.variables['HWCU'][1,:,:,:]-his_no  #+nc_fid.variables['HWCE'][1,:,:,:]+nc_fid.variables['HWCE'][1,:,:,:]
	nc_fid.close()
	#fixa
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_fixa.nc',mode='r')
	fix_no = nc_fid.variables['HWCV'][1,:,:,:]+nc_fid.variables['HWCP'][1,:,:,:]+nc_fid.variables['HWCU'][1,:,:,:]-his_no  #+nc_fid.variables['HWCE'][1,:,:,:]+nc_fid.variables['HWCD'][1,:,:,:]
	
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	## aerosol
	Aer_no =rcp_no-fix_no
	GHG_4160,GHG_8100 = total_no_slice(fix_no)
	Aer_4160,Aer_8100 = total_no_slice(Aer_no)
	nc_fid.close()
	return lon,lat,GHG_4160,GHG_8100,Aer_4160,Aer_8100
	
	
def CSI_gradient_pp():
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	#his
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_his.nc',mode='r')
	his_no = np.nansum((nc_fid.variables['CSCE'][1,:,:,:]+nc_fid.variables['CSCV'][1,:,:,:]+nc_fid.variables['CSCP'][1,:,:,:]+nc_fid.variables['CSCS'][1,:,:,:])[66:86,:,:],axis=0)
	nc_fid.close()	 
	#rcp85
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_rcp85.nc',mode='r')
	rcp_no = nc_fid.variables['CSCE'][1,:,:,:]+nc_fid.variables['CSCV'][1,:,:,:]+nc_fid.variables['CSCP'][1,:,:,:]+nc_fid.variables['CSCS'][1,:,:,:]-his_no#+nc_fid.variables['CSCD'][1,:,:,:]
	nc_fid.close()
	#fixa
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_fixa.nc',mode='r')
	fix_no = nc_fid.variables['CSCE'][1,:,:,:]+nc_fid.variables['CSCV'][1,:,:,:]+nc_fid.variables['CSCP'][1,:,:,:]+nc_fid.variables['CSCS'][1,:,:,:]-his_no#+nc_fid.variables['CSCD'][1,:,:,:]
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	## aerosol
	Aer_no =rcp_no-fix_no
	GHG_4160,GHG_8100 = total_no_slice(fix_no)
	Aer_4160,Aer_8100 = total_no_slice(Aer_no)
	nc_fid.close()
	return lon,lat,GHG_4160,GHG_8100,Aer_4160,Aer_8100
	
	
lon,lat,GHG_4160_H,GHG_8100_H,Aer_4160_H,Aer_8100_H =HWI_gradient_pp()
lon,lat,GHG_4160_C,GHG_8100_C,Aer_4160_C,Aer_8100_C =CSI_gradient_pp()

fig = plt.figure(facecolor='White',figsize=[12,7]);plot_setup();pad= 5;
colormap='Reds';#colormap = reverse_colourmap(colormap)
ax = plt.subplot(2,2,1);colorbar_min=0;colorbar_max=50;p_value = np.zeros((np.shape(GHG_8100_C)))
ax.annotate('GHG',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('Heat Wave',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(a)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,GHG_8100_H,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)

ax = plt.subplot(2,2,3);
ax.annotate('Aerosol',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
colormesh1 = spatial_figure(ax,Aer_8100_H,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.1, 0.05, 0.40, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

colormap='RdBu';colormap = reverse_colourmap(colormap)
colorbar_min=-5;colorbar_max=5;
ax = plt.subplot(2,2,2);
ax.annotate('Cold Spell',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(b)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
spatial_figure(ax,GHG_8100_C,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)

ax = plt.subplot(2,2,4);
ax.annotate('(d)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
colormesh1=spatial_figure(ax,Aer_8100_C,lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.56, 0.05, 0.40, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))


plt.subplots_adjust(left=0.05, bottom=0.10, right=0.96, top=0.92, wspace=0.04, hspace=0.05); 
plt.savefig('NoExceeding_CalDayThr.png', format='png', dpi=1000)
