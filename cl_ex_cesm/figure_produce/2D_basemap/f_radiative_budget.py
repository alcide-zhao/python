# -*- coding: utf-8 -*-
"""
This is to show the variation in  shortwave solar flux at surface from both clear-sky and all-sky 
and the surface latent and sensible flux
"""
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
    # 'lib'
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio



rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}

region = rergion_dic['ASIA'][:]
oceanmask=sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
nc_fid.close()
_,_,ocean_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_TS_PSL_Q_UV(Q)_LHFLX_FSNS(C)_FSNTOA(C).mat')
time_B = data['time'][0,:];lon = data['lon'][0,:];lat = data['lat'][0,:];
FSNTOA_rcp85 = data['FSNTOA_rcp85'];FSNTOA_fixa = data['FSNTOA_fixa'];
FSNSC_rcp85 = data['FSNSC_rcp85'];FSNSC_fixa = data['FSNSC_fixa'];
FSNS_rcp85 = data['FSNS_rcp85'];FSNS_fixa = data['FSNS_fixa'];
LHFLX_rcp85 = data['LHFLX_rcp85'];LHFLX_fixa = data['LHFLX_fixa'];
# TS_rcp85 = data['TS_rcp85'];TS_fixa = data['TS_fixa'];
del data

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
	# p_value=np.multiply(ocean_mask,p_value)
	return p_value

def spatial_diff_sig(region,lon,lat,variable_rcp85,variable_fixa):
	from scipy.stats import ttest_ind as test
	lons,lats,r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable_rcp85)
	lons,lats,fixa = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable_fixa)
	
	att_0615_R =stats.nanmean(r85[0:10,:,:],axis=0)
	att_0615_F =stats.nanmean(fixa[0:10,:,:],axis=0)
	att_0615_D = att_0615_R -att_0615_F
	att_3150_R = stats.nanmean(r85[25:45,:,:],axis=0)
	att_3150_F = stats.nanmean(fixa[25:45,:,:],axis=0)
	att_3150_D = att_3150_R -att_3150_F
	att_8100_R = stats.nanmean(r85[75:95,:,:],axis=0)
	att_8100_F = stats.nanmean(fixa[75:95,:,:],axis=0)
	att_8100_D = att_8100_R -att_8100_F
	p_value_thres = 0.05
	att_0615_P=mannwhitneyu_test(r85[0:10,:,:], fixa[0:10,:,:])
	att_3150_P=mannwhitneyu_test(r85[25:45,:,:], fixa[25:45,:,:])
	att_8100_P=mannwhitneyu_test(r85[75:95,:,:], fixa[75:95,:,:])
	
	return lons,lats,att_0615_D,att_3150_D,att_8100_D,att_0615_P,att_3150_P,att_8100_P

lons,lats,FSNSTOA_0615_D,FSNSTOA_3150_D,FSNSTOA_8100_D,FSNSTOA_0615_P,FSNSTOA_3150_P,FSNSTOA_8100_P = spatial_diff_sig(region,lon,lat,FSNTOA_rcp85,FSNTOA_fixa)
lons,lats,FSNSC_0615_D,FSNSC_3150_D,FSNSC_8100_D,FSNSC_0615_P,FSNSC_3150_P,FSNSC_8100_P = spatial_diff_sig(region,lon,lat,FSNSC_rcp85,FSNSC_fixa)
lons,lats,FSNS_0615_D,FSNS_3150_D,FSNS_8100_D,FSNS_0615_P,FSNS_3150_P,FSNS_8100_P = spatial_diff_sig(region,lon,lat,FSNS_rcp85,FSNS_fixa)
lons,lats,LHFLX_0615_D,LHFLX_3150_D,LHFLX_8100_D,LHFLX_0615_P,LHFLX_3150_P,LHFLX_8100_P = spatial_diff_sig(region,lon,lat,LHFLX_rcp85,LHFLX_fixa)
# lons,lats,TS_0615_D,TS_3150_D,TS_8100_D,TS_0615_P,TS_3150_P,TS_8100_P = spatial_diff_sig(region,lon,lat,TS_rcp85,TS_fixa)

# now plotting
# fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(5.9, 7), facecolor='White');plot_setup();pad= 5;
fig = plt.figure(facecolor='White',figsize=(5.9, 5));plot_setup();	pad= 5 

colorbar_min=-22;  colorbar_max = 22; colormap ='seismic';#colormap= reverse_colourmap(colormap)
# ax = plt.subplot(4,2,1)
# spatial_figure_norm(ax,FSNSC_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSC_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
# ax.annotate('2006-2015',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')

ax = plt.subplot(3,2,1)
spatial_figure_norm(ax,FSNSC_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSC_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('FSNSC',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('2031-2050',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
ax = plt.subplot(3,2,2)
ax.annotate('2081-2100',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
colormesh1 = spatial_figure_norm(ax,FSNSC_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSC_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

# ax=plt.subplot(4,2,4)
# spatial_figure_norm(ax,FSNS_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNS_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(3,2,3)
spatial_figure_norm(ax,FSNS_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNS_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('FSNS',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(3,2,4)
colormesh2=spatial_figure_norm(ax,FSNS_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNS_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.88, 0.37, 0.01, 0.55])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.05,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2));
char.set_label(r'$\mathrm{\mathsf{W\/m^{-2}}}$')
				
# ax=plt.subplot(4,2,7)
# spatial_figure_norm(ax,FSNSTOA_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSTOA_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
"""	
ax=plt.subplot(4,2,5)
spatial_figure_norm(ax,FSNSTOA_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSTOA_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('FSNTOA',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)	
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(4,2,6)
spatial_figure_norm(ax,FSNSTOA_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,FSNSTOA_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

"""

colorbar_min=-0.5;  colorbar_max = 0.5;	#colormap= reverse_colourmap(colormap)
# ax=plt.subplot(4,2,10)
# spatial_figure_norm(ax,LHFLX_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,LHFLX_0615_P, tb_lef=True,tb_bot=True )
# ax.annotate('(j)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
			
ax=plt.subplot(3,2,5)
spatial_figure_norm(ax,LHFLX_3150_D*0.03456,lons,lats,colormap,colorbar_min,colorbar_max,LHFLX_3150_P, tb_lef=True,tb_bot=True )
ax.annotate('Evaporation',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(3,2,6)
colormesh4 =spatial_figure_norm(ax,LHFLX_8100_D*0.03456,lons,lats,colormap,colorbar_min,colorbar_max,LHFLX_8100_P, tb_lef=False,tb_bot=True )
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.88, 0.05, 0.01, 0.27])

char = fig.colorbar(colormesh4, cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}$')

plt.subplots_adjust(left=0.1, bottom=0.03, right=0.86, top=0.95, wspace=0.05, hspace=0.05);
plt.savefig('Fig13.png', format='png', dpi=1000)
plt.show()
