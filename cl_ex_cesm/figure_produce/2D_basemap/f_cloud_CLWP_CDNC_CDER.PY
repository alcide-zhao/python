# -*- coding: utf-8 -*-
"""
This is to show the variation in  CLOUD FRACTION, the cloud droplet concentration and cloud effective rdroplet radius
"""
import numpy as np
from scipy import stats
import netCDF4 as nc4
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
import scipy.io as sio
oceanmask=sio.loadmat('//home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
nc_fid.close()
_,_,ocean_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_CLOUD_CDNC_CEDR_CLWP.mat')
time_B = data['time'][0,:];lon = data['lon'][0,:];lat = data['lat'][0,:];
CDNC_rcp85 = data['CDNC_rcp85']*10**(-10);CDNC_fixa = data['CDNC_fixa']*10**(-10);
CLWP_rcp85 =data['CLWP_rcp85']*1000;CLWP_fixa = data['CLWP_fixa']*1000;
CLOUD_rcp85 = data['CLOUD_rcp85'][:,23:29,:,:]*100;CLOUD_fixa =data['CLOUD_fixa'][:,23:29,:,:]*100; #850-912hpa
CLOUD_rcp85 = np.nanmax(CLOUD_rcp85,axis=1);CLOUD_fixa = np.nanmax(CLOUD_fixa,axis=1)

CDER_rcp85 =data['CDER_rcp85'][:,23:29,:,:];CDER_fixa = data['CDER_fixa'][:,23:29,:,:]; #print np.shape(CDER_fixa)
CDER_rcp85 =  stats.nanmean(CDER_rcp85,axis=1);CDER_fixa = stats.nanmean(CDER_fixa,axis=1)
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
	# att_3150_D = np.divide(att_3150_R -att_3150_F,att_3150_F)*100
	att_3150_D = att_3150_R -att_3150_F
	att_8100_R = stats.nanmean(r85[75:95,:,:],axis=0)
	att_8100_F = stats.nanmean(fixa[75:95,:,:],axis=0)
	# att_8100_D = np.divide(att_8100_R -att_8100_F,att_8100_F)*100
	att_8100_D = att_8100_R -att_8100_F
	p_value_thres = 0.05
	att_0615_P=mannwhitneyu_test(r85[25:45,:,:], fixa[25:45,:,:])
	att_3150_P=mannwhitneyu_test(r85[25:45,:,:], fixa[25:45,:,:])
	att_8100_P=mannwhitneyu_test(r85[75:95,:,:], fixa[75:95,:,:])	
	return lons,lats,att_0615_D,att_3150_D,att_8100_D,att_0615_P,att_3150_P,att_8100_P

lons,lats,CLOUD_0615_D,CLOUD_3150_D,CLOUD_8100_D,CLOUD_0615_P,CLOUD_3150_P,CLOUD_8100_P = spatial_diff_sig(region,lon,lat,CLOUD_rcp85,CLOUD_fixa)
lons,lats,CDNC_0615_D,CDNC_3150_D,CDNC_8100_D,CDNC_0615_P,CDNC_3150_P,CDNC_8100_P = spatial_diff_sig(region,lon,lat,CDNC_rcp85,CDNC_fixa)
lons,lats,CDER_0615_D,CDER_3150_D,CDER_8100_D,CDER_0615_P,CDER_3150_P,CDER_8100_P = spatial_diff_sig(region,lon,lat,CDER_rcp85,CDER_fixa)
lons,lats,CLWP_0615_D,CLWP_3150_D,CLWP_8100_D,CLWP_0615_P,CLWP_3150_P,CLWP_8100_P = spatial_diff_sig(region,lon,lat,CLWP_rcp85,CLWP_fixa)


# now plotting
# fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(5.9, 7), facecolor='White');plot_setup(); pad=5
fig = plt.figure(facecolor='White',figsize=(5.9, 7));plot_setup();	pad= 5;colormap ='seismic'; 
#CLOUD
# colorbar_min=-10; colorbar_max = 10; 
colorbar_min=-3.0; colorbar_max = 3.0; 
# ax = plt.subplot(4,3,1)
# spatial_figure_norm(ax,CLOUD_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,CLOUD_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
# ax.annotate('2006-2015',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')

ax = plt.subplot(4,2,1)
spatial_figure_norm(ax,CLOUD_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,CLOUD_3150_P, tb_lef=True,tb_bot=False )
ax.annotate('Cloud (%)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('2031-2050',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
ax = plt.subplot(4,2,2)
ax.annotate('2081-2100',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
colormesh1 = spatial_figure_norm(ax,CLOUD_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,CLOUD_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.92, 0.73, 0.01, 0.21])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2));
# char.set_label('%')

#CLWP
# colorbar_min=-25; colorbar_max =25;  
colorbar_min=-40; colorbar_max =40; 
# ax=plt.subplot(4,2,4)
# spatial_figure_norm(ax,CLWP_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,CLWP_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(4,2,5)
spatial_figure_norm(ax,CLWP_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,CLWP_3150_P, tb_lef=True,tb_bot=False )
ax.annotate(r'CLWP $\mathrm{\mathsf{g\/m^{-2}}}$)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(4,2,6)
colormesh2=spatial_figure_norm(ax,CLWP_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,CLWP_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.92, 0.27, 0.01, 0.21])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
# char.set_label(r'$\mathrm{\mathsf{g\/m^{-2}}}$')

#CDNC
# colorbar_min=-60; colorbar_max =60; 
colorbar_min=-6; colorbar_max =6; 
# ax=plt.subplot(4,2,7)
# spatial_figure_norm(ax,CDNC_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,CDNC_0615_P, tb_lef=True,tb_bot=False )
# ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax=plt.subplot(4,2,3)
spatial_figure_norm(ax,CDNC_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,CDNC_3150_P, tb_lef=True,tb_bot=False )
ax.annotate(r'CDNC ($\mathrm{\mathsf{10^{10}\/m^{-2}}}$)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(4,2,4)
colormesh2=spatial_figure_norm(ax,CDNC_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,CDNC_8100_P, tb_lef=False,tb_bot=False )
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.92, 0.50, 0.01, 0.21])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
# char.set_label()

# CDER
# colorbar_min=-10; colorbar_max =10;	 #colormap= reverse_colourmap(colormap)
colorbar_min=-0.3; colorbar_max =0.3;	
# ax=plt.subplot(4,2,10)
# spatial_figure_norm(ax,CDER_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,CDER_0615_P, tb_lef=True,tb_bot=True )
# ax.annotate('(j)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
		
ax=plt.subplot(4,2,7)
spatial_figure_norm(ax,CDER_3150_D,lons,lats,colormap,colorbar_min,colorbar_max,CDER_3150_P, tb_lef=True,tb_bot=True )
ax.annotate('CDER (micron)',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)	
ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax=plt.subplot(4,2,8)
colormesh3 =spatial_figure_norm(ax,CDER_8100_D,lons,lats,colormap,colorbar_min,colorbar_max,CDER_8100_P, tb_lef=False,tb_bot=True )
ax.annotate('(h)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.92, 0.04, 0.01, 0.21])
char = fig.colorbar(colormesh3, cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))
# char.set_label('micron')

plt.subplots_adjust(left=0.1, bottom=0.03, right=0.9, top=0.95, wspace=0.02, hspace=0.07);
plt.savefig('Fig12.png', format='png', dpi=1000)
plt.show()
