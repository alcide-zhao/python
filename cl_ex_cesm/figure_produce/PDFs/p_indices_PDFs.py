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
from numpy import linspace
from scipy.stats.kde import gaussian_kde
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


oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
nc_fid.close()
# _,_,ocean_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)


#############################################
### aerosol differebces
#############################################
	
def kde_bins(region, variable):
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
	
	att_3150_R = att_clipped_r85[25:45,:,:]
	att_3150_F = att_clipped_r85F[25:45,:,:]
	dis_min =min(np.nanmin(att_3150_R[:]),np.nanmin(att_3150_F[:]))
	dis_max =max(np.nanmax(att_3150_R[:]),np.nanmax(att_3150_F[:]))
	size =np.shape(att_3150_R)
	hist_CESM = np.reshape(att_3150_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds3150 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds3150)
	size =np.shape(att_3150_F)
	hist_CESM = np.reshape(att_3150_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM ); kde_F = kde(ds3150)
	kde_3150 = (kde_R -kde_F)*100
	
	
	att_8100_R = att_clipped_r85[75:95,:,:]
	att_8100_F = att_clipped_r85F[75:95,:,:]
	dis_min = min(np.nanmin(att_8100_R[:]),np.nanmin(att_8100_F[:]))
	dis_max = max(np.nanmax(att_8100_R[:]),np.nanmax(att_8100_F[:]))
	size =np.shape(att_8100_R)
	hist_CESM = np.reshape(att_8100_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds8100 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds8100)
	size =np.shape(att_8100_F)
	hist_CESM = np.reshape(att_8100_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM );kde_F = kde(ds8100)
	kde_8100 = (kde_R -kde_F)*100

	return ds3150,ds8100,kde_3150,kde_8100
"""
def kde_bins_r95pr(region):
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	#fixa
	file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = 100*np.divide(nc_fid.variables['r95p'][:],nc_fid.variables['total_precip'][:])
	att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	#RCP8.5
	file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	att_value = 100*np.divide(nc_fid.variables['r95p'][:],nc_fid.variables['total_precip'][:])
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	
	att_3150_R = att_clipped_r85[25:45,:,:]
	att_3150_F = att_clipped_r85F[25:45,:,:]
	dis_min =min(np.nanmin(att_3150_R[:]),np.nanmin(att_3150_F[:]))
	dis_max =max(np.nanmax(att_3150_R[:]),np.nanmax(att_3150_F[:]))
	size =np.shape(att_3150_R)
	hist_CESM = np.reshape(att_3150_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds3150 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds3150)
	size =np.shape(att_3150_F)
	hist_CESM = np.reshape(att_3150_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM ); kde_F = kde(ds3150)
	kde_3150 = (kde_R -kde_F)*100
	
	att_8100_R = att_clipped_r85[75:95,:,:]
	att_8100_F = att_clipped_r85F[75:95,:,:]
	dis_min = min(np.nanmin(att_8100_R[:]),np.nanmin(att_8100_F[:]))
	dis_max = max(np.nanmax(att_8100_R[:]),np.nanmax(att_8100_F[:]))
	size =np.shape(att_8100_R)
	hist_CESM = np.reshape(att_8100_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds8100 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds8100)
	size =np.shape(att_8100_F)
	hist_CESM = np.reshape(att_8100_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM );kde_F = kde(ds8100)
	kde_8100 = (kde_R -kde_F)*100

	return ds3150,ds8100,kde_3150,kde_8100
"""

region = rergion_dic['SA'][:]
# r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins_r95pr(region)
cwd_ds3150,cwd_ds8100,cwd_kde_3150,cwd_kde_8100 = kde_bins(region, variable = 'cwd')
sdii_ds3150,sdii_ds8100,sdii_kde_3150,sdii_kde_8100 = kde_bins(region, variable = 'sdii')
rx5day_ds3150,rx5day_ds8100,rx5day_kde_3150,rx5day_kde_8100= kde_bins(region, variable = 'rx5day')
ptot_ds3150,ptot_ds8100,ptot_kde_3150,ptot_kde_8100= kde_bins(region, variable = 'total_precip')
r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins(region,variable = 'r95p')	

fig = plt.figure(figsize=(8, 8), facecolor='White');plot_setup()

ax1 = plt.subplot(5,2,1);ax1.axhline(0, color='k',alpha=1);pad= 5;
# plt.plot(ptot_ds3150/92, ptot_kde_3150,'k--',alpha=1,label='AA_3150')
ax1.plot(ptot_ds8100/92, ptot_kde_8100,'r--',alpha=1,label='Rcp8.5-Rcp8.5_FixA')
ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax1.set_xlim([0,18]); ax1.set_xticks(np.arange(0,18.1,6));
ax1.set_ylim([-0.08,0.08]); ax1.set_yticks(np.arange(-0.08,0.081,0.04))
ax1.annotate(r'(a) PMEAN ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.annotate('South Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax1.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)
		
ax3 = plt.subplot(5,2,3);ax3.axhline(0, color='k',alpha=1);
ax3.plot(rx5day_ds8100, rx5day_kde_8100,'r--',alpha=1)
ax3.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax3.set_xlim([0,450]); ax3.set_xticks(np.arange(0,450.1,150));
ax3.set_ylim([-0.6,0.60]);ax3.set_yticks(np.arange(-0.6,0.61,0.3))
ax3.annotate('(b) RX5DAY (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
ax3.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)	
			
ax5 = plt.subplot(5,2,5);ax5.axhline(0, color='k',alpha=1);
ax5.plot(sdii_ds8100, sdii_kde_8100,'r--',alpha=1)
ax5.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax5.set_xlim([0,18]); ax5.set_xticks(np.arange(0,18.1,6));
ax5.set_ylim([-10,10]);ax5.set_yticks(np.arange(-10,10.1,5))
ax5.annotate(r'(c) SDII ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)							
ax5.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)	
				
ax7 = plt.subplot(5,2,7);ax7.axhline(0, color='k',alpha=1);
ax7.plot(cwd_ds8100, cwd_kde_8100,'r--',alpha=1)
ax7.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax7.set_xlim([0,90]); ax7.set_xticks(np.arange(0,90.1,30));
ax7.set_ylim([-0.8,0.8]);ax7.set_yticks(np.arange(-0.8,0.81,0.4))
ax7.annotate('(d) CWD (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
ax7.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)	
				
				
ax9 = plt.subplot(5,2,9);ax9.axhline(0, color='k',alpha=1);
ax9.plot(r95pr_ds8100, r95pr_kde_8100,'r--',alpha=1)
ax9.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
# ax9.set_xlim([0,60]); ax9.set_xticks(np.arange(0,60.1,20));
# ax9.set_ylim([-4,4]);ax9.set_yticks(np.arange(-4,4.1,2))
ax9.annotate('(e) R95P (%)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)		
ax9.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)	
	
	
	
region = rergion_dic['EA'][:]	
cwd_ds3150,cwd_ds8100,cwd_kde_3150,cwd_kde_8100 = kde_bins(region, variable = 'cwd')
sdii_ds3150,sdii_ds8100,sdii_kde_3150,sdii_kde_8100 = kde_bins(region, variable = 'sdii')
rx5day_ds3150,rx5day_ds8100,rx5day_kde_3150,rx5day_kde_8100= kde_bins(region, variable = 'rx5day')
ptot_ds3150,ptot_ds8100,ptot_kde_3150,ptot_kde_8100= kde_bins(region, variable = 'total_precip')	
r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins(region,variable = 'r95p')

		
ax2 = plt.subplot(5,2,2);ax2.axhline(0, color='k',alpha=1);
ax2.plot(ptot_ds8100/92, ptot_kde_8100,'r--',alpha=1)
ax2.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax2.set_xlim([0,15]); ax2.set_xticks(np.arange(0,15.1,5));
ax2.set_ylim([-0.08,0.08]); ax2.set_yticks(np.arange(-0.1,0.11,0.05))
ax2.annotate(r'(g) PMEAN ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax2.annotate('East Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')

							
ax4 = plt.subplot(5,2,4);ax4.axhline(0, color='k',alpha=1);
ax4.plot(rx5day_ds8100, rx5day_kde_8100,'r--',alpha=1)
ax4.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax4.set_xlim([0,240]);ax4.set_xticks(np.arange(0.0,250.1,80))
ax4.set_ylim([-0.8,0.8]);ax4.set_yticks(np.arange(-0.8,0.81,0.4));
ax4.annotate('(g) RX5DAY (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
					
ax6 = plt.subplot(5,2,6);ax6.axhline(0, color='k',alpha=1);
ax6.plot(sdii_ds8100, sdii_kde_8100,'r--',alpha=1)
ax6.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax6.set_xlim([0,15]);ax6.set_xticks(np.arange(0,15.1,5));
ax6.set_ylim([-15,15]);ax6.set_yticks(np.arange(-15,15.1,7.5))
ax6.annotate(r'(h) SDII ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)				
				
ax8 = plt.subplot(5,2,8);ax8.axhline(0, color='k',alpha=1)	;
ax8.plot(cwd_ds8100, cwd_kde_8100,'r--',alpha=1)
ax8.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax8.set_xlim([0,24]);ax8.set_xticks(np.arange(0,24.1,8));
ax8.set_ylim([-6,6]);ax8.set_yticks(np.arange(-6,6.1,3))
ax8.annotate('(i) CWD (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)					
				
				
ax10 = plt.subplot(5,2,10);ax10.axhline(0, color='k',alpha=1);	
ax10.plot(r95pr_ds8100, r95pr_kde_8100,'r--',alpha=1)
ax10.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
# ax10.set_xlim([0,60]);ax10.set_xticks(np.arange(0,60.1,20));
# ax10.set_ylim([-3,3]);ax10.set_yticks(np.arange(-3,3.1,1.5))
ax10.annotate('(j) R95P (%)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)			

######################################
#future-past
####################################
def kde_bins(region, variable):
	from scipy.stats import ttest_ind as test	
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	#fixa
	file_name = 'historical_ensumble_mean_PEI_global_1920_2005.nc'
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
	
	att_3150_R = att_clipped_r85[25:45,:,:]
	att_3150_F = att_clipped_r85F[66:86,:,:]
	dis_min =min(np.nanmin(att_3150_R[:]),np.nanmin(att_3150_F[:]))
	dis_max =max(np.nanmax(att_3150_R[:]),np.nanmax(att_3150_F[:]))
	size =np.shape(att_3150_R)
	hist_CESM = np.reshape(att_3150_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds3150 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds3150)
	size =np.shape(att_3150_F)
	hist_CESM = np.reshape(att_3150_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM ); kde_F = kde(ds3150)
	kde_3150 =(kde_R -kde_F)*100
	
	att_8100_R = att_clipped_r85[75:95,:,:]
	att_8100_F = att_clipped_r85F[66:86,:,:]
	dis_min = min(np.nanmin(att_8100_R[:]),np.nanmin(att_8100_F[:]))
	dis_max = max(np.nanmax(att_8100_R[:]),np.nanmax(att_8100_F[:]))
	size =np.shape(att_8100_R)
	hist_CESM = np.reshape(att_8100_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds8100 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds8100)
	size =np.shape(att_8100_F)
	hist_CESM = np.reshape(att_8100_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM );kde_F = kde(ds8100)
	kde_8100 = (kde_R -kde_F)*100

	return ds3150,ds8100,kde_3150,kde_8100
"""
def kde_bins_r95pr(region):
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	#fixa
	file_name = 'historical_ensumble_mean_PEI_global_1920_2005.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = 100*np.divide(nc_fid.variables['r95p'][:],nc_fid.variables['total_precip'][:])
	att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	#RCP8.5
	file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	att_value = 100*np.divide(nc_fid.variables['r95p'][:],nc_fid.variables['total_precip'][:])
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	
	att_3150_R = att_clipped_r85[25:45,:,:]
	att_3150_F = att_clipped_r85F[25:45,:,:]
	dis_min =min(np.nanmin(att_3150_R[:]),np.nanmin(att_3150_F[:]))
	dis_max =max(np.nanmax(att_3150_R[:]),np.nanmax(att_3150_F[:]))
	size =np.shape(att_3150_R)
	hist_CESM = np.reshape(att_3150_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds3150 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds3150)
	size =np.shape(att_3150_F)
	hist_CESM = np.reshape(att_3150_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM ); kde_F = kde(ds3150)
	kde_3150 = (kde_R -kde_F)*100
	
	
	att_8100_R = att_clipped_r85[75:95,:,:]
	att_8100_F = att_clipped_r85F[75:95,:,:]
	dis_min = min(np.nanmin(att_8100_R[:]),np.nanmin(att_8100_F[:]))
	dis_max = max(np.nanmax(att_8100_R[:]),np.nanmax(att_8100_F[:]))
	size =np.shape(att_8100_R)
	hist_CESM = np.reshape(att_8100_R,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds8100 = linspace( dis_min,dis_max, 100)
	kde = gaussian_kde( hist_CESM );kde_R = kde(ds8100)
	size =np.shape(att_8100_F)
	hist_CESM = np.reshape(att_8100_F,(size[0]*size[1]*size[2],1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	kde = gaussian_kde( hist_CESM );kde_F = kde(ds8100)
	kde_8100 = (kde_R -kde_F)*100

	return ds3150,ds8100,kde_3150,kde_8100
"""
# INDICES REPRODUCING
region = rergion_dic['EA'][:]	



cwd_ds3150,cwd_ds8100,cwd_kde_3150,cwd_kde_8100 = kde_bins(region, variable = 'cwd')
sdii_ds3150,sdii_ds8100,sdii_kde_3150,sdii_kde_8100 = kde_bins(region, variable = 'sdii')
rx5day_ds3150,rx5day_ds8100,rx5day_kde_3150,rx5day_kde_8100= kde_bins(region, variable = 'rx5day')
ptot_ds3150,ptot_ds8100,ptot_kde_3150,ptot_kde_8100= kde_bins(region, variable = 'total_precip')	
r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins(region,variable = 'r95p')

ax2.plot(ptot_ds8100/92, ptot_kde_8100,'r-',alpha=1)
ax4.plot(rx5day_ds8100, rx5day_kde_8100,'r-',alpha=1)
ax6.plot(sdii_ds8100, sdii_kde_8100,'r-',alpha=1)						
ax8.plot(cwd_ds8100, cwd_kde_8100,'r-',alpha=1)			
ax10.plot(r95pr_ds8100, r95pr_kde_8100,'r-',alpha=1)
# r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins_r95pr(region)
# plt.plot(ptot_ds3150/92, ptot_kde_3150,'k-',alpha=1)		
# plt.plot(rx5day_ds3150, rx5day_kde_3150,'k-',alpha=1)
# plt.plot(sdii_ds3150, sdii_kde_3150,'k-',alpha=1)			
# plt.plot(cwd_ds3150, cwd_kde_3150,'k-',alpha=1)
# plt.plot(r95pr_ds3150, r95pr_kde_3150,'k-',alpha=1)			

region = rergion_dic['SA'][:]

cwd_ds3150,cwd_ds8100,cwd_kde_3150,cwd_kde_8100 = kde_bins(region, variable = 'cwd')
sdii_ds3150,sdii_ds8100,sdii_kde_3150,sdii_kde_8100 = kde_bins(region, variable = 'sdii')
rx5day_ds3150,rx5day_ds8100,rx5day_kde_3150,rx5day_kde_8100= kde_bins(region, variable = 'rx5day')
ptot_ds3150,ptot_ds8100,ptot_kde_3150,ptot_kde_8100= kde_bins(region, variable = 'total_precip')	
r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins(region,variable = 'r95p')

ax1.plot(ptot_ds8100/92, ptot_kde_8100,'r-',alpha=1,label='RCP8.5')
ax3.plot(rx5day_ds8100, rx5day_kde_8100,'r-',alpha=1)
ax5.plot(sdii_ds8100, sdii_kde_8100,'r-',alpha=1)
ax7.plot(cwd_ds8100, cwd_kde_8100,'r-',alpha=1)
ax9.plot(r95pr_ds8100, r95pr_kde_8100,'r-',alpha=1)

# r95pr_ds3150,r95pr_ds8100,r95pr_kde_3150,r95pr_kde_8100= kde_bins_r95pr(region)
# plt.plot(ptot_ds3150/92, ptot_kde_3150,'k-',alpha=1,label='ALL_3150')
# plt.plot(rx5day_ds3150, rx5day_kde_3150,'k-',alpha=1)
# plt.plot(sdii_ds3150, sdii_kde_3150,'k-',alpha=1)
# plt.plot(cwd_ds3150, cwd_kde_3150,'k-',alpha=1)
# plt.plot(r95pr_ds3150, r95pr_kde_3150,'k-',alpha=1)


############################
# Bselines 1986-2005

############################
from scipy.interpolate import interp2d  as interp2d
import scipy.io as spio
oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720')['landocean']
ocean_mask = np.flipud(oceanmask)
cache = np.empty((720,1440))
cache[:,0:720] = ocean_mask[:,720:1440]; cache[:,720:1440] = ocean_mask[:,0:720]
ocean_mask_720_1440 = cache
OLM_lon = np.arange(0,360,0.25); OLM_lat = np.arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440, kind='linear')
ocean_mask_192_288 = f(lon, lat)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1



region = rergion_dic['SA'][:]
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical2.0/mean/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][66:86]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:86,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][66:86,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][66:86,:,:],axis=0);
sdii = stats.nanmean(nc_fid.variables['sdii'][66:86,:,:],axis=0);
ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
r95p= 100*np.divide(r95p,ptot)


RX5DAY_3D_masked = np.multiply(RX5DAY,ocean_mask_192_288)
lons,lats,RX5DAY_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,RX5DAY_3D_masked)
sdii_3D_masked = np.multiply(sdii,ocean_mask_192_288)
lons,lats,sdii_masked_clipped= range_clip(region[0],region[1],region[2],region[3],lon,lat,sdii_3D_masked)
cwd_3D_masked = np.multiply(cwd,ocean_mask_192_288)
lons,lats,cwd_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cwd_3D_masked)
r95p_3D_masked = np.multiply(r95p,ocean_mask_192_288)
lons,lats,r95p_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r95p_3D_masked)
ptot_3D_masked = np.multiply(ptot/92,ocean_mask_192_288)
lons,lats,ptot_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,ptot_3D_masked)


ax11 = ax1.twinx();
size =np.shape(ptot_masked_clipped)
hist_CESM = np.reshape(ptot_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 18, 100)
kde = gaussian_kde( hist_CESM )
ax11.plot(dist_space, 100*kde(dist_space),'k-',alpha=1, label = 'Baseline' )
ax11.set_xlim([0,18]);ax11.set_xticks(np.arange(0,18.1,6));
ax11.set_ylim([0,16]);ax11.set_yticks(np.arange(0,16.1,4))
align_yaxis(ax1,-0.08,ax11,0)


ax33 = ax3.twinx();
size =np.shape(RX5DAY_masked_clipped)
hist_CESM = np.reshape(RX5DAY_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 450, 100)
kde = gaussian_kde( hist_CESM )
ax33.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
ax33.set_xlim([0,450]);ax33.set_xticks(np.arange(0,450.1,150));
ax33.set_ylim([0,0.88]);ax33.set_yticks(np.arange(0,0.89,0.22))

align_yaxis(ax3,-0.6,ax33,0)

ax55 = ax5.twinx();
size =np.shape(sdii_masked_clipped)
hist_CESM = np.reshape(sdii_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 20, 100)
kde = gaussian_kde( hist_CESM )
ax55.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
ax55.set_xlim([0,18]);ax55.set_xticks(np.arange(0,18.1,6));
ax55.set_ylim([0,16]);ax55.set_yticks(np.arange(0,16.1,4));
align_yaxis(ax5,-10,ax55,0)

ax77 = ax7.twinx();
size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 90, 100)
kde = gaussian_kde( hist_CESM )
ax77.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
ax77.set_xlim([0,90]);ax77.set_xticks(np.arange(0,90.1,30));
ax77.set_ylim([0,2.4]);ax77.set_yticks(np.arange(0,2.41,0.6));
align_yaxis(ax7,-0.8,ax77,0)

ax99 = ax9.twinx();
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 80, 100)
kde = gaussian_kde( hist_CESM )
ax99.plot(dist_space, 100*kde(dist_space),'k-',alpha=1 )
ax99.set_xlim([0,75]);ax99.set_xticks(np.arange(0,75.1,25));
ax99.set_ylim([0,8]);ax99.set_yticks(np.arange(0,8.1,2));
align_yaxis(ax9,-4,ax99,0)


oceanmask=spio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720')['landocean']
ocean_mask = np.flipud(oceanmask)
cache = np.empty((720,1440))
cache[:,0:720] = ocean_mask[:,720:1440]; cache[:,720:1440] = ocean_mask[:,0:720]
ocean_mask_720_1440 = cache
OLM_lon = np.arange(0,360,0.25); OLM_lat = np.arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440, kind='linear')
ocean_mask_192_288 = f(lon, lat)
ocean_mask_192_288[ocean_mask_192_288<125] = np.nan; ocean_mask_192_288[ocean_mask_192_288>=125] =1;ocean_mask_192_288[0,:]=1

region = rergion_dic['EA'][:]
CESM_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical2.0/mean/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(CESM_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][66:86]
RX5DAY = stats.nanmean(nc_fid.variables['rx5day'][66:86,:,:],axis=0);
r95p =stats.nanmean( nc_fid.variables['r95p'][66:86,:,:],axis=0);
cwd = stats.nanmean(nc_fid.variables['cwd'][66:86,:,:],axis=0);
sdii = stats.nanmean(nc_fid.variables['sdii'][66:86,:,:],axis=0);
ptot = stats.nanmean(nc_fid.variables['total_precip'][66:86,:,:],axis=0);
r95p= 100*np.divide(r95p,ptot)


RX5DAY_3D_masked = np.multiply(RX5DAY,ocean_mask_192_288)
lons,lats,RX5DAY_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,RX5DAY_3D_masked)
sdii_3D_masked = np.multiply(sdii,ocean_mask_192_288)
lons,lats,sdii_masked_clipped= range_clip(region[0],region[1],region[2],region[3],lon,lat,sdii_3D_masked)
cwd_3D_masked = np.multiply(cwd,ocean_mask_192_288)
lons,lats,cwd_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,cwd_3D_masked)
r95p_3D_masked = np.multiply(r95p,ocean_mask_192_288)
lons,lats,r95p_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,r95p_3D_masked)
ptot_3D_masked = np.multiply(ptot/92,ocean_mask_192_288)
lons,lats,ptot_masked_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,ptot_3D_masked)


ax22 = ax2.twinx();
size =np.shape(ptot_masked_clipped)
hist_CESM = np.reshape(ptot_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 15, 100)
kde = gaussian_kde( hist_CESM )
ax22.plot(dist_space, 100*kde(dist_space),'k-',alpha=1 )
ax22.set_xlim([0,15]);ax22.set_xticks(np.arange(0,15.1,5));
ax22.set_ylim([0,24]);ax22.set_yticks(np.arange(0,24.1,6))
align_yaxis(ax2,-0.10,ax22,0)


ax44 = ax4.twinx();
size =np.shape(RX5DAY_masked_clipped)
hist_CESM = np.reshape(RX5DAY_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 300, 100)
kde = gaussian_kde( hist_CESM )
ax44.plot(dist_space, 100*kde(dist_space),'k-',alpha=1 )
ax44.set_xlim([0,270]);ax44.set_xticks(np.arange(0,270.1,90));
ax44.set_ylim([0,1.2]);ax44.set_yticks(np.arange(0,1.21,0.3));
align_yaxis(ax4,-0.8,ax44,0)

ax66 = ax6.twinx();
size =np.shape(sdii_masked_clipped)
hist_CESM = np.reshape(sdii_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 21, 100)
kde = gaussian_kde( hist_CESM )
ax66.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
ax66.set_xlim([0,18]);ax66.set_xticks(np.arange(0,18.1,6));
ax66.set_ylim([0,24]);ax66.set_yticks(np.arange(0,24.1,6));
align_yaxis(ax6,-15,ax66,0)

ax88 = ax8.twinx();
size =np.shape(cwd_masked_clipped)
hist_CESM = np.reshape(cwd_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 42, 100)
kde = gaussian_kde( hist_CESM )
ax88.plot(dist_space, 100*kde(dist_space),'k-',alpha=1 )
ax88.set_xlim([0,30]);ax88.set_xticks(np.arange(0,30.1,10));
ax88.set_ylim([0,10]);ax88.set_yticks(np.arange(0,10.1,2.5));
align_yaxis(ax8,-6,ax88,0)

ax00 = ax10.twinx();
size =np.shape(r95p_masked_clipped)
hist_CESM = np.reshape(r95p_masked_clipped,(size[0]*size[1],1))
hist_CESM [hist_CESM <=0]= np.nan
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]								
dist_space = linspace( 0, 80, 100)
kde = gaussian_kde( hist_CESM )
ax00.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
ax00.set_xlim([0,60]);ax00.set_xticks(np.arange(0,60.1,20));
ax00.set_ylim([0,8]);ax00.set_yticks(np.arange(0,8.1,2));
align_yaxis(ax10,-3,ax00,0)

legend1 = ax1.legend(shadow=False,ncol=1,bbox_to_anchor=(0.53, 0.93))	 
legend1.get_frame().set_facecolor('None');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.3)
legend11 = ax11.legend(shadow=False,ncol=1,bbox_to_anchor=(0.32, 1.05))	 
legend11.get_frame().set_facecolor('None');legend11.get_frame().set_edgecolor('None');legend11.get_frame().set_alpha(0.1)


							
plt.subplots_adjust(left=0.1, bottom=0.03, right=0.95, top=0.95, wspace=0.25, hspace=0.2);
plt.savefig('Fig14.png', format='png', dpi=1000)
plt.show()