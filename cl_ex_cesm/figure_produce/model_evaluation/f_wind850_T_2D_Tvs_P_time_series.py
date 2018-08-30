# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This is to plot the following fields
(a) The NCEP reanalysis 850hpa winds + surface temperature
(b) the CESM model produced 850hpa winds + surface temperature
(c) Time evelution of Mean surface temperature and pecipitation amount over the Asian Monsoon Region
(d) Time evelution of Mean surface temperature and pecipitation amount over the EASM
(e) Time evelution of Mean surface temperature and pecipitation amount over the SASM
"""
import scipy

import glob
import netCDF4 as nc4
import numpy as np
from scipy import stats
from scipy.stats import norm
from mpl_toolkits.basemap import Basemap
import math
import os; import site

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

# ocean land mask
import scipy.io as sio
oceanmask=sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720')['landocean']
ocean_mask = np.flipud(oceanmask)
cache = np.empty((720,1440))
cache[:,0:720] = ocean_mask[:,720:1440]; cache[:,720:1440] = ocean_mask[:,0:720]
ocean_mask_720_1440 = cache
rergion_dic={'ASIA':[60,150,-5,55],'EA':[110,135,22.5,45],'SA':[70,105,10,30]}
	
				
#######################################
# DATA INPUT and preprocessing 1971-2000
#######################################

file_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
file = file_path+'TS_P_850UV_model_eval.mat'
data = sio.loadmat(file)
# print data.keys()
year_series = data['time'][0,:];
lon_CESM = data['lon_CESM'][0,:];lat_CESM = data['lat_CESM'][0,:];
prep_JJA_mean_CESM = data['prep_mean_CESM'][:];  prep_JJA_std_CESM = data['prep_std_CESM'][:];
TS_JJA_mean_CESM = data['TS_JJA_mean_CESM'][:]-273.15;  TS_JJA_std_CESM = data['TS_JJA_std_CESM'][:];  
U_JJA_mean_CESM = data['U_JJA_mean_CESM'][:];  V_JJA_mean_CESM = data['V_JJA_mean_CESM'][:];  

lon_APHRO = data['lon_APHRO'][0,:]; lat_APHRO = data['lat_APHRO'][0,:]; 
prep_mean_APHRO = data['prep_mean_APHRO'][:]; temp_mean_APHRO = data['temp_mean_APHRO'][:]; 
lon_NCEP = data['lon_NCEP'][0,:];  lat_NCEP =data['lat_NCEP'][0,:]; 
TS_JJA_mean_NCEP = data['TS_JJA_mean_NCEP'][:];
U_JJA_mean_NCEP = data['U_JJA_mean_NCEP'][:];V_JJA_mean_NCEP = data['V_JJA_mean_NCEP'][:];


## Clip all the data over Asia Monsoon Region
region = rergion_dic['ASIA'][:]
lon_CESM_AM,lat_CESM_AM,prep_JJA_mean_CESM_AM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,prep_JJA_mean_CESM)
_,_,prep_JJA_std_CESM_AM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,prep_JJA_std_CESM)
_,_,TS_JJA_mean_CESM_AM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TS_JJA_mean_CESM)
_,_,TS_JJA_std_CESM_AM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TS_JJA_std_CESM)
_,_,U_JJA_mean_CESM_AM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,U_JJA_mean_CESM)
_,_,V_JJA_mean_CESM_AM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,V_JJA_mean_CESM)

lon_APHRO_AM,lat_APHRO_AM,prep_mean_APHRO_AM = range_clip(region[0],region[1],region[2],region[3],lon_APHRO,lat_APHRO,prep_mean_APHRO)
lon_NCEP_AM,lat_NCEP_AM,TS_JJA_mean_NCEP_AM = range_clip(region[0],region[1],region[2],region[3],lon_NCEP,lat_NCEP,TS_JJA_mean_NCEP)
# _,_,temp_mean_APHRO_AM = range_clip(region[0],region[1],region[2],region[3],lon_APHRO,lat_APHRO,temp_mean_APHRO)
_,_,U_JJA_mean_NCEP_AM = range_clip(region[0],region[1],region[2],region[3],lon_NCEP,lat_NCEP,U_JJA_mean_NCEP)
_,_,V_JJA_mean_NCEP_AM = range_clip(region[0],region[1],region[2],region[3],lon_NCEP,lat_NCEP,V_JJA_mean_NCEP)

# plt.imshow(temp_mean_APHRO_AM[10,:,:],origin ='lower');plt.show()
# print temp_mean_APHRO_AM[10,:,:]
################
# PLOTTING
################
import matplotlib.pyplot as plt
fig = plt.figure(facecolor='White',figsize=[7,3]);plot_setup();pad= 5;

###############################
# left panel : 2 plots
###############################

### 30year (1971-2000) mean
TS_JJA_mean_NCEP_AM_2D = stats.nanmean(TS_JJA_mean_NCEP_AM,axis =0);
# plt.imshow(temp_mean_APHRO_AM,origin ='lower');plt.show()
U_NCEP_2D_MA = stats.nanmean(U_JJA_mean_NCEP_AM,axis =0);#print np.shape(U_NCEP_2D_MA)
V_NCEP_2D_MA = stats.nanmean(V_JJA_mean_NCEP_AM,axis =0)

TS_CESM_2D_MA = stats.nanmean(TS_JJA_mean_CESM_AM,axis =0)
U_CESM_2D_MA = stats.nanmean(U_JJA_mean_CESM_AM,axis =0)
V_CESM_2D_MA = stats.nanmean(V_JJA_mean_CESM_AM,axis =0)

### interpolate  into NCEP grids
from scipy.interpolate import interp2d  as interp2d
lon_interp_AM = np.arange(60,152.5,2.5);lat_interp_AM = np.arange(-5,57.5,2.5);
f = interp2d(lon_CESM_AM, lat_CESM_AM, TS_CESM_2D_MA,kind='linear')
TS_CESM_2D_MA_interp = f(lon_interp_AM, lat_interp_AM)
f = interp2d(lon_CESM_AM, lat_CESM_AM, U_CESM_2D_MA,kind='linear')
U_CESM_2D_MA_interp = f(lon_interp_AM, lat_interp_AM)
f = interp2d(lon_CESM_AM, lat_CESM_AM, V_CESM_2D_MA,kind='linear')
V_CESM_2D_MA_interp = f(lon_interp_AM, lat_interp_AM)
# f = interp2d(lon_APHRO_AM, lat_APHRO_AM, temp_mean_APHRO_AM,kind='linear')
# temp_mean_APHRO_AM_interp = f(lon_interp_AM, lat_interp_AM);#print np.shape(temp_mean_APHRO_AM_interp)
# print temp_mean_APHRO_AM_interp

## subplot_(a) surface temperature and 850hpa winds from Observations

qk_caption = r'$\mathrm{\mathsf{10\/m\/s^{-1}}}$'; qk_scale =10
cb_min = -5;cb_max =40;colormap ='rainbow';
ax = plt.subplot2grid((1, 2), (0, 1), rowspan=3,colspan=1);
ax.annotate('(b) NCEP SAT+Winds850',xy=(0.37, 1.005), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
spatial_scaler_vector(ax,lon_NCEP_AM,lat_NCEP_AM,colormap,cb_min,cb_max,TS_JJA_mean_NCEP_AM_2D,U_NCEP_2D_MA,V_NCEP_2D_MA,qk_scale,qk_caption,qk_is = True, tb_lef=False,tb_bot=True)
## subplot_(b) surface temperature and 850hpa winds from CESM model 
ax = plt.subplot2grid((1, 2), (0, 0), rowspan=3,colspan=1);
ax.annotate('(a) CESM1 SAT+Winds850',xy=(0.37, 1.005), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
colormesh1 =spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,cb_min,cb_max,TS_CESM_2D_MA_interp,U_CESM_2D_MA_interp,V_CESM_2D_MA_interp,qk_scale,qk_caption,qk_is = False, tb_lef=True,tb_bot=True)

cbar_ax = fig.add_axes([0.25, 0.15, 0.55, 0.02])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.arange(-0.20,1.20,0.10)*(cb_max-cb_min)+cb_min)
char.set_label('Surface air temperature (degree)')

"""	
###############################
# right panel : 3 plots
###############################

# interpolate the ocean_mask to data gridcells
from scipy import arange
region = rergion_dic['ASIA'][:]
OLM_lon = arange(0,360,0.25); OLM_lat = arange(-90,90,0.25)
f = interp2d(OLM_lon, OLM_lat, ocean_mask_720_1440,kind='linear')
# the CESM model grid 
mask_cache = f(lon_CESM, lat_CESM)
mask_cache[mask_cache<125] = np.nan; mask_cache[mask_cache>=125] =1;mask_cache[0,:]=1
_,_,ocean_mask_CESM = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,mask_cache)	
## The APHRODIATE GRID
lon_global_APHRO = arange(0,360,0.5); lat_global_APHRO = arange(-90,90,0.5)
mask_cache = f(lon_global_APHRO, lat_global_APHRO)
mask_cache[mask_cache<125] = np.nan; mask_cache[mask_cache>=125] =1;mask_cache[0,:]=1
_,_,ocean_mask_APHRO = range_clip(region[0],region[1],region[2],region[3],lon_global_APHRO,lat_global_APHRO,mask_cache)	
## The NCEP GRID note the NCEP lettide is from 90 to -90 not -90 to 90
OLM_lon = arange(0,360,0.25); OLM_lat = arange(90,-90,-0.25)
f = interp2d(OLM_lon, OLM_lat, np.flipud(ocean_mask_720_1440),kind='linear')
mask_cache = f(lon_NCEP, lat_NCEP)
mask_cache[mask_cache<125] = np.nan; mask_cache[mask_cache>=125] =1;mask_cache[0,:]=1
_,_,ocean_mask_NCEP = range_clip(region[0],region[1],region[2],region[3],lon_NCEP,-lat_NCEP,mask_cache)
ocean_mask_NCEP = np.flipud(ocean_mask_NCEP)

# output_PATH ='//home/s1667168/coding/python/climate_extremes_cesm/external_data/'
# sio.savemat(output_PATH+'landoceanmask_NCEP.mat', {'landocean':ocean_mask_NCEP})



## mask off the ocean
prep_JJA_mean_CESM_AM_masked = np.multiply(prep_JJA_mean_CESM_AM,ocean_mask_CESM)
prep_JJA_std_CESM_AM_masked = np.multiply(prep_JJA_std_CESM_AM,ocean_mask_CESM)
TS_JJA_mean_CESM_AM_masked = np.multiply(TS_JJA_mean_CESM_AM,ocean_mask_CESM)
TS_JJA_std_CESM_AM_masked = np.multiply(TS_JJA_std_CESM_AM,ocean_mask_CESM)
TS_JJA_mean_NCEP_AM_masked = np.multiply(TS_JJA_mean_NCEP_AM,ocean_mask_NCEP)
prep_mean_APHRO_AM_masked = np.multiply(prep_mean_APHRO_AM,ocean_mask_APHRO[:-1,:-1])

## subplot_(c) Time evelution of surface temperature and precipitation over Asi
T_NCEP_mean = stats.nanmean(stats.nanmean(TS_JJA_mean_NCEP_AM_masked,axis=2),axis=1)
T_CESM_mean = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_AM_masked,axis=2),axis=1)
T_CESM_upper = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_AM_masked+1.96*TS_JJA_std_CESM_AM_masked,axis=2),axis=1)
T_CESM_lower = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_AM_masked-1.96*TS_JJA_std_CESM_AM_masked,axis=2),axis=1)
P_APHRO_mean = stats.nanmean(stats.nanmean(prep_mean_APHRO_AM_masked,axis=2),axis=1)
P_CESM_mean = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_AM_masked,axis=2),axis=1)
P_CESM_upper = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_AM_masked+1.96*prep_JJA_std_CESM_AM_masked,axis=2),axis=1)
P_CESM_lower = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_AM_masked-1.96*prep_JJA_std_CESM_AM_masked,axis=2),axis=1)

	
ax3 = plt.subplot2grid((6, 2), (0, 1), rowspan=3,colspan=1);plot_setup();pad =5
ax3.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax3.set_xlim([1970,2000]);ax3.set_xticklabels([]);ax3.locator_params(axis='y',nbins=4)
# ax3.set_ylim([-1.5,0.25]);
ax3.annotate('(c) Time Series',xy=(0.43, 1.005), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal')
## right y-axis the burrdens
ax3r = ax3.twinx();ax3r.locator_params(axis='y',nbins=4);ax3r.set_xlim([1971,2000])
# ax3r.set_ylim([-2.0,0.25]);
lins1 = ax3.plot(year_series,T_NCEP_mean,'-.',color="r",label='T_NCEP',linewidth=1.5)
lins2 = ax3.plot(year_series,T_CESM_mean,'-',color="r",label='T_CESM',linewidth=1.5)
lins3 = ax3.fill_between(year_series,T_CESM_upper,T_CESM_lower, where=T_CESM_upper>=T_CESM_lower, color="red", alpha = 0.2) 

lins4 = ax3r.plot(year_series,P_APHRO_mean,'-.',color="b",label='BC*10',linewidth=1.5)
lins5 = ax3r.plot(year_series,P_CESM_mean,'-',color="b",label='POA+SOA',linewidth=1.5)
lins6 = ax3r.fill_between(year_series,P_CESM_upper,P_CESM_lower, where=P_CESM_upper>=P_CESM_lower, color="blue", alpha = 0.2) 

ax3.set_ylabel('Surface ait temperature (degree)', color='k')
ax3r.set_ylabel('Daily precipitation mean 'r'$\mathrm{\mathsf{(mm\/day^{-1})}}$', color='k')


ax4 = plt.subplot2grid((6, 2), (3, 1), rowspan=3,colspan=1);plot_setup();pad =5
ax4.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
# ax3.set_xlim([1970,2000]);ax3.set_xticklabels([]);ax3.locator_params(axis='y',nbins=4)
# ax3.set_ylim([-1.5,0.25]);
ax4.annotate('(d) Precipitation against warming',xy=(0.43, 1.005), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal')
# size = np.shape(TS_JJA_mean_CESM_AM_masked)[0]* np.shape(TS_JJA_mean_CESM_AM_masked)[1]* np.shape(TS_JJA_mean_CESM_AM_masked)[2]
# T_CESM_mean=np.reshape(TS_JJA_mean_CESM_AM_masked,(size,1));
# P_CESM_mean=np.reshape(prep_JJA_mean_CESM_AM_masked,(size,1));
	
T_NCEP_mean=T_NCEP_mean-np.min(T_NCEP_mean);T_CESM_mean=T_CESM_mean-np.min(T_CESM_mean);		
P_APHRO_mean = (P_APHRO_mean-np.min(P_APHRO_mean))/np.min(P_APHRO_mean);
P_CESM_mean = (P_CESM_mean-np.min(P_CESM_mean))/np.min(P_CESM_mean)
				

lins1 = ax4.plot(T_NCEP_mean,P_APHRO_mean,'o',color="K",label='CESM',linewidth=1.5)
lins2 = ax4.plot(T_CESM_mean,P_CESM_mean,'o',color="R",label='NCEP',linewidth=1.5)


slope, intercept, r_value, p_value, std_err = stats.linregress(T_NCEP_mean,P_APHRO_mean)
ax4.annotate('Obser P= '+str(round(slope,2))+'T'+str(round(intercept,1))+' (p:'+str(round(p_value/2,2))+')',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=8)
slope, intercept, r_value, p_value, std_err = stats.linregress(T_CESM_mean,P_CESM_mean)
ax4.annotate('Model P= '+str(round(slope,2))+'T'+str(round(intercept,1))+' (p:'+str(round(p_value/2,2))+')',xy=(0.02, 0.70), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=8)

		
# align_yaxis(ax1,0,ax2,0)

## subplot_(d) Time evelution of surface temperature and precipitation over EA 
ax4 = plt.subplot2grid((6, 2), (2, 1), rowspan=2,colspan=1);pad =5
region = rergion_dic['EA'][:]
_,_,prep_JJA_mean_CESM_EA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,prep_JJA_mean_CESM_AM_masked)	
_,_,prep_JJA_std_CESM_EA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,prep_JJA_std_CESM_AM_masked)	
_,_,TS_JJA_mean_CESM_EA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,TS_JJA_mean_CESM_AM_masked)	
_,_,TS_JJA_std_CESM_EA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,TS_JJA_std_CESM_AM)	

_,_,TS_JJA_mean_NCEP_EA = range_clip(region[0],region[1],region[2],region[3],lon_NCEP_AM,lat_NCEP_AM,TS_JJA_mean_NCEP_AM_masked)	
_,_,prep_JJA_std_APHRO_EA = range_clip(region[0],region[1],region[2],region[3],lon_APHRO_AM,lat_APHRO_AM,prep_mean_APHRO_AM_masked)	


T_NCEP_mean = stats.nanmean(stats.nanmean(TS_JJA_mean_NCEP_EA,axis=2),axis=1)
T_CESM_mean = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_EA,axis=2),axis=1)
T_CESM_upper = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_EA+1.96*TS_JJA_std_CESM_EA,axis=2),axis=1)
T_CESM_lower = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_EA-1.96*TS_JJA_std_CESM_EA,axis=2),axis=1)
P_APHRO_mean = stats.nanmean(stats.nanmean(prep_JJA_std_APHRO_EA,axis=2),axis=1)
P_CESM_mean = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_EA,axis=2),axis=1)
P_CESM_upper = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_EA+1.96*prep_JJA_std_CESM_EA,axis=2),axis=1)
P_CESM_lower = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_EA-1.96*prep_JJA_std_CESM_EA,axis=2),axis=1)

# ax3 = plt.subplot2grid((6, 2), (0, 1), rowspan=3,colspan=1);
# plt.imshow(stats.nanmean(prep_JJA_mean_CESM_EA,axis=0))
# ax3 = plt.subplot2grid((6, 2), (3, 1), rowspan=3,colspan=1);
# plt.imshow(stats.nanmean(prep_JJA_std_APHRO_EA,axis=0))
# plt.show()

ax4.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax4.set_xlim([1971,2000]);ax4.set_xticklabels([]);ax4.locator_params(axis='y',nbins=4)
# ax4.set_ylim([-1.5,0.25]);
ax4.annotate('(d) East Asia (land)',xy=(0.26, 1.005), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal')
## right y-axis the burrdens
ax4r = ax4.twinx();ax4r.locator_params(axis='y',nbins=4);ax4r.set_xlim([1971,2000])
ax4.set_ylabel('Summertime mean surface tmperature (degree)', color='k')
ax4r.set_ylabel('Summertime mean precipitation amount 'r'$\mathrm{\mathsf{(mm\/day^{-1})}}$', color='k')
# ax4r.set_ylim([-2.0,0.25]);
lins1 = ax4.plot(year_series,T_NCEP_mean,'-.',color="r",label='T_NCEP',linewidth=1.5)
lins2 = ax4.plot(year_series,T_CESM_mean,'-',color="r",label='T_CESM',linewidth=1.5)
lins3 = ax4.fill_between(year_series,T_CESM_upper,T_CESM_lower, where=T_CESM_upper>=T_CESM_lower, color="red", alpha = 0.2) 

lins4 = ax4r.plot(year_series,P_APHRO_mean,'-.',color="b",label='P_APHRODITE',linewidth=1.5)
lins5 = ax4r.plot(year_series,P_CESM_mean,'-',color="b",label='P_CESM',linewidth=1.5)
lins6 = ax4r.fill_between(year_series,P_CESM_upper,P_CESM_lower, where=P_CESM_upper>=P_CESM_lower, color="blue", alpha = 0.2) 

slope, intercept, r_value, p_value, std_err = stats.linregress(T_NCEP_mean,P_APHRO_mean)
ax4.annotate('Obser P= '+str(round(slope,2))+'T'+str(round(intercept,1))+' (p:'+str(round(p_value/2,2))+')',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=8)
slope, intercept, r_value, p_value, std_err = stats.linregress(T_CESM_mean,P_CESM_mean)
ax4.annotate('Model P= '+str(round(slope,2))+'T'+str(round(intercept,1))+' (p:'+str(round(p_value/2,2))+')',xy=(0.02, 0.70), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=8)

## subplot_(e) Time evelution of surface temperature and precipitation over SA
ax5 = plt.subplot2grid((6, 2), (4, 1), rowspan=2,colspan=1)
region = rergion_dic['SA'][:]
_,_,prep_JJA_mean_CESM_SA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,prep_JJA_mean_CESM_AM_masked)	
_,_,prep_JJA_std_CESM_SA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,prep_JJA_std_CESM_AM_masked)	
_,_,TS_JJA_mean_CESM_SA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,TS_JJA_mean_CESM_AM_masked)	
_,_,TS_JJA_std_CESM_SA = range_clip(region[0],region[1],region[2],region[3],lon_CESM_AM,lat_CESM_AM,TS_JJA_std_CESM_AM)	

_,_,TS_JJA_mean_NCEP_SA = range_clip(region[0],region[1],region[2],region[3],lon_NCEP_AM,lat_NCEP_AM,TS_JJA_mean_NCEP_AM_masked)	
_,_,prep_JJA_std_APHRO_SA = range_clip(region[0],region[1],region[2],region[3],lon_APHRO_AM,lat_APHRO_AM,prep_mean_APHRO_AM_masked)	

T_NCEP_mean = stats.nanmean(stats.nanmean(TS_JJA_mean_NCEP_SA,axis=2),axis=1)
T_CESM_mean = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_SA,axis=2),axis=1)
T_CESM_upper = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_SA+1.96*TS_JJA_std_CESM_SA,axis=2),axis=1)
T_CESM_lower = stats.nanmean(stats.nanmean(TS_JJA_mean_CESM_SA-1.96*TS_JJA_std_CESM_SA,axis=2),axis=1)
P_APHRO_mean = stats.nanmean(stats.nanmean(prep_JJA_std_APHRO_SA,axis=2),axis=1)
P_CESM_mean = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_SA,axis=2),axis=1)
P_CESM_upper = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_SA+1.96*prep_JJA_std_CESM_SA,axis=2),axis=1)
P_CESM_lower = stats.nanmean(stats.nanmean(prep_JJA_mean_CESM_SA-1.96*prep_JJA_std_CESM_SA,axis=2),axis=1)

ax5.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax5.set_xlim([1971,2000]);ax5.locator_params(axis='y',nbins=4);
# ax4.set_ylim([-1.5,0.25]);
ax5.annotate('(e) South Assia (land)',xy=(0.30, 1.005), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal')
ax5.set_xlabel('Year', color='k');
## right y-axis the burrdens
ax5r = ax5.twinx();ax5r.locator_params(axis='y',nbins=4);ax5r.set_xlim([1971,2000])
# ax5r.set_ylim([-2.0,0.25]);
lins1 = ax5.plot(year_series,T_NCEP_mean,'-.',color="r",label='T_Obser',linewidth=1.5)
lins2 = ax5.plot(year_series,T_CESM_mean,'-',color="r",label='T_Model',linewidth=1.5)
lins3 = ax5.fill_between(year_series,T_CESM_upper,T_CESM_lower, where=T_CESM_upper>=T_CESM_lower, color="red", alpha = 0.2) 

lins4 = ax5r.plot(year_series,P_APHRO_mean,'-.',color="b",label='P_Obser',linewidth=1.5)
lins5 = ax5r.plot(year_series,P_CESM_mean,'-',color="b",label='P_Model',linewidth=1.5)
lins6 = ax5r.fill_between(year_series,P_CESM_upper,P_CESM_lower, where=P_CESM_upper>=P_CESM_lower, color="blue", alpha = 0.2) 

slope, intercept, r_value, p_value, std_err = stats.linregress(T_NCEP_mean,P_APHRO_mean)
ax5.annotate('Obser P='+str(round(slope,2))+'T+'+str(round(intercept,1))+' (p:'+str(round(p_value/2,2))+')',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=8)
slope, intercept, r_value, p_value, std_err = stats.linregress(T_CESM_mean,P_CESM_mean)
ax5.annotate('Model P= '+str(round(slope,2))+'T'+str(round(intercept,1))+'     (p:'+str(round(p_value/2,2))+')',xy=(0.02, 0.70), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=8)
##LEGEND
legend1 = ax5.legend(shadow=False,ncol=1,bbox_to_anchor=(0.52, 0.29))	 
legend1.get_frame().set_facecolor('gray');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)
legend2 = ax5r.legend(shadow=False,ncol=1,bbox_to_anchor=(0.8175, 0.29))	 
legend2.get_frame().set_facecolor('gray');legend2.get_frame().set_edgecolor('None');legend2.get_frame().set_alpha(0.1)
"""
plt.subplots_adjust(left=0.08, bottom=0.18, right=0.95, top=0.95, wspace=0.10, hspace=0.05);
plt.savefig('Fig1.png', format='png', dpi=1000)
plt.show()




