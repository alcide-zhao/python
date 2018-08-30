# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to show the time evelution of monsoon index caculated 
from CESM model outputs under the rcp85 and fixa scenario
"""

import scipy

import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import math
import site
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



def movingaverage (values, window=5):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result


file_path= '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'

rcp85_file = file_path + 'CESM_rcp85_Monsoon_index.nc'
fixa_file = file_path + 'CESM_fixa_Monsoon_index.nc'

nc_fid = nc4.Dataset(rcp85_file,mode='r')

lon = nc_fid.variables['lon'][:]; lat = nc_fid.variables['lat'][:];time = nc_fid.variables['time'][:];
WYI_rcp85 = nc_fid.variables['WYI'][:];
DIMI_rcp85  = nc_fid.variables['DIMI'][:]; 
WNPMI_rcp85  = nc_fid.variables['WNPMI'][:];
SAMI_rcp85  = nc_fid.variables['SAMI'][:]; 
UMI_rcp85  = nc_fid.variables['UMI'][:];
ASMI_rcp85  = nc_fid.variables['ASMI'][:]; 
EASI_rcp85  = nc_fid.variables['EASI'][:];
EASJ_rcp85  = nc_fid.variables['ASMJ'][:];

nc_fid = nc4.Dataset(fixa_file,mode='r')

lon = nc_fid.variables['lon'][:]; lat = nc_fid.variables['lat'][:];time = nc_fid.variables['time'][:]; #lev = nc_fid.variables['lev'][:]; 
WYI_fixa = nc_fid.variables['WYI'][:];
DIMI_fixa  = nc_fid.variables['DIMI'][:]; 
WNPMI_fixa  = nc_fid.variables['WNPMI'][:];
SAMI_fixa  = nc_fid.variables['SAMI'][:]; 
UMI_fixa  = nc_fid.variables['UMI'][:];
ASMI_fixa  = nc_fid.variables['ASMI'][:]; 
EASI_fixa  = nc_fid.variables['EASI'][:];
EASJ_fixa  = nc_fid.variables['ASMJ'][:];
plt.figure(1, facecolor='White',figsize=[8,8])

spst =[2,1]

# ax1=plt.subplot(spst[0],spst[1],1)
# ax1.plot(time,movingaverage(ASMI_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
# ax1.plot(time,movingaverage(ASMI_fixa),'-',color="blue",linewidth=2,label='FixA')
# ax1.set_xlabel('Year',fontsize=10); ax1.set_ylabel('ASMI/(m*s-1)',fontsize=10)
# ax1.set_xlim([2005,2100])
# ax1.set_title('Asia Summer monsoon index\n(U850 (5-20N, 60-120E)) - (U850 (20-30N, 110-140E))')
# legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(0.1, 0.9),fontsize=15)	

ax1=plt.subplot(spst[0],spst[1],1)
ax1.plot(time,movingaverage(EASI_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
ax1.plot(time,movingaverage(EASI_fixa),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_ylabel('WYI/(m*s-1)',fontsize=10)
ax1.set_xlim([2005,2100])
ax1.set_title('East Asian Summer Monsoon Index \n(U85 averaged over 10-N40 110-140E)')
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(0.1, 0.9),fontsize=15)

# ax1=plt.subplot(spst[0],spst[1],3)
# ax1.plot(time,movingaverage(WYI_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
# ax1.plot(time,movingaverage(WYI_fixa),'-',color="blue",linewidth=2,label='FixA')
# ax1.set_xlabel('Year',fontsize=10); ax1.set_ylabel('ASMI/(m*s-1)',fontsize=10)
# ax1.set_xlim([2005,2100])
# ax1.set_title('Webster Yang monsoon index\nU850-U200 averaged over 0-20N, 40-110E)')

# ax1=plt.subplot(spst[0],spst[1],4)
# ax1.plot(time,movingaverage(EASJ_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
# ax1.plot(time,movingaverage(EASJ_fixa),'-',color="blue",linewidth=2,label='FixA')
# ax1.set_xlabel('Year',fontsize=10); ax1.set_ylabel('WYI/(m*s-1)',fontsize=10)
# ax1.set_xlim([2005,2100])
# ax1.set_title('East Asian Summer Jet\nU200(30-50N, 110-140E)')

ax1=plt.subplot(spst[0],spst[1],2)
ax1.plot(time,movingaverage(DIMI_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
ax1.plot(time,movingaverage(DIMI_fixa),'-',color="blue",linewidth=2,label='FixA')
ax1.set_xlabel('Year',fontsize=10); ax1.set_ylabel('DMI/(m*s-1)',fontsize=10)
ax1.set_xlim([2005,2100])
ax1.set_title('Dynamic Indian monsoon index \n(U850(5-15N,40-80E)-(U850 20-30N,70-90E))')

# ax3=plt.subplot(spst[0],spst[1],6)
# ax3.plot(time,movingaverage(WNPMI_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
# ax3.plot(time,movingaverage(WNPMI_fixa),'-',color="blue",linewidth=2,label='FixA')
# ax3.set_xlabel('Year',fontsize=10); ax1.set_ylabel('EAWNPMI/(m*s-1)',fontsize=10)
# ax3.set_xlim([2005,2100])
# ax3.set_title('Western North Pacific Monsoon Index \n(U850 (5-15N, 100-130E) - U850 (20-30N, 110-140E))')




# ax4=plt.subplot(spst[0],spst[1],4)
# ax4.plot(time,movingaverage(SAMI_rcp85),'-',color="red",linewidth=2,label='RCP8.5')
# ax4.plot(time,movingaverage(SAMI_fixa),'-',color="blue",linewidth=2,label='FixA')
# ax4.set_xlabel('Year',fontsize=10); ax1.set_ylabel('SAMI/(m*s-1)',fontsize=10)
# ax4.set_xlim([2005,2100])
# ax4.set_title('South Asian Monsoon Index \n(U850 averaged over 5-22.5N, 35-97.5E)')


# ax5=plt.subplot(spst[0],spst[1],5)
# region=[110,140,10,40]
# UMI_rcp85_850 = UMI_rcp85[:,:,:]; UMI_fixa_850 = UMI_fixa[:,:,:]; 
# _,_,UMI_rcp85_850 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UMI_rcp85_850)
# _,_,UMI_fixa_850 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UMI_fixa_850)
# ax5.plot(time,stats.nanmean(stats.nanmean(UMI_rcp85_850,axis =2),axis =1),'-',color="red",linewidth=2,label='RCP8.5')
# ax5.plot(time,stats.nanmean(stats.nanmean(UMI_fixa_850,axis =2),axis =1),'-',color="blue",linewidth=2,label='FixA')
# ax5.set_xlabel('Year',fontsize=10); ax1.set_ylabel('UMI/(m*s-1)',fontsize=10)
# ax5.set_xlim([2005,2100])
# ax5.set_title('East Asian monsoon index over 10-40N,110-140E')

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.03, hspace=0.3);
plt.show()