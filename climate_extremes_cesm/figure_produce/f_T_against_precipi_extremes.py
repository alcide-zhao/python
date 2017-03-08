# -*- coding: utf-8 -*-
"""
Created on Feb 9 2017
This is to show the relationships betwen India surface Temperature and precipitation evelutions
@author: Alcide.Zhao
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as spio
import glob  
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale as scale

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)

site.addsitedir(lib_path)
# print lib_path
# from lib import ncdump

from lib import range_clip




#######################################
# 0.1 varibale definition
#######################################
time_s = 2006
time_e = 2100
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}

region = rergion_dic['ASIA'][:]


oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

rcp85_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/ensumble_mean_PEI_global_2006_2100_rcp85.nc'
fixa_file = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/ensumble_mean_PEI_global_2006_2100_fixa.nc'
TS_file = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/ensumble_mean_TS_200602_210101.nc'

# Rcp85 mean and extreme precipitations
nc_fid = nc4.Dataset(rcp85_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
time = nc_fid.variables['time'][:]
rx5day_rcp85 = nc_fid.variables['rx5day'][:]
sdii_rcp85 = nc_fid.variables['sdii'][:]
cwd_rcp85 = nc_fid.variables['cwd'][:]
cdd_rcp85 = nc_fid.variables['cdd'][:]
r99p_rcp85 = nc_fid.variables['r99p'][:]
r10_rcp85 = nc_fid.variables['r10'][:]
mean_precip_rcp85 = nc_fid.variables['mean_precip'][:]

# FixA mean and extreme precipitations
nc_fid = nc4.Dataset(fixa_file,mode='r')
rx5day_fixa = nc_fid.variables['rx5day'][:]
sdii_fixa = nc_fid.variables['sdii'][:]
cwd_fixa = nc_fid.variables['cwd'][:]
cdd_fixa = nc_fid.variables['cdd'][:]
r99p_fixa = nc_fid.variables['r99p'][:]
r10_fixa = nc_fid.variables['r10'][:]
mean_precip_fixa = nc_fid.variables['mean_precip'][:]

# CESM surface temperature
nc_fid = nc4.Dataset(TS_file,mode='r')
TS_time = nc_fid.variables['time'][:]
TS_rcp85 = nc_fid.variables['rcp85'][:]
RCP85_MV = nc_fid.variables['rcp85'].missing_value
TS_rcp85[TS_rcp85 == RCP85_MV] = np.nan
TS_fixa = nc_fid.variables['rcp85_fixA'][:]
fixA_MV = nc_fid.variables['rcp85_fixA'].missing_value
TS_fixa[TS_fixa == fixA_MV] = np.nan
units =  nc_fid.variables['rcp85_fixA'].units
# long_name =  nc_fid.variables['rcp85_fixA'].long_name
long_name = 'surface temperature'
# Extract the JJA mean from monthly data

_,_,oceanmasks = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)
size = np.shape(TS_fixa)
year= [ value/10000 for value in map(int,TS_time)]
year_series = [value for value in np.unique(year) if value <=2100]

TS_rcp85_JJAmean = np.empty((95,size[1],size[2])); TS_rcp85_JJAmean[:] =np.nan
TS_fixa_JJAmean = np.empty((95,size[1],size[2])); TS_fixa_JJAmean[:] =np.nan
layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
# print year_series
for iyear in year_series:
	layer_b = layer_e + 10
	layer_e = layer_b + 2
	cache = TS_rcp85[layer_b:layer_e+1,:,:]
	TS_rcp85_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	cache = TS_fixa[layer_b:layer_e+1,:,:]
	TS_fixa_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)

def mask_clip(oceanmask,range,data):
	data_masked = np.multiply(data,oceanmask)
	lons,lats,cache = range_clip(region[0],region[1],region[2],region[3],lon,lat,data_masked)
	cache = stats.nanmean(stats.nanmean(cache,axis =2),axis=1)
	data_pp= np.reshape(cache[:],(1,np.size(cache)))[0,:]
	return lons, lats, data_pp
	
# all the data to be pped before plotting the figures
OG_dic = {'mean_precip_rcp85':mean_precip_rcp85,'mean_precip_fixa':mean_precip_fixa,'TS_rcp85_JJAmean':TS_rcp85_JJAmean,'TS_fixa_JJAmean':TS_fixa_JJAmean,\
'rx5day_rcp85':rx5day_rcp85,'rx5day_fixa':rx5day_fixa,'cwd_rcp85':cwd_rcp85,'cwd_fixa':cwd_fixa,'cdd_rcp85':cdd_rcp85,'cdd_fixa':cdd_fixa,\
'r99p_rcp85':r99p_rcp85,'r99p_fixa':r99p_fixa,'sdii_rcp85':sdii_rcp85,'sdii_fixa':sdii_fixa}

pp_dic = {'mean_precip_rcp85':[],'mean_precip_fixa':[],'TS_rcp85_JJAmean':[],'TS_fixa_JJAmean':[],\
'rx5day_rcp85':[],'rx5day_fixa':[],'cwd_rcp85':[],'cwd_fixa':[],'cdd_rcp85':[],'cdd_fixa':[],'r99p_rcp85':[],'r99p_fixa':[],'sdii_rcp85':[],'sdii_fixa':[]}

for key in OG_dic.keys():
	_,_,pp_dic[key] = mask_clip(oceanmask,range,OG_dic[key][:])
	

plt.figure(1, facecolor='White',figsize=[12.5,5]); spst=[2,3]
# surface temperature v.s. mean precipitaton
sub_plot = 1 
ax=plt.subplot(spst[0],spst[1],sub_plot)
x = pp_dic['TS_rcp85_JJAmean']
y = pp_dic['mean_precip_rcp85']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='v',color="green",label='Rcp85 '+str(round(cc[0,1],2)))
x = pp_dic['TS_fixa_JJAmean']
y = pp_dic['mean_precip_fixa']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='o',color="red",label='FixA '+str(round(cc[0,1],2)))
ax.set_xlabel('Mean Surface Temperature (K)',fontsize=15)
ax.set_ylabel('Mean Precipitation (mm/day)',fontsize=15)
ax.set_title('Surf Temp V.S. Mean Precip',fontsize=15)
legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.9),fontsize=15)	

sub_plot = 2 
ax=plt.subplot(spst[0],spst[1],sub_plot)
x = pp_dic['TS_rcp85_JJAmean']
y = pp_dic['rx5day_rcp85'] 
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='v',color="green",label='Rcp85  '+str(round(cc[0,1],2)))
x = pp_dic['TS_fixa_JJAmean'] 
y = pp_dic['rx5day_fixa']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='o',color="red",label='FixA '+str(round(cc[0,1],2)))
ax.set_xlabel('Mean Surface Temperature (K)',fontsize=15)
ax.set_ylabel('RX5DAY (mm)',fontsize=15)
ax.set_title('Surf Temp V.S. RX5DAY',fontsize=15)
legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.9),fontsize=15)	

sub_plot = 3 
ax=plt.subplot(spst[0],spst[1],sub_plot)
x = pp_dic['TS_rcp85_JJAmean'] 
y = pp_dic['sdii_rcp85']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='v',color="green",label='Rcp85  '+str(round(cc[0,1],2)))
x = pp_dic['TS_fixa_JJAmean']
y = pp_dic['sdii_fixa']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='o',color="red",label='FixA '+str(round(cc[0,1],2)))
ax.set_xlabel('Mean Surface Temperature (K)',fontsize=15)
ax.set_ylabel('sdii (mm)',fontsize=15)
ax.set_title('Surf Temp V.S. SDII',fontsize=15)
legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.9),fontsize=15)	

sub_plot = 4 
ax=plt.subplot(spst[0],spst[1],sub_plot)
x = pp_dic['TS_rcp85_JJAmean']
y = pp_dic['cdd_rcp85']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='v',color="green",label='Rcp85  '+str(round(cc[0,1],2)))
x = pp_dic['TS_fixa_JJAmean']
y = pp_dic['cdd_fixa']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='o',color="red",label='FixA '+str(round(cc[0,1],2)))
ax.set_xlabel('Mean Surface Temperature (K)',fontsize=15)
ax.set_ylabel('CDD (days)',fontsize=15)
ax.set_title('Surf Temp V.S. CDD',fontsize=15)
legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.9),fontsize=15)	

sub_plot = 5 
ax=plt.subplot(spst[0],spst[1],sub_plot)
x = pp_dic['TS_rcp85_JJAmean']
y = pp_dic['cwd_rcp85']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='v',color="green",label='Rcp85  '+str(round(cc[0,1],2)))
x = pp_dic['TS_fixa_JJAmean']
y = pp_dic['cwd_fixa']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='o',color="red",label='FixA '+str(round(cc[0,1],2)))
ax.set_xlabel('Mean Surface Temperature (K)',fontsize=15)
ax.set_ylabel('CWD (days)',fontsize=15)
ax.set_title('Surf Temp V.S. CWD',fontsize=15)
legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.9),fontsize=15)	

sub_plot = 6 
ax=plt.subplot(spst[0],spst[1],sub_plot)
x = pp_dic['TS_rcp85_JJAmean']
y = pp_dic['r99p_rcp85']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='v',color="green",label='Rcp85  '+str(round(cc[0,1],2)))
x = pp_dic['TS_fixa_JJAmean']
y = pp_dic['r99p_fixa']
cc = np.ma.corrcoef(x,y)
ax.scatter(x,y,marker='o',color="red",label='FixA '+str(round(cc[0,1],2)))
ax.set_xlabel('Mean Surface Temperature (K)',fontsize=15)
ax.set_ylabel('r99p (days)',fontsize=15)
ax.set_title('Surf Temp V.S. R99P',fontsize=15)
legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.2, 0.9),fontsize=15)

plt.show()