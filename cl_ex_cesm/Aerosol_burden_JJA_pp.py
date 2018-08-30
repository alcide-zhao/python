
# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to prprocess the CESM model produced aerosol burdens data into JJA means from monthly mean 
resolution 192*288
"""
import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
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
    # 'lib'
)

site.addsitedir(lib_path)

from lib import *

#######################################
# 0.0 data input
#######################################
oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan


def JJA_mean(variable_name):
	file_name = 'ensumble_mean_'+variable_name+'_200602_210101.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]

	#RCP8.5
	rcp85 = nc_fid.variables['rcp85'][:]
	# rcp85 = np.multiply(rcp85,oceanmask)
	# time_series_att_rcp85 = stats.nanmean(stats.nanmean(att_rcp85_clipped,axis=2),axis=1)
	#FixA
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:]
	# rcp85_fixA = np.multiply(rcp85_fixA,oceanmask)
	nc_fid.close()
	
	###JJA mean
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	rcp85_jja = np.empty((95,192,288))
	rcp85_jja[:] =np.nan
	fixA_jja = np.empty((95,192,288))
	fixA_jja[:] =np.nan
	layer_e = -6  
	# print year_series
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = rcp85[layer_b:layer_e+1]
		rcp85_jja[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
		cache = rcp85_fixA[layer_b:layer_e+1]
		fixA_jja[iyear-2006] = stats.nanmean(cache,axis=0)
	return lon,lat,year_series,rcp85_jja,fixA_jja
	
	
################################################
# Aerosol burden and AOD from CESM model outputs
################################################
input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'

lon,lat,year_series,AOD_RCP85,AOD_FIXA = JJA_mean(variable_name='AODVIS')
_,_,_,BC_RCP85,BC_FIXA = JJA_mean(variable_name='BURDENBC')
_,_,_,POM_RCP85,POM_FIXA = JJA_mean('BURDENPOM')
_,_,_,SOA_RCP85,SOA_FIXA = JJA_mean('BURDENSOA')
_,_,_,SO4_RCP85,SO4_FIXA = JJA_mean('BURDENSO4')

OC_RCP85 = SOA_RCP85+POM_RCP85
OC_FIXA  = SOA_FIXA +POM_FIXA 

import scipy.io as sio
sio.savemat(input_path+'Aerosol_Burden_AOD_JJA.mat', {'time':year_series,'lon':lon,'lat':lat,\
'BC_RCP85':BC_RCP85,'BC_FIXA':BC_FIXA,'OC_RCP85':OC_RCP85,'OC_FIXA':OC_FIXA,\
'AOD_RCP85':AOD_RCP85,'AOD_FIXA':AOD_FIXA,'SO4_RCP85':SO4_RCP85,'SO4_FIXA':SO4_FIXA})
data = sio.loadmat(input_path+'Aerosol_Burden_AOD_JJA.mat')
print sio.whosmat(input_path+"Aerosol_Burden_AOD_JJA.mat")

