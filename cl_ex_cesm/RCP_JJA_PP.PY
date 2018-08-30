# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
tHIS IS TO PREPROCESS THE RCP EMISSIONS INTO JJA MEANS IN THREE SECTOES: so2,bc and oc
@author: Alcide.Zhao
"""

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
from scipy import stats
import math
# from scipy.interpolate import interp2d



input_path = '/exports/csce/datastore/geos/users/s1667168/RCP/'
def JJA_RCP_EMISSION():
	RCP_jja_emissions = {'BC':np.zeros((96,360,720)),'OC':np.zeros((96,360,720)),'SO2':np.zeros((96,360,720))}
	oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask_360_720.mat')['landoceanmask']	
	oceanmask[oceanmask==0]=np.nan
	for key in RCP_jja_emissions.keys():
		file_name ='accmip_interpolated_emissions_RCP85_'+key+'_2005_2100_0.5x0.5.nc'
		file = input_path+file_name
		nc_fid = nc4.Dataset(input_path+file_name,mode='r')
		lon = nc_fid.variables['lon'][:]
		lat = nc_fid.variables['lat'][:]
		time =range(2005,2101)
		ships = nc_fid.variables['ships'][:]
		anthropogenic = nc_fid.variables['anthropogenic'][:]
		biomass_burning = nc_fid.variables['biomass_burning'][:]
		# value = np.multiply(ships+anthropogenic+biomass_burning,oceanmask)
		value =ships+anthropogenic+biomass_burning
			###JJA mean
		
		layer_e = -5   # RCP hae the record for 2005 and begins from Jan
		# print year_series
		# yeaar_se_cache = np.empty((96))
		for iyear in range(2005,2101):
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = value[layer_b:layer_e+1,:,:]
			# print stats.nanmean(cache,axis=0)
			RCP_jja_emissions[key][iyear-2005,:,:] = stats.nanmean(cache,axis=0)
	return lon,lat,time,RCP_jja_emissions

lon,lat,time,RCP_jja_emissions = JJA_RCP_EMISSION()

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
import scipy.io as sio
sio.savemat(input_path+'RCP_jja_emissions.mat', {'time':time,'lon':lon,'lat':lat,'BC':RCP_jja_emissions['BC'],\
'OC':RCP_jja_emissions['OC'],'SO2':RCP_jja_emissions['SO2']})
data = sio.loadmat('RCP_jja_emissions.mat')
print sio.whosmat("RCP_jja_emissions.mat")









