# -*- coding: utf-8 -*-
'''
This function is to process thw diagnostic field into JJA mean from the CESM
monthly data
The input is the file which contains both rcp85 and fixA simulations
The output is should be normaly lons, lats, time_series and range clipped data
the lev information willl also be returned for 4D field
'''

import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as sio

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.

def Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name):
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	rcp85 = nc_fid.variables['rcp85'][:] 
	RCP85_MV = nc_fid.variables['rcp85'].missing_value
	rcp85[rcp85 == RCP85_MV] = np.nan
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:] 
	RCP85_fixA_MV = nc_fid.variables['rcp85_fixA'].missing_value
	rcp85_fixA[rcp85_fixA == RCP85_fixA_MV] = np.nan
	# lons,lats,rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85)
	# lons,lats,fixA_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85_fixA)	
	size = np.shape(rcp85)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	###JJA mean
	if (np.rank(rcp85) == 4):
		
		rcp85_JJAmean = np.empty((95,30,size[2],size[3])); rcp85_JJAmean[:] =np.nan
		fixa_JJAmean = np.empty((95,30,size[2],size[3])); fixa_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = rcp85[layer_b:layer_e+1,:,:,:]
			rcp85_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
			cache = rcp85_fixA[layer_b:layer_e+1,:,:,:]
			fixa_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
	else:
		rcp85_JJAmean = np.empty((95,size[1],size[2])); rcp85_JJAmean[:] =np.nan
		fixa_JJAmean = np.empty((95,size[1],size[2])); fixa_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = rcp85[layer_b:layer_e+1,:,:]
			rcp85_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
			cache = rcp85_fixA[layer_b:layer_e+1,:,:]
			fixa_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	nc_fid.close()
	return lon, lat, year_series,rcp85_JJAmean,fixa_JJAmean
	

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'

file_name = input_path+'ensumble_mean_PS_200602_210101.nc' 
lon, lat, year_series,PS_rcp85,PS_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# import scipy.io as sio
# lev = np.array([3.64,7.59,14.36,24.61,39.27,54.60,72.01,87.82,103.31,121.55,142.89,168.23,197.91,232.83,273.91,322.24,379.10,445.99,524.68,609.77,691.39,763.40,820.86,859.53,887.02,912.64,936.20,957.49,976.33,992.56])
output_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
sio.savemat(output_path+'CESM_2006_2100_RCP_FIXA_P.mat',\
{'time':year_series,'lon':lon,'lat':lat,\
'PS_rcp85':PS_rcp85,'PS_fixa':PS_fixa})

# file_name = input_path+'ensumble_mean_CLOUD_200602_210101.nc' 
# lon, lat, year_series,CLOUD_rcp85,CLOUD_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_CDNUMC_200602_210101.nc' 
# _, _, _,CDNC_rcp85,CDNC_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_AREL_200602_210101.nc'
# _, _, _,CDER_rcp85,CDER_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_TGCLDLWP_200602_210101.nc'
# _, _, _,CLWP_rcp85,CLWP_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)


# output_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
# sio.savemat(output_path+'CESM_2006_2100_RCP_FIXA_CLOUD_CDNC_CEDR_CLWP.mat',\
# {'time':year_series,'lon':lon,'lat':lat,\
# 'CLOUD_rcp85':CLOUD_rcp85,'CLOUD_fixa':CLOUD_fixa,\
# 'CDNC_rcp85':CDNC_rcp85,'CDNC_fixa':CDNC_fixa,\
# 'CLWP_rcp85':CLWP_rcp85,'CLWP_fixa':CLWP_fixa,\
# 'CDER_rcp85':CDER_rcp85,'CDER_fixa':CDER_fixa})


# file_name = input_path+'ensumble_mean_TS_200602_210101.nc' 
# lon, lat, year_series,TS_rcp85,TS_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_PSL_200602_210101.nc' 
# lon, lat, year_series,PSL_rcp85,PSL_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_U_200602_210101.nc'
# _, _, _,U850_rcp85,U850_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_V_200602_210101.nc'
# _, _, _,V850_rcp85,V850_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_UQ_200602_210101.nc'
# lon, lat, year_series,UQ_rcp85,UQ_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_VQ_200602_210101.nc'
# _, _, _,VQ_rcp85,VQ_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
# sio.savemat(output_path+'CESM_2006_2100_RCP_FIXA_UQ_VQ.mat',\
# {'time':year_series,'lon':lon,'lat':lat,'lev':lev,\
# 'UQ_rcp85':UQ_rcp85,'UQ_fixa':UQ_fixa,'VQ_rcp85':VQ_rcp85,'VQ_fixa':VQ_fixa})

# file_name = input_path+'ensumble_mean_Q_200602_210101.nc'
# _, _, _,Q850_rcp85,Q850_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_FSNS_200602_210101.nc'
# _, _, _,FSNS_rcp85,FSNS_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_FSNSC_200602_210101.nc'
# _, _, _,FSNSC_rcp85,FSNSC_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_FSNTOA_200602_210101.nc'
# _, _, _,FSNTOA_rcp85,FSNTOA_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_FSNTOAC_200602_210101.nc'
# _, _, _,FSNTOAC_rcp85,FSNTOAC_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# file_name = input_path+'ensumble_mean_LHFLX_200602_210101.nc'
# _, _, _,LHFLX_rcp85,LHFLX_fixa = Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# import scipy.io as sio
# output_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
# sio.savemat(output_path+'CESM_2006_2100_RCP_FIXA_TS_PSL_Q_UV(Q)_LHFLX_FSNS(C)_FSNTOA(C).mat',\
# {'time':year_series,'lon':lon,'lat':lat,\
# 'TS_rcp85':TS_rcp85,'TS_fixa':TS_fixa,'PSL_rcp85':PSL_rcp85,'PSL_fixa':PSL_fixa,\
# 'U850_rcp85':U850_rcp85,'U850_fixa':U850_fixa,'V850_rcp85':V850_rcp85,'V850_fixa':V850_fixa,\
# 'UQ850_rcp85':UQ850_rcp85,'UQ850_fixa':UQ850_fixa,'VQ850_rcp85':VQ850_rcp85,'VQ850_fixa':VQ850_fixa,\
# 'Q850_rcp85':Q850_rcp85,'Q850_fixa':Q850_fixa,\
# 'FSNS_rcp85':FSNS_rcp85,'FSNS_fixa':FSNS_fixa,'FSNSC_rcp85':FSNSC_rcp85,'FSNSC_fixa':FSNSC_fixa,\
# 'FSNTOA_rcp85':FSNTOA_rcp85,'FSNTOA_fixa':FSNTOA_fixa,'FSNTOAC_rcp85':FSNTOAC_rcp85,'FSNTOAC_fixa':FSNTOAC_fixa,\
# 'LHFLX_rcp85':LHFLX_rcp85,'LHFLX_fixa':LHFLX_fixa})

