
region_dic = {'NorthIndia':[70,90,22,31],'SouthIndia':[70,90,8,22],'SouthEastChina':[105,125,20,30],'CentralEastChina':[105,125,30,40]}
import netCDF4 as nc4
import numpy as np
from scipy import stats

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio
ocean_mask_NCEP = sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720_0_360_-90_90')['landoceanmask']
ocean_mask_NCEP[ocean_mask_NCEP==0]=np.nan
# import matplotlib.pyplot as plt

# OLM_lon = np.arange(0,360,0.25); OLM_lat = np.arange(90,-90,-0.25)
# f = interp2d(OLM_lon, OLM_lat, np.flipud(ocean_mask_720_1440),kind='linear')
# mask_cache = f(lon, lat)
# mask_cache[mask_cache<125] = np.nan; mask_cache[mask_cache>=125] =1;mask_cache[0,:]=1
# ocean_mask_NCEP = np.flipud(mask_cache)
# output_PATH ='//home/s1667168/coding/python/climate_extremes_cesm/external_data/'
# sio.savemat(output_PATH+'landoceanmask_NCEP.mat', {'landocean':ocean_mask_NCEP})
# plt.imshow(ocean_mask_NCEP);plt.show()

def mask_match(country_key,lon,lat):
	"""
	Readin the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	from scipy.interpolate import interp2d  as interp2d
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/countries_mask_lon_lat_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask= np.flipud(ocean_mask[country_key][:])
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
	mask = f(lon, lat); 
	mask[mask > 0] = 1;
	mask[mask < 1] = np.nan
	return mask

def percentile_JJAS_temp(variable,country_key,region_box):
	region = region_dic[region_box][:]
	file_name =  '/exports/csce/datastore/geos/users/s1667168/OBSERVATIONS/NCEP_'+variable+'.2m.gauss.1971_2001.nc'
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	mask = mask_match(country_key,lon,lat);
	# plt.imshow(mask);plt.show()
	value = nc_fid.variables[variable][:,151:273,:,:];
	variable_MV = nc_fid.variables[variable].missing_value; value[value == variable_MV] = np.nan;#print variable_MV
	variable_valid_range = nc_fid.variables[variable].valid_range; #print variable_valid_range
	value[value>variable_valid_range[1] ] = np.nan; value[value<variable_valid_range[0]] = np.nan;
	value = np.multiply(value,mask)-273.15
	# plt.imshow(value[1,1,:,:]);plt.show()
	_,_,value_range = range_clip(region[0],region[1],region[2],region[3],lon,lat,value);#print np.shape(value_range)
	value_range_TiSe = stats.nanmean(stats.nanmean(value_range,axis=3),axis=2);#print value_range_TiSe
	length = np.shape(value_range_TiSe)[0]*np.shape(value_range_TiSe)[1]
	obkject = np.reshape(value_range_TiSe,(length,1))
	p50,bins,pdfs=percentile(obkject,percent=50);
	return p50,bins,pdfs
SEC_bins_pdfs = np.zeros((4,100));CEC_bins_pdfs = np.zeros((4,100));
SID_bins_pdfs = np.zeros((4,100));NID_bins_pdfs = np.zeros((4,100));

variable = 'tmin'; _,SEC_bins_pdfs[0,:],SEC_bins_pdfs[1,:]=percentile_JJAS_temp(variable,country_key = 'China',region_box = 'SouthEastChina');
variable = 'tmin'; _,CEC_bins_pdfs[0,:],CEC_bins_pdfs[1,:]=percentile_JJAS_temp(variable,country_key = 'China',region_box = 'CentralEastChina');
variable = 'tmin'; _,NID_bins_pdfs[0,:],NID_bins_pdfs[1,:]=percentile_JJAS_temp(variable,country_key = 'India',region_box = 'NorthIndia');
variable = 'tmin'; _,SID_bins_pdfs[0,:],SID_bins_pdfs[1,:]=percentile_JJAS_temp(variable,country_key = 'India',region_box = 'SouthIndia');

variable = 'tmax'; _,SEC_bins_pdfs[2,:],SEC_bins_pdfs[3,:]=percentile_JJAS_temp(variable,country_key = 'China',region_box = 'SouthEastChina');
variable = 'tmax'; _,CEC_bins_pdfs[2,:],CEC_bins_pdfs[3,:]=percentile_JJAS_temp(variable,country_key = 'China',region_box = 'CentralEastChina');
variable = 'tmax'; _,NID_bins_pdfs[2,:],NID_bins_pdfs[3,:]=percentile_JJAS_temp(variable,country_key = 'India',region_box = 'NorthIndia');
variable = 'tmax'; _,SID_bins_pdfs[2,:],SID_bins_pdfs[3,:]=percentile_JJAS_temp(variable,country_key = 'India',region_box = 'SouthIndia');



data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'
sio.savemat(data_path+'TNX_BINS_PDF_1971_2000_JJAS_NCEP.mat',{'SEC_bins_pdfs':SEC_bins_pdfs,'CEC_bins_pdfs':CEC_bins_pdfs,'SID_bins_pdfs':SID_bins_pdfs,'NID_bins_pdfs':NID_bins_pdfs})
