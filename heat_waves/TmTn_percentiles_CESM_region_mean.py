
"""
tHIS IS TO CALCULATE THE DESIRED PERCENTILES OF THE JJAS temperature (MX and MN)
input:
	the folders include enseble members of TX and TN
	The ountries_mask_lon_lat_720_360.mat mask 
output:
	the desired percentiles for the regions defined in the region _dic in .mat
	would be results for all the ensembles
	The PDFs of the TX and TN stored into .mat
"""

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as sio

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


def mask_match(country_key,lon,lat):
	"""
	Readin the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	from scipy.interpolate import interp2d  as interp2d
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/countries_mask_lon_lat_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask= ocean_mask[country_key][:]
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
	mask = f(lon, lat); 
	mask[mask > 0] = 1;
	mask[mask < 1] = np.nan
	return mask

######################################################
# The climatology reference 1961-1990  percentiles
######################################################
def PDF_percentile_JJAS_TmTn(file_name,country_key,region_box):
	"""
	readout the Mask file and interpolate it into the climatic variable resolution 
	read the climate variables from .nc files
	Make the domain average for each time points
	"""
	region = region_dic[region_box][:]
	nc_fid = nc4.Dataset(file_name,mode='r')
	year = nc_fid.variables['year'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	mask = mask_match(country_key,lon,lat)
	# mask off ocean and then box clip
	TX = np.multiply(nc_fid.variables['TX'][50:81,:,:,:],mask) # 1971-2000
	TN = np.multiply(nc_fid.variables['TN'][50:81,:,:,:],mask)
	_,_,TX_range = range_clip(region[0],region[1],region[2],region[3],lon,lat,TX)
	_,_,TN_range = range_clip(region[0],region[1],region[2],region[3],lon,lat,TN)
	# spatial average into 30yrs*122days
	TX_TiSe = stats.nanmean(stats.nanmean(TX_range,axis=3),axis=2)
	TN_TiSe = stats.nanmean(stats.nanmean(TN_range,axis=3),axis=2)
	# Reshape for PDF computing
	length = np.shape(TN_TiSe)[0]*np.shape(TN_TiSe)[1]
	TX_TiSe = np.reshape(TX_TiSe,(length,1));TN_TiSe = np.reshape(TN_TiSe,(length,1))
	# bins PDFs and percentiles
	TX50,TX_bins,TX_pdfs=percentile(TX_TiSe,percent=50);TX98,_,_=percentile(TX_TiSe,percent=98);
	TN50,TN_bins,TN_pdfs=percentile(TN_TiSe,percent=50);TN98,_,_=percentile(TN_TiSe,percent=98)
	return TN50,TN98,TX50,TX98,TN_bins,TN_pdfs,TX_bins,TX_pdfs,

######################
###     Main       ###
######################
region_dic = {'NorthIndia':[70,90,22,33],'SouthIndia':[70,90,8,22],'SouthEastChina':[105,125,20,30],'CentralEastChina':[105,125,30,40]}

input_pat = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/historical/'
os.system('find ' + os.path.abspath(input_pat) + ' -name "*' + '.nc' + '" -print | sort > ' + input_pat + '/file_list.txt')
text_file = open(input_pat + '/file_list.txt', "r")
text_content = text_file.readlines()	
LEN = len(text_content) # NUMBER OF ENSEMBLE WHIHC CAN BE OBTAINED



# the first dimensiion for ensemble no, 2nd for [TN_bins,TN_pdfs,TX_bins,TX_pdfs]
SEC_bins_pdfs = np.zeros((30,4,100));CEC_bins_pdfs = np.zeros((30,4,100));
SID_bins_pdfs = np.zeros((30,4,100));NID_bins_pdfs = np.zeros((30,4,100));
# print np.shape(SEC_bins_pdfs[:])
for en_no in range(LEN): 
	file_name = text_content[en_no][:-1]
	_,_,_,_,TN_bins,TN_pdfs,TX_bins,TX_pdfs= PDF_percentile_JJAS_TmTn(file_name,country_key = 'China',region_box = 'SouthEastChina')
	print np.shape(TN_bins)
	SEC_bins_pdfs[en_no,0,:] = TN_bins;SEC_bins_pdfs[en_no,1,:] = TN_pdfs;
	SEC_bins_pdfs[en_no,2,:] = TX_bins;SEC_bins_pdfs[en_no,3,:] = TX_pdfs;
	_,_,_,_,TN_bins,TN_pdfs,TX_bins,TX_pdfs= PDF_percentile_JJAS_TmTn(file_name,country_key = 'China',region_box = 'CentralEastChina')
	CEC_bins_pdfs[en_no,0,:] = TN_bins;CEC_bins_pdfs[en_no,1,:] = TN_pdfs;
	CEC_bins_pdfs[en_no,2,:] = TX_bins;CEC_bins_pdfs[en_no,3,:] = TX_pdfs;
	_,_,_,_,TN_bins,TN_pdfs,TX_bins,TX_pdfs= PDF_percentile_JJAS_TmTn(file_name,country_key = 'India',region_box = 'SouthIndia')
	SID_bins_pdfs[en_no,0,:] = TN_bins;SID_bins_pdfs[en_no,1,:] = TN_pdfs;
	SID_bins_pdfs[en_no,2,:] = TX_bins;SID_bins_pdfs[en_no,3,:] = TX_pdfs;	
	_,_,_,_,TN_bins,TN_pdfs,TX_bins,TX_pdfs= PDF_percentile_JJAS_TmTn(file_name,country_key = 'India',region_box = 'NorthIndia')
	NID_bins_pdfs[en_no,0,:] = TN_bins;NID_bins_pdfs[en_no,1,:] = TN_pdfs;
	NID_bins_pdfs[en_no,2,:] = TX_bins;NID_bins_pdfs[en_no,3,:] = TX_pdfs;

data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'
sio.savemat(data_path+'TNX_BINS_PDF_1971_2000_JJAS_CESM30.mat',{'SEC_bins_pdfs':SEC_bins_pdfs,'CEC_bins_pdfs':CEC_bins_pdfs,'SID_bins_pdfs':SID_bins_pdfs,'NID_bins_pdfs':NID_bins_pdfs})


	



