import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats


import os; import site
lib_path = os.path.join(
	os.path.realpath(
       os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio

from scipy.stats.kde import gaussian_kde
from numpy import linspace

ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask_CESM.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;ocean_mask_CESM[0:27,:]=np.nan


# linuc path
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/his/'
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()


TX90 = np.empty((1));TX90[:] = np.nan   
TN90 = np.empty((1));TN90[:] = np.nan
TX5 = np.empty((1));TX5[:] = np.nan   
TN5 = np.empty((1));TN5[:] = np.nan  
TX_cache = np.empty((1,30*365));TX_cache[:] = np.nan 
TN_cache = np.empty((1,30*365));TN_cache[:] = np.nan 

fig = plt.figure(facecolor='White',figsize=(20, 5));plot_setup();pad= 5 
for ensumble_member in range(0,1): 
	print ensumble_member
	nc_f = text_content[ensumble_member][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	#  RUSSIA 168 80  South America 92 236   South China 121,88
	TX = np.multiply(nc_fid.variables['TX'],ocean_mask_CESM)[41:71,:,121,88]
	TN = np.multiply(nc_fid.variables['TN'],ocean_mask_CESM)[41:71,:,121,88]
	#  varibale initialization 
	TX_cache[ensumble_member,:] = np.reshape(TX,(30*365))
	ax1 = plt.subplot(2,1,1);
	# time_series = range(0,365*30)
	# ax1.plot(time_series,TX_cache[ensumble_member,:],'k')
	dist_space = linspace( -50, 50, 100)
	kde = gaussian_kde( TX_cache[ensumble_member,:] )
	ax1.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
	
	TN_cache[ensumble_member,:]  = np.reshape(TN,(30*365))
	ax2 = plt.subplot(2,1,2);
	# time_series = range(0,365*30)
	# ax2.plot(time_series,TN_cache[ensumble_member,:],'k')
	kde = gaussian_kde( TN_cache[ensumble_member,:] )
	ax2.plot(dist_space, 100*kde(dist_space),'k-',alpha=1)
	
threshold_data ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/TmTn_percentile_365.nc'
nc_fid = nc4.Dataset(threshold_data,mode='r')
lat = nc_fid.variables['lat']
lon = nc_fid.variables['lon']
TX90P = np.multiply(nc_fid.variables['TX90P'],ocean_mask_CESM)[121,88]
TN90P = np.multiply(nc_fid.variables['TN90P'],ocean_mask_CESM)[121,88]
TX10P = np.multiply(nc_fid.variables['TX10P'],ocean_mask_CESM)[121,88]
TN10P = np.multiply(nc_fid.variables['TN10P'],ocean_mask_CESM)[121,88]

ax1 = plt.subplot(2,1,1);
# time_series = range(0,365*30)
# ax1.plot(time_series,TX90P*np.ones([365*30]),'r')
# ax1.plot(time_series,TX5P*np.ones([365*30]),'b')
ax1.axvline(TX90P, color='r',alpha=1);
ax1.axvline(TX10P, color='b',alpha=1);

ax2 = plt.subplot(2,1,2);
# ax2.plot(time_series,TN90P*np.ones([365*30]),'r')
# ax2.plot(time_series,TN5P*np.ones([365*30]),'b')
ax2.axvline(TN90P, color='r',alpha=1);
ax2.axvline(TN10P, color='b',alpha=1);


plt.subplots_adjust(left=0.1, bottom=0.03, right=0.95, top=0.95, wspace=0.25, hspace=0.2);
plt.savefig('China_PDF.png', format='png', dpi=200)
plt.show()





