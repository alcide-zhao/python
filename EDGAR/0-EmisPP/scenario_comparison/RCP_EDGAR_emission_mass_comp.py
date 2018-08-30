# -*- coding: utf-8 -*-


"""
This script is created to verify aerosol molecule and particle numbers produced using my methods are correct.
The emissions in mass are compared to see the difference between EDGAR and RCP8.5 for 2010
EDGAR emissions are 1800*3600 while RCP85 is in 360*720
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
import scipy.io as sio
from mpl_toolkits.basemap import interp

# from scipy.interpolate import interp2d

#######################################
# 0.0 data input
#######################################

#######################################################
# 0.2 functions and subroutines                       #
#######################################################
def movingaverage (values, window=3):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result
def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

	
def range_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,data):
	"""
	clip the data based on given range
	"""
	lon = np.array(lon)
	lat = np.array(lat)
	colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
	# print colum_s
	colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
	row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
	row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
	if (colum_s> colum_e):
		cache = colum_e; colum_e = colum_s; colum_s = cache;
	if (row_s> row_e):
		cache = row_e; row_e = row_s; row_s = cache;
	lon_clipped = lon[colum_s:colum_e+1]
	lat_clipped = lat[row_s:row_e+1]
	if (np.rank(data) == 2):
		data_clipped = data[row_s:row_e+1,colum_s:colum_e+1]
	elif (np.rank(data) == 3):
		data_clipped = data[:,row_s:row_e+1,colum_s:colum_e+1]
	elif (np.rank(data) == 4):
		data_clipped = data[:,:,row_s:row_e+1,colum_s:colum_e+1]
	return lon_clipped, lat_clipped, data_clipped
	
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
region_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,60],'EA':[100,145,20,50],'SA':[65,100,5,30],'EUR':[0,40,48,75],'NAM':[210,300,0,60]}

## read the CESM grids to interpolate the RCP emissionsn firstly before caculations
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
lons,lats = np.meshgrid (lon_CESM,lat_CESM)
lat_res = lat_CESM[1] - lat_CESM[0];lon_res = lon_CESM[1] - lon_CESM[0];
CESM_gc_area = getArea(lons,lons+lon_res,lats,lats+lat_res)

# CESM_fid.close()



"""
###############################################################
# RCP8.5 SO2 emissions for the year 2005, 2010 ,2050 , 2100
# annual total emissions with gridcell area considered
###############################################################
print '####################################################'
print '                       RCP85                       '
print '####################################################'
# emissions
file_path = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
rcp85 = file_path+'accmip_interpolated_emissions_RCP85_SO2_2005_2100_0.5x0.5.nc'
nc_fid = nc4.Dataset(rcp85,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]

rcp_anth=nc_fid.variables['anthropogenic'][:]+nc_fid.variables['ships'][:]
years=[2005,2010,2030,2050,2080,2100]
for iyear in years:
	print iyear
	layer_s = (iyear-2005)*12;layer_e = (iyear-2005+1)*12;
	# annual total emissions, from kg/m2/s to kg/m2/year
	rcp_annual_sum = 10**(-9)*365*24*60*60*stats.nanmean(rcp_anth[layer_s:layer_e,:,:],axis=0)#/12 #
	# from kg/m2/year to kg/(gridcell)/year
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid (lon,lat)  #meshgrid
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	rcp_gc_sum = np.multiply(rcp_annual_sum,area)   # area of the grid cell in m2 * emission/m2/year
	for key in region_dic:
		region = region_dic[key]
		_,_,rcp_region_gc = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp_gc_sum)
		region_sum = np.nansum(np.nansum(rcp_region_gc,axis=1),axis=0)
		print key, round(region_sum,1)

###############################################################
# EDGAR SO2 emissions for the year 2010
# annual total (all anthropogenic sectors) emissions with gridcell area considered 
###############################################################
"""
print '####################################################'
print '                       EDGAR                        '
print '####################################################'
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
region_dic={'GLOBE':[0,360,-90,90]}#,'ASIA':[60,150,-15,60],'EA':[100,145,20,50],'SA':[65,100,5,30],'EUR':[0,40,48,75],'NAM':[210,300,0,60]}
calendar_day = {'1':31,'2':28,'3':31,'4':30,'5':31,'6':30,'7':31,'8':31,'9':30,'10':31,'11':30,'12':31}
species=['SO2'] 
emission='2010' #2010 or stat_tech
file_path ="/exports/csce/datastore/geos/users/s1667168/EDGAR/"
# Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH','CDS','CRS','LTO','SPS']
Sectors = ['AGR','SHIP','SWD','OTH','RCO','TRO','TNG','CDS','CRS','LTO','SPS']
# Sectors = ['ENE','IND','TRF','REF','PRO','PPA']
# Sectors = ['ENE','IND','TRO']

def regrid_edgar2cesm(lon,lat,emission,CESM_gc_area):
	x=lon;y=np.linspace(lat[0],lat[-1],1824);
	# print y
	xout,yout=np.meshgrid(x,y)
	cesm_emission = np.zeros((96,144));
	emission_1824 = interp(emission, lon, lat, xout, yout, checkbounds=False, masked=False, order=1)
	lons,lats = np.meshgrid (x,y)
	lat_res = y[1] - y[0];lon_res = x[1] - x[0];
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	emission_gc_mass = np.multiply(emission_1824,area)
	# plt.imshow(emission_gc_mass,origin='lower');plt.show()
	x_cesm=range(0,144);y_cesm=range(0,96)
	for row in range(0,96):
		for column in range(0,144):	
			cesm_emission[row,column] =np.nansum(np.nansum(emission_gc_mass[row*19:(row+1)*19,column*25:(column+1)*25],axis=1),axis=0)
	cesm_emission = np.divide(cesm_emission,CESM_gc_area)
	# plt.imshow(cesm_emission,origin='lower');plt.show()
	return cesm_emission
	
months=range(0,12) 

for specie in species:
	SO2_mass=np.zeros((12,96,144));
	for month in months:
		day=calendar_day[str(month+1)]
		for sector in Sectors:
			if emission == 'stag_tech':
				if sector in ['ENE','IND','TRO']:
					file = file_path+"stag_tech/"+specie+'/'+'JRC_PEGASOS_V2_STAG_TECH_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
			else:
				# print sector
				file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
			# print file
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			SO2=nc_fid.variables['emi_'+specie.lower()][:]
			SO2_mass[month,:,:] =SO2_mass[month,:,:]+10**(-9)*day*24*60*60*regrid_edgar2cesm(lon,lat,SO2,CESM_gc_area) #interp(SO2, lon, lat, lons, lats, checkbounds=False, masked=False, order=1)
	# edgar_annual_sum = 10**(-9)*365*24*60*60*stats.nanmean(SO2_mass,axis=0) #
	# f = interp2d(lon,lat,edgar_annual_sum,kind='linear')
	# edgar_annual_sum = f(lon_CESM, lat_CESM);
	# lon_res = lon_CESM[1] - lon_CESM[0];lat_res = lat_CESM[1] - lat_CESM[0];
	# lons,lats = np.meshgrid (lon_CESM,lat_CESM)  #meshgrid
	# for key in region_dic:
		# region = region_dic[key]
		# _,_,edgar_region_gc = range_clip(region[0],region[1],region[2],region[3],lon,lat,edgar_gc_sum)
	 	# region_sum = np.nansum(np.nansum(edgar_region_gc,axis=1),axis=0)
		# print key, round(region_sum,3)
	edgar_gc_sum = np.multiply(SO2_mass,CESM_gc_area)   # area of the grid cell in m2 * emission/m2/year

	#time series of 12 months
	edgar_ts=  np.nansum(np.nansum(edgar_gc_sum,axis=2),axis=1)
	print np.sum(edgar_ts)
	

###############################################################
# RCP8.5 by sectors
###############################################################
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]

print '####################################################'
print '                       RCP85                       '
print '####################################################'
file_path ="/exports/csce/datastore/geos/users/s1667168/RCP/"
rcp_file = file_path+'accmip_interpolated_emissions_RCP85_SO2_ships_2010_0.5x0.5.nc'
nc_fid = nc4.Dataset(rcp_file,mode='r')
SHP =  nc_fid.variables['emiss_shp'][:]; # 0.1 from kg/m2/s to g/cm2/s
nc_fid.close()
rcp_file = file_path+'accmip_interpolated_emissions_RCP85_SO2_anthropogenic_2010_0.5x0.5.nc'
nc_fid = nc4.Dataset(rcp_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
IND =  nc_fid.variables['emiss_ind'][:];
ENE =  nc_fid.variables['emiss_ene'][:];
WST =  nc_fid.variables['emiss_wst'][:];
AWB =  nc_fid.variables['emiss_awb'][:];
DOM =  nc_fid.variables['emiss_dom'][:];
TRA =  nc_fid.variables['emiss_tra'][:];
nc_fid.close()


rcp_anth=TRA+AWB+DOM+WST+SHP#+ENE+IND#+
months=range(0,12)
rcp_annual_sum=np.zeros((12,96,144))
for imonth in months:
	day=calendar_day[str(month+1)]
	SO2= rcp_anth[imonth,:,:]
	rcp_annual_sum [imonth,:,:]= 10**(-9)*day*24*60*60*interp(SO2, lon, lat, lons, lats, checkbounds=False, masked=False, order=1)
rcp_gc_sum = np.multiply(rcp_annual_sum,CESM_gc_area)  # area of the grid cell in m2 * emission/m2/year
rcp_ts=np.nansum(np.nansum(rcp_gc_sum,axis=2),axis=1)
	# for key in region_dic:
		# region = region_dic[key]
		# _,_,rcp_region_gc = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp_gc_sum)
		# region_sum = np.nansum(np.nansum(rcp_region_gc,axis=1),axis=0)
		# print key, round(region_sum,1)
	
print np.sum(rcp_ts)
x=range(0,12)
plt.plot(x,edgar_ts,'r')
plt.plot(x,rcp_ts)
plt.show()

# fig = plt.figure(facecolor='White',figsize=(8.5, 3.5));	
# plt.subplot(1, 2, 1)
# rcp_region_gc[rcp_region_gc==0]=np.nan
# plt.imshow(rcp_region_gc,cmap='rainbow',origin='lower');
# # plt.clim(0,0.00001)
# plt.subplot(1, 2, 2)
# edgar_region_gc[edgar_region_gc==0]=np.nan
# plt.imshow(edgar_region_gc,cmap='rainbow',origin='lower');
# # plt.clim(0,0.00001)
# # plt.colorbar()

