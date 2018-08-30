# -*- coding: utf-8 -*-
"""
Produce the annual cycle of BC, OC and SO4 aerosols emissions
plot out and compare between scenarios
"""

import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import math
import site
import scipy.io as sio

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

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
	calculate the earth area in m2
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
	
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.region_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,60],'EA':[100,145,20,50],'SA':[65,100,5,30],'EUR':[0,40,48,75],'NAM':[210,300,0,60]}

def MonthWeight(annual_data,area):
	"""
	weight the days to each month and sum the mass up throughout the year
	"""
	calendar_day = {'1':31,'2':28,'3':31,'4':30,'5':31,'6':30,'7':31,'8':31,'9':30,'10':31,'11':30,'12':31}
	month_sum= np.empty((np.shape(annual_data)))
	for month in range(0,12):
		day_no = calendar_day[str(month+1)]
		month_sum[month,:,:] =annual_data[month,:,:]*day_no*24*60*60/(10**(9))
	month_total = np.nansum(np.nansum(np.multiply(month_sum,area),axis=2),axis=1)
	return month_total

def get_RCP85_emissions(species):
	"""
	get the sum of emissions over 2010 from all the three sectors
	"""
	rcp85_path = '/exports/csce/datastore/geos/users/s1667168/RCP85/'
	rcp85 = rcp85_path+'accmip_interpolated_emissions_RCP85_'+species.upper()+'_2005_2100_0.5x0.5.nc'
	nc_fid = nc4.Dataset(rcp85,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	emission = (nc_fid.variables['anthropogenic'][:]+nc_fid.variables['ships'][:]+nc_fid.variables['biomass_burning'][:])[60:72,:,:]   	
	lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
	lons,lats = np.meshgrid (lon,lat)  #meshgrid
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
		# for key in region_dic:
		# region = region_dic[key]
		# _,_,rcp_region_gc = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp_gc_sum)
		# region_sum = np.nansum(np.nansum(rcp_region_gc,axis=1),axis=0)
		# print key, round(region_sum,1)
	AreaWeighMass = MonthWeight(emission,area)
	# AreaWeighMass = 10**(-9)*24*60*60*np.nansum(np.nansum(np.multiply(np.multiply(emission,area),axis=2),axis=1),calendar_day)
	return AreaWeighMass

	
def get_EDGAR_emissions(scenario, specie):
	file_path ="/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/"
	if scenario == '1970':
		if specie=='OC':
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH']
		elif specie=='NH3':
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PPA']
		else:
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH']
	else:
		if specie=='OC':
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH','CDS','CRS','LTO','SPS'] #
		elif specie=='NH3':
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PPA']
		elif species == "NOx":
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH']
		else:
			Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH','CDS','CRS','LTO','SPS'] #
	emiss=np.zeros((12,1800,3600));
	print scenario, specie
	for month in range(0,12):
		for sector in Sectors:
			if  scenario == 'STAG_TECH':
				if sector in ['ENE','IND','TRO']:
					file = file_path+scenario+"/"+specie+'/'+'JRC_PEGASOS_V2_STAG_TECH_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
			elif scenario == 'STAG_ENE':
				if sector in ['ENE','IND','TRO']:
					file = file_path+scenario+"/"+specie+'/'+'JRC_PEGASOS_V2_STAG_ENERGY_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
			elif scenario == '2010':
				file = file_path+scenario+"/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_'+scenario+'_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
			else:
				file = file_path+scenario+"/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_'+scenario+'_'+str(month+1)+'_'+sector+'.nc'
			# print file
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
			lons,lats = np.meshgrid (lon,lat)  #meshgrid
			area = getArea(lons,lons+lon_res,lats,lats+lat_res)		
			emis_month=nc_fid.variables['emi_'+specie.lower()][:]
			emiss[month,:,:] = emiss[month,:,:]+emis_month
	# edgar_ts = MonthWeight(emiss,area)
	# edgar_ts=  np.multiply(np.nansum(np.nansum(np.multiply(emiss,area),axis=2),axis=1),calendar_day)*10**(-9)*24*60*60
	# return edgar_ts
	return lon,lat,emiss,area
	
def get_EDGAR_NMVOC(scenario,specie):
	print scenario,specie
	file_path ="/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/"
	bio=np.zeros((12,1800,3600));fossil=np.zeros((12,1800,3600));
	#1970
	if scenario == '1970':
		Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH']
		for month in range(0,12):
			for sector in Sectors:
				file = file_path+scenario+"/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_'+scenario+'_'+str(month+1)+'_'+sector+'.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')
				lat = nc_fid.variables['lat'][:]
				lon = nc_fid.variables['lon'][:]
				lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
				lons,lats = np.meshgrid (lon,lat)  #meshgrid
				area = getArea(lons,lons+lon_res,lats,lats+lat_res)		
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				bio[month,:,:] = bio[month,:,:]+emis_month 	
		# edgar_ts = MonthWeight(bio,area)
		emis = bio
		# edgar_ts=  np.multiply(np.nansum(np.nansum(np.multiply(bio,area),axis=2),axis=1),calendar_day)*10**(-9)*24*60*60
	#2010
	elif scenario == '2010':
		for month in range(0,12):
			Sectors = ['ENE','IND','AGR','RCO','REF','SWD','TRO','TRF']
			for sector in Sectors:
				file = file_path+scenario+"/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_bio_'+scenario+'_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')	
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				bio[month,:,:] = bio[month,:,:]+emis_month 	
			Sectors = ['ENE','IND','SHIP','SWD','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH','CDS','CRS','LTO','SPS']
			for sector in Sectors:
				file = file_path+scenario+"/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_fossil_'+scenario+'_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')
				lat = nc_fid.variables['lat'][:]
				lon = nc_fid.variables['lon'][:]
				lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
				lons,lats = np.meshgrid (lon,lat)  #meshgrid
				area = getArea(lons,lons+lon_res,lats,lats+lat_res)		
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				fossil[month,:,:] = fossil[month,:,:]+emis_month 			
		# edgar_ts = MonthWeight(bio+fossil,area)
		emis = bio+fossil
		# edgar_ts=  np.multiply(np.nansum(np.nansum(np.multiply(bio+fossil,area),axis=2),axis=1),calendar_day)*10**(-9)*24*60*60	
		
	#STAG_TECH
	elif scenario == 'STAG_TECH':
		for month in range(0,12):
			Sectors = ['ENE','IND','AGR','RCO','REF','SWD','TRO','TRF']
			for sector in Sectors:
				if sector in ['ENE','IND','TRO']:
					file = file_path+scenario+"/"+specie+'/'+'JRC_PEGASOS_V2_STAG_TECH_'+specie+'_bio_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_bio_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')	
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				bio[month,:,:] = bio[month,:,:]+emis_month 	
			Sectors = ['ENE','IND','SHIP','SWD','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH','CDS','CRS','LTO','SPS']
			for sector in Sectors:
				if sector in ['ENE','IND','TRO']:
					file = file_path+scenario+"/"+specie+'/'+'JRC_PEGASOS_V2_STAG_TECH_'+specie+'_fossil_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_fossil_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')
				lat = nc_fid.variables['lat'][:]
				lon = nc_fid.variables['lon'][:]
				lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
				lons,lats = np.meshgrid (lon,lat)  #meshgrid
				area = getArea(lons,lons+lon_res,lats,lats+lat_res)		
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				fossil[month,:,:] = fossil[month,:,:]+emis_month 			
		# edgar_ts = MonthWeight(bio+fossil,area)
		emis = bio+fossil
		# edgar_ts=  np.multiply(np.nansum(np.nansum(np.multiply(bio+fossil,area),axis=2),axis=1),calendar_day)*10**(-9)*24*60*60	
	#stag_ene
	elif scenario == 'STAG_ENE':
		for month in range(0,12):
			Sectors = ['ENE','IND','AGR','RCO','REF','SWD','TRO','TRF']
			for sector in Sectors:
				if sector in ['ENE','IND','TRO']:
					file = file_path+scenario+"/"+specie+'/'+'JRC_PEGASOS_V2_STAG_ENERGY_'+specie+'_bio_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_bio_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')	
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				bio[month,:,:] = bio[month,:,:]+emis_month 	
			Sectors = ['ENE','IND','SHIP','SWD','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH','CDS','CRS','LTO','SPS']
			for sector in Sectors:
				if sector in ['ENE','IND','TRO']:
					file = file_path+scenario+"/"+specie+'/'+'JRC_PEGASOS_V2_STAG_ENERGY_'+specie+'_fossil_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_fossil_2010_'+str(month+1)+'_'+sector+'.0.1x0.1.nc'
				# print file
				nc_fid = nc4.Dataset(file,mode='r')
				lat = nc_fid.variables['lat'][:]
				lon = nc_fid.variables['lon'][:]
				lat_res=lat[1]-lat[0];lon_res=lon[1]-lon[0];
				lons,lats = np.meshgrid (lon,lat)  #meshgrid
				area = getArea(lons,lons+lon_res,lats,lats+lat_res)		
				emis_month=nc_fid.variables['emi_'+specie.lower()][:]
				fossil[month,:,:] = fossil[month,:,:]+emis_month 			
		# edgar_ts = MonthWeight(bio+fossil,area)
		emis = bio+fossil
		# edgar_ts=  np.multiply(np.nansum(np.nansum(np.multiply(bio+fossil,area),axis=2),axis=1),calendar_day)*10**(-9)*24*60*60			
	return emis
	
output_path ="/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/"


scenario = "1970"  # 2010  STAG_TECH  STAG_ENE 1970
species = "BC";lon,lat,BC,area = get_EDGAR_emissions(scenario,species);
species = "OC";lon,lat,OC,area = get_EDGAR_emissions(scenario,species);
species = "SO2";lon,lat,SO2,area = get_EDGAR_emissions(scenario,species);
species = "CO";lon,lat,CO,area = get_EDGAR_emissions(scenario,species);
species = "NH3";lon,lat,NH3,area = get_EDGAR_emissions(scenario,species);
species = "NOx";lon,lat,NOx,area = get_EDGAR_emissions(scenario,species);
species = "NMVOC";NMVOC=get_EDGAR_NMVOC(scenario,species)

file_name = output_path+'EDGAR4.3.1_Sector_merged.'+scenario+'.nc'
f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('month', 12)

months = f.createVariable('month',np.float64, ('month'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))

area = f.createVariable('gridarea',np.float32,('lat','lon'))
BCs = f.createVariable('BC',np.float32,('month','lat','lon'))
OCs = f.createVariable('OC',np.float32,('month','lat','lon'))
SO2s = f.createVariable('SO2',np.float32,('month','lat','lon'))
NMVOCs = f.createVariable('NMVOC',np.float32,('month','lat','lon'))
COs = f.createVariable('CO',np.float32,('month','lat','lon'))
NH3s = f.createVariable('NH3',np.float32,('month','lat','lon'))
NOxs = f.createVariable('NOx',np.float32,('month','lat','lon'))	

months[:] = range(1,13)
latitudes[:] = lat
longitudes[:] = lon
BCs[:]= BC;OCs[:]= OC;SO2s[:]= SO2;NMVOCs[:]= NMVOC;
COs[:]= CO;NOxs[:]= NOx;NH3s[:]= NH3;
area.units="m2"
f.description = 'EDHAR aerosol emissions: all sector merged, units = kg/m2/s'

f.close()
"""
BC_1970=get_EDGAR_emissions('1970','BC'); OC_1970=get_EDGAR_emissions('1970','OC'); SO2_1970=get_EDGAR_emissions('1970','SO2'); 
# BC_2010=get_EDGAR_emissions('2010','BC'); OC_2010=get_EDGAR_emissions('2010','OC'); SO2_2010=get_EDGAR_emissions('2010','SO2'); 
# BC_TECH=get_EDGAR_emissions('STAG_TECH','BC'); OC_TECH=get_EDGAR_emissions('STAG_TECH','OC'); SO2_TECH=get_EDGAR_emissions('STAG_TECH','SO2'); 
# BC_ENE=get_EDGAR_emissions('STAG_ENE','BC'); OC_ENE=get_EDGAR_emissions('STAG_ENE','OC'); SO2_ENE=get_EDGAR_emissions('STAG_ENE','SO2'); 
# NMVOC_1970=get_EDGAR_NMVOC('1970','NMVOC'); NMVOC_2010=get_EDGAR_NMVOC('2010','NMVOC'); NMVOC_ENE=get_EDGAR_NMVOC('STAG_ENE','NMVOC'); NMVOC_TECH=get_EDGAR_NMVOC('STAG_TECH','NMVOC'); 

# sio.savemat(file_path+'BC_OC_SO2_NMVOC_EDGAR_RCP8.5.mat',{'BC_RCP85':BC_RCP85,'BC_1970':BC_1970,'BC_2010':BC_2010,'BC_TECH':BC_TECH,'BC_ENE':BC_ENE,\
# 'OC_RCP85':OC_RCP85,'OC_1970':OC_1970,'OC_2010':OC_2010,'OC_TECH':OC_TECH,'OC_ENE':OC_ENE,\
# 'SO2_RCP85':SO2_RCP85,'SO2_1970':SO2_1970,'SO2_2010':SO2_2010,'SO2_TECH':SO2_TECH,'SO2_ENE':SO2_ENE,\
# 'NMVOC_1970':NMVOC_1970,'NMVOC_2010':NMVOC_2010,'NMVOC_TECH':NMVOC_TECH,'NMVOC_ENE':NMVOC_ENE})

data = sio.loadmat(file_path+'BC_OC_SO2_NMVOC_EDGAR_RCP8.5.mat')
BC_RCP85 = data['BC_RCP85'][0,:];BC_1970 = data['BC_1970'][0,:];BC_2010 = data['BC_2010'][0,:];BC_TECH = data['BC_TECH'][0,:];BC_ENE = data['BC_ENE'][0,:];
OC_RCP85 = data['OC_RCP85'][0,:];OC_1970 = data['OC_1970'][0,:];OC_2010 = data['OC_2010'][0,:];OC_TECH = data['OC_TECH'][0,:];OC_ENE = data['OC_ENE'][0,:];
SO2_RCP85 = data['SO2_RCP85'][0,:];SO2_1970 = data['SO2_1970'][0,:];SO2_2010 = data['SO2_2010'][0,:];SO2_TECH = data['SO2_TECH'][0,:];SO2_ENE = data['SO2_ENE'][0,:];
NMVOC_1970 = data['NMVOC_1970'][0,:];NMVOC_2010 = data['NMVOC_2010'][0,:];NMVOC_TECH = data['NMVOC_TECH'][0,:];NMVOC_ENE = data['NMVOC_ENE'][0,:];


print np.sum(SO2_2010),np.sum(SO2_TECH),np.sum(SO2_RCP85)
time_series=range(1,13)
fig1 = plt.figure(facecolor='White',figsize=[8,6]);plot_setup();
ax1 = plt.subplot(2,2,1);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.2,0.6]);ax1.set_yticks(np.arange(0.20,0.61,0.1))
ax1.set_xticks(np.array([1,4,7,10,12]));ax1.set_xticklabels(['Jan','Apr','July','Oct','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(a) BC',xy=(0.05, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
lins1 = ax1.plot(time_series,np.divide(BC_1970[:],1),'-',color="k",label='REF_1970',linewidth=1.5)
lins2 = ax1.plot(time_series,np.divide(BC_2010[:],1),'-',color="b",label='REF_2010',linewidth=1.5)
# lins3 = ax1.plot(time_series,BC_RCP85[:],'-',color="g",label='RCP85_2010',linewidth=1.5)
lins4 = ax1.plot(time_series,np.divide(BC_TECH[:],1),'-',color="r",label='TECH_2010',linewidth=1.5)
lins5 = ax1.plot(time_series,np.divide(BC_ENE[:],1),'-',color="g",label='ENE_2010',linewidth=1.5)

ax1 = plt.subplot(2,2,2);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.5,1.5]);ax1.set_yticks(np.arange(0.5,1.51,0.25))
ax1.set_xticks(np.array([1,4,7,10,12]));ax1.set_xticklabels(['Jan','Apr','July','Oct','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(b) OC',xy=(0.05, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
lins1 = ax1.plot(time_series,np.divide(OC_1970[:],1),'-',color="k",label='REF_1970',linewidth=1.5)
lins2 = ax1.plot(time_series,np.divide(OC_2010[:],1),'-',color="b",label='REF_2010',linewidth=1.5)
# lins3 = ax1.plot(time_series,OC_RCP85[:],'-',color="g",label='RCP85_2010',linewidth=1.5)
lins4 = ax1.plot(time_series,np.divide(OC_TECH[:],1),'-',color="r",label='TECH_2010',linewidth=1.5)
lins5 = ax1.plot(time_series,np.divide(OC_ENE[:],1),'-',color="g",label='ENE_2010',linewidth=1.5)


ax1 = plt.subplot(2,2,3);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([3,13]);ax1.set_yticks(np.arange(3,13.1,2.5))
ax1.set_xticks(np.array([1,4,7,10,12]));ax1.set_xticklabels(['Jan','Apr','Jul','Oct','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(c) SO2',xy=(0.05, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
lins1 = ax1.plot(time_series,np.divide(SO2_1970[:],1),'-',color="k",label='REF_1970',linewidth=1.5)
lins2 = ax1.plot(time_series,np.divide(SO2_2010[:],1),'-',color="b",label='REF_2010',linewidth=1.5)
# lins3 = ax1.plot(time_series,SO2_RCP85[:],'-',color="g",label='RCP85_2010',linewidth=1.5)
lins4 = ax1.plot(time_series,np.divide(SO2_TECH[:],1),'-',color="r",label='TECH_2010',linewidth=1.5)
lins5 = ax1.plot(time_series,np.divide(SO2_ENE[:],1),'-',color="g",label='ENE_2010',linewidth=1.5)

legend = ax1.legend(shadow=False,ncol=5,bbox_to_anchor=(1.8, -0.1))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.1)


ax1 = plt.subplot(2,2,4);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([6,14]);ax1.set_yticks(np.arange(6,14.1,2))
ax1.set_xticks(np.array([1,4,7,10,12]));ax1.set_xticklabels(['Jan','Apr','July','Oct','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(d) NMVOC',xy=(0.05, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
lins1 = ax1.plot(time_series,np.divide(NMVOC_1970[:],1),'-',color="k",label='REF_1970',linewidth=1.5)
lins2 = ax1.plot(time_series,np.divide(NMVOC_2010[:],1),'-',color="b",label='REF_2010',linewidth=1.5)
# lins3 = ax1.plot(time_series,SO2_RCP85[:],'-',color="g",label='RCP85_2010',linewidth=1.5)
lins4 = ax1.plot(time_series,np.divide(NMVOC_TECH[:],1),'-',color="r",label='TECH_2010',linewidth=1.5)
lins5 = ax1.plot(time_series,np.divide(NMVOC_ENE[:],1),'-',color="g",label='ENE_2010',linewidth=1.5)


plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3);
plt.savefig('f_BC_OC_SO2_anuual_cycle_EDGAR_RCP8.5.png', format='png', dpi=1000)
# plt.show()
"""