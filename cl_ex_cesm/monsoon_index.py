# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 2017
This is to produce the monsoon index for the Asia major monsoon regions

Webster-yang monsoon index : (U850-U200 averaged over 0-20oN, 40oE-110oE)
Australian monsoon index (U850 averaged over 2.5oS-15oS, 110oE-150oE)
South Asian monsoon index (V850-V200 averaged over 10oN-30oN, 70oE-110oE)
Dynamic Indian monsoon index (U850 (5oN-15oN, 40oE-80oE) - (U850 20oN-30oN, 70oE-90oE))
East Asian - Western North Pacific monsoon index (U850 (5oN-15oN, 90oE-130oE) - U850 (20oN-30oN, 110oE-140oE))

The result would be written to a txt file
"""

import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import numpy as np
import os
import site
lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)

site.addsitedir(lib_path)
from lib import range_clip

#######################################
#  Monsoon Index caculation function
#######################################

# region definitions and data input
region_dic = {'WYI':[40,110,0,20],'AUSI':[110,150,-15,-2.5],'SAMI':[35,97.5,5,22.5],'DIMI1':[40,80,5,15],\
'DIMI2':[70,90,25,35],'EAWNPMI1':[100,130,5,15],'EAWNPMI2':[110,140,20,35]}


#Webster_Yang_monsoon_index		
def Webster_Yang_monsoon_index(lon,lat,u850,u200):
	region = region_dic['WYI']
	u_shear = u850-u200
	_,_,u_shear_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,u_shear)
	WYI = stats.nanmean(stats.nanmean(u_shear_clipped,axis = 2),axis =1)
	return WYI

# Australian monsoon index
def Australian_monsoon_index(lon,lat,u850):
	region = region_dic['AUSI']
	_,_,u850_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	AMI =stats.nanmean(stats.nanmean(u850_clipped,axis = 2),axis =1)
	return AMI

#South Asian monsoon index	
def South_Asian_monsoon_index(lon,lat,u850):
	region = region_dic['SAMI']
	_,_,u850_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	SAMI = stats.nanmean(stats.nanmean(u850_clipped,axis = 2),axis =1)
	return SAMI
	
#Dynamic Indian monsoon index
def Dynamic_Indian_monsoon_index(lon,lat,u850):
	region = region_dic['DIMI1']
	_,_,u850_DIMI1 = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	region = region_dic['DIMI2']
	_,_,u850_DIMI2 = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	DIMI = stats.nanmean(stats.nanmean(u850_DIMI1,axis = 2),axis =1)-stats.nanmean(stats.nanmean(u850_DIMI2,axis = 2),axis =1)
	return DIMI

#DEast Asian - Western North Pacific monsoon index 
def WNP_monsoon_index(lon,lat,u850):
	region = region_dic['EAWNPMI1']
	_,_,u850_EAWNPMI1 = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	region = region_dic['EAWNPMI2']
	_,_,u850_EAWNPMI2 = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	WNPMI = stats.nanmean(stats.nanmean(u850_EAWNPMI1,axis = 2),axis =1)-stats.nanmean(stats.nanmean(u850_EAWNPMI2,axis = 2),axis =1)
	return WNPMI
	
def East_Asia_Monsoon_index(lon,lat,u850):
	region = [110,140,10,40]
	_,_,u850_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	EASI = stats.nanmean(stats.nanmean(u850_clipped,axis = 2),axis =1)
	return EASI
# East Asian Summer Jet	
def East_Asian_Summer_Jet(lon,lat,u200):
	region = [110,140,30,50]
	_,_,u200_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,u200)
	EASJ = stats.nanmean(stats.nanmean(u200_clipped,axis = 2),axis =1)
	return EASJ

def Asia_Summer_Monsson_Index(lon,lat,u850):
	region = [60,120,5,20]
	_,_,u850_ASMI1 = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	region = [110,140,20,30]
	_,_,u850_ASMI2 = range_clip(region[0],region[1],region[2],region[3],lon,lat,u850)
	ASMI = stats.nanmean(stats.nanmean(u850_ASMI1,axis = 2),axis =1)-stats.nanmean(stats.nanmean(u850_ASMI2,axis = 2),axis =1)
	return ASMI
	
def Uniform_monsoon_index(time,lat,lon,u,v):
	'''
	refer to Li, J., and Q. Zeng, 2002: A unified monsoon index.
	Geophys. Res. Lett.,29, NO. 8, 1274, 10.1029/2001GL013874.
	input:monthly u and v at the same pressure level 
	oytput: the monsoon index at each point
	'''
	region = [59,151,-15,56] #  Asia Monsoon Region with 1 degree margins
	lons,lats,u_clipped =  range_clip(region[0],region[1],region[2],region[3],lon,lat,u)
	_,_,v_clipped =  range_clip(region[0],region[1],region[2],region[3],lon,lat,v)
	
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	
	U_Jan_Clima = stats.nanmean([u_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==101],axis=0)
	V_Jan_Clima = stats.nanmean([v_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==101],axis=0)
	U_June = [u_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==601]
	V_June = [v_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==601]
	U_July = [u_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==701]
	V_July = [v_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==701]
	U_Aug = [u_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==801]
	V_Aug = [v_clipped[lev] for lev in range(len(time)) if time[lev]%10000 ==801]
	U_Jul_Clima = stats.nanmean(U_July,axis =0)
	V_Jul_Clima = stats.nanmean(V_July,axis =0)
	
	# the member of the fraction
	def norm_value(year_series,lon,lat,u,v):
		lon_res= np.deg2rad(lon[1]-lon[0])
		lat_res= np.deg2rad(lat[1]-lat[0])
		lat = np.deg2rad(lat)
		NormValue = np.empty((np.shape(u))); NormValue[:]=np.nan
		ws2 = np.power(u,2)+np.power(v,2) # wind speed aquare
		for iyear in year_series:
			for row in range(1,len(lat)-1):
				for column in range(1,len(lon)-1):
					NormValue[iyear-year_series[0],row,column] = 6371.22/4*lon_res*lat_res*\
					np.sqrt(math.cos(lat[row])*(ws2[iyear-year_series[0],row,column-1]+4*ws2[iyear-year_series[0],row,column]+ws2[iyear-year_series[0],row,column+1])\
					+math.cos(lat[row-1])*ws2[iyear-year_series[0],row-1,column]\
					+math.cos(lat[row+1])*ws2[iyear-year_series[0],row+1,column])
		return NormValue
	# the dominator of the fraction
	def norm_value2D(lon,lat,u,v):
		lon_res= np.deg2rad(lon[1]-lon[0])
		lat_res= np.deg2rad(lat[1]-lat[0])
		lat = np.deg2rad(lat)
		NormValue = np.empty((np.shape(u))); NormValue[:]=np.nan
		ws2 = np.power(u,2)+np.power(v,2) # wind speed aquare
		for row in range(1,len(lat)-1):
			for column in range(1,len(lon)-1):
				NormValue[row,column] =  6371.22/4*lon_res*lat_res*\
				np.sqrt(math.cos(lat[row])*(ws2[row,column-1]+4*ws2[row,column]+ws2[row,column+1])\
				+math.cos(lat[row-1])*ws2[row-1,column]\
				+math.cos(lat[row+1])*ws2[row+1,column])
		return NormValue
		
	# print np.shape(member_u)
	member_u = -(U_June-U_Jan_Clima); member_v = -(V_June-V_Jan_Clima)
	June_member = norm_value(year_series,lons,lats,member_u,member_v)
	member_u = -(U_July-U_Jan_Clima); member_v = -(V_July-V_Jan_Clima)
	July_member = norm_value(year_series,lons,lats,member_u,member_v)
	member_u = -(U_Aug-U_Jan_Clima); member_v = -(V_Aug-V_Jan_Clima)
	Aug_member = norm_value(year_series,lons,lats,member_u,member_v)
	NV_member =(June_member+July_member+Aug_member)/3
	
	denominator_u = (U_Jan_Clima+U_Jul_Clima)/2; denominator_v = (V_Jan_Clima+V_Jul_Clima)/2;
	NV_denominator = norm_value2D(lons,lats,denominator_u,denominator_v)
	# print np.shape(NV_member)
	UMI= np.divide(NV_member,NV_denominator)-2
	return lons, lats, UMI
		 
######################################
# Index computing
######################################
file_path= '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'

u_file = 'ensumble_mean_U_200602_210101.nc'
nc_fid = nc4.Dataset(file_path+u_file,mode='r')
time = nc_fid.variables['time'][:]; lat = nc_fid.variables['lat'][:]; lev = nc_fid.variables['lev'][:]; lon = nc_fid.variables['lon'][:]
# u = nc_fid.variables['rcp85'][:]
u = nc_fid.variables['rcp85'][:]   #_fixA
nc_fid.close()

v_file = 'ensumble_mean_V_200602_210101.nc'
nc_fid = nc4.Dataset(file_path+v_file,mode='r')
# v_850 = nc_fid.variables['rcp85'][:,23,:,:]
v_850 = nc_fid.variables['rcp85'][:,23,:,:] #_fixA
nc_fid.close()


year= [ value/10000 for value in map(int,time)]
year_series = [value for value in np.unique(year) if value <=2100]

def JJA_mean(data):
	'''
	here data must be a 3D arrqay
	'''
	result = np.empty((95,192,288))
	layer_e = -6  #The Cesm data begins from Feb of 2006
	for iyear in year_series:
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = data[layer_b:layer_e+1,:,:]
		result[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	return result
u_850 = u[:,23,:,:]; u_200 = u[:,12,:,:]
u_JJA_850= JJA_mean(u_850); u_JJA_200= JJA_mean(u_200);v_JJA_850= JJA_mean(v_850);

WYI = Webster_Yang_monsoon_index(lon,lat,u_JJA_850,u_JJA_200)
SAMI=South_Asian_monsoon_index(lon,lat,u_JJA_850)
DIMI = Dynamic_Indian_monsoon_index(lon,lat,u_JJA_850)
WNPMI = WNP_monsoon_index(lon,lat,u_JJA_850)
EASI = East_Asia_Monsoon_index(lon,lat,u_JJA_850)
EASJ = East_Asian_Summer_Jet(lon,lat,u_JJA_200)
ASMI = Asia_Summer_Monsson_Index(lon,lat,u_JJA_850)
lon_UMI, lat_UMI, UMI= Uniform_monsoon_index(time,lat,lon,u_850,v_850)


#################################
# writting to nc file
#################################
file_name_out = file_path + 'CESM_rcp85_Monsoon_index.nc'
f = nc4.Dataset(file_name_out,'w', format='NETCDF4') #'w' stands for write
f.createDimension('lat', len(lat_UMI))
f.createDimension('lon', len(lon_UMI))
f.createDimension('time',len(year_series))
# f.createDimension('lev', len(lev))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
# levs = f.createVariable('lev',np.float32, ('lev'))
WYIs = f.createVariable('WYI',np.float32, ('time'))
SAMIs = f.createVariable('SAMI',np.float32, ('time'))
DIMIs = f.createVariable('DIMI',np.float32, ('time'))
WNPMIs = f.createVariable('WNPMI',np.float32, ('time'))
EASIs = f.createVariable('EASI',np.float32, ('time'))
ASMIs = f.createVariable('ASMI',np.float32, ('time'))
UMIs = f.createVariable('UMI',np.float32, ('time','lat','lon'))
EASJs = f.createVariable('ASMJ',np.float32, ('time'))

times[:] = year_series; latitudes[:] = lat_UMI; longitudes[:] = lon_UMI; 
WYIs[:] = WYI; SAMIs[:] = SAMI; DIMIs[:] = DIMI; WNPMIs[:] = WNPMI; 
EASIs[:] = EASI; ASMIs[:] = ASMI; UMIs[:] = UMI; EASJs[:] = EASJ

times.long_name = 'Year'
times.units = 'Year'
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'
longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'


WYIs.units='m/s'
WYIs.standard_name ='Webster_Yang_monsoon_index'
WYIs.long_name ='U850-U200 averaged over 0-20N, 40-110E)'
SAMIs.units='m/s'
SAMIs.standard_name ='South_Asian_monsoon_index'
SAMIs.long_name ='V850-V200 averaged over 10-30N, 70-110E)'
DIMIs.units='m/s'
DIMIs.standard_name =' Dynamic_Indian_monsoon_index'
DIMIs.long_name ='(U850 (5-15N, 40-80E)) - (U850 (20-30N, 70-90E))'
WNPMIs.units='m/s'
WNPMIs.standard_name ='WNP_monsoon_index'
WNPMIs.long_name ='(U850 (5-15N, 100-130E)) - (U850 (20-35N, 110-140E))'
EASIs.units='m/s'
EASIs.standard_name ='East_Asian_Summer_monsoon_index'
EASIs.long_name ='V850(20-40N, 120-140E)'
EASJs.units = 'm/s'
EASJs.standard_name ='East_Asian_Summer_Jet'
EASJs.long_name ='U200(30-50N, 110-140E)'
ASMIs.units='m/s'
ASMIs.standard_name ='Asia_Summer_monsoon_index'
ASMIs.long_name ='(U850 (5-20N, 60-120E)) - (U850 (20-30N, 110-140E))'
UMIs.units='m/s'
UMIs.standard_name ='A unified East-Asia_Summer_monsoon_index'
UMIs.long_name ='Refer to Jianping li and Zeng(2002)'

import time as clock
f.description = '2006-2100 monsoon index from CESM model'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()
