"""
This is to process EDGAR NMVOC_bio and NMVOC_Fossil into number of molecules
			
"""
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import glob 
import os
import matplotlib.pyplot as plt


emission='1970' #2010 or stag_tech stag_energy
Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH'] #,'CDS','CRS','LTO','SPS'
NA=6.022*10**(23)		#NA This factor (= 1.0e3 * avogadro_number) is needed to work with the CAM3 surface emissions code

## INPUT AND OUTPUT FILE DIRECTORY
file_path ="/exports/csce/datastore/geos/users/s1667168/EDGAR/"
fn_dic = {'1970':'REFERENCE','2010':'REFERENCE','stag_tech':'STAG_TECH','stag_energy':'STAG_ENE'}
fn_op = fn_dic[emission]

months=range(1,13) 

for month in months:
	bio_share=np.zeros((1800,3600));fossil_share=np.zeros((1800,3600));
	bio_scenario=np.zeros((1800,3600));fossil_scenario=np.zeros((1800,3600));
	input_path = "/exports/csce/datastore/geos/users/s1667168/EDGAR/"+emission+'/NMVOC/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*_NMVOC*1970_' +str(month)+'_*.nc' + '" -print | sort > ' + input_path +'file_list'+str(month)+'.txt')
	text_file = open(input_path + 'file_list'+str(month)+'.txt', "r")
	text_content = text_file.readlines()
	# print 'bio'
	for slot in range(0,13):
		nc_f = text_content[slot][:-1]
		# print nc_f
		nc_fid = nc4.Dataset(nc_f,mode='r')
		lat = nc_fid.variables['lat']
		lon = nc_fid.variables['lon']
		bio_scenario =bio_scenario + nc_fid.variables['emi_nmvoc'][:]
	# print 'fossil'
	bio = (bio_share+bio_scenario)*0.1   # 0.1 from kg/m2/s in EDGAR into molecules/cm2/s in CESM  1000/100/100/
	fossil = (fossil_share+fossil_scenario)*0.1   # 0.1 from kg/m2/s in EDGAR into molecules/cm2/s in CESM  1000/100/100/
	BIGALK = (fossil*0.712+bio*0.491)*0.05*1.5/1.4/12.0*NA   # SOAG yield 5% for BIGALK
	BIGENE = (fossil*0.136+bio*0.352)*0.05*1.5/1.4/12.0*NA   # SOAG yield 5% for BIGENE
	TOLUEN = (fossil*0.153+bio*0.158)*0.15*1.5/1.4/12.0*NA   # SOAG yield 15% for TOLUEN
	# output results into .nc file 
	file_out = "/exports/csce/datastore/geos/users/s1667168/EDGAR/"+emission+'/'+fn_op+'_NMVOC_1970_'+str(month)+'_molecules.0.1x0.1.nc'
	f = nc4.Dataset(file_out,'w', format='NETCDF4') 
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	latitudes[:] = lat
	longitudes[:] = lon		
	BIGALKs=f.createVariable('BIGALK',np.float32,('lat','lon'))	
	BIGENEs=f.createVariable('BIGENE',np.float32,('lat','lon'))
	TOLUENs=f.createVariable('TOLUEN',np.float32,('lat','lon'))

	BIGALKs[:]=BIGALK; BIGENEs[:]=BIGENE; TOLUENs[:]=TOLUEN;
	
	
	BIGALKs.units='molecules/cm2/s'
	BIGALKs.long_name ='SOAG converted from BIGALK of NMVOC_fossil and NMVOC_bio of EDGAR 4.3.1, mass yield 5%'
	
	BIGENEs.units='molecules/cm2/s'
	BIGENEs.long_name ='SOAG converted from BIGENE of NMVOC_fossil and NMVOC_bio of EDGAR 4.3.1, mass yield 5%'

	TOLUENs.units='molecules/cm2/s'
	TOLUENs.long_name ='SOAG converted from TOLUEN of NMVOC_fossil and NMVOC_bio of EDGAR 4.3.1, mass yield 15%'

	
	f.title='SOAG, mass yield is increased by a factor of 1.5'
	f.close()
		


