"""
This is to process BC eissions from EDGAR into BC and OC into accumulation mode
surface and elevated files
1.4*OC mass is taken as POM.
					
"""

import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import glob  
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d

# # AGR: Agriculture
# # ENERGY: Power generation
# # INDUSTRY: Manufacturing industry
# # SHIP: Shipping emissions
# # RCO: Residential combustion
# # SWD: Solid waste disposal
# # TNG: Non road transport
# # TRO: Road transport
# # AVI: Aviation transport

emission='1970' #2010 or stag_tech stag_energy
specie='OC' 
Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH'] 

NA=6.022*10**(23)

## mass of BC particles
densities={'SO4':1.77,'BC':1.7,'POM':1.0}   # conversion to number emissions uses dey density of 1.77 1.7 1.0 g/cm3 for BC, BC and POM   g/cm3
density=densities['POM'];
diameter=0.134*10**(-4);   #convert from micrometer into centimeter
volume=np.pi*np.power(diameter,3)/6.0;
m_p = volume*density

# input and output file directories
file_path ="/exports/csce/datastore/geos/users/s1667168/EDGAR/"
fn_dic = {'1970':'REFERENCE','2010':'REFERENCE','stag_tech':'STAG_TECH','stag_energy':'STAG_ENE'}
fn_op = fn_dic[emission]

## read the CESM grids to interpolate the RCP emissionsn firstly before caculations
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
# CESM_fid.close()

months=range(1,13) 
for month in months:
	POM_a1_srf_mol=np.zeros((1800,3600))
	POM_a1_srf_num=np.zeros((1800,3600))
	for sector in Sectors:
		# print sector
		if emission in ['stag_tech','stag_energy']:
			if sector in ['ENE','IND','TRO']:
				file = file_path+emission+"/"+specie+'/'+'JRC_PEGASOS_V2_'+emission.upper()+'_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
			else:
				file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc' #
		else:
			file = file_path+"1970/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_1970_'+str(month)+'_'+sector+'.nc' #
		# print file
		nc_fid = nc4.Dataset(file,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		# 1.4 is to convert OC into POM 
		POM_mass = nc_fid.variables['emi_'+specie.lower()][:]*0.1*1.4    # kg/m2/s into g/cm2/s in CESM  1000/100/100/12
		# f = interp2d(lon,lat,POM,kind='linear')
		# POM_mass = f(lon_CESM, lat_CESM);
		POM_a1_srf_mol=POM_a1_srf_mol+POM_mass*NA/12.0
		POM_a1_srf_num=POM_a1_srf_num+POM_mass*10**(26)/m_p
		nc_fid.close()
	# output results into .nc file 
	file_out = file_path+emission+'/'+fn_op+'_'+specie+'_1970_'+str(month)+'_NUM_MOL.0.1x0.1.nc'
	f = nc4.Dataset(file_out,'w', format='NETCDF4') #'w' stands for write
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	latitudes[:] = lat
	longitudes[:] = lon		
	POM_a1_srf_mols=f.createVariable('POM_a1_srf_mol',np.float32,('lat','lon'))	
	POM_a1_srf_nums=f.createVariable('POM_a1_srf_num',np.float32,('lat','lon'))

	POM_a1_srf_mols[:]=POM_a1_srf_mol
	POM_a1_srf_nums[:]=POM_a1_srf_num	
	
	POM_a1_srf_mols.units='particles/cm2/s'
	POM_a1_srf_mols.long_name ='POM surface emission in molecule numbers, with all sectors added up,1.4 factor applied'


	POM_a1_srf_nums.units='particles/cm2/s*6.0e26'
	POM_a1_srf_nums.long_name ='POM surface emission in particle numbers, with all sectors added up, 1.4 factor applied'

	f.title='POM emissions molecules/particles numbers for a1 surface level'
	f.close()

