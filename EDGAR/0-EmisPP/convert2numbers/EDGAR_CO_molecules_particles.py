"""
This is to map sum CO emissions from all the sectors in EDGAR 4.3.1			
"""
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import glob  
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d


species=['CO'] 
emission='STAG_TECH' #2010 or STAG_TECH STAG_ENERGY
NA=6.022*10**(23)		#NA This factor (= 1.0e3 * avogadro_number) is needed to work with the CAM3 surface emissions code
Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH','CDS','CRS','LTO','SPS']

## INPUT AND OUTPUT FILE DIRECTORY
file_path ="/exports/csce/datastore/geos/users/s1667168/EDGAR/"
fn_dic = {'2010':'REFERENCE','STAG_TECH':'STAG_TECH','STAG_ENERGY':'STAG_ENE'}
fn_op = fn_dic[emission]


months=range(1,13) 
for specie in species:
	for month in months:
		mass=np.zeros((1800,3600));
		for sector in Sectors:
			if emission in ['STAG_TECH','STAG_ENERGY']:
				if sector in ['ENE','IND','TRO']:
					file = file_path+emission+"/"+specie+'/'+'JRC_PEGASOS_V2_'+emission.upper()+'_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
			else:
				file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
			# print file
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			mass = mass+nc_fid.variables['emi_'+specie.lower()][:]		
			# nc_fid.close()
		# output results into .nc file 
		file_out = file_path+emission+'/'+fn_op+'_'+specie+'_2010_'+str(month)+'_MOL.0.1x0.1.nc'
		f = nc4.Dataset(file_out,'w', format='NETCDF4') 
		f.createDimension('lat', len(lat))
		f.createDimension('lon', len(lon))
		latitudes = f.createVariable('lat',np.float32, ('lat'))
		longitudes = f.createVariable('lon',np.float32, ('lon'))
		latitudes[:] = lat
		longitudes[:] = lon		
		CO_srfs=f.createVariable('CO_srf',np.float32,('lat','lon'))

		CO_srfs[:]=mass*NA*0.1/28.01          # kg/m2/s in EDGAR into molecules/cm2/s in CESM  1000/100/100/
			
		CO_srfs.units='molecules/cm2/s'
		CO_srfs.long_name ='CO surface emission in molecules'

		f.title='Co from all emissionss'
		f.close()
		
		


