import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  


# # AGR: Agriculture
# # ENE: Power generation
# # IND: Manufacturing industry
# # OTH: Other sources
# # PPA: Process emissions during production and application
# # PRO: Fuel production/transmission
# # RCO: Residential combustion
# # REF: Oil production and refining
# # SWD: Solid waste disposal
# # TNR: Non road transport
# # TRF: Transformation industry
# # TRO: Road transport

#'Merge EDGAR emission setors into one 
months=range(1,13) 
species=['OC'] #,'SO2','BC',
# sectors= ['CDS','CRS','LTO','SPS']  # Aviation
sectors=['IND','TRF','REF','PRO']   # INDUSTRY  ,'PPA'  no PPA for OC
# sectors=['ENE','OTH']     # ENERGY     
file_path ="/exports/csce/datastore/geos/users/s1667168/EDGAR/2010/"
for specie in species:
	for month in months:
		emis_total=np.zeros((1800,3600));
		for sector in sectors:
			file = file_path+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
			# print file
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]   
			emis_total =emis_total+ nc_fid.variables['emi_'+specie.lower()][:]
			sn=nc_fid.variables['emi_'+specie.lower()].standard_name
			ln=nc_fid.variables['emi_'+specie.lower()].long_name
			unit=nc_fid.variables['emi_'+specie.lower()].units
			cm=nc_fid.variables['emi_'+specie.lower()].cell_method
			title=nc_fid.title
			nc_fid.close()
		# output results into .nc file 
		file_out = file_path+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month)+'_INDUSTRY.0.1x0.1.nc'
		f = nc4.Dataset(file_out,'w', format='NETCDF4') #'w' stands for write
		f.createDimension('lat', len(lat))
		f.createDimension('lon', len(lon))
		latitudes = f.createVariable('lat',np.float32, ('lat'))
		longitudes = f.createVariable('lon',np.float32, ('lon'))
		emis = f.createVariable('emi_'+specie.lower(),np.float32,('lat','lon'))
		latitudes[:] = lat
		longitudes[:] = lon
		emis[:]= emis_total
		emis.units=unit
		emis.long_name =ln
		emis.standard_name =ln
		emis.cell_method =cm
		f.title=title
		f.close()