# -*- coding: utf-8 -*-
"""
Sum up all aerosol species (BC, POM,SOA,SEA-SALT,DUST,SO4) to produce PM2.5
using deaily model output from the CAM-CEHM runs
For sea-salt and  dust, only a poertion (45.25%, based on the Gaussian distribution (mean =2.225, dtd=1.80))
are included
"""
import site
import os
import time as clock
import numpy as np
import netCDF4 as nc4

def netcdf4_write(file_name,time,lon,lat,data,data_units,time_units):

	f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	
	f.createDimension('time', len(time))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	
	times = f.createVariable('time',np.float64, ('time'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	datas = f.createVariable('PM2p5',np.float32,('time','lat','lon'))
	times[:] = time
	latitudes[:] = lat
	longitudes[:] = lon
	datas[:] = data
	
	
	times.long_name = 'time'
	times.units = time_units
	times.calender='noleap'
	
	latitudes.long_name = 'latitude'
	latitudes.units = 'degree_north'
	longitudes.long_name = 'longitude'
	longitudes.units = 'degree_east'
	
	txxs.standard_name = 'PM2.5'
	txxs.long_name = 'PM2.5 in kg/kg from all aerosol species from MAM3'
	datas.units =data_units

	f.description = 'PM2.5 from the CAM-CHEM runs,all species ll aerosol species (BC, POM,SOA,SEA-SALT,DUST,SO4) to produce PM2.5 using deaily model output from the CAM-CEHM runs; For sea-salt and  dust, only a poertion (45.25%, based on the Gaussian distribution (mean =2.225, dtd=1.80)) are included'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'	
	f.close()
	
	
	
def get_all_species(scenario):
	def data_netcdf(scenario,variable):
		var_path = '/nerc/n02/n02/alcide/cesm122_AC_pp/'+scenario+'day/atm/'+scenario+'.atm.day.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		time = nc_fid.variables['time'][:];
		time_units = nc_fid.variables['time'].units
		data = nc_fid.variables[variable][:]
		data_units =  nc_fid.variables[variable][:].units
		return time,lon,lat,data,data_units,time_units
	variable ='bc_a1_SRF';	time,lon,lat,bc_a1,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='pom_a1_SRF';	time,lon,lat,pom_a1,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='soa_a1_SRF';	time,lon,lat,soa_a1,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='soa_a2_SRF';	time,lon,lat,soa_a2,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='so4_a1_SRF';	time,lon,lat,so4_a1,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='so4_a2_SRF';	time,lon,lat,so4_a2,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='so4_a3_SRF';	time,lon,lat,so4_a3,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='dst_a1_SRF';	time,lon,lat,dst_a1,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='dst_a3_SRF';	time,lon,lat,dst_a3,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='ncl_a1_SRF';	time,lon,lat,ncl_a1,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='ncl_a2_SRF';	time,lon,lat,ncl_a2,data_units,time_units =  data_netcdf(scenario,variable)
	variable ='ncl_a3_SRF';	time,lon,lat,ncl_a3,data_units,time_units =  data_netcdf(scenario,variable)
	sum = bc_a1+pom_a1+soa_a1+soa_a2+so4_a1+so4_a2+so4_a3+dst_a1+ncl_a1+ncl_a2+(dst_a3+ncl_a3)*0.4525
	return time,lon,lat,sum,data_units,time_units

for scenario in ["EdgTechP","EdgEneP","EdgRef70AP","EdgRefP","Edg1970"]
	time,lon,lat,sum,data_units,time_units = get_all_species(scenario)
	file_nam = '/nerc/n02/n02/alcide/cesm122_AC_pp/'+scenario+'day/atm/'+scenario+'.atm.day.PM2.5_summed.nc'
	netcdf4_write(file_name,time,lon,lat,sum,data_units,time_units)
