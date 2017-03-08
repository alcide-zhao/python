# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to priduce rhe 	RCP 85 EMISSIONS INTO TIME EVELUTION SERIES
@author: Alcide.Zhao
"""
import site
import os
import time as clock
import numpy as np
import netCDF4 as nc4
from scipy import stats

########################################################
#0. setting variables
########################################################

input_path = '/exports/csce/datastore/geos/users/s1667168/RCP/'
species_list = ('OC','SO2')
sector_list = ('anthropogenic','ships','biomassburning')
SCENARIO='RCP85'
RES='0.5x0.5'

########################################################
#2. funtions
########################################################
def netcdf4_write(file_name,value_dic,units_dic,molecular_weight_dic):
	# print np.shape(value_dic['time'])
	f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	
	f.createDimension('time', 1152)
	f.createDimension('lat', 360)
	f.createDimension('lon', 720)
	
	times = f.createVariable('time',np.float64, ('time'))
	lats = f.createVariable('lat',np.float32, ('lat'))
	lons = f.createVariable('lon',np.float32, ('lon'))
	emiss_anthropogenic=f.createVariable('anthropogenic',np.float32,('time','lat','lon'))
	emiss_ships=f.createVariable('ships',np.float32,('time','lat','lon'))
	emiss_biomass_burnings=f.createVariable('biomass_burning',np.float32,('time','lat','lon'))

	times[:] = value_dic['time']
	lats[:] = value_dic['lat']
	lons[:] = value_dic['lon']
	
	emiss_anthropogenic[:]=value_dic['anthropogenic']
	emiss_ships[:]=value_dic['ships']
	emiss_biomass_burnings[:]=value_dic['biomass_burning']

	times.units = units_dic['time']
	lats.units = units_dic['lat']
	lons.units = units_dic['lon']
	
	emiss_anthropogenic.units=units_dic['anthropogenic']
	emiss_ships.units=units_dic['ships']
	emiss_biomass_burnings.units=units_dic['biomass_burning']

	emiss_anthropogenic.molecular_value=molecular_weight_dic['anthropogenic']
	emiss_ships.molecular_value=molecular_weight_dic['ships']
	emiss_biomass_burnings.molecular_value=molecular_weight_dic['biomass_burning']
	
	f.description = 'RCP85 emissions'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	f.close()
	return





########################################################
#2. mian program
########################################################
value_dic = {'time':[],'lon':[],'lat':[],'anthropogenic':[],'ships':[],'biomass_burning':[]}
units_dic = {'time':[],'lon':[],'lat':[],'anthropogenic':[],'ships':[],'biomass_burning':[]}
molecular_weight_dic = {'time':[],'lon':[],'lat':[],'anthropogenic':[],'ships':[],'biomass_burning':[]}

# units_value =np.empty((1152,360,720))
for SPECIES in species_list:
	time = np.empty((1152))
	emiss_awb = np.empty((1152,360,720));emiss_dom = np.empty((1152,360,720));emiss_ene = np.empty((1152,360,720))
	emiss_ind = np.empty((1152,360,720));emiss_tra = np.empty((1152,360,720));emiss_wst = np.empty((1152,360,720))
	emiss_shp = np.empty((1152,360,720));emiss_gra = np.empty((1152,360,720));emiss_for = np.empty((1152,360,720))
	layer_e = -1
	for YEAR in range(2005,2101):
		layer_b = layer_e+1
		layer_e =layer_b+11
		file_name = input_path+'accmip_interpolated_emissions_'+SCENARIO+'_'+SPECIES+'_anthropogenic_'+str(YEAR)+'_'+RES+'.nc'
		nc_fid = nc4.Dataset(file_name,mode='r')
		
		value_dic['lat'] = nc_fid.variables['lat'][:]
		units_dic['lat'] = nc_fid.variables['lat'].units
		
		value_dic['lon'] = nc_fid.variables['lon'][:]
		units_dic['lon'] = nc_fid.variables['lon'].units
		
		time[layer_b:layer_e+1] = nc_fid.variables['time'][:]
		units_dic['time'] = nc_fid.variables['time'].units
		
		emiss_awb[layer_b:layer_e+1] = nc_fid.variables['emiss_awb'][:]
		units_dic['anthropogenic'] = nc_fid.variables['emiss_awb'].units
		molecular_weight_dic['anthropogenic'] = nc_fid.variables['emiss_awb'].molecular_weight
		
		emiss_dom[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_dom'][:]	
		emiss_ene[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_ene'][:]
		emiss_ind[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_ind'][:]
		emiss_tra[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_tra'][:]
		emiss_wst[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_wst'][:]
		
		file_name = input_path+'accmip_interpolated_emissions_'+SCENARIO+'_'+SPECIES+'_ships_'+str(YEAR)+'_'+RES+'.nc'
		nc_fid = nc4.Dataset(file_name,mode='r')

		emiss_shp[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_shp'][:]
		units_dic['ships'] = nc_fid.variables['emiss_shp'].units
		molecular_weight_dic['ships'] = nc_fid.variables['emiss_shp'].molecular_weight				
		nc_fid.close()
		
		file_name = input_path+'accmip_interpolated_emissions_'+SCENARIO+'_'+SPECIES+'_biomassburning_'+str(YEAR)+'_'+RES+'.nc'
		nc_fid = nc4.Dataset(file_name,mode='r')

		emiss_gra[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_gra'][:]
		units_dic['biomass_burning'] = nc_fid.variables['emiss_gra'].units
		molecular_weight_dic['biomass_burning'] = nc_fid.variables['emiss_gra'].molecular_weight	
		
		emiss_for[layer_b:layer_e+1,:,:] = nc_fid.variables['emiss_for'][:]		
		nc_fid.close()
	
	value_dic['time'] = time;
	value_dic['anthropogenic'] = emiss_awb+emiss_dom+emiss_ene+emiss_ind+emiss_tra+emiss_wst
	value_dic['ships'] = emiss_shp
	value_dic['biomass_burning'] = emiss_gra+emiss_for
	# print np.shape(value_dic['anthropogenic'])
	file_name = input_path+'accmip_interpolated_emissions_'+SCENARIO+'_'+SPECIES+'_2005_2100_'+RES+'.nc'
	netcdf4_write(file_name,value_dic,units_dic,molecular_weight_dic)
		
		
		
		
		
