"""
This is to compare the grid-area weighted aerosol 
in number of molecules and particles between the four scenario
	ref207,ref2010,stag_tech,stag_ene
	
"""

import netCDF4 as nc4
import numpy as np
import math
import sys
#from mpl_toolkits.basemap import interp
import scipy.io as sio

data_path ='/work/n02/n02/alcide/cesm1_2_2/inputdata/atm/cam/chem/trop_mozart_aero/emis/EDGAR/'
scenarios =['ref1970','ref2010','stag_tech','stag_ene']
layer_dic ={'ref1970':156,'ref2010':24,'stag_tech':24,'stag_ene':24}
region_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,60],'EA':[100,145,20,50],'SA':[65,100,5,30],'EUR':[0,40,48,75],'NAM':[210,300,0,60]}
region = region_dic['GLOBE']

def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area


so2_surf_file=data_path+'ref2010/RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
lons,lats = np.meshgrid (lon_CESM,lat_CESM)
lat_res = lat_CESM[1] - lat_CESM[0];lon_res = lon_CESM[1] - lon_CESM[0];
CESM_gc_area = getArea(lons,lons+lon_res,lats,lats+lat_res)

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
	

###SO2 and SO4
secotrs_so24={'so2_surf':'emiss_awb','so2_elev':'emiss_ene','so4_a1_surf':'emiss_awb','so4_a1_elev':'emiss_ene','so4_a2_surf':'emiss_dom'}
substance_so24={'so2_surf':np.zeros((4,12)),'so2_elev':np.zeros((4,12)),'so4_a1_surf':np.zeros((4,12)),'so4_a1_elev':np.zeros((4,12)),'so4_a2_surf':np.zeros((4,12))}

for sector in secotrs_so24:
	nc_variable = secotrs_so24[sector]
	scenario_layer = 0
	for scenario in scenarios:
		layer_b = layer_dic[scenario]
		if (scenario == 'ref1970'):
			file = data_path+scenario+'/'+'ar5_mam3_'+sector+'_1850-2005_c090804.nc'			
		else:
			file = data_path+scenario+'/'+'RCP85_mam3_'+sector+'_2000-2100_c20110913.nc'
		# print (file)
		if sector in ['so2_elev', 'so4_a1_elev']:
			nc_fid=nc4.Dataset(file,mode='r')
			variable = nc_fid.variables[nc_variable][layer_b:layer_b+12,0,:,:]*12600/0.13
			nc_fid.close()
		else:
			nc_fid=nc4.Dataset(file,mode='r')
			variable = nc_fid.variables[nc_variable][layer_b:layer_b+12,:,:]
			nc_fid.close()			
		area_weighted = np.multiply(variable,CESM_gc_area)
		_,_,area_weight_clip = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,area_weighted)
		monthly_variable = np.nansum(np.nansum(area_weight_clip,axis=2),axis=1)
		substance_so24[sector][scenario_layer,:]=monthly_variable; scenario_layer=scenario_layer+1;


sio.savemat(data_path+'sulphate_molecules.mat',\
{'so2_surf':substance_so24['so2_surf'][:],'so2_elev':substance_so24['so2_elev'][:],'so4_a1_surf':substance_so24['so4_a1_surf'][:],'so4_a1_elev':substance_so24['so4_a1_elev'][:],'so4_a2_surf':substance_so24['so4_a2_surf'][:]})


### BC & OC
secotrs_BOC={'bc_surf':'emiss_awb','oc_surf':'emiss_awb'}
substance_BOC={'bc_surf':np.zeros((4,12)),'oc_surf':np.zeros((4,12))}

for sector in secotrs_BOC:
	nc_variable = secotrs_BOC[sector]
	scenario_layer = 0
	for scenario in scenarios:
		layer_b = layer_dic[scenario]
		if (scenario == 'ref1970'):
			file = data_path+scenario+'/'+'ar5_mam3_'+sector+'_1850-2005_c090804.nc'			
		else:
			file = data_path+scenario+'/'+'RCP85_mam3_'+sector+'_2000-2100_c20110913.nc'
		# print (file)
		nc_fid=nc4.Dataset(file,mode='r')
		variable = nc_fid.variables[nc_variable][layer_b:layer_b+12,:,:]
		nc_fid.close()			
		area_weighted = np.multiply(variable,CESM_gc_area)
		_,_,area_weight_clip = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,area_weighted)
		monthly_variable = np.nansum(np.nansum(area_weight_clip,axis=2),axis=1)
		substance_BOC[sector][scenario_layer,:]=monthly_variable; scenario_layer=scenario_layer+1;

## soag

soag=np.zeros((4,12))
scenario_layer = 0
for scenario in scenarios:
	layer_b = layer_dic[scenario]
	if (scenario == 'ref1970'):
		file = data_path+scenario+'/ar5_mam3_soag_1.5_surf_1850-2005_c100429.nc'			
	else:
		file = data_path+scenario+'/RCP85_soag_1.5_surf_2000-2100_c20110913.nc'
	# print file
	nc_fid=nc4.Dataset(file,mode='r')
	variable = nc_fid.variables['SOAG_BIGALK'][layer_b:layer_b+12,:,:]+\
			nc_fid.variables['SOAG_BIGENE'][layer_b:layer_b+12,:,:]+\
			nc_fid.variables['SOAG_ISOPRENE'][layer_b:layer_b+12,:,:]+\
			nc_fid.variables['SOAG_TERPENE'][layer_b:layer_b+12,:,:]+\
			nc_fid.variables['SOAG_TOLUENE'][layer_b:layer_b+12,:,:]
	nc_fid.close()			
	area_weighted = np.multiply(variable,CESM_gc_area)
	_,_,area_weight_clip = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,area_weighted)
	monthly_variable = np.nansum(np.nansum(area_weight_clip,axis=2),axis=1)
	soag[scenario_layer,:]=monthly_variable; scenario_layer=scenario_layer+1;

sio.savemat(data_path+'BC_OC_soag_molecules.mat',\
{'bc_surf':substance_BOC['bc_surf'][:],'oc_surf':substance_BOC['oc_surf'][:],'soag':soag[:]})

### no of particles
fle_particles={'SO4_a1_surf':'num_a1_surf','SO4_a1_elev':'num_a1_elev','SO4_a2_surf':'num_a2_surf','OC_a1_surf':'num_a1_surf','BC_a1_surf':'num_a1_surf'}
secotrs_particles={'SO4_a1_surf':'SO4_emiss_awb','SO4_a1_elev':'SO4_emiss_ene','SO4_a2_surf':'SO4_emiss_dom','OC_a1_surf':'OC_emiss_awb','BC_a1_surf':'BC_emiss_awb'}
substance_particles={'SO4_a1_surf':np.zeros((4,12)),'SO4_a1_elev':np.zeros((4,12)),'SO4_a2_surf':np.zeros((4,12)),'OC_a1_surf':np.zeros((4,12)),'BC_a1_surf':np.zeros((4,12))}

####
for sector in secotrs_particles:
	nc_variable= secotrs_particles[sector]
	scenario_layer = 0
	for scenario in scenarios:
		layer_b = layer_dic[scenario]
		if (scenario == 'ref1970'):
			file = data_path+scenario+'/'+'ar5_mam3_'+fle_particles[sector]+'_1850-2005_c090804.nc'			
		else:
			file = data_path+scenario+'/'+'RCP85_mam3_'+fle_particles[sector]+'_2000-2100_c20110913.nc'
		if sector in ['SO4_a1_elev']:
			nc_fid=nc4.Dataset(file,mode='r')
			variable = nc_fid.variables[nc_variable][layer_b:layer_b+12,0,:,:]*12600/0.13
			nc_fid.close()
		else:
			nc_fid=nc4.Dataset(file,mode='r')
			variable = nc_fid.variables[nc_variable][layer_b:layer_b+12,:,:]
			nc_fid.close()			
		area_weighted = np.multiply(variable,CESM_gc_area)
		_,_,area_weight_clip = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,area_weighted)
		monthly_variable = np.nansum(np.nansum(area_weight_clip,axis=2),axis=1)
		substance_particles[sector][scenario_layer,:]=monthly_variable;
		scenario_layer=scenario_layer+1;


sio.savemat(data_path+'aerosol_particles.mat',\
{'SO4_a1_surf':substance_particles['SO4_a1_surf'][:],'SO4_a1_elev':substance_particles['SO4_a1_elev'][:],'SO4_a2_surf':substance_particles['SO4_a2_surf'][:],'OC_a1_surf':substance_particles['OC_a1_surf'][:],'BC_a1_surf':substance_particles['BC_a1_surf'][:]})
	
