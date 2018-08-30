"""
Replace RCP8.5 2010 emissions with that from EDGAR4.3.1.
	1. replace the monthly data one after the other
	2. For one mode, all emissions are put in just one sector, all other sectors should be as zero,
		except for Volcanic, grassfire and forest fire
	3. Oxidants are kept unchanged, i.e. use that in RCP8.5 emussions
	4. SOAG is kept unchanged, i.e. the RCP8.5 one
"""

import netCDF4 as nc4
import numpy as np
import math
import sys
from mpl_toolkits.basemap import interp

# print sys.argv
# print len(sys.argv)
# emission =str(sys.argv); 
# print emission
emission = 'stag_energy'
months=range(0,12)
path_EDGAR = '/exports/csce/datastore/geos/users/s1667168/EDGAR/'+emission+'/'
path_rcp85 = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
fn_dic = {'1970':'REFERENCE','2010':'REFERENCE','stag_tech':'STAG_TECH','stag_energy':'STAG_ENE'}
year_dic={'1970':'1970','2010':'2010','stag_tech':'2010','stag_energy':'2010'}
layer_dic ={'1970':156,'2010':24,'stag_tech':24,'stag_energy':24}
layer_b = layer_dic[emission]
fn_op = fn_dic[emission]
year=year_dic[emission]


def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_num_a1_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
lons,lats = np.meshgrid (lon_CESM,lat_CESM)
lat_res = lat_CESM[1] - lat_CESM[0];lon_res = lon_CESM[1] - lon_CESM[0];
CESM_gc_area = getArea(lons,lons+lon_res,lats,lats+lat_res)


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
	

####################################
# molecule numbers
####################################

####BC_surf
# file=path_rcp85+'RCP85_mam3_bc_surf_2000-2100_c20110913.nc'
file=path_rcp85+'RCP85_mam3_bc_surf_2000-2100_c20110913.nc'
nc_surf = nc4.Dataset(file,mode='r+')
for imonth in months:
	edgar_file = path_EDGAR+fn_op+'_BC_'+str(year)+'_'+str(imonth+1)+'_NUM_MOL.0.1x0.1.nc'  #filename month starts from 1
	edgar_fid=nc4.Dataset(edgar_file,mode='r')
	BC_a1_surf_mol = edgar_fid.variables['BC_a1_srf_mol'][:] 
	lon = edgar_fid.variables['lon'][:]
	lat = edgar_fid.variables['lat'][:] 
	BC_a1_surf_mol=regrid_edgar2cesm(lon,lat,BC_a1_surf_mol,CESM_gc_area)
	nc_surf.variables['emiss_awb'][layer_b+imonth,:,:]=BC_a1_surf_mol
	nc_surf.variables['emiss_dom'][layer_b+imonth,:,:]=0
	nc_surf.variables['emiss_tra'][layer_b+imonth,:,:]=0
	nc_surf.variables['emiss_wst'][layer_b+imonth,:,:]=0
	nc_surf.variables['emiss_shp'][layer_b+imonth,:,:]=0
	nc_surf.variables['emiss_ene'][layer_b+imonth,:,:]=0
	nc_surf.variables['emiss_ind'][layer_b+imonth,:,:]=0	
	edgar_fid.close()
nc_surf.close()

####################################
# particle numbers
####################################

#### BC_a1_srf_num
# num_a1_surf_file=path_rcp85+'RCP85_mam3_num_a1_surf_2000-2100_c20110913.nc'
num_a1_surf_file=path_rcp85+'RCP85_mam3_num_a1_surf_2000-2100_c20110913.nc'
num_a1_surf = nc4.Dataset(num_a1_surf_file,mode='r+')
for imonth in months:
	edgar_file = path_EDGAR+fn_op+'_BC_'+str(year)+'_'+str(imonth+1)+'_NUM_MOL.0.1x0.1.nc'  #filename month starts from 1
	edgar_fid=nc4.Dataset(edgar_file,mode='r')
	BC_a1_srf_num= edgar_fid.variables['BC_a1_srf_num'][:] 
	lon = edgar_fid.variables['lon'][:]
	lat = edgar_fid.variables['lat'][:] 
	BC_a1_srf_num=regrid_edgar2cesm(lon,lat,BC_a1_srf_num,CESM_gc_area)
	num_a1_surf.variables['BC_emiss_awb'][layer_b+imonth,:,:]=BC_a1_srf_num
	num_a1_surf.variables['BC_emiss_wst'][layer_b+imonth,:,:]=0
	num_a1_surf.variables['BC_emiss_shp'][layer_b+imonth,:,:]=0
	num_a1_surf.variables['BC_emiss_dom'][layer_b+imonth,:,:]=0
	num_a1_surf.variables['BC_emiss_tra'][layer_b+imonth,:,:]=0
	num_a1_surf.variables['BC_emiss_ind'][layer_b+imonth,:,:]=0
	num_a1_surf.variables['BC_emiss_ene'][layer_b+imonth,:,:]=0
	edgar_fid.close()
num_a1_surf.close()

#### num_a2_elev are from volcanic



	





