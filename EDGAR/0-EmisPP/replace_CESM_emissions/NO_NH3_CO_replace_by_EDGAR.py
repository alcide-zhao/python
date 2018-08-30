"""
Replace RCP8.5 2010 emissions with that from EDGAR4.3.1.
	1. replace the monthly data one after the other
	2. For one mode, all emissions are put in just the anthropogenic sector, all other sectors (ocean,soil,bb) are fixed the same as RCP85,
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
emission = '1970' #'STAG_ENERGY','STAG_TECH' '2010'
months=range(0,12)
path_EDGAR = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/'+emission+'/'
path_rcp85 = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/CESM_rcp85/1970/'
fn_dic = {'1970':'REFERENCE','2010':'REFERENCE','STAG_TECH':'STAG_TECH','STAG_ENERGY':'STAG_ENE'}
year_dic={'1970':'1970','2010':'2010','STAG_TECH':'2010','STAG_ENERGY':'2010'}
# layer_b = 97  #2010
layer_b = 216   #1970
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

path_CESM = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/CESM_rcp85/'
so2_surf_file=path_CESM+'rcp85_finn_2002-2015_NO_1.9x2.5_mol_c20160128.nc'
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
	

####CO_surf
file=path_rcp85+'maccity_maccity_corrdates_CO_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_surf = nc4.Dataset(file,mode='r+')
for imonth in months:
	edgar_file = path_EDGAR+fn_op+'_CO_'+str(year)+'_'+str(imonth+1)+'_MOL.0.1x0.1.nc'  #filename month starts from 1
	edgar_fid=nc4.Dataset(edgar_file,mode='r')
	CO_srf = edgar_fid.variables['CO_srf'][:] 
	lon = edgar_fid.variables['lon'][:]
	lat = edgar_fid.variables['lat'][:] 
	CO_srf=regrid_edgar2cesm(lon,lat,CO_srf,CESM_gc_area)
	nc_surf.variables['emiss_anthro'][layer_b+imonth,:,:]=CO_srf	
	edgar_fid.close()
nc_surf.close()

####NO_surf
file=path_rcp85+'maccity_maccity_corrdates_NO_1960-2008_1.9x2.5_mol_c130314.nc'
nc_surf = nc4.Dataset(file,mode='r+')
for imonth in months:
	edgar_file = path_EDGAR+fn_op+'_NOx_'+str(year)+'_'+str(imonth+1)+'_MOL.0.1x0.1.nc'  #filename month starts from 1
	edgar_fid=nc4.Dataset(edgar_file,mode='r')
	NO_srf = edgar_fid.variables['NO_srf'][:] 
	lon = edgar_fid.variables['lon'][:]
	lat = edgar_fid.variables['lat'][:] 
	NO_srf=regrid_edgar2cesm(lon,lat,NO_srf,CESM_gc_area)
	nc_surf.variables['emiss_anthro'][layer_b+imonth,:,:]=NO_srf	
	edgar_fid.close()
nc_surf.close()


####NH3_surf
file=path_rcp85+'maccity_maccity_corrdates_NH3_1960-2008_1.9x2.5_mol_c130314.nc'
nc_surf = nc4.Dataset(file,mode='r+')
for imonth in months:
	edgar_file = path_EDGAR+fn_op+'_NH3_'+str(year)+'_'+str(imonth+1)+'_MOL.0.1x0.1.nc'  #filename month starts from 1
	edgar_fid=nc4.Dataset(edgar_file,mode='r')
	NH3_srf = edgar_fid.variables['NH3_srf'][:] 
	lon = edgar_fid.variables['lon'][:]
	lat = edgar_fid.variables['lat'][:] 
	NH3_srf=regrid_edgar2cesm(lon,lat,NH3_srf,CESM_gc_area)
	nc_surf.variables['emiss_anthro'][layer_b+imonth,:,:]=NH3_srf	
	edgar_fid.close()
nc_surf.close()
