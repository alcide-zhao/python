"""
Replace RCP8.5 2010 emissions with that from EDGAR4.3.1.
"""

import netCDF4 as nc4
import numpy as np
import math
import sys
from mpl_toolkits.basemap import interp

# emission =str(sys.argv[:])
emission = '1970'
months=range(0,12)
path_EDGAR = '/exports/csce/datastore/geos/users/s1667168/EDGAR/'+emission+'/'
path_rcp85 = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
fn_dic = {'1970':'REFERENCE','2010':'REFERENCE','stag_tech':'STAG_TECH','stag_energy':'STAG_ENE'}
year_dic={'1970':'1970','2010':'2010','stag_tech':'2010','stag_energy':'2010'}
layer_dic ={'1970':156,'2010':24,'stag_tech':24,'stag_energy':24}
layer_b = layer_dic[emission]
year=year_dic[emission]
fn_op = fn_dic[emission]

def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
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
	emission_18layer_b = interp(emission, lon, lat, xout, yout, checkbounds=False, masked=False, order=1)
	lons,lats = np.meshgrid (x,y)
	lat_res = y[1] - y[0];lon_res = x[1] - x[0];
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	emission_gc_mass = np.multiply(emission_18layer_b,area)
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

# soag_file=path_rcp85+'RCP85_soag_1.5_surf_2000-2100_c20110913.nc'
soag_file=path_rcp85+'ar5_mam3_soag_1.5_surf_1850-2005_c100429.nc'
nc_surf = nc4.Dataset(soag_file,mode='r+')
for imonth in months:
	edgar_file = path_EDGAR+fn_op+'_NMVOC_'+str(year)+'_'+str(imonth+1)+'_molecules.0.1x0.1.nc'  #filename month starts from 1
	edgar_fid=nc4.Dataset(edgar_file,mode='r')
	BIGALK = edgar_fid.variables['BIGALK'][:]
	BIGENE = edgar_fid.variables['BIGENE'][:]
	TOLUEN = edgar_fid.variables['TOLUEN'][:]
	lon = edgar_fid.variables['lon'][:]
	lat = edgar_fid.variables['lat'][:] 
	BIGALK=regrid_edgar2cesm(lon,lat,BIGALK,CESM_gc_area)
	BIGENE=regrid_edgar2cesm(lon,lat,BIGENE,CESM_gc_area)
	TOLUEN=regrid_edgar2cesm(lon,lat,TOLUEN,CESM_gc_area)
	nc_surf.variables['SOAG_BIGALK'][layer_b+imonth,:,:]=BIGALK
	nc_surf.variables['SOAG_TOLUENE'][layer_b+imonth,:,:]=BIGENE
	nc_surf.variables['SOAG_BIGENE'][layer_b+imonth,:,:]=TOLUEN
	edgar_fid.close()
nc_surf.close()
	
	
	
	
	
	
	
	

	





