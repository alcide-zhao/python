import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import interp

def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

	
"""
############
#compare the time series (2006-2100) of no of molecules between CESM and EDGAR
############

######### CESM ORIGINAL files
file= '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/Original/RCP85_mam3_so2_elev_2000-2100_c20110913.nc'
nc_fid = nc4.Dataset(file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:] 
lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
lons,lats = np.meshgrid (lon,lat)  #meshgrid
area = getArea(lons,lons+lon_res,lats,lats+lat_res)
CESM = nc_fid.variables['emiss_ene'][:,2,:,:]+ nc_fid.variables['emiss_ind'][:,2,:,:]
# CESM = nc_fid.variables['emiss_awb'][:]+ nc_fid.variables['emiss_dom'][:]\
	# + nc_fid.variables['emiss_tra'][:]+ nc_fid.variables['emiss_shp'][:]+ nc_fid.variables['emiss_wst'][:]
nc_fid.close()

######### CESM replaced with EDGAR
file= '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/RCP85_mam3_so2_elev_2000-2100_c20110913.nc'
nc_fid = nc4.Dataset(file,mode='r')
EDGAR = nc_fid.variables['emiss_ene'][:,2,:,:]+ nc_fid.variables['emiss_ind'][:,2,:,:]
# EDGAR = nc_fid.variables['emiss_awb'][:]+ nc_fid.variables['emiss_dom'][:]\
	# + nc_fid.variables['emiss_tra'][:]+ nc_fid.variables['emiss_shp'][:]+ nc_fid.variables['emiss_wst'][:]

edgar_ts= np.nansum(np.nansum(np.multiply(EDGAR,area),axis=2),axis=1)
cesm_ts= np.nansum(np.nansum(np.multiply(CESM,area),axis=2),axis=1)

x=range(0,144)
plt.plot(x,edgar_ts,'r')
plt.plot(x,cesm_ts)
plt.show()

"""

############
#compare the monthly time series of no of molecules between RCP85 and EDGAR
# using the same methods
############
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

#########  RCP85
file_path ="/exports/csce/datastore/geos/users/s1667168/RCP/"
rcp_file = file_path+'RCP85_2010_SO2_NUM_MOL.0.1x0.1.nc'
nc_fid = nc4.Dataset(rcp_file,mode='r')
lon =  nc_fid.variables['lon'][:]; 
lat =  nc_fid.variables['lat'][:];
So2_RCP85 =  nc_fid.variables['SO2_ele_mol'][:]; 
nc_fid.close()
lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
lons,lats = np.meshgrid (lon,lat)  #meshgrid	
area = getArea(lons,lons+lon_res,lats,lats+lat_res)
rcp_gc_sum = np.multiply(So2_RCP85,area)   # area of the grid cell in m2 * emission/m2/year
rcp_ts=np.nansum(np.nansum(rcp_gc_sum,axis=2),axis=1)
print np.sum(rcp_ts)


## read the CESM grids to interpolate the RCP emissionsn firstly before caculations
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
lons,lats = np.meshgrid (lon_CESM,lat_CESM)
# CESM_fid.close()

#########  EDGAR
file_path = "/exports/csce/datastore/geos/users/s1667168/EDGAR/2010/"
SO2_EDGAR=np.zeros((12,96,144))
months=range(0,12) 
for imonth in months:
	file = file_path+'REFERENCE_SO2_2010_'+str(imonth+1)+'_NUM_MOL.0.1x0.1.nc'
	nc_fid = nc4.Dataset(file,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	SO2=nc_fid.variables['SO2_ele_mol'][:]
	SO2_EDGAR[imonth,:,:] =regrid_edgar2cesm(lon,lat,SO2,CESM_gc_area)
	nc_fid.close()
lons,lats = np.meshgrid (lon_CESM,lat_CESM);
lat_res = lat_CESM[1] - lat_CESM[0];lon_res = lon_CESM[1] - lon_CESM[0];
area = getArea(lons,lons+lon_res,lats,lats+lat_res)
edgar_gc_sum = np.multiply(SO2_EDGAR,area)   # area of the grid cell in m2 * emission/m2/year
	
#time series of 12 months
edgar_ts=  np.nansum(np.nansum(edgar_gc_sum,axis=2),axis=1)
print np.sum(edgar_ts)

x=range(0,12)
plt.plot(x,edgar_ts,'r')
plt.plot(x,rcp_ts)
plt.show()






