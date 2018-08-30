
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import scipy.io as sio
ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan#;ocean_mask_CESM[ocean_mask_CESM==0]=1

def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# area=np.multiply(area,ocean_mask_CESM)
	# plt.imshow(area);plt.show()
	area = area/np.nansum(np.nansum(area,axis=1),axis=0)
	return area

	
data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp'
osci_data =data_path+'/TS_oscilation_rhis_1960_2005.nc'
nc_fid = nc4.Dataset(osci_data)

lon=nc_fid.variables['lon'][:];lat=nc_fid.variables['lat'][:];
lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
lons,lats = np.meshgrid (lon,lat)
area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)

# np.nansum(np.nansum(annual_data,axis=2),axis=1))
his=np.nansum(np.nansum(np.multiply(nc_fid.variables['his'][0,:,:,:],area),axis=2),axis=1)-273.15
osci_data = data_path+'/TS_oscilation_rcp85_fixa_2006_2100.nc'
nc_fid = nc4.Dataset(osci_data)
rcp=np.nansum(np.nansum(np.multiply(nc_fid.variables['rcp85'][0,:,:,:],area),axis=2),axis=1)-273.15
fix=np.nansum(np.nansum(np.multiply(nc_fid.variables['fixa'][0,:,:,:],area),axis=2),axis=1)-273.15

print np.shape(his)
year = range(1961,2006);
fig = plt.figure(facecolor='White',figsize=[8,8])
ax1 = plt.subplot(1,1,1);
ax1.plot(year,his,'-',color="k",linewidth=2)
year = range(2006,2101);
ax1.plot(year,rcp,'-',color="r",linewidth=2);#ax1.set_ylim([16.65,18.5]);#ax1.set_ylim([16.65,18.8]);
ax1.plot(year,fix,'-',color="b",linewidth=2);#ax1.set_ylim([16.65,18.5]);#ax1.set_ylim([16.65,18.8]);

plt.show()

