"""
This is to convert the 3D EDGAR data in the .txt file into meshgrid netcdf format
"""
import netCDF4 as nc4
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from numpy import loadtxt as loadtxt

def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the gridcell area in m2
	'''
	radius = 6371000;
	area = (np.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

lon = np.arange(0.05, 360.05, .1)
lat = np.arange(-89.95, 90.05, .1)	
lons,lats = np.meshgrid (lon,lat)
lat_res = lat[1] - lat[0];lon_res = lon[1] - lon[0];
EDGAR_gc_area = getArea(lons,lons+lon_res,lats,lats+lat_res)


def read_EDGAR_TXT(f_name,no_day):	
	lat,lon,data= loadtxt(f_name, delimiter=';', skiprows=3, usecols=(0,1,2), unpack=True)
	lon = lon%360  # from -180-180 to 0-360
	data= data *1000/no_day/24/60/60 # from tons/gridbox/month to kg/gridbox/s  
	map = np.zeros((1800,3600))
	for i in range(len(lon)):
		# print (lat[i]+90)*10-1,(lon[i]*10)-1
		map[int((lat[i]+90)*10-1),int(lon[i]*10-1)] = data[i]
	map = np.divide(map,EDGAR_gc_area)
	return map

# Species =['BC','SO2','NMVOC','OC']
Species =['NH3'] #,'NOx','NH3'
# Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH','PPA']
calendar_day = {'1':31,'2':28,'3':31,'4':30,'5':31,'6':30,'7':31,'8':31,'9':30,'10':31,'11':30,'12':31}
months=range(1,13)
filepath ="/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/1970/"

for species in Species:
	if (species == 'CO' or species == 'NOx'):
		Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH','PPA','CDS','CRS','LTO']
	elif species == 'NH3':
		Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PPA']
	else:
		Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','OTH','PPA']
	for sector in Sectors:
		for imonth in months:
			no_day = calendar_day[str(imonth)]
			f_in = filepath+species+'/'+'v431_v2_REFERENCE_'+species+'_1970_'+str(imonth)+'_'+sector+'.txt'
			map_2D = read_EDGAR_TXT(f_in,no_day)
			lon = np.arange(0.05, 360.05, .1)
			lat = np.arange(-89.95, 90.05, .1)
			f_out = filepath+species+'/'+'v431_v2_REFERENCE_'+species+'_1970_'+str(imonth)+'_'+sector+'.nc'
			f = nc4.Dataset(f_out,'w', format='NETCDF4') 
			f.createDimension('lat', len(lat))
			f.createDimension('lon', len(lon))
			latitudes = f.createVariable('lat',np.float32, ('lat'))
			longitudes = f.createVariable('lon',np.float32, ('lon'))
			latitudes[:] = lat
			longitudes[:] = lon		
			substances=f.createVariable('emi_'+species.lower(),np.float32,('lat','lon'))	
			substances[:]=map_2D
			substances.units='kg/m2/s'
			f.title= species +'emissions: 1970 from EDGAR_v431'
			f.close()

"""
X, Y = np.meshgrid(x,y)
map[map==0]=np.nan
plt.imshow(map,origin ='lower'); plt.show()
# plt.pcolor(X,Y, map, cmap='rainbow');plt.clim(10**(-5),10**(0));plt.show()
"""
