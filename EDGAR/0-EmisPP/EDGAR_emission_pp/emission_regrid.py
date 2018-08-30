import numpy as np 


def getArea(lon1,lon2,lat1,lat2):
	'''
	Get the gridcell areas of the lon-lat meshgrids (rank == 2)
	input: lon1,lon2,lat1,lat2 must be in the format of meshgrid 
	       lon2-lon1 denotes the zonal grid resolution at each point
		   lon2-lon1 denotes the meridional grid resolution at each point
	Output:
	       gc_area: a matrix the same size as input matix, stores the gridcell areas
	!!! Trust me this algorithm is better than that used by the IRIS package
	'''
	radius = 6371000;  # radius of the earth (m)
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return gc_area


def regrid_emission(lon_o,lat_o,emis_o,lon_n,lat_n):
    """
    input:
	    lon_o,lat_o: original emission file lat/lon (rank == 1)
		emis_o: original emission file (2D) in kg/m2/time
		lon_n,lat_n: target emission file grids (rank == 1)
    output:
		emis : regridded emission file (2D) in kg/m2/time
	Note The new grid resolution must be of integer number of the old grid resolution otherwise you will lose mass 
	at boundry regions (either arctic/antarctic or 0/180 longitude band)
	To get rid of this, some preprocessing must be aplied: i.e interpolate the original emission file to teh closest
	meshgrid where the new grid resolution is integer number of the interpolated grid resolution
	"""
	emis = np.zeros((len(lat_n),len(lon_n)));
	
	# Derive the original grids $area_o
	lon_m,lat_m=np.meshgrid(lon_o,lat_o)   # meshgrid both lon and lat into 2D fields
	lat_res = lat_o[1] - lat_o[0];lon_res = lon_o[1] - lon_o[0];  # This applies only for regular regrid, otherwise change it to array
	area_o = getArea(lon_m,lon_m+lon_res,lat_m,lat_m+lat_res)
	
	# Derive the target grids $area_n
	lon_m,lat_m=np.meshgrid(lon_n,lat_n)   # meshgrid both lon and lat into 2D fields
	lat_res = lat_n[1] - lat_n[0];lon_res = lon_n[1] - lon_n[0];  # This applies only for regular regrid, otherwise change it to array
	area_n = getArea(lon_m,lon_m+lon_res,lat_m,lat_m+lat_res)
	
	# Weighted emission at each gridcell $emis_o_w
	emis_o_w =  np.multiply(emis_o,area_n)
	
	# regridding
	No_box_lat=int(1.0*(lat_n[1] - lat_n[0])/(lat_o[1] - lat_o[0]))  # how many rows is putting into the new gridbox
	No_box_lon=int(1.0*(lon_n[1] - lon_n[0])/(lon_o[1] - lon_o[0]))  # how many columns is putting into the new gridbox
	for row in range(0,len(lat_n)):
		for column in range(0,len(lon_n)):	
			emis[row,column] =np.nansum(np.nansum(emis_o_w[row*No_box_lat:(row+1)*No_box_lat,column*No_box_lon:(column+1)*No_box_lon],axis=1),axis=0)
	emis = np.divide(emis,CESM_gc_area)
	return emis
	

	
##############
#    Main    #
##############
lon_o,lat_o,emis_o =...
lon_n,lat_n =...
emis = regrid_emission(lon_o,lat_o,emis_o,lon_n,lat_n)