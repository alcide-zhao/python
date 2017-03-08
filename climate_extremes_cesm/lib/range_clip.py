import numpy as np


def range_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,data):
	"""
	clip the data based on given range
	"""
	lon = np.array(lon)
	lat = np.array(lat)
	colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
	colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
	row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
	row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
	lon_clipped = lon[colum_s:colum_e+1]
	lat_clipped = lat[row_s:row_e+1]
	if (np.rank(data) == 2):
		data_clipped = data[row_s:row_e+1,colum_s:colum_e+1]
	elif (np.rank(data) == 3):
		data_clipped = data[:,row_s:row_e+1,colum_s:colum_e+1]
	elif (np.rank(data) == 4):
		data_clipped = data[:,:,row_s:row_e+1,colum_s:colum_e+1]
	return lon_clipped, lat_clipped, data_clipped