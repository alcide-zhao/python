import glob   
path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/PREP/6080//*.nc'   
files=glob.glob(path)   
for file in files: 
	print file