# -*- coding: utf-8 -*-
"""
Created on Dec 10 2016
Thisis to import cf-nc data and plot the global precipitation of a particular day
@author: Alcide.Zhao
"""
import site
import os
import time as clock
import numpy as np
import netCDF4 as nc4
from scipy import stats


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,
    'lib'
)

site.addsitedir(lib_path)

from lib.climate_extreme_indeciess_calculation import temperature_extreme_indeces as teic
 


########################################################
# nc data writing
########################################################
def netcdf4_write(file_name,year_series,lon,lat,txx, txn, tnx, tnn, dtr, fd0, su25, id0, tr20, tn10p, tn90p, tx10p, tx90p):

	"""
	Here is how you would normally create and store data in a netCDF file:
    1) Open/create a netCDF dataset.
    2) Define the dimensions of the data.
    3) Construct netCDF variables using the defined dimensions.
    4) Pass data into the netCDF variables.
    5) Add attributes to the variables and dataset (optional but recommended).
    6) Close the netCDF dataset.
	"""
	# 1) Open/create a netCDF dataset.
	f = nc4.Dataset(file_name,'w', format='NETCDF4') #'w' stands for write	
	"""
	The above line creates a netCDF file called "file_name" in the filepath folder.
	f is a netCDF Dataset object that provides methods for storing data to the file. 
	f also doubles as the root group. A netCDF group is basically a directory or 
	folder within the netCDF dataset. This allows you to organize data as you would 
	in a unix file system. Let's create a group for the heck of it
	"""
	
	# 2) Define the dimensions of the data.
	"""
	netCDF defines the sizes of all variables in terms of dimensions, so before any 
	variables can be created the dimen sions they use must be created first. A special 
	case, not often used in practice, is that of a scalar variable, which has no dimensions.
	A dimension is created using the Dataset.createDimension method of a Dataset. A Python 
	string is used to set the name of the dimension, and an integer value is used to set 
	the size. To create an unlimited dimension (a dimension that can be appended to), t
	he size value is set to None or 0. In this example, the time dimension is unlimited. 
	In netCDF4 files you can have more than one unlimited dimension, in netCDF 3 files 
	there may be only one, and it must be the first (leftmost) dimension of the variable
	"""
	f.createDimension('time', len(year_series))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))
	
	
	#3) Construct netCDF variables using the defined dimensions.
	"""
	netCDF variables behave much like python multidimensional array objects supplied by the 
	numpy module.  However, unlike numpy arrays, netCDF4 variables can be appended to along 
	one or more  ‘unlimited’ dimensions. To create a netCDF variable, use the Dataset.createVariable()
	method. The Dataset.createVariable() method has two mandatory arguments, the variable name 
	(a Python string), and the variable datatype. The variable’s dimensions are given by a tuple
	containing the dimension names (defined previously with Dataset.createDimension()). To create
	a scalar variable, simply leave out the dimensions keyword. The variable primitive datatypes
	correspond to the dtype attribute of a numpy array. You can specify the datatype as a numpy 
	dtype object, or anything that can be converted to a numpy dtype object. The unsigned integer 
	types and the 64-bit integer type can only be used if the file format is NETCDF4. The dimensions 
	themselves are usually also defined as variables, called coordinate variables. The Dataset.createVariable()
	method returns an instance of the Variable class whose methods can be used later to access 
	and set variable data and attributes.
	"""
	times = f.createVariable('time',np.float64, ('time'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	txxs = f.createVariable('txx',np.float32,('time','lat','lon'))
	txns = f.createVariable('txn',np.float32,('time','lat','lon'))
	tnxs = f.createVariable('tnx',np.float32,('time','lat','lon'))
	tnns = f.createVariable('tnn',np.float32,('time','lat','lon'))
	dtrs = f.createVariable('dtr',np.float32,('time','lat','lon'))
	fd0s = f.createVariable('fd0',np.float32,('time','lat','lon'))
	su25s = f.createVariable('su25',np.float32,('time','lat','lon'))
	id0s = f.createVariable('id0',np.float32,('time','lat','lon'))
	tr20s = f.createVariable('tr20',np.float32,('time','lat','lon'))
	tn10ps = f.createVariable('tn10p',np.float32,('time','lat','lon'))
	tn90ps = f.createVariable('tn90p',np.float32,('time','lat','lon'))
	tx10ps = f.createVariable('tx10p',np.float32,('time','lat','lon'))
	tx90ps = f.createVariable('tx90p',np.float32,('time','lat','lon'))
	# min_precips =f.createVariable('min_precip',np.float32,('time','lat','lon'))
	# print np.shape(rx1days)
	# print np.shape(rx1day)
	#4) Passing data into variables
	times[:] = year_series
	latitudes[:] = lat
	longitudes[:] = lon
	txxs[:] = txx
	txns[:] = txn
	tnxs[:] = tnx
	tnns[:] = tnn
	dtrs[:] = dtr
	fd0s[:] = fd0
	su25s[:] = su25
	id0s[:] = id0
	tr20s[:] = tr20
	tn10ps[:] =  tn10p
	tn90ps[:] =  tn90p
	tx10ps[:] = tx10p 
	tx90ps[:] = tx90p
	
	# 5) Add attributes to the variables and dataset (optional but recommended).
	"""
	There are two types of attributes in a netCDF file, global and variable. Global attributes provide information
	about the entire dataset, as a whole. Variable attributes provide information about one of the variables. Global
	attributes are set by assigning values to Dataset instance variables. Variable attributes are set by assigning
	values to Variable instances variables. Attributes can be strings, numbers or sequences. 
	"""
	
	# Global Attributes
	"""
	title       : What’s in the file
	institution : Where it was produced
	source      : How it was produced e.g. model version, instrument type
	history     : Audit trail of processing operations
	references  : Pointers to publications or web documentation
	comment     : Miscellaneous
	"""
	f.description = 'Temperature extrme indeces calculated using the CESM model outputs'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	
	# Variable Attributes
	"""
	units               : mandatory for all variables containing data other than dimensionless numbers. 
	standard_name       : identifies the quantity. Units must be consistent with standard name and any 
						  statistical processing.Standard name does not include coordinate or processing information.
	long_name           : not standardised.
	ancillary variables : a pointer to variables providing metadata about the individual data values e.g. standard error 
						  or data quality information. Numeric data variables may have 
						  ##############################################################################
						  # FillValue, missing_value, valid_max, valid_min, valid_range. missing_value #
						  ##############################################################################
	"""
	times.long_name = 'Year'
	times.units = '1'
	
	latitudes.long_name = 'latitude'
	latitudes.units = 'degree_north'
	longitudes.long_name = 'longitude'
	longitudes.units = 'degree_east'
	
	txxs.standard_name = 'Max Tmax'
	txxs.long_name = 'JJA maximum of daily temperature maximam'
	txxs.units = 'Degrees'
	
	txns.standard_name = 'Min Tmax'
	txns.long_name = 'JJA minimum of daily temperature maxmam'
	txns.units = 'Degrees'
	
	tnxs.standard_name = 'Max Tmin'
	tnxs.long_name = 'JJA maximum of daily temperature minimam'
	tnxs.units = 'Degrees'
	
	tnns.standard_name = 'Min Tmin'
	tnns.long_name = 'JJA minimum of daily temperature minimam'
	tnns.units = 'Degrees'
	
	dtrs.standard_name = 'Diurnal temperature range'
	dtrs.long_name = 'JJA mean between daily TX and TM'
	dtrs.units = 'Degrees'	

	fd0s.standard_name = 'Frost day '
	fd0s.long_name = 'JJA count when daily minimum temperature <0 degree'
	fd0s.units = 'days'	
	
	su25s.standard_name = 'Summer days'
	su25s.long_name = 'JJA count when daily minimum temperature >25 degree'
	su25s.units = 'days'	
	
	id0s.standard_name = 'Ice day'
	id0s.long_name = 'JJA count when daily maximum temperature <0 degree'
	id0s.units = 'days'	
	
	tr20s.standard_name = 'Tropical night'
	tr20s.long_name = 'JJA count when daily minimum temperature >25 degree'
	tr20s.units = 'days'	
	
	tn10ps.standard_name = 'Cool night'
	tn10ps.long_name = 'Percent of days when TN <10th percentile'
	tn10ps.units = 'Percntile'	

	tn90ps.standard_name = 'warm night'
	tn90ps.long_name = 'Percent of days when TN >90th percentile'
	tn90ps.units = 'Percntile'		
	
	tx10ps.standard_name = 'cool day'
	tx10ps.long_name = 'Percent of days when TX <10th percentile'
	tx10ps.units = 'Percntile'
	
	tx90ps.standard_name = 'Warm days'
	tx90ps.long_name = 'Percent of days when TX >90th percentile'
	tx90ps.units = 'Percntile'
	
	f.close()
	
	return
	
	
########################################################
#0. setting variables
########################################################
# linux path
input_path_mn = '/exports/csce/datastore/geos/users/s1667168/CESM/TEMP/medium_unemble_rcp8.5_fixA2005_2006_2010_trefhtmn'
input_path_mx = '/exports/csce/datastore/geos/users/s1667168/CESM/TEMP/medium_unemble_rcp8.5_fixA2005_2006_2010_trefhtmx'

########################################################
#1. find all .nc files that need to be extracted to be read
########################################################
os.system('find ' + os.path.abspath(input_path_mn) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path_mn + '/file_list.txt')
text_file_mn = open(input_path_mn + '/file_list.txt', "r")
text_content_mn = text_file_mn.readlines()

os.system('find ' + os.path.abspath(input_path_mx) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path_mx + '/file_list.txt')
text_file_mx = open(input_path_mx + '/file_list.txt', "r")
text_content_mx = text_file_mx.readlines()

# rntprint text_file 
    
#######################################################
#2. read in the .nc file and then calculate indexies, 
#######################################################

for ensumble_member in range(8,10): 
	#######################
	# 2.1 file loading
	#######################
	#    name  = os.path.basename(line)[-25:-17]    
	#    day   = name[6:8]
	#    month = name[4:6]
	#    year  = name[0:4]
	nc_f_mn = text_content_mn[ensumble_member][:-1]
	# print nc_f_mn
	nc_fid_mn = nc4.Dataset(nc_f_mn,mode='r')
	# nc_attrs, nc_dims, nc_vars = ncdump(nc_fid, False)
    #print nc_attrs
    #print nc_dims
    #print nc_vars
	lat = nc_fid_mn.variables['lat']
	lon = nc_fid_mn.variables['lon']
	time = nc_fid_mn.variables['date']
	temp_mn = nc_fid_mn.variables['TREFHTMN'][:,:,:]-273.15
    #sm = sm_all_time[0,:,:]#only one time slice in this dataset
    #sm = sm_all_time[0,::-1,:]#without ::-1, the image is upside down.
	nc_f_mx = text_content_mx[ensumble_member][:-1]
	nc_fid_mx = nc4.Dataset(nc_f_mx,mode='r')	
	temp_mx = nc_fid_mx.variables['TREFHTMX'][:,:,:]-273.15
	#######################
	# 2.2 calculation 
	#######################
	size_data = np.shape(temp_mn)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	# year_series = np.array([2006,2007])
	# print year_series
	#  varibale initialization 
	txx = np.empty((len(year_series),size_data[1],size_data[2]))
	txx[:] = np.nan
	txn = np.empty((len(year_series),size_data[1],size_data[2]))
	txn[:] = np.nan
	tnx = np.empty((len(year_series),size_data[1],size_data[2]))
	tnx[:] = np.nan
	tnn = np.empty((len(year_series),size_data[1],size_data[2]))
	tnn[:] = np.nan
	dtr = np.empty((len(year_series),size_data[1],size_data[2]))
	dtr[:] = np.nan
	fd0 = np.empty((len(year_series),size_data[1],size_data[2]))
	fd0[:] = np.nan
	su25 = np.empty((len(year_series),size_data[1],size_data[2]))
	su25[:] = np.nan
	id0 = np.empty((len(year_series),size_data[1],size_data[2]))
	id0[:] = np.nan
	tr20 = np.empty((len(year_series),size_data[1],size_data[2]))
	tr20[:] = np.nan
	tn10p = np.empty((len(year_series),size_data[1],size_data[2]))
	tn10p[:] = np.nan
	tn90p = np.empty((len(year_series),size_data[1],size_data[2]))
	tn90p[:] = np.nan
	tx10p = np.empty((len(year_series),size_data[1],size_data[2]))
	tx10p[:] = np.nan	
	tx90p = np.empty((len(year_series),size_data[1],size_data[2]))
	tx90p[:] = np.nan	
	
	
	
	# produce the reference data for calculating tn10p and so on
	temp_mn_refer = np.empty((5,92,size_data[1],size_data[2]))
	temp_mn_refer[:] = np.nan	
	temp_mx_refer = np.empty((5,92,size_data[1],size_data[2]))
	temp_mx_refer[:] = np.nan	
	tn10p_min_threshold = np.empty((size_data[1],size_data[2]))
	tn10p_min_threshold[:] = np.nan	
	tn90p_min_threshold = np.empty((size_data[1],size_data[2]))
	tn90p_min_threshold[:] = np.nan	
	tx10p_max_threshold = np.empty((size_data[1],size_data[2]))
	tx10p_max_threshold[:] = np.nan	
	tx90p_max_threshold = np.empty((size_data[1],size_data[2]))
	tx90p_max_threshold[:] = np.nan	
	
	layer_e = -123
	for layer in range(5):
		layer_b = layer_e + 274
		layer_e = layer_b + 91
		temp_mn_refer[layer,:,:,:] = temp_mn[layer_b:layer_e+1,:,:]
		temp_mx_refer[layer,:,:,:] = temp_mx[layer_b:layer_e+1,:,:]
	temp_mn_refer_3d =stats.nanmean(temp_mn_refer,axis=0)
	temp_mx_refer_3d =stats.nanmean(temp_mx_refer,axis=0)
	tn10p_min_threshold[:] = np.percentile(temp_mn_refer_3d,10,axis=0) 
	tn90p_min_threshold[:] = np.percentile(temp_mn_refer_3d,90,axis=0) 
	tx10p_max_threshold[:] = np.percentile(temp_mx_refer_3d,10,axis=0) 
	tx90p_max_threshold[:] = np.percentile(temp_mx_refer_3d,90,axis=0) 
	layer_output = 0 # the time dimenssion of the output variable
	layer_e = -123   # In order to let the first year behins from the 151 day 274-213=151
	
	for iyear in year_series:	
		layer_b = layer_e + 274
		layer_e = layer_b + 91 # this should depend on if the year is leapyear or not,but here 364 because all CESM year are with 364 days
		# layer_index = [layer for layer in range(len(year)) if year[layer] == iyear]
		# layer_b = np.min(layer_index)
		# layer_e = np.max(layer_index)
		temp_jja_mn = temp_mn[layer_b:layer_e+1,:,:]
		temp_jja_mx = temp_mx[layer_b:layer_e+1,:,:]
		# r10, r20, rnm, sdii, precptot, rx5day, rx1day, r95p, r99p, cdd, cwd
		
		txx[layer_output,:,:], txn[layer_output,:,:], tnx[layer_output,:,:], tnn[layer_output,:,:], dtr[layer_output,:,:], fd0[layer_output,:,:], su25[layer_output,:,:],\
		id0[layer_output,:,:], tr20[layer_output,:,:], tn10p[layer_output,:,:], tn90p[layer_output,:,:], tx10p[layer_output,:,:], tx90p[layer_output,:,:] = teic(temp_jja_mn,temp_jja_mx,tn10p_min_threshold,\
		tn90p_min_threshold,tx10p_max_threshold,tx90p_max_threshold)
		
		print str(iyear)+" passed to layer " +str(layer_output) 
		layer_output = layer_output+1
		
	#######################
	# 2.3 write results into nc file
	#######################
	file_name = input_path_mn + '/'+os.path.basename(input_path_mn)+str(['_0'+str(ensumble_member+1) if ensumble_member<9 else '_'+str(ensumble_member+1)])[2:5] + '_EI_JJA.nc'
	# time = year_series.astype('int')
	# print "r99p"
	# print r99p
	
	netcdf4_write(file_name, year_series, lon, lat, txx, txn, tnx, tnn, dtr, fd0, su25, id0, tr20, tn10p, tn90p, tx10p, tx90p)
	print "Data stored for unsemble "+ str(ensumble_member+1)
	nc_fid_mx.close()
	nc_fid_mn.close()
	
	
