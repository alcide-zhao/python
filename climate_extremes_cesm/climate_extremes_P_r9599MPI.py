# -*- coding: utf-8 -*-
"""
Created on Mon WED 16  2016
Thisis to import cf-nc data and plot the global precipitation of a particular day
@author: Alcide.Zhao
"""
import site
import os
import time as clock
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    'lib'
)

site.addsitedir(lib_path)
from lib.climate_extreme_indeciess_calculation import precip_extreme_indeces as peic
 
#leapyear judgement
def leapyear(iyear):
	leapyear = False
	if (iyear % 400 == 0):
		leapyear = True
	else:
		if (iyear % 100 == 0):
			leapyear = False
		else:
			if (iyear % 4 == 0):
				leapyear = True
			else:
				leapyear = False



########################################################
# nc data writing
########################################################
def netcdf4_write(file_name,year_series,lon,lat,rx1day, rx5day, sdii, r10, r20, rnm, cdd, cwd, r95p, r99p, precptot, total_precip, mean_prep, std_prep,MPI):
	
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
	rx1days = f.createVariable('rx1day',np.float32,('time','lat','lon'))
	rx5days = f.createVariable('rx5day',np.float32,('time','lat','lon'))
	sdiis = f.createVariable('sdii',np.float32,('time','lat','lon'))
	r10s = f.createVariable('r10',np.float32,('time','lat','lon'))
	r20s = f.createVariable('r20',np.float32,('time','lat','lon'))
	rnms = f.createVariable('rnm',np.float32,('time','lat','lon'))
	cdds = f.createVariable('cdd',np.float32,('time','lat','lon'))
	cwds = f.createVariable('cwd',np.float32,('time','lat','lon'))
	r95ps = f.createVariable('r95p',np.float32,('time','lat','lon'))
	r99ps = f.createVariable('r99p',np.float32,('time','lat','lon'))
	prcptots = f.createVariable('precptot',np.float32,('time','lat','lon'))
	total_precips = f.createVariable('total_precip',np.float32,('time','lat','lon'))
	mean_preps = f.createVariable('mean_precip',np.float32,('time','lat','lon'))
	std_preps = f.createVariable('std_precip',np.float32,('time','lat','lon'))
	MPIs = f.createVariable('MPI',np.float32,('time','lat','lon'))
	
	#4) Passing data into variables
	times[:] = year_series
	latitudes[:] = lat
	longitudes[:] = lon
	rx1days[:,:,:] = rx1day[:,:,:]
	rx5days[:] = rx5day
	sdiis[:] = sdii
	r10s[:] = r10
	r20s[:] = r20
	rnms[:] = rnm
	cdds[:] = cdd
	cwds[:] = cwd
	r95ps[:] = r95p
	r99ps[:] =  r99p
	prcptots[:] =  precptot
	total_precips[:] = total_precip 
	mean_preps[:] = mean_prep
	std_preps[:] = std_prep
	MPIs[:] = MPI
	
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
	f.description = 'Precipitation extrme indeces calculated using the CESM model outputs'
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
	
	rx1days.standard_name = 'maximum 1-day precip amount'
	rx1days.long_name = 'maximum 1-day precipitation amount'
	rx1days.units = 'mm'
	
	rx5days.standard_name = 'maximum 5-day precip amount'
	rx5days.long_name = 'maximum 5-day precipitation amount'
	rx5days.units = 'mm'
	
	sdiis.standard_name = 'simple daily intensity index'
	sdiis.long_name = 'annual precip amount devided by the number of wet days (> 1mm)'
	sdiis.units = 'm/s/day'
	
	r10s.standard_name = 'number of heavy precip days'
	r10s.long_name = 'nnual count of days woth precipitation >= 10mm'
	r10s.units = 'day'
	
	r20s.standard_name = 'number of very heavy precip days'
	r20s.long_name = 'nnual count of days woth precipitation >= 20mm'
	r20s.units = 'day'	

	rnms.standard_name = 'number of days above nn mm here defines as 30 mm'
	rnms.long_name = 'annual count of days woth precipitation >= 30mm'
	rnms.units = 'day'	
	
	cdds.standard_name = 'consective dry days'
	cdds.long_name = 'maximum consecive days with precipitation<1mm'
	cdds.units = 'day'	
	
	cwds.standard_name = 'consective wet days'
	cwds.long_name = 'maximum consecive days with precipitation>=1mm'
	cwds.units = 'day'	
	
	r95ps.standard_name = 'very wet days'
	r95ps.long_name = 'annual total precipitation >95th percetile'
	r95ps.units = 'day'	
	
	r99ps.standard_name = 'extremely wet days'
	r99ps.long_name = 'annual total precipitation >99th percetile'
	r99ps.units = 'day'	

	prcptots.standard_name = 'annual total wet day prrecip'
	prcptots.long_name = 'annual total precipitation in wet days (rr>1mm)'
	prcptots.units = 'mm'		
	
	total_precips.standard_name = 'annual total prrecip'
	total_precips.long_name = 'annual total precipitation'
	total_precips.units = 'mm'
	
	MPIs.standard_name = 'Monsoon Precipitation Index'
	MPIs.Long_name = 'Ratio of local summer(MJJAS)-winter(NDJFM) to the annual total for the North Hemisphere'
	f.close()
	
	return
	
	
########################################################
#0. setting variables
########################################################
# linuc path
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/PRECT/rcp85/'
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
text_file = open(input_path + '/file_list.txt', "r")
text_content = text_file.readlines()
# rntprint text_file 
#######################################################
#2. read in the .nc file and then calculate indexies, 
#######################################################
threshold_file = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/precip_9599_thre.nc'
nc_fid = nc4.Dataset(threshold_file,mode='r')
r95p_threshold = nc_fid.variables['r95pt']
r99p_threshold = nc_fid.variables['r99pt']


for ensumble_member in range(0,30): 
	#######################
	# 2.1 file loading
	#######################
	#    name  = os.path.basename(line)[-25:-17]    
	#    day   = name[6:8]
	#    month = name[4:6]
	#    year  = name[0:4]
	nc_f = text_content[ensumble_member][:-1]
	nc_fid = nc4.Dataset(nc_f,mode='r')
	# nc_attrs, nc_dims, nc_vars = ncdump(nc_fid, False)
	lat = nc_fid.variables['lat']
	lon = nc_fid.variables['lon']
	time = nc_fid.variables['date']
	precip = nc_fid.variables['PRECT'][:,:,:]
    #sm = sm_all_time[0,:,:]#only one time slice in this dataset
    #sm = sm_all_time[0,::-1,:]#without ::-1, the image is upside down.
	precip[precip <0]=np.nan
	
	#######################
	# 2.2 calculation 
	#######################
	size_data = np.shape(precip)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	
	#  varibale initialization 
	r10 = np.empty((len(year_series),size_data[1],size_data[2]));	r10[:] = np.nan
	r20 = np.empty((len(year_series),size_data[1],size_data[2]));	r20[:] = np.nan
	rnm = np.empty((len(year_series),size_data[1],size_data[2]));	rnm[:] = np.nan
	sdii = np.empty((len(year_series),size_data[1],size_data[2]));	sdii[:] = np.nan
	precptot = np.empty((len(year_series),size_data[1],size_data[2]));	precptot[:] = np.nan
	rx5day = np.empty((len(year_series),size_data[1],size_data[2]));	rx5day[:] = np.nan
	rx1day = np.empty((len(year_series),size_data[1],size_data[2]));	rx1day[:] = np.nan
	cdd = np.empty((len(year_series),size_data[1],size_data[2]));	cdd[:] = np.nan
	cwd = np.empty((len(year_series),size_data[1],size_data[2]));	cwd[:] = np.nan
	total_precip = np.empty((len(year_series),size_data[1],size_data[2])); total_precip[:] = np.nan	
	r95p = np.empty((len(year_series),size_data[1],size_data[2])); r95p[:] = np.nan
	r99p = np.empty((len(year_series),size_data[1],size_data[2])); r99p[:] = np.nan
	mean_prep = np.empty((len(year_series),size_data[1],size_data[2])); mean_prep[:] = np.nan
	std_prep = np.empty((len(year_series),size_data[1],size_data[2]));std_prep[:] = np.nan	
	MPI = np.empty((len(year_series),size_data[1],size_data[2])); MPI[:] = np.nan
	
	layer_output = 0 # the time dimenssion of the output variable
	layer_JJA_e = -123   # In order to let the first day behins from the 152 day 274-123=151
	layer_annual_e = -1 # In order to let the first day behins from the 0 day  1-1=0
	layer_summer_e = -93 # In order to let the first day behins from the 0 day 213-93=120
	layer_winter_e = -62 # In order to let the first day behins from the 152-62=90 day 
	for iyear in year_series:	
		""" 
		CESM DATA
		"""
		# layer_range = range((iyear-year_series[0])*365,(iyear-year_series[0]+1)*365)
		# summer_pa = np.sum([precip[layer,:,:] for layer in layer_range if (time[layer] >= iyear*10000+501 and time[layer] <= iyear*10000+930)],axis=0)
		# winter_pa = np.sum([precip[layer,:,:] for layer in layer_range if ((time[layer] >= iyear*10000+1101 and time[layer] <= iyear*10000+1231) or (time[layer] >= iyear*10000+101 and time[layer] <= iyear*10000+331) ) ],axis=0)			

		# The Monsoon Precipitation Index
		layer_annual_b = layer_annual_e+1
		layer_annual_e = layer_annual_b+364
		annual_pa = np.nansum(precip[layer_annual_b:layer_annual_e+1,:,:],axis=0)
		
		layer_summer_b = layer_summer_e+213
		layer_summer_e = layer_summer_b+152
		summer_pa = np.nansum(precip[layer_summer_b:layer_summer_e+1,:,:],axis=0)		
		 
		layer_winter_b = layer_winter_e+152
		layer_winter_e = layer_winter_b+213
		winter_pa = annual_pa - np.nansum(precip[layer_winter_b:layer_winter_e+1,:,:],axis=0)		
		
		MPI[layer_output,:,:] = np.divide((summer_pa-winter_pa),annual_pa)
		
		# layer_b = [layer for layer in layer_range if time[layer] == iyear*10000+601][0]
		# layer_e = [layer for layer in layer_range if time[layer] == iyear*10000+831][0]
		layer_JJA_b = layer_JJA_e+274
		layer_JJA_e = layer_JJA_b+91
		precp_year_data = precip[layer_JJA_b:layer_JJA_e+1,:,:]*24*60*60*1000  # from m/s to mmper day 
		"""
		APHRO Observations
		"""	
		# if leapyear(iyear):
			# layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*1000+153][0]
			# layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*1000+243][0]
		# else:
			# layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*1000+152][0]
			# layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*1000+242][0]
		# precp_year_data = precip[layer_b:layer_e+1,:,:] 		
		# r10, r20, rnm, sdii, precptot, rx5day, rx1day, r95p, r99p, cdd, cwd
		_, _, _, _, _, _, _,\
		r95p[layer_output,:,:], r99p[layer_output,:,:], _, _, _,mean_prep[layer_output,:,:], std_prep[layer_output,:,:]= peic(precp_year_data,r95p_threshold, r99p_threshold,rnt = 30)
		print str(iyear)+" passed to layer " +str(layer_output) 
		layer_output = layer_output+1
	#######################
	# 2.3 write results into nc file
	#######################
	# note here ensumble_member begins from 1 while the loop begins from 0
	# this is the ensumble_member+1 is performed for the output nc file name
	file_name = input_path + 'CESM_PrepEI_2006_2100'+str(['_0'+str(ensumble_member+1) if ensumble_member<9 else '_'+str(ensumble_member+1)])[2:5] + '_JJA.nc'
	# time = year_series.astype('int')
	# print "r99p"
	# print r99p
	netcdf4_write(file_name, year_series, lon, lat, rx1day, rx5day, sdii, r10, r20, rnm, cdd, cwd, r95p, r99p, precptot, total_precip,mean_prep,std_prep,MPI)
	print "Data stored for unsemble "+ str(ensumble_member+1)
	nc_fid.close()