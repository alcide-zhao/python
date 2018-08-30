import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  

#############################################
# 0.0 data input
#############################################
variable='TS'
output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name_out = output_path+'ensumble_mean_'+variable+'_192001-200512.nc'
# b.e11.B20TRC5CNBDRD.f09_g16.001.cam.h1.PRECT.18500101-20051231

#############################################
# 1 ensumble mean 3D data
#############################################
def ensumble_mean_of_CESM_TP(files):
	att_ensumble_mean = np.empty((31390,192,288))
	for time in range(31390):
		att_value = np.empty((len(files),192,288))
		ensumble_no = 0
		for file in files:
			# print file
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			lev = nc_fid.variables['lev'][:]
			lev_sn= nc_fid.variables['lev'].standard_name
			time_series = nc_fid.variables['date'][:]
			att_value[ensumble_no,:,:] = nc_fid.variables[variable][time,:,:]
			units = nc_fid.variables[variable].units
			long_name = nc_fid.variables[variable].long_name
			# missing_value = nc_fid.variables[variable].missing_values
			missing_value = np.nan
			ensumble_no=ensumble_no+1
			nc_fid.close()
		att_value[att_value == missing_value] =np.nan
		att_ensumble_mean[time,:,:] = stats.nanmean(att_value,axis = 0)
	return lon,lat,lev,lev_sn,time_series,att_ensumble_mean,units,missing_value,long_name
	
	
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'+variable+'/rcp85/*.nc'
files=sorted(glob.glob(input_path))
lon,lat,time_series,att_ensumble_mean_rcp85,rcp85_units,rcp85_mv,rcp85_ln = ensumble_mean_of_CESM_TP(files)
# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'+variable+'/fixa/*.nc'
# files=sorted(glob.glob(input_path))
# lon,lat,time_series,att_ensumble_mean_rcp85_fixA,fixA_units,fixA_mv,fixA_ln = ensumble_mean_of_CESM_TP(files)


######################################
#1.1 writting the results into nc files
######################################

f = nc4.Dataset(file_name_out,'w', format='NETCDF4') #'w' stands for write
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))
f.createDimension('time', len(time_series))

times = f.createVariable('time',np.float64, ('time'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))

PRECTs = f.createVariable('PRECT',np.float32,('time','lat','lon'))
# rcp85_fixAs = f.createVariable('rcp85_fixA',np.float32,('time','lat','lon'))
	
times[:] = time_series
latitudes[:] = lat
longitudes[:] = lon
PRECTs[:]= att_ensumble_mean_rcp85
# rcp85_fixAs[:]= att_ensumble_mean_rcp85_fixA

times.long_name = 'Month'
times.units = 'Month'
latitudes.long_name = 'latitude'
latitudes.units = 'degree_north'
longitudes.long_name = 'longitude'
longitudes.units = 'degree_east'

PRECTs.units=rcp85_units
PRECTs.missing_value=rcp85_mv
PRECTs.long_name =rcp85_ln
# rcp85_fixAs.units=fixA_units
# rcp85_fixAs.missing_value=fixA_mv
# rcp85_fixAs.long_name =fixA_ln

f.description = 'Time series of nsumble mean of '+variable+' calculated using the CESM model outputs'
f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
f.institution = 'Alcide Zhao at the university of Edinburgh'
f.close()