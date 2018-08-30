'''
This is to process the caculated HWI and CSI indices into ensemble
mean, median, max, min, 75th and 25th

	rcp8.5   rcp8.5_FixA  His (1961-1990) rcp4.5 (-2080)
'''

import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import os
import glob 


#############################################
# 2 ensumble mean 4D data
#############################################
def ensemble_pp(text_content,variable,YearNo):
	en_mean = np.empty((YearNo,192,288,8));
	en_median = np.empty((YearNo,192,288,8));
	en_max = np.empty((YearNo,192,288,8));
	en_min = np.empty((YearNo,192,288,8));
	en_75 = np.empty((YearNo,192,288,8));
	en_25 = np.empty((YearNo,192,288,8));
	for iyear in range(YearNo):
		att_value = np.empty((ensemble,192,288,8))
		for en_no in range(0,ensemble):
			# print en_no
			nc_f = text_content[en_no][:-1]
			nc_fid = nc4.Dataset(nc_f,mode='r')
			att_value[en_no,:,:,:] = nc_fid.variables[variable][iyear,:,:,:]
			nc_fid.close()
		en_mean[iyear,:,:,:] = stats.nanmean(att_value,axis = 0)
		en_median[iyear,:,:,:] = stats.nanmedian(att_value,axis = 0)
		en_max[iyear,:,:,:] = np.nanpercentile(att_value,95,axis = 0)
		en_min[iyear,:,:,:] = np.nanpercentile(att_value,5,axis = 0)
		en_75[iyear,:,:,:] = np.nanpercentile(att_value,75,axis = 0)
		en_25[iyear,:,:,:] = np.nanpercentile(att_value,25,axis = 0)
	return en_mean,en_median,en_max,en_min,en_75,en_25


# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/abs_v2/'+scenario+'/*.nc'
# files=sorted(glob.glob(input_path))


scenario_dic ={'rcp85':[21,0,95]}  #,'fixa':[14,0,95],'rcp45':[15,0,75],'his':[22,0,86],'rcp85':[20,0,95]
directory = 'CalDayThr_season'  #CalDayThr_InterOsi
for scenario in scenario_dic:
	print scenario
	ensemble = scenario_dic[scenario][0]; YearNo=scenario_dic[scenario][2];
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'+directory+'/'+scenario+'/'  #_InterOsi
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + 'file_list.txt', "r")
	text_content = text_file.readlines()

	variable='HWI_MJJASO';MJJASO_mean,MJJASO_median,MJJASO_max,MJJASO_min,MJJASO_75,MJJASO_25 = ensemble_pp(text_content,variable,YearNo)
	variable='HWI_NDJFMA';NDJFMA_mean,NDJFMA_median,NDJFMA_max,NDJFMA_min,NDJFMA_75,NDJFMA_25 = ensemble_pp(text_content,variable,YearNo)
	threshold_data ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'+directory+'/'+scenario+'/HWI_CSI_01_'+scenario+'.nc'
	print threshold_data
	nc_fid = nc4.Dataset(threshold_data,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	year= nc_fid.variables['time'][:]
	
	output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'+directory+'/'
	file_name_out = output_path+'HWI_CSI_EnPP_'+scenario+'.nc'
	f = nc4.Dataset(file_name_out,'w', format='NETCDF4')
	f.createDimension('stats', 6)
	f.createDimension('year', len(year))
	f.createDimension('lat', len(lat))
	f.createDimension('lon', len(lon))

	statss = f.createVariable('stats',np.float32, ('stats'))
	years = f.createVariable('year',np.float32, ('year'))
	latitudes = f.createVariable('lat',np.float32, ('lat'))
	longitudes = f.createVariable('lon',np.float32, ('lon'))
	statss[:]=range(6)
	years[:] = year
	latitudes[:] = lat
	longitudes[:] = lon
	HW_MJJASO_Ms= f.createVariable('HW_MJJASO_M',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_Xs= f.createVariable('HW_MJJASO_X',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_NOs=f.createVariable('HW_MJJASO_NO',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_DCs= f.createVariable('HW_MJJASO_DC',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_DuXs= f.createVariable('HW_MJJASO_DuX',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_DuMs= f.createVariable('HW_MJJASO_DuM',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_InXs= f.createVariable('HW_MJJASO_InX',np.float32,('stats','year','lat','lon'))
	HW_MJJASO_InMs= f.createVariable('HW_MJJASO_InM',np.float32,('stats','year','lat','lon'))

	HW_MJJASO_Ms[:] = np.array([MJJASO_mean[:,:,:,0],MJJASO_median[:,:,:,0],MJJASO_max[:,:,:,0],MJJASO_min[:,:,:,0],MJJASO_75[:,:,:,0],MJJASO_25[:,:,:,0]])
	HW_MJJASO_Xs[:] = np.array([MJJASO_mean[:,:,:,1],MJJASO_median[:,:,:,1],MJJASO_max[:,:,:,1],MJJASO_min[:,:,:,1],MJJASO_75[:,:,:,1],MJJASO_25[:,:,:,1]])
	HW_MJJASO_NOs[:] = np.array([MJJASO_mean[:,:,:,2],MJJASO_median[:,:,:,2],MJJASO_max[:,:,:,2],MJJASO_min[:,:,:,2],MJJASO_75[:,:,:,2],MJJASO_25[:,:,:,2]])
	HW_MJJASO_DCs[:] = np.array([MJJASO_mean[:,:,:,3],MJJASO_median[:,:,:,3],MJJASO_max[:,:,:,3],MJJASO_min[:,:,:,3],MJJASO_75[:,:,:,3],MJJASO_25[:,:,:,3]])
	HW_MJJASO_DuXs[:] = np.array([MJJASO_mean[:,:,:,4],MJJASO_median[:,:,:,4],MJJASO_max[:,:,:,4],MJJASO_min[:,:,:,4],MJJASO_75[:,:,:,4],MJJASO_25[:,:,:,4]])
	HW_MJJASO_DuMs[:] = np.array([MJJASO_mean[:,:,:,5],MJJASO_median[:,:,:,5],MJJASO_max[:,:,:,5],MJJASO_min[:,:,:,5],MJJASO_75[:,:,:,5],MJJASO_25[:,:,:,5]])
	HW_MJJASO_InXs[:] = np.array([MJJASO_mean[:,:,:,6],MJJASO_median[:,:,:,6],MJJASO_max[:,:,:,6],MJJASO_min[:,:,:,6],MJJASO_75[:,:,:,6],MJJASO_25[:,:,:,6]])
	HW_MJJASO_InMs[:] = np.array([MJJASO_mean[:,:,:,7],MJJASO_median[:,:,:,7],MJJASO_max[:,:,:,7],MJJASO_min[:,:,:,7],MJJASO_75[:,:,:,7],MJJASO_25[:,:,:,7]])
	
	statss.long_name='0-ensemble mean,1-ensemble median,2-ensemble max,3-ensemble min,4-ensemble 75th percentile,5-ensemble 25th percentile'
	HW_MJJASO_Ms.long_name='Heatwvae magnitude: MJJASO mean'
	HW_MJJASO_Xs.long_name='Heatwvae magnitude: MJJASO max'
	HW_MJJASO_NOs.long_name='Heatwvae : MJJASO no of evnents'
	HW_MJJASO_DCs.long_name='total dayso of HW in MJJASO'
	HW_MJJASO_DuXs.long_name='maximum Duration of MJJASO heatwaves '
	HW_MJJASO_DuMs.long_name='mean Duration of MJJASO heatwaves '
	HW_MJJASO_InXs.long_name='maximum intensity of MJJASO heatwaves'
	HW_MJJASO_InMs.long_name='mean intensity of MJJASO heatwaves'
	
	HW_NDJFMA_Ms= f.createVariable('HW_NDJFMA_M',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_Xs= f.createVariable('HW_NDJFMA_X',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_NOs=f.createVariable('HW_NDJFMA_NO',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_DCs= f.createVariable('HW_NDJFMA_DC',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_DuXs= f.createVariable('HW_NDJFMA_DuX',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_DuMs= f.createVariable('HW_NDJFMA_DuM',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_InXs= f.createVariable('HW_NDJFMA_InX',np.float32,('stats','year','lat','lon'))
	HW_NDJFMA_InMs= f.createVariable('HW_NDJFMA_InM',np.float32,('stats','year','lat','lon'))

	HW_NDJFMA_Ms[:] = np.array([NDJFMA_mean[:,:,:,0],NDJFMA_median[:,:,:,0],NDJFMA_max[:,:,:,0],NDJFMA_min[:,:,:,0],NDJFMA_75[:,:,:,0],NDJFMA_25[:,:,:,0]])
	HW_NDJFMA_Xs[:] = np.array([NDJFMA_mean[:,:,:,1],NDJFMA_median[:,:,:,1],NDJFMA_max[:,:,:,1],NDJFMA_min[:,:,:,1],NDJFMA_75[:,:,:,1],NDJFMA_25[:,:,:,1]])
	HW_NDJFMA_NOs[:] = np.array([NDJFMA_mean[:,:,:,2],NDJFMA_median[:,:,:,2],NDJFMA_max[:,:,:,2],NDJFMA_min[:,:,:,2],NDJFMA_75[:,:,:,2],NDJFMA_25[:,:,:,2]])
	HW_NDJFMA_DCs[:] = np.array([NDJFMA_mean[:,:,:,3],NDJFMA_median[:,:,:,3],NDJFMA_max[:,:,:,3],NDJFMA_min[:,:,:,3],NDJFMA_75[:,:,:,3],NDJFMA_25[:,:,:,3]])
	HW_NDJFMA_DuXs[:] = np.array([NDJFMA_mean[:,:,:,4],NDJFMA_median[:,:,:,4],NDJFMA_max[:,:,:,4],NDJFMA_min[:,:,:,4],NDJFMA_75[:,:,:,4],NDJFMA_25[:,:,:,4]])
	HW_NDJFMA_DuMs[:] = np.array([NDJFMA_mean[:,:,:,5],NDJFMA_median[:,:,:,5],NDJFMA_max[:,:,:,5],NDJFMA_min[:,:,:,5],NDJFMA_75[:,:,:,5],NDJFMA_25[:,:,:,5]])
	HW_NDJFMA_InXs[:] = np.array([NDJFMA_mean[:,:,:,6],NDJFMA_median[:,:,:,6],NDJFMA_max[:,:,:,6],NDJFMA_min[:,:,:,6],NDJFMA_75[:,:,:,6],NDJFMA_25[:,:,:,6]])
	HW_NDJFMA_InMs[:] = np.array([NDJFMA_mean[:,:,:,7],NDJFMA_median[:,:,:,7],NDJFMA_max[:,:,:,7],NDJFMA_min[:,:,:,7],NDJFMA_75[:,:,:,7],NDJFMA_25[:,:,:,7]])
	
	statss.long_name='0-ensemble mean,1-ensemble median,2-ensemble max,3-ensemble min,4-ensemble 75th percentile,5-ensemble 25th percentile'
	HW_NDJFMA_Ms.long_name='Heatwvae magnitude: NDJFMA mean'
	HW_NDJFMA_Xs.long_name='Heatwvae magnitude: NDJFMA max'
	HW_NDJFMA_NOs.long_name='Heatwvae : NDJFMA no of evnents'
	HW_NDJFMA_DCs.long_name='total dayso of HW in NDJFMA'
	HW_NDJFMA_DuXs.long_name='maximum Duration of NDJFMA heatwaves '
	HW_NDJFMA_DuMs.long_name='mean Duration of NDJFMA heatwaves '
	HW_NDJFMA_InXs.long_name='maximum intensity of NDJFMA heatwaves'
	HW_NDJFMA_InMs.long_name='mean intensity of NDJFMA heatwaves'

	f.description = 'Heat wave and Cold Spell indexs'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	f.close()
	nc_fid.close()

