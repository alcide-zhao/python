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
	en_mean = np.empty((YearNo,192,288,21));
	en_median = np.empty((YearNo,192,288,21));
	en_max = np.empty((YearNo,192,288,21));
	en_min = np.empty((YearNo,192,288,21));
	en_75 = np.empty((YearNo,192,288,21));
	en_25 = np.empty((YearNo,192,288,21));
	for iyear in range(YearNo):
		att_value = np.empty((ensemble,192,288,21))
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


scenario_dic ={'rcp85':[20,0,95]}  #,'fixa':[14,0,95],'rcp45':[15,0,75],'his':[22,0,86],'rcp85':[21,0,95]
directory = 'CalDayThr_InterOsi'  #CalDayThr_InterOsi
for scenario in scenario_dic:
	print scenario
	ensemble = scenario_dic[scenario][0]; YearNo=scenario_dic[scenario][2];
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'+directory+'/'+scenario+'/'  #_InterOsi
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + 'file_list.txt', "r")
	text_content = text_file.readlines()

	variable='HWI';HWI_mean,HWI_median,HWI_max,HWI_min,HWI_75,HWI_25 = ensemble_pp(text_content,variable,YearNo)
	variable='CSI';CSI_mean,CSI_median,CSI_max,CSI_min,CSI_75,CSI_25 = ensemble_pp(text_content,variable,YearNo)
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
	
	HWI_Ms= f.createVariable('HWI_M',np.float32,('stats','year','lat','lon'))
	HWI_Ns= f.createVariable('HWI_N',np.float32,('stats','year','lat','lon'))
	HWI_Xs= f.createVariable('HWI_X',np.float32,('stats','year','lat','lon'))
	HWI_NOs=f.createVariable('HWI_NO',np.float32,('stats','year','lat','lon'))
	HWCMs= f.createVariable('HWCM',np.float32,('stats','year','lat','lon'))
	HWCNs= f.createVariable('HWCN',np.float32,('stats','year','lat','lon'))
	HWCDs= f.createVariable('HWCD',np.float32,('stats','year','lat','lon'))
	HWCSs= f.createVariable('HWCS',np.float32,('stats','year','lat','lon'))
	HWCEs= f.createVariable('HWCE',np.float32,('stats','year','lat','lon'))
	HWCVs= f.createVariable('HWCV',np.float32,('stats','year','lat','lon'))
	HWCPs= f.createVariable('HWCP',np.float32,('stats','year','lat','lon'))
	HWCUs= f.createVariable('HWCU',np.float32,('stats','year','lat','lon'))
	HWRFMs= f.createVariable('HWRFM',np.float32,('stats','year','lat','lon'))
	HWRFNs= f.createVariable('HWRFN',np.float32,('stats','year','lat','lon'))
	HWRFXs= f.createVariable('HWRFX',np.float32,('stats','year','lat','lon'))
	HWDCs= f.createVariable('HWDC',np.float32,('stats','year','lat','lon'))
	RFSMs= f.createVariable('RFSM',np.float32,('stats','year','lat','lon'))
	HWDuXs= f.createVariable('HWDuX',np.float32,('stats','year','lat','lon'))
	HWDuMs= f.createVariable('HWDuM',np.float32,('stats','year','lat','lon'))
	HWInXs= f.createVariable('HWInX',np.float32,('stats','year','lat','lon'))
	HWInMs= f.createVariable('HWInM',np.float32,('stats','year','lat','lon'))

	HWI_Ms[:] = np.array([HWI_mean[:,:,:,0],HWI_median[:,:,:,0],HWI_max[:,:,:,0],HWI_min[:,:,:,0],HWI_75[:,:,:,0],HWI_25[:,:,:,0]])
	HWI_Ns[:] = np.array([HWI_mean[:,:,:,1],HWI_median[:,:,:,1],HWI_max[:,:,:,1],HWI_min[:,:,:,1],HWI_75[:,:,:,1],HWI_25[:,:,:,1]])
	HWI_Xs[:] = np.array([HWI_mean[:,:,:,2],HWI_median[:,:,:,2],HWI_max[:,:,:,2],HWI_min[:,:,:,2],HWI_75[:,:,:,2],HWI_25[:,:,:,2]])
	HWI_NOs[:]= np.array([HWI_mean[:,:,:,3],HWI_median[:,:,:,3],HWI_max[:,:,:,3],HWI_min[:,:,:,3],HWI_75[:,:,:,3],HWI_25[:,:,:,3]])
	HWCMs[:] = np.array([HWI_mean[:,:,:,4],HWI_median[:,:,:,4],HWI_max[:,:,:,4],HWI_min[:,:,:,4],HWI_75[:,:,:,4],HWI_25[:,:,:,4]])
	HWCNs[:] = np.array([HWI_mean[:,:,:,5],HWI_median[:,:,:,5],HWI_max[:,:,:,5],HWI_min[:,:,:,5],HWI_75[:,:,:,5],HWI_25[:,:,:,5]])
	HWCDs[:] = np.array([HWI_mean[:,:,:,6],HWI_median[:,:,:,6],HWI_max[:,:,:,6],HWI_min[:,:,:,6],HWI_75[:,:,:,6],HWI_25[:,:,:,6]])
	HWCSs[:] = np.array([HWI_mean[:,:,:,7],HWI_median[:,:,:,7],HWI_max[:,:,:,7],HWI_min[:,:,:,7],HWI_75[:,:,:,7],HWI_25[:,:,:,7]])
	HWCEs[:] = np.array([HWI_mean[:,:,:,8],HWI_median[:,:,:,8],HWI_max[:,:,:,8],HWI_min[:,:,:,8],HWI_75[:,:,:,8],HWI_25[:,:,:,8]])
	HWCVs[:] = np.array([HWI_mean[:,:,:,9],HWI_median[:,:,:,9],HWI_max[:,:,:,9],HWI_min[:,:,:,9],HWI_75[:,:,:,9],HWI_25[:,:,:,9]])
	HWCPs[:] = np.array([HWI_mean[:,:,:,10],HWI_median[:,:,:,10],HWI_max[:,:,:,10],HWI_min[:,:,:,10],HWI_75[:,:,:,10],HWI_25[:,:,:,10]])
	HWCUs[:] = np.array([HWI_mean[:,:,:,11],HWI_median[:,:,:,11],HWI_max[:,:,:,11],HWI_min[:,:,:,11],HWI_75[:,:,:,11],HWI_25[:,:,:,11]])
	HWRFMs[:] = np.array([HWI_mean[:,:,:,12],HWI_median[:,:,:,12],HWI_max[:,:,:,12],HWI_min[:,:,:,12],HWI_75[:,:,:,12],HWI_25[:,:,:,12]])
	HWRFNs[:] = np.array([HWI_mean[:,:,:,13],HWI_median[:,:,:,13],HWI_max[:,:,:,13],HWI_min[:,:,:,13],HWI_75[:,:,:,13],HWI_25[:,:,:,13]])
	HWRFXs[:] = np.array([HWI_mean[:,:,:,14],HWI_median[:,:,:,14],HWI_max[:,:,:,14],HWI_min[:,:,:,14],HWI_75[:,:,:,14],HWI_25[:,:,:,14]])
	HWDCs[:] = np.array([HWI_mean[:,:,:,15],HWI_median[:,:,:,15],HWI_max[:,:,:,15],HWI_min[:,:,:,15],HWI_75[:,:,:,15],HWI_25[:,:,:,15]])
	RFSMs[:] = np.array([HWI_mean[:,:,:,16],HWI_median[:,:,:,16],HWI_max[:,:,:,16],HWI_min[:,:,:,16],HWI_75[:,:,:,16],HWI_25[:,:,:,16]])
	HWDuXs[:] = np.array([HWI_mean[:,:,:,17],HWI_median[:,:,:,17],HWI_max[:,:,:,17],HWI_min[:,:,:,17],HWI_75[:,:,:,17],HWI_25[:,:,:,17]])
	HWDuMs[:] = np.array([HWI_mean[:,:,:,18],HWI_median[:,:,:,18],HWI_max[:,:,:,18],HWI_min[:,:,:,18],HWI_75[:,:,:,18],HWI_25[:,:,:,18]])
	HWInXs[:] = np.array([HWI_mean[:,:,:,19],HWI_median[:,:,:,19],HWI_max[:,:,:,19],HWI_min[:,:,:,19],HWI_75[:,:,:,19],HWI_25[:,:,:,19]])
	HWInMs[:] = np.array([HWI_mean[:,:,:,20],HWI_median[:,:,:,20],HWI_max[:,:,:,20],HWI_min[:,:,:,20],HWI_75[:,:,:,20],HWI_25[:,:,:,20]])

	statss.long_name='0-ensemble mean,1-ensemble median,2-ensemble max,3-ensemble min,4-ensemble 75th percentile,5-ensemble 25th percentile'
	HWI_Ms.long_name='Heatwvae magnitude: annual mean'
	HWI_Ns.long_name='Heatwvae magnitude: annual min'
	HWI_Xs.long_name='Heatwvae magnitude: annual max'
	HWI_NOs.long_name='Heatwvae : annual no of evnents'
	HWCMs.long_name='no of annual heatwaves with normalized heatwave magitude lies between [-inf,0]'
	HWCNs.long_name='no of annual heatwaves with normalized heatwave magitude lies between (0,2]'
	HWCDs.long_name='no of annual heatwaves with normalized heatwave magitude lies between (2,4]'
	HWCSs.long_name='no of annual heatwaves with normalized heatwave magitude lies between (4,8]'
	HWCEs.long_name='no of annual heatwaves with normalized heatwave magitude lies between (8,16]'
	HWCVs.long_name='no of annual heatwaves with normalized heatwave magitude lies between (16,32]'
	HWCPs.long_name='no of annual heatwaves with normalized heatwave magitude lies between (32,64]'
	HWCUs.long_name='no of annual heatwaves with normalized heatwave magitude lies between [-64,inf)'
	HWRFMs.long_name='recover factor of heatwaves: event mean'
	HWRFNs.long_name='recover factor of heatwaves: event min'
	HWRFXs.long_name='recover factor of heatwaves: event max'
	HWDCs.long_name='annual total number of HeatWave days '
	RFSMs.long_name='recover factor of summer: JJAS'
	HWDuXs.long_name='maximum Duration of annual heatwaves '
	HWDuMs.long_name='mean Duration of annual heatwaves '
	HWInXs.long_name='maximum intensity of annual heatwaves'
	HWInMs.long_name='mean intensity of annual heatwaves'
	
	CSI_Ms= f.createVariable('CSI_M',np.float32,('stats','year','lat','lon'))
	CSI_Ns= f.createVariable('CSI_N',np.float32,('stats','year','lat','lon'))
	CSI_Xs= f.createVariable('CSI_X',np.float32,('stats','year','lat','lon'))
	CSI_NOs=f.createVariable('CSI_NO',np.float32,('stats','year','lat','lon'))
	CSCMs= f.createVariable('CSCM',np.float32,('stats','year','lat','lon'))
	CSCNs= f.createVariable('CSCN',np.float32,('stats','year','lat','lon'))
	CSCDs= f.createVariable('CSCD',np.float32,('stats','year','lat','lon'))
	CSCSs= f.createVariable('CSCS',np.float32,('stats','year','lat','lon'))
	CSCEs= f.createVariable('CSCE',np.float32,('stats','year','lat','lon'))
	CSCVs= f.createVariable('CSCV',np.float32,('stats','year','lat','lon'))
	CSCPs= f.createVariable('CSCP',np.float32,('stats','year','lat','lon'))
	CSCUs= f.createVariable('CSCU',np.float32,('stats','year','lat','lon'))
	CSRFMs= f.createVariable('CSRFM',np.float32,('stats','year','lat','lon'))
	CSRFNs= f.createVariable('CSRFN',np.float32,('stats','year','lat','lon'))
	CSRFXs= f.createVariable('CSRFX',np.float32,('stats','year','lat','lon'))
	CSDCs= f.createVariable('CSDC',np.float32,('stats','year','lat','lon'))
	RFWTs= f.createVariable('RFWT',np.float32,('stats','year','lat','lon'))
	CSDuXs= f.createVariable('CSDuX',np.float32,('stats','year','lat','lon'))
	CSDuMs= f.createVariable('CSDuM',np.float32,('stats','year','lat','lon'))
	CSInXs= f.createVariable('CSInX',np.float32,('stats','year','lat','lon'))
	CSInMs= f.createVariable('CSInM',np.float32,('stats','year','lat','lon'))

	CSI_Ms[:] = np.array([CSI_mean[:,:,:,0],CSI_median[:,:,:,0],CSI_max[:,:,:,0],CSI_min[:,:,:,0],CSI_75[:,:,:,0],CSI_25[:,:,:,0]])
	CSI_Ns[:] = np.array([CSI_mean[:,:,:,1],CSI_median[:,:,:,1],CSI_max[:,:,:,1],CSI_min[:,:,:,1],CSI_75[:,:,:,1],CSI_25[:,:,:,1]])
	CSI_Xs[:] = np.array([CSI_mean[:,:,:,2],CSI_median[:,:,:,2],CSI_max[:,:,:,2],CSI_min[:,:,:,2],CSI_75[:,:,:,2],CSI_25[:,:,:,2]])
	CSI_NOs[:]= np.array([CSI_mean[:,:,:,3],CSI_median[:,:,:,3],CSI_max[:,:,:,3],CSI_min[:,:,:,3],CSI_75[:,:,:,3],CSI_25[:,:,:,3]])
	CSCMs[:] = np.array([CSI_mean[:,:,:,4],CSI_median[:,:,:,4],CSI_max[:,:,:,4],CSI_min[:,:,:,4],CSI_75[:,:,:,4],CSI_25[:,:,:,4]])
	CSCNs[:] = np.array([CSI_mean[:,:,:,5],CSI_median[:,:,:,5],CSI_max[:,:,:,5],CSI_min[:,:,:,5],CSI_75[:,:,:,5],CSI_25[:,:,:,5]])
	CSCDs[:] = np.array([CSI_mean[:,:,:,6],CSI_median[:,:,:,6],CSI_max[:,:,:,6],CSI_min[:,:,:,6],CSI_75[:,:,:,6],CSI_25[:,:,:,6]])
	CSCSs[:] = np.array([CSI_mean[:,:,:,7],CSI_median[:,:,:,7],CSI_max[:,:,:,7],CSI_min[:,:,:,7],CSI_75[:,:,:,7],CSI_25[:,:,:,7]])
	CSCEs[:] = np.array([CSI_mean[:,:,:,8],CSI_median[:,:,:,8],CSI_max[:,:,:,8],CSI_min[:,:,:,8],CSI_75[:,:,:,8],CSI_25[:,:,:,8]])
	CSCVs[:] = np.array([CSI_mean[:,:,:,9],CSI_median[:,:,:,9],CSI_max[:,:,:,9],CSI_min[:,:,:,9],CSI_75[:,:,:,9],CSI_25[:,:,:,9]])
	CSCPs[:] = np.array([CSI_mean[:,:,:,10],CSI_median[:,:,:,10],CSI_max[:,:,:,10],CSI_min[:,:,:,10],CSI_75[:,:,:,10],CSI_25[:,:,:,10]])
	CSCUs[:] = np.array([CSI_mean[:,:,:,11],CSI_median[:,:,:,11],CSI_max[:,:,:,11],CSI_min[:,:,:,11],CSI_75[:,:,:,11],CSI_25[:,:,:,11]])
	CSRFMs[:] = np.array([CSI_mean[:,:,:,12],CSI_median[:,:,:,12],CSI_max[:,:,:,12],CSI_min[:,:,:,12],CSI_75[:,:,:,12],CSI_25[:,:,:,12]])
	CSRFNs[:] = np.array([CSI_mean[:,:,:,13],CSI_median[:,:,:,13],CSI_max[:,:,:,13],CSI_min[:,:,:,13],CSI_75[:,:,:,13],CSI_25[:,:,:,13]])
	CSRFXs[:] = np.array([CSI_mean[:,:,:,14],CSI_median[:,:,:,14],CSI_max[:,:,:,14],CSI_min[:,:,:,14],CSI_75[:,:,:,14],CSI_25[:,:,:,14]])
	CSDCs[:] = np.array([CSI_mean[:,:,:,15],CSI_median[:,:,:,15],CSI_max[:,:,:,15],CSI_min[:,:,:,15],CSI_75[:,:,:,15],CSI_25[:,:,:,15]])
	RFWTs[:] = np.array([CSI_mean[:,:,:,16],CSI_median[:,:,:,16],CSI_max[:,:,:,16],CSI_min[:,:,:,16],CSI_75[:,:,:,16],CSI_25[:,:,:,16]])
	CSDuXs[:] = np.array([CSI_mean[:,:,:,17],CSI_median[:,:,:,17],CSI_max[:,:,:,17],CSI_min[:,:,:,17],CSI_75[:,:,:,17],CSI_25[:,:,:,17]])
	CSDuMs[:] = np.array([CSI_mean[:,:,:,18],CSI_median[:,:,:,18],CSI_max[:,:,:,18],CSI_min[:,:,:,18],CSI_75[:,:,:,18],CSI_25[:,:,:,18]])
	CSInXs[:] = np.array([CSI_mean[:,:,:,19],CSI_median[:,:,:,19],CSI_max[:,:,:,19],CSI_min[:,:,:,19],CSI_75[:,:,:,19],CSI_25[:,:,:,19]])
	CSInMs[:] = np.array([CSI_mean[:,:,:,20],CSI_median[:,:,:,20],CSI_max[:,:,:,20],CSI_min[:,:,:,20],CSI_75[:,:,:,20],CSI_25[:,:,:,20]])

	CSI_Ms.long_name='ColdSpell magnitude: annual mean'
	CSI_Ns.long_name='ColdSpell magnitude: annual min'
	CSI_Xs.long_name='ColdSpell magnitude: annual max'
	CSI_NOs.long_name='ColdSpell : annual no of evnents'
	CSCMs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between [-inf,0]'
	CSCNs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between (0,2]'
	CSCDs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between (2,4]'
	CSCSs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between (4,8]'
	CSCEs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between (8,16]'
	CSCVs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between (16,32]'
	CSCPs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between (32,64]'
	CSCUs.long_name='no of annual ColdSpell with normalized heatwave magitude lies between [-64,inf)'
	CSRFMs.long_name='recover factor of ColdSpell: event mean'
	CSRFNs.long_name='recover factor of ColdSpell: event min'
	CSRFXs.long_name='recover factor of ColdSpell: event max'
	CSDCs.long_name='annual total number of coldSpell days'
	RFWTs.long_name='recover factor of winter: JJAS'
	CSDuXs.long_name='maximum Duration of annual ColdSpell '
	CSDuMs.long_name='mean Duration of annual ColdSpell '
	CSInXs.long_name='maximum intensity of annual ColdSpell'
	CSInMs.long_name='mean intensity of annual ColdSpell'

	f.description = 'Heat wave and Cold Spell indexs'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	f.close()
	nc_fid.close()

