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



# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/abs_v2/'+scenario+'/*.nc'
# files=sorted(glob.glob(input_path))

scenario='NCEP'

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/abs_v2/'+scenario+'/'
nc_f = input_path+'HWI_CSI_NCEP_abs.nc'
nc_fid = nc4.Dataset(nc_f,mode='r')
HWI = nc_fid.variables['HWI'][:]
CSI = nc_fid.variables['CSI'][:]
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
year = nc_fid.variables['time'][:]
			

output_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/abs_v2/'
file_name_out = output_path+'HWI_CSI_EnPP_'+scenario+'.nc'
f = nc4.Dataset(file_name_out,'w', format='NETCDF4')
f.createDimension('year', len(year))
f.createDimension('lat', len(lat))
f.createDimension('lon', len(lon))

times = f.createVariable('year',np.float32, ('year'))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
times[:] = year
latitudes[:] = lat
longitudes[:] = lon

HWI_Ms= f.createVariable('HWI_M',np.float32,('year','lat','lon'))
HWI_Ns= f.createVariable('HWI_N',np.float32,('year','lat','lon'))
HWI_Xs= f.createVariable('HWI_X',np.float32,('year','lat','lon'))
HWI_NOs=f.createVariable('HWI_NO',np.float32,('year','lat','lon'))
HWCMs= f.createVariable('HWCM',np.float32,('year','lat','lon'))
HWCNs= f.createVariable('HWCN',np.float32,('year','lat','lon'))
HWCDs= f.createVariable('HWCD',np.float32,('year','lat','lon'))
HWCSs= f.createVariable('HWCS',np.float32,('year','lat','lon'))
HWCEs= f.createVariable('HWCE',np.float32,('year','lat','lon'))
HWCVs= f.createVariable('HWCV',np.float32,('year','lat','lon'))
HWCPs= f.createVariable('HWCP',np.float32,('year','lat','lon'))
HWCUs= f.createVariable('HWCU',np.float32,('year','lat','lon'))
HWRFMs= f.createVariable('HWRFM',np.float32,('year','lat','lon'))
HWRFNs= f.createVariable('HWRFN',np.float32,('year','lat','lon'))
HWRFXs= f.createVariable('HWRFX',np.float32,('year','lat','lon'))
HWDCs= f.createVariable('HWDC',np.float32,('year','lat','lon'))
RFSMs= f.createVariable('RFSM',np.float32,('year','lat','lon'))
HWDuXs= f.createVariable('HWDuX',np.float32,('year','lat','lon'))
HWDuMs= f.createVariable('HWDuM',np.float32,('year','lat','lon'))
HWInXs= f.createVariable('HWInX',np.float32,('year','lat','lon'))
HWInMs= f.createVariable('HWInM',np.float32,('year','lat','lon'))

HWI_Ms[:] =HWI[:,:,:,0]
HWI_Ns[:] = HWI[:,:,:,1]
HWI_Xs[:] = HWI[:,:,:,2]
HWI_NOs[:]= HWI[:,:,:,3]
HWCMs[:] = HWI[:,:,:,4]
HWCNs[:] = HWI[:,:,:,5]
HWCDs[:] = HWI[:,:,:,6]
HWCSs[:] = HWI[:,:,:,7]
HWCEs[:] = HWI[:,:,:,8]
HWCVs[:] = HWI[:,:,:,9]
HWCPs[:] = HWI[:,:,:,10]
HWCUs[:] = HWI[:,:,:,11]
HWRFMs[:] = HWI[:,:,:,12]
HWRFNs[:] = HWI[:,:,:,13]
HWRFXs[:] = HWI[:,:,:,14]
HWDCs[:] = HWI[:,:,:,15]
RFSMs[:] = HWI[:,:,:,16]
HWDuXs[:] = HWI[:,:,:,17]
HWDuMs[:] = HWI[:,:,:,18]
HWInXs[:] = HWI[:,:,:,19]
HWInMs[:] = HWI[:,:,:,20]

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

CSI_Ms= f.createVariable('CSI_M',np.float32,('year','lat','lon'))
CSI_Ns= f.createVariable('CSI_N',np.float32,('year','lat','lon'))
CSI_Xs= f.createVariable('CSI_X',np.float32,('year','lat','lon'))
CSI_NOs=f.createVariable('CSI_NO',np.float32,('year','lat','lon'))
CSCMs= f.createVariable('CSCM',np.float32,('year','lat','lon'))
CSCNs= f.createVariable('CSCN',np.float32,('year','lat','lon'))
CSCDs= f.createVariable('CSCD',np.float32,('year','lat','lon'))
CSCSs= f.createVariable('CSCS',np.float32,('year','lat','lon'))
CSCEs= f.createVariable('CSCE',np.float32,('year','lat','lon'))
CSCVs= f.createVariable('CSCV',np.float32,('year','lat','lon'))
CSCPs= f.createVariable('CSCP',np.float32,('year','lat','lon'))
CSCUs= f.createVariable('CSCU',np.float32,('year','lat','lon'))
CSRFMs= f.createVariable('CSRFM',np.float32,('year','lat','lon'))
CSRFNs= f.createVariable('CSRFN',np.float32,('year','lat','lon'))
CSRFXs= f.createVariable('CSRFX',np.float32,('year','lat','lon'))
CSDCs= f.createVariable('CSDC',np.float32,('year','lat','lon'))
RFWTs= f.createVariable('RFWT',np.float32,('year','lat','lon'))
CSDuXs= f.createVariable('CSDuX',np.float32,('year','lat','lon'))
CSDuMs= f.createVariable('CSDuM',np.float32,('year','lat','lon'))
CSInXs= f.createVariable('CSInX',np.float32,('year','lat','lon'))
CSInMs= f.createVariable('CSInM',np.float32,('year','lat','lon'))

CSI_Ms[:] = CSI[:,:,:,0]
CSI_Ns[:] = CSI[:,:,:,1]
CSI_Xs[:] = CSI[:,:,:,2]
CSI_NOs[:]= CSI[:,:,:,3]
CSCMs[:] = CSI[:,:,:,4]
CSCNs[:] = CSI[:,:,:,5]
CSCDs[:] = CSI[:,:,:,6]
CSCSs[:] = CSI[:,:,:,7]
CSCEs[:] = CSI[:,:,:,8]
CSCVs[:] = CSI[:,:,:,9]
CSCPs[:] = CSI[:,:,:,10]
CSCUs[:] = CSI[:,:,:,11]
CSRFMs[:] = CSI[:,:,:,12]
CSRFNs[:] =CSI[:,:,:,13]
CSRFXs[:] = CSI[:,:,:,14]
CSDCs[:] =CSI[:,:,:,15]
RFWTs[:] = CSI[:,:,:,16]
CSDuXs[:] = CSI[:,:,:,17]
CSDuMs[:] = CSI[:,:,:,18]
CSInXs[:] =CSI[:,:,:,19]
CSInMs[:] = CSI[:,:,:,20]

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

