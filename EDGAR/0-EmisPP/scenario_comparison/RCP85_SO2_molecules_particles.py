"""
This is to process SO2 eissions into SO4 (0.25%) and so2 (99.75%) for the Aikten mode (a2) and accumulation mode (a1)
The output would be SO2_suf, SO2_ele, SO4_a1_suf,SO4_a1_ele, SO4_a2_suf, SO4_a2_ele				
"""
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d

## read the CESM grids to interpolate the RCP emissionsn firstly before caculations
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
CESM_fid.close()
# rcp emissions
file_path ="/exports/csce/datastore/geos/users/s1667168/RCP/"
rcp_file = file_path+'accmip_interpolated_emissions_RCP85_SO2_ships_2010_0.5x0.5.nc'
nc_fid = nc4.Dataset(rcp_file,mode='r')
SHP =  nc_fid.variables['emiss_shp'][:]*0.1; # 0.1 from kg/m2/s to g/cm2/s
nc_fid.close()
rcp_file = file_path+'accmip_interpolated_emissions_RCP85_SO2_anthropogenic_2010_0.5x0.5.nc'
nc_fid = nc4.Dataset(rcp_file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
IND =  nc_fid.variables['emiss_ind'][:]*0.1;
ENE =  nc_fid.variables['emiss_ene'][:]*0.1;
WST =  nc_fid.variables['emiss_wst'][:]*0.1;
AWB =  nc_fid.variables['emiss_awb'][:]*0.1;
DOM =  nc_fid.variables['emiss_dom'][:]*0.1;
TRA =  nc_fid.variables['emiss_tra'][:]*0.1;
# nc_fid.close()
## arrays to store caculated no. of particles and molecules 
densities={'SO4':1.77,'BC':1.7,'POM':1.0}    # conversion to number emissions uses dey density of 1.77 1.7 1.0 g/cm3 for SO4, BC and POM   g/cm3
density=1.77
NA=6.022*10**(23)		#NA This factor (= 1.0e3 * avogadro_number) is needed to work with the CAM3 surface emissions code
Sectors = {'ENE':ENE,'IND':IND,'SHP':SHP,'WST':WST,'AWB':AWB,'DOM':DOM,'TRA':TRA}
diameters={'ENE':0.261,'IND':0.261,'SHP':0.261,'WST':0.134,'AWB':0.134,'DOM':0.0504 ,'TRA':0.0504} 

SO2_surf_mol=np.zeros((12,96,144));SO2_ele_mol=np.zeros((12,96,144));
S04_a1_srf_mol=np.zeros((12,96,144));S04_a1_ele_mol=np.zeros((12,96,144));
S04_a2_srf_mol=np.zeros((12,96,144));
S04_a1_srf_num=np.zeros((12,96,144));S04_a1_ele_num=np.zeros((12,96,144));
S04_a2_srf_num=np.zeros((12,96,144));	
months= range(0,12)
for imonth in months:
	for sector in Sectors:
		diameter=diameters[sector]*10**(-4);   #convert from micrometer into centimeter
		volume=np.pi*np.power(diameter,3)/6.0;  #pi/6*D3
		m_p = volume*density
		SO2 = Sectors[sector][imonth,:,:]    # kg/m2/s in EDGAR into molecules/cm2/s in CESM  1000/100/100/
		f = interp2d(lon,lat,SO2,kind='linear')
		SO2_mass = f(lon_CESM, lat_CESM);
		# print stats.nanmean(stats.nanmean(SO2_mass,axis=0),axis=0)
		if sector in ['AWB','DOM','TRA','WST','SHP']:
			SO2_surf_mol[imonth,:,:] =SO2_surf_mol[imonth,:,:]+SO2_mass*0.975*NA/64  
			if sector in ['AWB','WST','SHP']:
				S04_a1_srf_mol[imonth,:,:]=S04_a1_srf_mol[imonth,:,:]+SO2_mass*0.025*1.8*NA/115
				S04_a1_srf_num[imonth,:,:]=S04_a1_srf_num[imonth,:,:]+6.022*10**(26)*SO2_mass*0.025*1.8/m_p
			else:
				S04_a2_srf_mol[imonth,:,:]=S04_a2_srf_mol[imonth,:,:]+SO2_mass*0.025*1.8*NA/115
				S04_a2_srf_num[imonth,:,:]=S04_a2_srf_num[imonth,:,:]+6.022*10**(26)*SO2_mass*0.025*1.8/m_p
		else:
			print sector
			SO2_ele_mol[imonth,:,:] = SO2_ele_mol[imonth,:,:]+SO2_mass*0.975*NA/64  
			S04_a1_ele_mol[imonth,:,:] = S04_a1_ele_mol[imonth,:,:]+SO2_mass*0.025*1.8*NA/115
			S04_a1_ele_num[imonth,:,:] = S04_a1_ele_num[imonth,:,:]+6.022*10**(26)*SO2_mass*0.025*1.8/m_p
					
# output results into .nc file 
file_out = file_path+'RCP85_2010_SO2_NUM_MOL.0.1x0.1.nc'
f = nc4.Dataset(file_out,'w', format='NETCDF4') #'w' stands for write
f.createDimension('time', 12)
f.createDimension('lat', len(lat_CESM))
f.createDimension('lon', len(lon_CESM))
latitudes = f.createVariable('lat',np.float32, ('lat'))
longitudes = f.createVariable('lon',np.float32, ('lon'))
times=np.arange(1,13)+2010*100
latitudes[:] =lat_CESM
longitudes[:] = lon_CESM		
SO2_surf_mols=f.createVariable('SO2_surf_mol',np.float32,('time','lat','lon'))
SO2_ele_mols=f.createVariable('SO2_ele_mol',np.float32,('time','lat','lon'))		
S04_a1_srf_mols=f.createVariable('S04_a1_srf_mol',np.float32,('time','lat','lon'))	
S04_a1_ele_mols=f.createVariable('S04_a1_ele_mol',np.float32,('time','lat','lon'))
S04_a2_srf_mols=f.createVariable('S04_a2_srf_mol',np.float32,('time','lat','lon'))
S04_a1_srf_nums=f.createVariable('S04_a1_srf_num',np.float32,('time','lat','lon'))
S04_a1_ele_nums=f.createVariable('S04_a1_ele_num',np.float32,('time','lat','lon'))
S04_a2_srf_nums=f.createVariable('S04_a2_srf_num',np.float32,('time','lat','lon'))

SO2_surf_mols[:]=SO2_surf_mol
SO2_ele_mols[:]=SO2_ele_mol
S04_a1_srf_mols[:]=S04_a1_srf_mol
S04_a2_srf_nums[:]=S04_a2_srf_num
S04_a1_ele_mols[:]=S04_a1_ele_mol		
S04_a2_srf_mols[:]=S04_a2_srf_mol
S04_a1_srf_nums[:]=S04_a1_srf_num
S04_a1_ele_nums[:]=S04_a1_ele_num		

SO2_surf_mols.units='molecules/cm2/s'
SO2_surf_mols.long_name ='SO2 surface emission in molecules,97.5 percent from the origiinal emissions with all sectors added up'

SO2_ele_mols.units='molecules/cm2/s'
SO2_ele_mols.long_name ='SO2 elevated emission in molecules,97.5 percent from the origiinal emissions with all sectors added up'

S04_a1_srf_mols.units='particles/cm2/s*6.0e26'
S04_a1_srf_mols.long_name ='SO4 surface emission in particle numbers, 2.5 percent from the origiinal emissions with all sectors added up'

S04_a1_ele_mols.units='particles/cm2/s*6.0e26'
S04_a1_ele_mols.long_name ='SO4 elevated emission in particle numbers, 2.5 percent from the origiinal emissions with all sectors added up'


S04_a2_srf_mols.units='particles/cm2/s*6.0e26'
S04_a2_srf_mols.long_name ='SO4 surface emission in particle numbers, 2.5 percent from the origiinal emissions with all sectors added up'


S04_a1_srf_nums.units='particles/cm2/s*6.0e26'
S04_a1_srf_nums.long_name ='SO4 surface emission in particle numbers, 2.5 from the origiinal emissions with all sectors added up'

S04_a1_ele_nums.units='particles/cm2/s*6.0e26'
S04_a1_ele_nums.long_name ='SO4 elevated emission in particle numbers, 2.5 percent from the origiinal emissions with all sectors added up'


S04_a2_srf_nums.units='particles/cm2/s*6.0e26'
S04_a2_srf_nums.long_name ='SO4 surface emission in particle numbers, 2.5 percent from the origiinal emissions with all sectors added up'	

f.title='SO2 emissions and sulphate aerosol molecules/particles numbers for a1 and a2 at both surface level and elevated levels'
f.close()
	
		
"""
file ="/exports/csce/datastore/geos/users/s1667168/EDGAR/2010/SO2/v431_v2_REFERENCE_SO2_2010_2_NUM_MOL.0.1x0.1.nc"
print file
nc_fid = nc4.Dataset(file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:] 
var=nc_fid.variables['SO2_ele_mol'][:] ; var[var==0]=np.nan
plt.imshow(var,origin='lower');plt.show()

"""

