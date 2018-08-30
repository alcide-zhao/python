"""
This is to process SO2 eissions into SO4 (0.25%) and so2 (99.75%) for the Aikten mode (a2) and accumulation mode (a1)
The output would be SO2_suf, SO2_ele, SO4_a1_suf,SO4_a1_ele, SO4_a2_suf, SO4_a2_ele				
"""
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import glob  
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d

# # AGR: Agriculture
# # ENERGY: Power generation
# # INDUSTRY: Manufacturing industry
# # SHIP: Shipping emissions
# # RCO: Residential combustion
# # SWD: Solid waste disposal
# # TNG: Non road transport
# # TRO: Road transport
# # AVI: Aviation transport

species=['SO2'] 
emission='1970' #2010 or stag_tech stag_energy

Sectors = ['ENE','IND','SHIP','SWD','AGR','RCO','TRO','TNG','TRF','REF','PRO','PPA','OTH']
diameters={'ENE':0.261,'IND':0.261,'SHIP':0.261,'SWD':0.134,'AGR':0.134,'RCO':0.0504,'TRO':0.0504,'TNG':0.0504,\
	'TRF':0.261,'REF':0.261,'PRO':0.261,'PPA':0.261,'OTH':0.134,'CDS':0.0504,'CRS':0.0504,'LTO':0.0504,'SPS':0.0504} #micormeter
densities={'SO4':1.77,'BC':1.7,'POM':1.0}    # conversion to number emissions uses dey density of 1.77 1.7 1.0 g/cm3 for SO4, BC and POM   g/cm3
density=densities['SO4'];
NA=6.022*10**(23)		#NA This factor (= 1.0e3 * avogadro_number) is needed to work with the CAM3 surface emissions code


## INPUT AND OUTPUT FILE DIRECTORY
file_path ="/exports/csce/datastore/geos/users/s1667168/EDGAR/"
fn_dic = {'1970':'REFERENCE','2010':'REFERENCE','stag_tech':'STAG_TECH','stag_energy':'STAG_ENE'}
fn_op = fn_dic[emission]

## read the CESM grids to interpolate the RCP emissionsn firstly before caculations
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
# CESM_fid.close()

months=range(1,13) 
for specie in species:
	for month in months:
		SO2_surf_mol=np.zeros((1800,3600));SO2_ele_mol=np.zeros((1800,3600));
		S04_a1_srf_mol=np.zeros((1800,3600));S04_a1_ele_mol=np.zeros((1800,3600));
		S04_a2_srf_mol=np.zeros((1800,3600));
		S04_a1_srf_num=np.zeros((1800,3600));S04_a1_ele_num=np.zeros((1800,3600));
		S04_a2_srf_num=np.zeros((1800,3600));
		for sector in Sectors:
			diameter=diameters[sector]*10**(-4);   #convert from micrometer into centimeter
			volume=np.pi*np.power(diameter,3)/6.0;
			m_p = volume*density
			# print m_p
			# for stag_tech
			if emission in ['stag_tech','stag_energy']:
				if sector in ['ENE','IND','TRO']:
					file = file_path+emission+"/"+specie+'/'+'JRC_PEGASOS_V2_'+emission.upper()+'_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
				else:
					file = file_path+"2010/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_2010_'+str(month)+'_'+sector+'.0.1x0.1.nc'
			else:
				file = file_path+"1970/"+specie+'/'+'v431_v2_REFERENCE_'+specie+'_1970_'+str(month)+'_'+sector+'.nc'
			# print file
			nc_fid = nc4.Dataset(file,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			SO2_mass = nc_fid.variables['emi_'+specie.lower()][:]*0.1          # kg/m2/s in EDGAR into molecules/cm2/s in CESM  1000/100/100/
			# f = interp2d(lon,lat,SO2,kind='linear')
			# SO2_mass = f(lon_CESM, lat_CESM);
			if sector in ['AGR','SHIP','SWD','OTH','RCO','TRO','TNG','CDS','CRS','LTO','SPS']:		
				SO2_surf_mol =SO2_surf_mol+SO2_mass*0.975*NA/64  
				if sector in ['AGR','SHIP','SWD','OTH']: #
					S04_a1_srf_mol=S04_a1_srf_mol+SO2_mass*0.025*1.8*NA/115
					S04_a1_srf_num=S04_a1_srf_num+6.022*10**(26)*SO2_mass*0.025*1.8/m_p
				else:
					S04_a2_srf_mol=S04_a2_srf_mol+SO2_mass*0.025*1.8*NA/115
					S04_a2_srf_num=S04_a2_srf_num+6.022*10**(26)*SO2_mass*0.025*1.8/m_p
			else:
				# print sector
				SO2_ele_mol = SO2_ele_mol+SO2_mass*0.975*NA/64  
				S04_a1_ele_mol = S04_a1_ele_mol+SO2_mass*0.025*1.8*NA/115
				S04_a1_ele_num = S04_a1_ele_num+6.022*10**(26)*SO2_mass*0.025*1.8/m_p			
			nc_fid.close()
		# output results into .nc file 
		file_out = file_path+emission+'/'+fn_op+'_'+specie+'_1970_'+str(month)+'_NUM_MOL.0.1x0.1.nc'
		f = nc4.Dataset(file_out,'w', format='NETCDF4') 
		f.createDimension('lat', len(lat))
		f.createDimension('lon', len(lon))
		latitudes = f.createVariable('lat',np.float32, ('lat'))
		longitudes = f.createVariable('lon',np.float32, ('lon'))
		latitudes[:] = lat
		longitudes[:] = lon		
		SO2_surf_mols=f.createVariable('SO2_surf_mol',np.float32,('lat','lon'))
		SO2_ele_mols=f.createVariable('SO2_ele_mol',np.float32,('lat','lon'))		
		S04_a1_srf_mols=f.createVariable('S04_a1_srf_mol',np.float32,('lat','lon'))	
		S04_a1_ele_mols=f.createVariable('S04_a1_ele_mol',np.float32,('lat','lon'))
		S04_a2_srf_mols=f.createVariable('S04_a2_srf_mol',np.float32,('lat','lon'))
		S04_a1_srf_nums=f.createVariable('S04_a1_srf_num',np.float32,('lat','lon'))
		S04_a1_ele_nums=f.createVariable('S04_a1_ele_num',np.float32,('lat','lon'))
		S04_a2_srf_nums=f.createVariable('S04_a2_srf_num',np.float32,('lat','lon'))

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

