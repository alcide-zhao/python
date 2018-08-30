"""
This is to rescale the NMVOCS from BIGENE BIGALK to the species that need to be emitted in CAM-CHEM
For this, assumption are made that 
"""


import netCDF4 as nc4
import numpy as np
import sys

emission = 'stag_tech' 
months=range(0,12)
layer_b = 97;       #2010 stag_ene  stag_tech
# layer_b = 132;    #1970

file_path = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/CESM_rcp85/'+emission+'/'


BIGENE_F = file_path+'rcp85_finn_2002-2015_BIGENE_1.9x2.5_mol_c20160128.nc' #filename month starts from 1
# BIGENE_F = file_path+'maccity_maccity_corrdates_BIGENE_1960-2008_1.9x2.5_mol_c130314.nc' #filename month starts from 1
BIGENE_fid=nc4.Dataset(BIGENE_F,mode='r')
BIGENE = BIGENE_fid.variables['anthro'][layer_b:layer_b+12,:,:]*56  # from molecules to gram assuming a molecular weight of 56
BIGENE_fid.close()

BIGALK_F = file_path+'rcp85_finn_2002-2015_BIGALK_1.9x2.5_mol_c20160128.nc' #filename month starts from 1
# BIGALK_F = file_path+'maccity_maccity_corrdates_BIGALK_1960-2008_1.9x2.5_mol_c130314.nc' #filename month starts from 1
BIGALK_fid=nc4.Dataset(BIGALK_F,mode='r')
BIGALK = BIGALK_fid.variables['anthro'][layer_b:layer_b+12,:,:]*72  # from molecules to gram assuming a molecular weight of 72
BIGALK_fid.close()


## C2H2
SPECIES_F = file_path +'rcp85_finn_2002-2015_C2H2_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_C2H2_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGENE*22.44/100/26
nc_fid.close()
	
## C2H4
SPECIES_F = file_path +'rcp85_finn_2002-2015_C2H4_woBiog_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_C2H4_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGENE*24.49/100/28
nc_fid.close()
	
## C3H6
SPECIES_F = file_path +'rcp85_finn_2002-2015_C3H6_woBiog_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_C3H6_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGENE*53.06/100/42
nc_fid.close()	
	
## C2H5OH
SPECIES_F = file_path +'rcp85_finn_2002-2015_C2H5OH_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_C2H5OH_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*5.37/100/46
nc_fid.close()		
	
## C2H6
SPECIES_F = file_path +'rcp85_finn_2002-2015_C2H6_woBiog_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_C2H6_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*7.89/100/30
nc_fid.close()		

## c3h8
SPECIES_F = file_path +'rcp85_finn_2002-2015_C3H8_woBiog_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_C3H8_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*9.23/100/44
nc_fid.close()	

## ch2o
SPECIES_F = file_path +'rcp85_finn_2002-2015_CH2O_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_CH2O_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*0.84/100/30
nc_fid.close()	

## ch3cho
SPECIES_F = file_path +'rcp85_finn_2002-2015_CH3CHO_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_CH3CHO_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*1.68/100/44
nc_fid.close()	

## ch3coch3
SPECIES_F = file_path +'rcp85_finn_2002-2015_CH3CHO_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_CH3COCH3_woBiog_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*1.8/100/58
nc_fid.close()

## ch3cooh
SPECIES_F = file_path +'rcp85_finn_2002-2015_CH3CHO_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_CH3COOH_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*3.0/100/60
nc_fid.close()

## hcooh
SPECIES_F = file_path +'rcp85_finn_2002-2015_CH3CHO_1.9x2.5_mol_c20160128.nc'
# SPECIES_F = file_path +'maccity_maccity_corrdates_HCOOH_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(SPECIES_F,mode='r+')
nc_fid.variables['anthro'][layer_b:layer_b+12,:,:]=BIGALK*3.0/100/46
nc_fid.close()