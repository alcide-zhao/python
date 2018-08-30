"""
For CAM-CHEM, VOCS species including BIGALK, BIGENE and TOLUNE have to be emitted seperately. 
This is to map the SOAG emissions from these VOCs back to VOCS species assuming molecular mass of 72,56 and 92 respectively
As such, the mass yield from VOC to SOAG of 5%, 5% and 15% have to be omitted, and also the POM/OC ratio (1.4)
"""


import netCDF4 as nc4
import numpy as np
import sys

emission = '1970' 
months=range(0,12)
path_SOAG = '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/CESM_rcp85/'+emission+'/'
path_VOC =  '/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/CESM_rcp85/VOCs/'
# SOAG_layer_b = 24; VOC_layer_b = 97;  #2010 stag_ene  stag_tech
SOAG_layer_b = 156; VOC_layer_b = 132;  #1970


edgar_SOAG = path_SOAG+'ar5_mam3_soag_1.5_surf_1850-2005_c100429.nc' #filename month starts from 1
print edgar_SOAG
edgar_fid=nc4.Dataset(edgar_SOAG,mode='r')
BIGALK = edgar_fid.variables['SOAG_BIGALK'][SOAG_layer_b:SOAG_layer_b+12,:,:]*12*1.4/0.05/1.5/72
BIGENE = edgar_fid.variables['SOAG_BIGENE'][SOAG_layer_b:SOAG_layer_b+12,:,:]*12*1.4/0.05/1.5/56
TOLUENE = edgar_fid.variables['SOAG_TOLUENE'][SOAG_layer_b:SOAG_layer_b+12,:,:]*12*1.4/0.15/1.5/92

# voc_BIGALK = path_VOC +'rcp85_finn_2002-2015_BIGALK_1.9x2.5_mol_c20160128.nc'
voc_BIGALK = path_VOC +'maccity_maccity_corrdates_BIGALK_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(voc_BIGALK,mode='r+')
nc_fid.variables['emiss_anthro'][VOC_layer_b:VOC_layer_b+12,:,:]=BIGALK
nc_fid.close()

# voc_BIGENE = path_VOC +'rcp85_finn_2002-2015_BIGENE_1.9x2.5_mol_c20160128.nc'
voc_BIGENE = path_VOC +'maccity_maccity_corrdates_BIGENE_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(voc_BIGENE,mode='r+')
nc_fid.variables['emiss_anthro'][VOC_layer_b:VOC_layer_b+12,:,:]=BIGENE
nc_fid.close()

# voc_TOLUEN = path_VOC +'rcp85_finn_2002-2015_TOLUENE_1.9x2.5_mol_c20160128.nc'
voc_TOLUEN = path_VOC +'maccity_maccity_corrdates_TOLUENE_SOA_1960-2008_1.9x2.5_mol_c130314.nc'
nc_fid = nc4.Dataset(voc_TOLUEN,mode='r+')
nc_fid.variables['emiss_anthro'][VOC_layer_b:VOC_layer_b+12,:,:]=TOLUENE
nc_fid.close()
	
	
	
	
	
	
	

	





