import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d

months=range(0,12)
path_EDGAR = '/exports/csce/datastore/geos/users/s1667168/EDGAR/2010/'
path_CESM = '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/'

def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area



	
#################################
# compare RCP85 with CESM
#################################
"""
####so2_surf
so2_surf_file=path_CESM+'RCP85_mam3_num_a1_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
print 'month', '        CESM ', '           EDGAR ', '        CESM-EDGAR ',  '   % '
for imonth in months:
	EDGAR = path_EDGAR+'REFERENCE_SO2_2010_'+str(imonth+1)+'_NUM_MOL.0.1x0.1.nc'  #filename month starts from 1
	# print EDGAR
	EDGAR_fid=nc4.Dataset(EDGAR,mode='r')
	SO2_EDGAR = EDGAR_fid.variables['S04_a1_srf_num'][:] 
	emiss_awb=CESM_fid.variables['SO4_emiss_awb'][24+imonth,:,:]
	# emiss_dom=CESM_fid.variables['emiss_dom'][24+imonth,:,:]
	# emiss_tra=CESM_fid.variables['emiss_tra'][24+imonth,:,:]
	emiss_shp=CESM_fid.variables['SO4_emiss_shp'][24+imonth,:,:]
	emiss_wst=CESM_fid.variables['SO4_emiss_wst'][24+imonth,:,:]
	SO2_CESM=emiss_awb+emiss_shp+emiss_wst#emiss_dom+emiss_tra+;
	# f = interp2d(lon_RCP_85,Lat_RCP_85,SO2_RCP_o,kind='linear')
	# SO2_RCP = f(lon_CESM, lat_CESM);
	lon_res = lon_CESM[1] - lon_CESM[0];lat_res = lat_CESM[1] - lat_CESM[0];
	lons,lats = np.meshgrid (lon_CESM,lat_CESM)  #meshgrid
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	EDGAR_gc_sum = np.multiply(SO2_EDGAR,area)   # area of the grid cell in m2 * emission/m2/year
	CESM_gc_sum = np.multiply(SO2_CESM,area)   # area of the grid cell in m2 * emission/m2/year
	EDGAR = np.nansum(np.nansum(EDGAR_gc_sum,axis=1),axis=0)
	CESM = np.nansum(np.nansum(CESM_gc_sum,axis=1),axis=0)
	print imonth, CESM, EDGAR, CESM-EDGAR,round(100*(CESM-EDGAR)/EDGAR,2)
EDGAR_fid.close()
CESM_fid.close()

############ S04_a1_srf_mol  

print '#################################'
print '       S04_a2_srf_mol             '
print '#################################'
so2_surf_file=path_CESM+'RCP85_mam3_so4_a2_surf_2000-2100_c20110913.nc'

CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
print 'month', '        CESM ', '           RCP ', '        CESM-RCP ',  '   % '
for imonth in months:
	RCP_85 = path_RCP85+'RCP85_2010_SO2_NUM_MOL.0.1x0.1.nc'  #filename month starts from 1
	RCP_85_fid=nc4.Dataset(RCP_85,mode='r')
	SO2_RCP = RCP_85_fid.variables['S04_a2_srf_mol'][imonth,:,:]
	lon_RCP_85 = RCP_85_fid.variables['lon'][:]
	Lat_RCP_85 = RCP_85_fid.variables['lat'][:]
	
	# ####S04_a1_srf_mol
	# emiss_awb=CESM_fid.variables['emiss_awb'][24+imonth,:,:]
	# emiss_wst=CESM_fid.variables['emiss_wst'][24+imonth,:,:]
	# emiss_shp=CESM_fid.variables['emiss_shp'][24+imonth,:,:]
	# SO2_CESM=emiss_awb+emiss_wst+emiss_shp
	
	###S04_a2_srf_mol
	emiss_tra=CESM_fid.variables['emiss_tra'][24+imonth,:,:]
	emiss_dom=CESM_fid.variables['emiss_dom'][24+imonth,:,:]
	SO2_CESM=emiss_dom+emiss_tra
	
	###num_a1_SO4 
	# SO4_emiss_awb=CESM_fid.variables['SO4_emiss_awb'][24+imonth,:,:]
	# SO4_emiss_wst=CESM_fid.variables['SO4_emiss_wst'][24+imonth,:,:]
	# SO4_emiss_shp=CESM_fid.variables['SO4_emiss_shp'][24+imonth,:,:]
	# SO2_CESM=SO4_emiss_awb+SO4_emiss_wst+SO4_emiss_shp
	
	# # ####num_a2_SO4
	# SO4_emiss_tra=CESM_fid.variables['SO4_emiss_tra'][24+imonth,:,:]
	# SO4_emiss_dom=CESM_fid.variables['SO4_emiss_dom'][24+imonth,:,:]
	# SO2_CESM=SO4_emiss_dom+SO4_emiss_tra
	
	####S04_a1_elev_mol
	# emiss_ene=CESM_fid.variables['emiss_ene'][24+imonth,:,:,:]
	# emiss_ind=CESM_fid.variables['emiss_ind'][24+imonth,:,:,:]
	# SO2_CESM=np.nansum(emiss_ene+emiss_ind,axis=0)
	# SO2_CESM=emiss_awb+emiss_wst+emiss_shp+emiss_dom+emiss_tra
	
	# f = interp2d(lon_RCP_85,Lat_RCP_85,SO2_RCP_o,kind='linear')
	# SO2_RCP = f(lon_CESM, lat_CESM);
	
	lon_res = lon_CESM[1] - lon_CESM[0];lat_res = lat_CESM[1] - lat_CESM[0];
	lons,lats = np.meshgrid (lon_CESM,lat_CESM)  #meshgrid
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	rcp_gc_sum = np.multiply(SO2_RCP,area)   # area of the grid cell in m2 * emission/m2/year
	CESM_gc_sum = np.multiply(SO2_CESM,area)   # area of the grid cell in m2 * emission/m2/year
	RCP = np.nansum(np.nansum(rcp_gc_sum,axis=1),axis=0)
	CESM = np.nansum(np.nansum(CESM_gc_sum,axis=1),axis=0)
	print imonth, CESM, RCP, CESM-RCP,round(100*(CESM-RCP)/RCP,2)
RCP_85_fid.close()
CESM_fid.close()
"""

#################################
# compare RCP85 with EDGAR
#################################

PATH_EDGAR = '/exports/csce/datastore/geos/users/s1667168/EDGAR/2010/'
path_RCP85 = '/exports/csce/datastore/geos/users/s1667168/RCP/'
####so2_surf
so2_surf_file=path_CESM+'RCP85_mam3_so2_surf_2000-2100_c20110913.nc'
CESM_fid = nc4.Dataset(so2_surf_file,mode='r')
lat_CESM = CESM_fid.variables['lat'][:]
lon_CESM = CESM_fid.variables['lon'][:]
print 'month', '        EDGAR ', '           RCP ', '        EDGAR-RCP ',  '   % '
fig = plt.figure(facecolor='White',figsize=(10, 3));
for imonth in months:
	RCP_85 = path_RCP85+'RCP85_2010_SO2_NUM_MOL.0.1x0.1.nc'  #filename month starts from 1
	RCP_85_fid=nc4.Dataset(RCP_85,mode='r')
	SO2_RCP = RCP_85_fid.variables['SO2_ele_mol'][imonth,:,:]
	lon_RCP_85 = RCP_85_fid.variables['lon'][:]
	Lat_RCP_85 = RCP_85_fid.variables['lat'][:]
	
    # RCP_85 = path_RCP85+'RCP85_2010_SO2_NUM_MOL.0.1x0.1.nc'  #filename month starts from 1
	#stag_tech
	# tech=PATH_stag+'STAG_TECH_OC_2010_'+str(imonth+1)+'_NUM_MOL.0.1x0.1.nc'
	# tech_fid=nc4.Dataset(tech,mode='r')
	# SO2_tech = tech_fid.variables['POM_a1_srf_mol'][:]
	
	EDGAR_FILE =PATH_EDGAR+'REFERENCE_SO2_2010_'+str(imonth+1)+'_NUM_MOL.0.1x0.1.nc'
	# print EDGAR
	EDGAR_fid=nc4.Dataset(EDGAR_FILE,mode='r')
	SO2_EDGAR = EDGAR_fid.variables['SO2_ele_mol'][:] 	
	lon_res = lon_CESM[1] - lon_CESM[0];lat_res = lat_CESM[1] - lat_CESM[0];
	lons,lats = np.meshgrid(lon_CESM,lat_CESM)  #meshgrid
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	rcp_gc_sum = np.multiply(SO2_RCP,area)   # area of the grid cell in m2 * emission/m2/year
	EDGAR_gc_sum = np.multiply(SO2_EDGAR,area)   # area of the grid cell in m2 * emission/m2/year
	RCP = np.nansum(np.nansum(rcp_gc_sum,axis=1),axis=0)
	EDGAR = np.nansum(np.nansum(EDGAR_gc_sum,axis=1),axis=0)
	print imonth, EDGAR, RCP, EDGAR-RCP,round(100*(EDGAR-RCP)/RCP,2)
	plt.subplot(1, 2, 1)
	rcp_gc_sum[rcp_gc_sum==0]=np.nan
	plt.imshow(rcp_gc_sum,cmap='rainbow',origin='lower');
	plt.clim(10**(10),10**(20))
	plt.subplot(1, 2, 2)
	EDGAR_gc_sum[EDGAR_gc_sum==0]=np.nan
	plt.imshow(EDGAR_gc_sum,cmap='rainbow',origin='lower');
	plt.clim(10**(10),10**(20))
	plt.colorbar()
	plt.show()
tech_fid.close()
EDGAR_fid.close()
