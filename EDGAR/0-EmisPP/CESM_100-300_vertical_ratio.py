import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import time as clock
import glob  
import matplotlib.pyplot as plt

file= '/exports/csce/datastore/geos/users/s1667168/EDGAR/CESM_rcp85/RCP85_mam3_so2_elev_2000-2100_c20110913.nc'
nc_fid = nc4.Dataset(file,mode='r')
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:] 
BC_mass = nc_fid.variables['emiss_ene'][0,:,:,:]
BC_l0=BC_mass[0,:,:]*126
BC_l1=BC_mass[1,:,:]*152
BC_l2=BC_mass[2,:,:]*176
sum= BC_l0+BC_l1+BC_l2
r00=100*np.divide(BC_l0,sum)
r10=100*np.divide(BC_l1,sum)#-r_10;r10[r10==0]=np.nan
r21=100*np.divide(BC_l2,sum)#-r_21;r21[r21==0]=np.nan

print stats.nanmean(stats.nanmean(r00,axis=1),axis=0)
print stats.nanmean(stats.nanmean(r10,axis=1),axis=0)
print stats.nanmean(stats.nanmean(r21,axis=1),axis=0)


