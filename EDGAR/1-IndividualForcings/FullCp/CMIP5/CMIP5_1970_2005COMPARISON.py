"""
This is to plot the mean SAT and Precip from the CMIP5 single forcing runs
we are interested in the all-forcing , GHG and aerosol runs
time clice 1965-1975, 
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import matplotlib.pyplot as plt
import math

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

def CMIP5_data_read_pp(var):
	datapath = '/exports/csce/datastore/geos/groups/TITAN/sabine/'
	his_path = datapath+'CMIP5_DATA--hist--'+var+'/*.nc'
	files=sorted(glob.glob(input_path))
	ensumble_no = 0
	for file in files:
		# print file
		nc_fid = nc4.Dataset(file,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		data = nc_fid.variables[var][time,:,:]
		ms_value = nc_fid.variables[var].missing_value
		fi_value = nc_fid.variables[var]._FillValue
		data[data == ms_value]=np.nan;data[data==fi_value]=np.nan;
		ensumble_no=ensumble_no+1
		nc_fid.close()
	att_value[att_value == missing_value] =np.nan
	att_ensumble_mean[time,:,:] = stats.nanmean(att_value,axis = 0)
	AAs_path = datapath+'CMIP5_DATA--histMisc--'+var
	GHG_path =  datapath+'CMIP5_DATA--histGHG--'+var
	
	
	