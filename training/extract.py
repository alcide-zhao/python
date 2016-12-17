# -*- coding: utf-8 -*-
"""
Created on Mon May 04 19:50:51 2015

this code is to extract NetCDF data out as Tiff, and reproject, clip, and resample the data into certain size:
(1)reproject the soil moisture raster, to ensure sm raster is the same projection with plot polygon.  
(2)clip sm raster to a smaller region, because now is global scale
(3)resample it to a smaller pixel size (25m by default). To successfully extract raster with polygon, polygon should cover more than 1 pixel  

@author: s1251670
"""
import site
import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import pdb

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    'lib'
)

site.addsitedir(lib_path)

import geoTiff
from osgeo import gdal, gdal_array


def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print "NetCDF Global Attributes:"
        for nc_attr in nc_attrs:
            print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
            if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars


########################################################
#0. setting variables
########################################################
#input_path = '/exports/csce/datastore/geos/groups/guasha/Yaqing/WCM/soil_moisture/passive'
#ref_tif = '/exports/csce/datastore/geos/groups/guasha/Yaqing/WCM/soil_moisture/ref_20070601.tif'#here we used a soil moisture tif created from ArcMap to get the information....because it's too difficult to code :-(
#value_layer = 'sm'
input_path = '/exports/csce/datastore/geos/groups/guasha/Yaqing/precipitation'
ref_tif = '/exports/csce/datastore/geos/groups/guasha/Yaqing/precipitation/20070101sample.tif'#here we used a soil moisture tif created from ArcMap to get the information....because it's too difficult to code :-(
value_layer = 'pcp'

proj_flag= True
zone = 'UTM36s'
if zone == 'UTM36s':
    EPSG_code = '32736'
elif zone == 'UTM37s':
    EPSG_code = '32737'

clip_flag= True
outline_shp = '/exports/csce/datastore/geos/groups/guasha/Yaqing/precipitation/district_dissolve.shp'

## read tif into an array which is later used to provide information for .nc to export as .tif
g = gdal.Open(ref_tif)
s0 = gdal_array.DatasetReadAsArray(g)


########################################################
#1. find all .nc files that need to be extracted to .tif
########################################################
os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print > ' + input_path + '/nc_to_tiff.txt')
text_file = open(input_path + '/nc_to_tiff.txt', "r")

    
#######################################################
#2. read in the .nc file, extract soil moisture layer #
#######################################################
for line in text_file: 
    print line
#    name  = os.path.basename(line)[-25:-17]    
#    day   = name[6:8]
#    month = name[4:6]
#    year  = name[0:4]
    name  = os.path.basename(line)[5:13]    
    day   = name[6:8]
    month = name[4:6]
    year  = name[0:4]
    
    nc_f = line[:-1]
    nc_fid = Dataset(nc_f,mode='r')
    nc_attrs, nc_dims, nc_vars = ncdump(nc_fid)
    #print nc_attrs
    #print nc_dims
    #print nc_vars
#    lats = nc_fid.variables['lat'][:]
#    lons = nc_fid.variables['lon'][:]
    lats = nc_fid.variables['latitude'][:]
    lons = nc_fid.variables['longitude'][:]
    sm_all_time = nc_fid.variables[value_layer][:]
 #   sm_all_time = nc_fid.variables['sm'][:]
    #sm = sm_all_time[0,:,:]#only one time slice in this dataset
    sm = sm_all_time[0,::-1,:]#without ::-1, the image is upside down.
  #  sm[sm<0]=np.nan    
    
    #print sm.shape
    
#    sm_un_all_time = nc_fid.variables['sm_uncertainty'][:] 
#    sm_uncertainty = sm_un_all_time[0,:,:]
#    band = nc_fid.variables['freqband'][:] 
#    sensor = nc_fid.variables['sensor'][:]
#    time = nc_fid.variables['t0']
#######################################################
# 3. plot the figure                                  #
#######################################################
#plt.figure(figsize = (100,100))
#plt.imshow(sm, origin = 'lower')
#plt.title(nc_fid.title)
#plt.show()

#fig = plt.figure()
## Setup the map. See http://matplotlib.org/basemap/users/mapsetup.html
## for other projections.
#m = Basemap(llcrnrlon=-35.,llcrnrlat=-30,urcrnrlon=80.,urcrnrlat=50.,resolution='l',area_thresh=1000.,projection='poly',lat_0=0.,lon_0=20.)#just show Africa
#m.drawcoastlines()
#m.drawmapboundary()
## Create 2D lat/lon arrays for Basemap
#lon2d, lat2d = np.meshgrid(lons, lats)
## Transforms lat/lon into plotting coordinates for projection
#x, y = m(lon2d, lat2d)
## Plot of air temperature with 11 contour intervals
#cs = m.contourf(x, y, sm, 11, cmap=plt.cm.Spectral_r)
#cbar = plt.colorbar(cs, orientation='horizontal', shrink=0.5)
##cbar.set_label("%s (%s)" % (nc_fid.variables['sm'].var_desc,nc_fid.variables['sm'].units))
#plt.imshow(sm[0,:,:], origin = 'lower')
#plt.title(nc_fid.title)
#plt.show()
########################################################
#4. export the NetCDF data into .tif                   #
########################################################
    out_tif = os.path.dirname(line) + '/' + name + '.tif'    
    geoTiff.create(out_tif, g, s0, sm, noData = np.nan)
    
#######################################################
#5. reproject ...........this is not working on windows machine, due to some problem on gdal installed
#######################################################
    os.system('gdalwarp -t_srs EPSG:'+EPSG_code+' -r bilinear -dstnodata None -of GTiff -overwrite ' + out_tif + ' ' + os.path.dirname(line) + '/' + name + zone+'_proj.tif') 


#####clip tif with outline of a shapefile    
    out_clip_tif = os.path.dirname(line) + '/' + name + '_clip.tif'   
    print 'cutting  by outline of ' + str(outline_shp)
    os.system('gdalwarp -cutline '+ outline_shp +' -crop_to_cutline '+ '-overwrite ' +os.path.dirname(line) + '/' + name + 'wgs84_36s_proj.tif' + ' '+ out_clip_tif)
    print 'saving to....' + out_clip_tif
    

    nc_fid.close()