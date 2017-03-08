import numpy as np 
from mpl_toolkits.basemap import Basemap
from lib import discrete_cmap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib as mpl


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
		
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r
		
def spatial_figure(data,lons,lats,colormap,colorbar_min,colorbar_max,p_value):
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	
	# calculate the origin of the map
	lon_0 = lons.mean()
	lat_0 = lats.mean()
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = int((lon_e -lon_b)/3)
	lat_bin = int((lat_e -lat_b)/3)
	
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, data)
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='x',color='k',zorder=10)
	
	# Add Grid Lines
	map.drawparallels(np.arange(round(lat_b,1), round(lat_e,1), lat_bin), labels=[1,0,0,0],linewidth=0.0, fontsize=10)
	map.drawmeridians(np.arange(round(lon_b,1), round(lon_e,1), lon_bin), labels=[0,0,0,1],linewidth=0.0, fontsize=10)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines()
	# map.drawstates()
	map.drawcountries()
	# Add Colorbar

	# cmap.set_over((1., 0., 0.))
	# cmap.set_under((0., 0., 1.))
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max)
	cmap= plt.cm.get_cmap(colormap)
	# cmap= reverse_colourmap(cmap)
	cmap.set_bad('w',alpha = 1.0)
	cmap.set_over('k')
	cmap.set_under('w')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max) #
	return colormesh

def spatial_figure_norm(data,lons,lats,colormap,colorbar_min,colorbar_max,p_value):
	"""
	This functiion is designed specifially for spatial divergent fields with 0 as the center.
	"""
	
	# calculate the origin of the map
	lon_0 = lons.mean()
	lat_0 = lats.mean()
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = int((lon_e -lon_b)/3)
	lat_bin = int((lat_e -lat_b)/3)
	
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, data)
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='x',color='k',zorder=10)
	
	# Add Grid Lines
	map.drawparallels(np.arange(round(lat_b,1), round(lat_e,1), lat_bin), labels=[1,0,0,0],linewidth=0.0, fontsize=10)
	map.drawmeridians(np.arange(round(lon_b,1), round(lon_e,1), lon_bin), labels=[0,0,0,1],linewidth=0.0, fontsize=10)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines()
	# map.drawstates()
	map.drawcountries()
	# Add Colorbar

	data[np.isnan(data)] = 0.
	# masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = np.ma.masked_where(np.isnan(data), data)
	norm = MidpointNormalize(midpoint=0.,vmin = colorbar_min, vmax = colorbar_max)
	# cmap = discrete_cmap(50,colormap)
	cmap= plt.cm.get_cmap(colormap)
	# cmap= reverse_colourmap(cmap)
	cmap.set_bad('w',alpha = 1.0);cmap.set_over('k'); cmap.set_under('w')
	colormesh = map.pcolormesh(xi, yi, data,cmap=cmap,norm=norm,vmin=colorbar_min,vmax=colorbar_max) #
	return colormesh
# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/'

def spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,scaler,vector_u,vector_v,projection='cyl'):
	# calculate the origin of the map
	lon_0 = lons.mean();lat_0 = lats.mean()
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = int((lon_e -lon_b)/3); lat_bin = int((lat_e -lat_b)/3)
	
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# plot SLP contours.
	# set desired contour levels.
	clevs = np.arange(np.nanmin(scaler[:]),np.nanmax(scaler[:]),(np.nanmax(scaler[:])-np.nanmin(scaler[:]))/20)
	cmap = discrete_cmap(20,colormap)
	# cmap= reverse_colourmap(cmap)
	cmap.set_bad('w',alpha = 1.0); cmap.set_over('k'); cmap.set_under('w')
	masked_obj = np.ma.masked_where(np.isnan(scaler), scaler)
	# norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max)
	CS1 = map.contourf(xi,yi,masked_obj,clevs,linewidths=0.5,colors='k',animated=True)
	CS2 = map.contourf(xi,yi,masked_obj,clevs,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)
	# # plot wind vectors on projection grid.
	# # first, shift grid so it goes from -180 to 180 (instead of 0 to 360
	# # in longitude).  Otherwise, interpolation is messed up.
	# ugrid,newlons = shiftgrid(180.,vector_u,lons,start=False)
	# vgrid,newlons = shiftgrid(180.,vector_v,lons,start=False)
	# # transform vectors to projection grid.
	# uproj,vproj,xx,yy = \
	# map.transform_vector(ugrid,vgrid,newlons,lats,31,31,returnxy=True,masked=True)
	# now plot.
	# urot,vrot,x,y = map.rotate_vector(vector_u,vector_v,lons,lats],returnxy=True)
	Q = map.quiver(xi,yi,vector_u,vector_v,scale=None)
	# make quiver key.
	qk = plt.quiverkey(Q, 0.9, 1.1, 2, '2m/s', labelpos='W') #2 m/s
	
	# Add Grid Lines
	map.drawparallels(np.arange(round(lat_b,1), round(lat_e,1), lat_bin), labels=[1,0,0,0],linewidth=0.0, fontsize=10)
	map.drawmeridians(np.arange(round(lon_b,1), round(lon_e,1), lon_bin), labels=[0,0,0,1],linewidth=0.0, fontsize=10)
	map.drawcoastlines(); map.drawcountries()
	return CS2