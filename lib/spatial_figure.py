import numpy as np 
from mpl_toolkits.basemap import Basemap, maskoceans
from colormap_modify import discrete_cmap,reverse_colourmap
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

		
def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,p_value,tb_lef=True,tb_bot=True ): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	# lons[lons>180]-=360; 
	# calculate the origin of the map
	lon_0 = lons.mean(); 
	lat_0 = lats.mean(); 
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	# lon_bin = int((lon_e -lon_b)/5)
	# lat_bin = int((lat_e -lat_b)/5)
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	s = map.pcolor(xi, yi, data)
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='.',color='k',zorder=10)
	
	# Add Grid Lines
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines();# map.drawstates(); map.drawcountries()
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); cmap.set_under('b')
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min-0.01, vmax=colorbar_max) #,latlon=True
	return colormesh

def spatial_figure_norm(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=True ):
	"""
	This functiion is designed specifially for spatial divergent fields with 0 as the center.
	"""
	# calculate the origin of the map
	lons[lons>180]-=360; 
	lon_0 = lons.mean()
	lat_0 = lats.mean()
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	# lon_bin = int((lon_e -lon_b)/5);lat_bin = int((lat_e -lat_b)/5)
	lon_bin = 20; lat_bin = 15
	
	
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# Add Grid Lines
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines()
	# map.drawstates()
	map.drawcountries()

	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(50,colormap)
	cmap.set_bad('w');cmap.set_over('r'); cmap.set_under('w');
	norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max) #
	# clevs = np.arange(np.nanmin(data[:]),np.nanmax(data[:]),(np.nanmax(data[:])-np.nanmin(data[:]))/20)
	# CS1 = map.contourf(xi,yi,data,clevs,linewidths=0.5,colors='k',animated=True)
	# colormesh = map.contourf(xi,yi,data,clevs,cmap=cmap,vmin=colorbar_min,norm=norm,vmax=colorbar_max) #
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min,norm=norm,vmax=colorbar_max,latlon=True) #
	sc = map.scatter(xi[p_value==1], yi[p_value==1],p_value[p_value==1],marker='.',color='k',zorder=10)
	return colormesh

def spatial_scaler_vector(axs,lons,lats,colormap,colorbar_min,colorbar_max,scaler,vector_u,vector_v,qk_scale=None,qk_caption=None,qk_is = False,tb_lef=True,tb_bot=True ):
	# calculate the origin of the map
	lon_0 = lons.mean();lat_0 = lats.mean()
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 20; lat_bin = 15
	# lon_bins = np.arange(60.0,150.0,20.0); lat_bins = np.arange(-10.0,60.0,20.0)
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# plot contours.
	# set desired contour levels.
	# clevs = np.arange(np.nanmin(scaler[:]),np.nanmax(scaler[:]),(np.nanmax(scaler[:])-np.nanmin(scaler[:]))/20)
	cmap = discrete_cmap(10,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over([0,0,0]); cmap.set_under('b')
	masked_obj = np.ma.masked_where(np.isnan(scaler), scaler)
	norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max)
	# CS1 = map.contourf(xi,yi,masked_obj,linewidths=0.5,colors='k',animated=True)
	# CS2 = map.contourf(xi,yi,masked_obj,cmap=cmap,vmin=colorbar_min,vmax=colorbar_max,norm=norm) #
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)  #norm=norm,
	Q = map.quiver(xi,yi,vector_u,vector_v,scale=None)
	def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax):
		xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
		ys = [latmin,latmin,latmax,latmax,latmin]
		bmap.plot(xs, ys,latlon = True,linewidth = 1,color='k')
	# plot_rectangle(map, 100,145,20,55); x,y = map(102,48);plt.text(x,y,'East Asia',color='k')
	# plot_rectangle(map, 65,100,5,30); x,y = map(67,7);plt.text(x,y,'South Asia',color='k')
	
	# make quiver key.
	if qk_is:
		qk = plt.quiverkey(Q,1.11,0.05,qk_scale, qk_caption, labelpos='E') #10 m/s ,coordinates='data'   1.12,0.1 1.16,0.05  0.80,-0.20
		qk.text.set_backgroundcolor('w');qk.text.set_fontsize(10);
		qk.text.set_horizontalalignment('left')
	# Add Grid Lines
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0), round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	else:
		map.drawparallels(np.arange(round(lat_b,0), round(lat_e,0), lat_bin), labels=[0,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	else:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,0],linewidth=0.0,fontsize=8)
	map.drawcoastlines(); map.drawcountries()
	return colormesh