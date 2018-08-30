"""
This is to plot the T and P extreme indices calculated from the EDGAR experiments
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


def spatial_figure(axs,data,p_value,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	lons[lons>180]-=360; 
	# calculate the origin of the map
	lon_0 = lons.mean(); 
	lat_0 = lats.mean(); 
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# s = map.pcolor(xi, yi, data)
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines(); #map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	ax.set_ylim([-60,90])
	return colormesh
		
def data_readin(variable,scale):	
	def data_netcdf(scenario,scale,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/analysis/'+scale+'/'
		var_path = input_path+scenario+'.TempPrecp.extremes.'+scale+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		data = nc_fid.variables[variable][0:30,:,:]
		nc_fid.close()
		return lon,lat,data
	
	lon,lat,Edg70GO = data_netcdf('Edg70GO',scale,variable)
	_,_,T1970 = data_netcdf('T1970RCP',scale,variable)
	_,_,EdgRef = data_netcdf('EdgRef',scale,variable)
	_,_,Edg70Oz = data_netcdf('Edg70Oz',scale,variable)
	_,_,Edg70T10SOZ = data_netcdf('Edg70T10SOZ',scale,variable)
	return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ
	
def print_domain_mean(variable,scale):
	print "-------------------------------" + variable+"-------------------------------"
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""
	def AreaWeight(lon1,lon2,lat1,lat2):
		'''
		calculate the earth radius in m2
		'''
		radius = 6371000;
		area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
		(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
		# print np.nansum(np.nansum(area,axis=1),axis=0)
		return area
		
	def box_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,mask):
		"""
		fill the range outside the box with 0
		"""
		lon = np.array(lon)
		lat = np.array(lat)
		colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
		colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
		row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
		row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
		if (colum_s> colum_e):
			cache = colum_e; colum_e = colum_s; colum_s = cache;
		if (row_s> row_e):
			cache = row_e; row_e = row_s; row_s = cache;
		mask[:,0:colum_s] =0; mask[:,colum_e:-1] =0
		# plt.imshow(mask,origin='lower');plt.show()
		mask[0:row_s,:] =0; mask[row_e:-1,:] =0
		# plt.imshow(mask,origin='lower');plt.show()
		return mask

	def mask_weight(region_key,lon,lat,reverse=False):
		"""
		Read in the country mask
		interpolate it to the required resolution grids with lon_interp,lat_interp 
		
		"""
		lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
		lons,lats = np.meshgrid(lon,lat)
		area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)

		##OCEAN_MASKS FOR COUNTRIES
		ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')
		lon_mask = ocean_mask['lon'][0,:];
		lat_mask = ocean_mask['lat'][0,:];
		box_region_dic={'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
		if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'GloLand'):
			mask= ocean_mask[region_key][:]
		elif  region_key in box_region_dic:
			mask= ocean_mask['All'][:]
			box = box_region_dic[region_key]
			mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
		else:
			print "error region name"
		
		# interpolate from 360*720 to 192*288
		mask[np.isnan(mask)]=0;	mask[mask>0]=1;
		f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon, lat);
		# plt.imshow(mask,origin='lower');plt.show()
		mask[mask >= 1] = 1;mask[mask < 1] = 0;
		# weight each grid cell by its area weight against the total area
		if reverse:
			mask=1-mask
		mask[mask==0] = np.nan
		mask=np.multiply(mask,area); 
		mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		return mask_weighted
		
	oceanmask=sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
	oceanmask[oceanmask==0]=np.nan;
	####calculate the difference and their significance
	def diff_sig_T(vairable1,variable2):
		def mannwhitneyu_test(vairable1,variable2):
			p_threshold=0.05
			size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
			p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
			from scipy.stats import mannwhitneyu as test
			for x in range(size[0]):
				for y in range(size[1]):
					cache1 = vairable1[:,x,y]
					cache2 = variable2[:,x,y]
					if np.array_equal(cache1,cache2): p_value[x,y] = np.nan;
					else: _,p_value[x,y] = test(cache1,cache2);
			p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
			p_value = np.multiply(oceanmask,p_value)
			# plt.imshow(p_value);plt.show()
			return p_value

		dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
		sig = mannwhitneyu_test(vairable1, variable2)
		# sig_dif = np.multiply(dif,sig)
		return dif,sig
		
	def diff_sig_P(vairable1,variable2,ref):
		def mannwhitneyu_test(vairable1,variable2):
			p_threshold=0.05
			size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
			p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
			from scipy.stats import mannwhitneyu as test
			for x in range(size[0]):
				for y in range(size[1]):
					cache1 = vairable1[:,x,y]
					cache2 = variable2[:,x,y]
					if np.array_equal(cache1,cache2): p_value[x,y] = np.nan;
					else: _,p_value[x,y] = test(cache1,cache2);
			p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
			p_value = np.multiply(oceanmask,p_value)
			return p_value

		dif = np.divide(np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0),np.nanmean(ref,axis=0))*100
		sig = mannwhitneyu_test(vairable1, variable2)
		return dif,sig

	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ = data_readin(variable,scale);
	
	if INDEX in ['TXX','TNN','TSM','TSSTD','CDD','R95','R99','CWD']:
		Total,Total_s = diff_sig_T(EdgRef,T1970)
		GHG,GHG_s = diff_sig_T(Edg70Oz,Edg70GO)
		AEROSOL,AEROSOL_s = diff_sig_T(Edg70GO,T1970)
		TrO3,TrO3_s = diff_sig_T(EdgRef, Edg70T10SOZ)
		StO3,StO3_s = diff_sig_T(Edg70T10SOZ, Edg70Oz)
	elif INDEX in ['RX5DAY','SDII','R10','RX1DAY']: 
		Total,Total_s = diff_sig_P(EdgRef,T1970,T1970)
		GHG,GHG_s = diff_sig_P(Edg70Oz,Edg70GO,T1970)
		AEROSOL,AEROSOL_s = diff_sig_P(Edg70GO,T1970,T1970)
		TrO3,TrO3_s = diff_sig_P(EdgRef, Edg70T10SOZ,T1970)
		StO3,StO3_s = diff_sig_P(Edg70T10SOZ, Edg70Oz,T1970)
	
	region_keys =['Globe','EA','India','Europe','USA','SESA']
	for region_key in region_keys:
		# print region_key
		mask = mask_weight(region_key,lon,lat); #print mask
		print "=========="+region_key+"=========="
		print 'Total2010 ',round(np.nansum(np.nansum(np.multiply(mask,Total),axis=1),axis=0),4)
		print 'GHG       ',round(np.nansum(np.nansum(np.multiply(mask,GHG),axis=1),axis=0),4)
		print 'AEROSOL   ',round(np.nansum(np.nansum(np.multiply(mask,AEROSOL),axis=1),axis=0),4)	
		print 'Trop.O3   ',round(np.nansum(np.nansum(np.multiply(mask,TrO3),axis=1),axis=0),4)
		print 'Stra.O3   ',round(np.nansum(np.nansum(np.multiply(mask,StO3),axis=1),axis=0),4)	
	return lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s
	
	
INDEX = 'R99'; AorS='annual'  #TXX TNN CDD RX5DAY R10 R99

lon,lat,Total,GHG,AEROSOL,TrO3,StO3,Total_s,GHG_s,AEROSOL_s,TrO3_s,StO3_s= print_domain_mean(variable=INDEX,scale=AorS);

if INDEX in ['TXX','TNN','TSM','TSSTD']:
	colormap='RdBu';colormap = reverse_colourmap(colormap);
	colorbar_min=-2;colorbar_max=2;
elif INDEX in ['RX5DAY','R10','CWD','CDD','SDII','RX1DAY','R95','R99']:
	colormap='BrBG';
	if INDEX == 'CDD': colormap = reverse_colourmap(colormap);
	if INDEX in ['RX5DAY','SDII','RX1DAY','R10','R20']:
		colorbar_min=-25;colorbar_max=25;
	elif INDEX in ['R95','R99']:
		colorbar_min=-2;colorbar_max=2;
	else:
		colorbar_min=-5;colorbar_max=5;
else:
	print "make sure of a right extreme index"

fig = plt.figure(facecolor='White',figsize=[9.5,7.5]);pad= 5
# colormap='RdBu_r';  #colormap = reverse_colourmap(colormap);
ax = plt.subplot(3,2,1);
ax.annotate('(a) Total',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
spatial_figure(ax,Total,Total_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,2);
ax.annotate('(b) GHGs',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
colormesh0 = spatial_figure(ax,GHG,GHG_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)


ax = plt.subplot(3,2,3);
ax.annotate('(c) AAs',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,4);
ax.annotate('(d) Trop. O3',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
colormesh1=spatial_figure(ax,TrO3,TrO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(3,2,5);
ax.annotate('(e) Strat. O3',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
colormesh1=spatial_figure(ax,StO3,StO3_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.10, 0.03, 0.75, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.02, bottom=0.08, right=0.98, top=0.95, wspace=0.04, hspace=0.15); 
plt.savefig(INDEX+'.png', format='png', dpi=1000)
