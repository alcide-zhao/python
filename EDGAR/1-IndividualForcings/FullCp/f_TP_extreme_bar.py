"""
This is to plot the TS shift using bars: mean and std of the effects from GHG,aerosol, ozone are splitted
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
	
def domain_stats(variable,scale):
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted) and std
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

	def mask_weight(region_key,lon,lat):
		"""
		Read in the country mask
		interpolate it to the required resolution grids with lon_interp,lat_interp 
		
		"""
		lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
		lons,lats = np.meshgrid (lon,lat)
		area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
		##OCEAN_MASKS FOR COUNTRIES
		ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
		lon_mask = ocean_mask['lon'][0,:];
		lat_mask = ocean_mask['lat'][0,:];
		box_region_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,60],'EA':[105,145,15,35],'SA':[65,100,5,30],'SESA':[295,315,-40,-25]}
		if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'Globe'):
			mask= ocean_mask[region_key][:]
		elif (region_key == 'ASIA' or region_key == 'EA' or region_key == 'SA' or region_key == 'SESA'):
			mask= ocean_mask['Globe'][:]
			box = box_region_dic[region_key]
			mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
		else:
			print "error region name"
		# interpolate from 360*720 to 192*288
		mask[np.isnan(mask)]=0;	mask[mask>0]=1;
		f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
		mask = f(lon, lat);
		mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
		# weight each grid cell by its area weight against the total area
		mask=np.multiply(mask,area);  
		mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		# print np.nansum(np.nansum(mask_weighted,axis=1),axis=0)
		return mask_weighted	
	
	def dif_std_for_region(var1,var2,region_key):
		"""
		for a rgional sif and std. the first step is to evaluate the area-weighted mean, and then the differences and average of difference
		"""
		mask = mask_weight(region_key,lon,lat);
		var1_domain_mean = np.nansum(np.nansum(np.multiply(mask,var1),axis=2),axis=1)
		var2_domain_mean =np.nansum(np.nansum(np.multiply(mask,var2),axis=2),axis=1)
		dif =  np.nanmean(var1_domain_mean - var2_domain_mean,axis=0)
		p25 =  dif - np.nanpercentile(var1_domain_mean - var2_domain_mean,25)
		p75 =  np.nanpercentile(var1_domain_mean - var2_domain_mean,75)- dif
		dif_p25_p75 =np.array([dif,p25/2,p75/2])
		return dif_p25_p75

	lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,Edg70T10SOZ = data_readin(variable,scale);	
		
	# region_stats ={'Globe':np.empty((5,3)),'EA':np.empty((5,3)),'India':np.empty((5,3)),'Europe':np.empty((5,3)),'USA':np.empty((5,3)),'SESA':np.empty((5,3))}
	region_stats ={'ASIA':np.empty((5,3)),'USA':np.empty((5,3)),'Europe':np.empty((5,3)),'SESA':np.empty((5,3))}
	for region_key in region_stats:
		dif_p25_p75_total = dif_std_for_region(EdgRef,T1970,region_key)
		dif_p25_p75_ghg = dif_std_for_region(Edg70Oz,Edg70GO,region_key)
		dif_p25_p75_aerosol = dif_std_for_region(Edg70GO,T1970,region_key)		
		dif_p25_p75_TrO3 = dif_std_for_region(EdgRef, Edg70T10SOZ,region_key)
		dif_p25_p75_StO3 = dif_std_for_region(Edg70T10SOZ, Edg70Oz,region_key)		
		region_stats[region_key][:] = np.array([dif_p25_p75_total,dif_p25_p75_ghg,dif_p25_p75_aerosol,dif_p25_p75_TrO3,dif_p25_p75_StO3]) 
	return region_stats
	
def plot_bar_err(var):
	print "------"+var+"--------"
	def autolabel(X_pos,values,height_lift):
		"""
		Attach a text label above each bar displaying its height
		"""
		height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
		for i in range(len(height)):
			ax.text(X_pos[i],y_pos[i],'%4.2f' % height[i], ha='center', va='bottom',size=8)
	
	var_stats = domain_stats(variable=var,scale='annual');
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
	# ax.axvline(x=3,color='k',linewidth = 2);
	ax.axhline(y=0,color='k',linewidth = 2);
	ax.set_xlim([-0.5,12]);ax.set_xticks(np.arange(0.75,12.5,2.5));ax.set_xticklabels(('Total','GHGs','AAs','Trop. O3','Str. O3'),fontsize=15);
	# ax.set_ylim([-1,2.0]);ax.set_yticks(np.arange(-1,2.01,1))
	ax.axvspan(-0.5, 2, alpha=0.1, color='r',edgecolor='none');
	ax.axvspan(4.5, 7, alpha=0.1, color='y',edgecolor='none');
	ax.axvspan(9.5, 12, alpha=0.1, color='g',edgecolor='none');
	
	X_pos=np.arange(0,12.5,2.5);dis=0.5
	
	# region_key = 'Globe'
	# y_value = var_stats[region_key][:,0];yerr=abs(np.array([var_stats[region_key][:,1],var_stats[region_key][:,2]])); 
	# rects1=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='r', ecolor='r',capsize=0,alpha=0.7,width=0.5,lw=0)
	# autolabel(X_pos,y_value,1.1)
	# print var_stats[region_key]

	region_key = 'ASIA'
	y_value = var_stats[region_key][:,0];yerr=abs(np.array([var_stats[region_key][:,1],var_stats[region_key][:,2]]));  
	rects1=ax.bar(X_pos+dis*0, y_value, yerr=yerr, align='center',color='r', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*0,y_value,1.1)
	
	region_key = 'Europe'
	y_value = var_stats[region_key][:,0];yerr=abs(np.array([var_stats[region_key][:,1],var_stats[region_key][:,2]])); 
	rects2=ax.bar(X_pos+dis*1, y_value, yerr=yerr, align='center',color='orange', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*1,y_value,1.1)

	# region_key = 'USA'
	# y_value = var_stats[region_key][:,0];yerr=abs(np.array([var_stats[region_key][:,1],var_stats[region_key][:,2]])); 
	# rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr, align='center',color='y', ecolor='y',capsize=0,alpha=0.7,width=0.5,lw=0)
	# autolabel(X_pos+dis*2,y_value,1.1)

	region_key = 'USA'
	y_value = var_stats[region_key][:,0];yerr=abs(np.array([var_stats[region_key][:,1],var_stats[region_key][:,2]]));  
	rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr, align='center',color='g', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*2,y_value,1.1)
	
	region_key = 'SESA'
	y_value = var_stats[region_key][:,0];yerr=abs(np.array([var_stats[region_key][:,1],var_stats[region_key][:,2]]));  
	rects4=ax.bar(X_pos+dis*3, y_value, yerr=yerr, align='center',color='m', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
	autolabel(X_pos+dis*3,y_value,1.1)
	return rects1,rects2,rects3,rects4


fig = plt.figure(facecolor='White',figsize=[16,9]);pad= 5;
colormap='RdBu';colormap = reverse_colourmap(colormap);
"""
ax = plt.subplot(3,2,1);
ax.annotate('(a) SAT (K) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('TSM') #   #,rects5[0],rects6[0]
legend1=ax.legend((rects1[0], rects2[0],rects3[0], rects4[0]),\
(  'Asia', 'Europe','USA', 'SES. America'), shadow=False,ncol=2,bbox_to_anchor=(0.9, 0.95))  #'Globe','USA',
legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)
"""

ax = plt.subplot(2,2,1);
ax.annotate('(a) TXX (K) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('TXX') #   #,rects5[0],rects6[0]
legend1=ax.legend((rects1[0], rects2[0],rects3[0], rects4[0]),\
(  'Asia', 'Europe','USA', 'SES. America'), shadow=False,ncol=2,bbox_to_anchor=(0.9, 0.95))  #'Globe','USA',
legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)

ax = plt.subplot(2,2,2);
ax.annotate('(b) TNN (K) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('TNN')
"""
ax = plt.subplot(3,2,2);
ax.annotate('(d) Precip (mm/day) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('PRECPM')
"""
ax = plt.subplot(2,2,3);
ax.annotate('(c) CDD (days) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('CDD')

ax = plt.subplot(2,2,4);
ax.annotate('(d) RX5DAY (mm) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('RX5DAY')

"""
ax = plt.subplot(3,2,5);
ax.annotate('(e) R99 (mm/day) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
rects1,rects2,rects3,rects4 = plot_bar_err('R99')


ax = plt.subplot(3,2,6);
ax.annotate('(f) R10P (mm) ',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
rects1,rects2,rects3,rects4 = plot_bar_err('R10')
"""
plt.subplots_adjust(left=0.05, bottom=0.08, right=0.98, top=0.95, wspace=0.20, hspace=0.30); 
plt.savefig('P_T_EXTREMES_fraction_bar.png', format='png', dpi=1000)
