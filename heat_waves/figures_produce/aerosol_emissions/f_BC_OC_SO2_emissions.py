"""
"""
from scipy import stats
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import site
import scipy.io as sio

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,
	os.path.pardir,
	os.path.pardir, 
)
site.addsitedir(lib_path)
from lib import *



#######country mask  not use sio.whosmat to see the valuables in the mat files
land_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
lon=land_mask['lon'][:];lat=land_mask['lat'][:]
# print lon, lat
# GLO=np.empty((360,720));GLO[:]=1
GLO=land_mask['Globe'][:];
AUS=land_mask['Australia'][:];EUR=land_mask['Europe'][:];CHA=land_mask['China'][:];
USA=land_mask['USA'][:];IND=land_mask['India'][:];BRZ=land_mask['Brazil'][:]; SAF=land_mask['Southern African'][:]; 
 
# plt.imshow(USA,origin='lower');plt.show()

aerosol_path = '/exports/csce/datastore/geos/users/s1667168/emissions/'
aerosol_total={'BC':np.empty((8,96)),'OC':np.empty((8,96)),'SO2':np.empty((8,96))} 
#######################################################
# 0.2 functions and subroutines                       #
#######################################################
def movingaverage (values, window=3):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result
	
def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area

def MonthWeightAnnualSum(annual_data):
	"""
	weight the days to each month and sum the mass up throughout the year
	"""
	calendar_day = {'1':31,'2':28,'3':31,'4':30,'5':31,'6':30,'7':31,'8':31,'9':30,'10':31,'11':30,'12':31}
	month_sum= np.empty((np.shape(annual_data)))
	for month in range(0,12):
		day_no = calendar_day[str(month+1)]
		month_sum[month,:,:] =annual_data[month,:,:]*day_no*24*60*60/(10**(9))
	annual_sum = np.nansum(month_sum,axis=0)
	return annual_sum
	
def get_AOD():
	AOD_PATH = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_AODVIS_200602_210101.nc'
	nc_fid = nc4.Dataset(AOD_PATH,mode='r');
	AOD_diff = nc_fid.variables['rcp85'][:]-nc_fid.variables['rcp85_fixA'][:];
	AOD=np.empty((95,192,288)); AOD[:]=np.nan
	for iyear in range(0,95):
		AOD[iyear,:,:] = stats.nanmean(AOD_diff[iyear*12:(iyear+1)*12,:,:],axis=0)
	lon = nc_fid.variables['lon'][:];lat = nc_fid.variables['lat'][:];
	# nc_fid.close()
	return lon,lat,AOD

def spatial_figure(axs,data,lons,lats,colormap,colorbar_min,colorbar_max,p_value,tb_lef=True,tb_bot=True ): #c_bad,c_under,c_over,c_number=20,
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
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(10,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_under('k') #cmap.set_over('r'); 
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min-0.0000001, vmax=colorbar_max,latlon=True) #
	return colormesh
	
for species in aerosol_total:
	file_name = aerosol_path+'accmip_interpolated_emissions_RCP85_'+species+'_2005_2100_0.5x0.5.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	data = nc_fid.variables['anthropogenic'][:] +nc_fid.variables['ships'][:]+nc_fid.variables['biomass_burning'][:]
	lon = nc_fid.variables['lon'];lat = nc_fid.variables['lat'];
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid (lon,lat)
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	nc_fid.close()
	annual_weighted=np.empty((96,360,720)); annual_weighted[:]=np.nan
	store=np.empty((8,96));store[:]=np.nan
	for iyear in range(0,96):
		annual_weighted[iyear,:,:] = np.multiply(MonthWeightAnnualSum(data[iyear*12:(iyear+1)*12,:,:]),area)
	store[0,:] = np.nansum(np.nansum(annual_weighted,axis=2),axis=1)
	store[1,:] = np.nansum(np.nansum(np.multiply(annual_weighted,IND),axis=2),axis=1) 	
	store[2,:] = np.nansum(np.nansum(np.multiply(annual_weighted,CHA),axis=2),axis=1) 
	store[3,:] = np.nansum(np.nansum(np.multiply(annual_weighted,EUR),axis=2),axis=1) 
	store[4,:] = np.nansum(np.nansum(np.multiply(annual_weighted,USA),axis=2),axis=1) 
	store[5,:] = np.nansum(np.nansum(np.multiply(annual_weighted,BRZ),axis=2),axis=1) 	
	store[6,:] = np.nansum(np.nansum(np.multiply(annual_weighted,AUS),axis=2),axis=1)
	store[7,:] = np.nansum(np.nansum(np.multiply(annual_weighted,SAF),axis=2),axis=1)

	aerosol_total[species.upper()][:] = store


	
year = range(2005,2101)	;
fig = plt.figure(facecolor='White',figsize=[11,6.5]);plot_setup();pad= 5;
# ax1 = plt.subplot(2,2,1);
ax1 =plt.subplot2grid((2, 4), (0, 0), colspan=1)
BC=aerosol_total['BC'][:]; 
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_xlim([2005,2100]);
ax1.set_ylim([0,8]);ax1.set_yticks(np.arange(0,8.1,2))
ax1.annotate(r'(a) $\mathrm{BC\/(Tg\/year^{-1})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal')		
ax1.plot(year,BC[0,:],'-',color="k",linewidth=0.3)
ax1.plot(year,BC[1,:],'-',color="k",linewidth=0.3)
ax1.plot(year,np.nansum(BC[1:3,:],axis=0),'-',color="k",linewidth=0.3)
ax1.plot(year,np.nansum(BC[1:4,:],axis=0),'-',color="k",linewidth=0.3)
ax1.plot(year,np.nansum(BC[1:5,:],axis=0),'-',color="k",linewidth=0.3)
ax1.plot(year,np.nansum(BC[1:6,:],axis=0),'-',color="k",linewidth=0.3)
ax1.plot(year,np.nansum(BC[1:7,:],axis=0),'-',color="k",linewidth=0.3)
ax1.plot(year,np.nansum(BC[1:8,:],axis=0),'-',color="k",linewidth=0.3)
ax1.fill_between(year,BC[1,:],0, where=BC[1,:]>=0,facecolor='m',label='India',alpha=0.7)
ax1.fill_between(year,np.nansum(BC[1:3,:],axis=0),np.nansum(BC[1:2,:],axis=0),facecolor="r",label='China',alpha=0.7)
ax1.fill_between(year,np.nansum(BC[1:4,:],axis=0),np.nansum(BC[1:3,:],axis=0),facecolor="y",label='Europe',alpha=0.7)
ax1.fill_between(year,np.nansum(BC[1:5,:],axis=0),np.nansum(BC[1:4,:],axis=0),facecolor="g",label='USA',alpha=0.7)
ax1.fill_between(year,np.nansum(BC[1:6,:],axis=0),np.nansum(BC[1:5,:],axis=0),facecolor="c",label='Brazil',alpha=0.7)
ax1.fill_between(year,np.nansum(BC[1:7,:],axis=0),np.nansum(BC[1:6,:],axis=0),facecolor="b",label='Australia',alpha=0.7)
ax1.fill_between(year,np.nansum(BC[1:8,:],axis=0),np.nansum(BC[1:7,:],axis=0),facecolor="orange",label='Southern Africa',alpha=0.7)
ax1.fill_between(year,BC[0,:],np.nansum(BC[1:8,:],axis=0),facecolor="wheat",label='Others',alpha=0.7)


ax2 =plt.subplot2grid((2, 4), (0, 1), colspan=1)
OC=aerosol_total['OC'][:]; 
ax2.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax2.set_xlim([2005,2100]);
ax2.set_ylim([0,36]);ax2.set_yticks(np.arange(0,36.1,9))
ax2.annotate(r'(b) $\mathrm{OC\/(Tg\/year^{-1})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal')		
ax2.plot(year,OC[0,:],'-',color="k",linewidth=0.3)
ax2.plot(year,OC[1,:],'-',color="k",linewidth=0.3)
ax2.plot(year,np.nansum(OC[1:3,:],axis=0),'-',color="k",linewidth=0.3)
ax2.plot(year,np.nansum(OC[1:4,:],axis=0),'-',color="k",linewidth=0.3)
ax2.plot(year,np.nansum(OC[1:5,:],axis=0),'-',color="k",linewidth=0.3)
ax2.plot(year,np.nansum(OC[1:6,:],axis=0),'-',color="k",linewidth=0.3)
ax2.plot(year,np.nansum(OC[1:7,:],axis=0),'-',color="k",linewidth=0.3)
ax2.plot(year,np.nansum(OC[1:8,:],axis=0),'-',color="k",linewidth=0.3)
ax2.fill_between(year,OC[1,:],0, where=OC[1,:]>=0,facecolor='m',label='India',alpha=0.7)
ax2.fill_between(year,np.nansum(OC[1:3,:],axis=0),np.nansum(OC[1:2,:],axis=0),facecolor="r",label='China',alpha=0.7)
ax2.fill_between(year,np.nansum(OC[1:4,:],axis=0),np.nansum(OC[1:3,:],axis=0),facecolor="y",label='Europe',alpha=0.7)
ax2.fill_between(year,np.nansum(OC[1:5,:],axis=0),np.nansum(OC[1:4,:],axis=0),facecolor="g",label='USA',alpha=0.7)
ax2.fill_between(year,np.nansum(OC[1:6,:],axis=0),np.nansum(OC[1:5,:],axis=0),facecolor="c",label='Brazil',alpha=0.7)
ax2.fill_between(year,np.nansum(OC[1:7,:],axis=0),np.nansum(OC[1:6,:],axis=0),facecolor="b",label='Australia',alpha=0.7)
ax2.fill_between(year,np.nansum(OC[1:8,:],axis=0),np.nansum(OC[1:7,:],axis=0),facecolor="orange",label='Southern Africa',alpha=0.7)
ax2.fill_between(year,OC[0,:],np.nansum(OC[1:8,:],axis=0),facecolor="wheat",label='Others',alpha=0.7)

ax3 =plt.subplot2grid((2, 4), (0, 2), colspan=1)
SO2=aerosol_total['SO2'][:]; print SO2[0,:]
ax3.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax3.set_xlim([2005,2100]);
ax3.set_ylim([0,120]);ax3.set_yticks(np.arange(0,120.1,30))
ax3.annotate(r'(c) $\mathrm{SO2\/(Tg\/year^{-1})}$',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal')		
ax3.plot(year,SO2[0,:],'-',color="k",linewidth=0.3)
ax3.plot(year,SO2[1,:],'-',color="k",linewidth=0.3)
ax3.plot(year,np.nansum(SO2[1:3,:],axis=0),'-',color="k",linewidth=0.3)
ax3.plot(year,np.nansum(SO2[1:4,:],axis=0),'-',color="k",linewidth=0.3)
ax3.plot(year,np.nansum(SO2[1:5,:],axis=0),'-',color="k",linewidth=0.3)
ax3.plot(year,np.nansum(SO2[1:6,:],axis=0),'-',color="k",linewidth=0.3)
ax3.plot(year,np.nansum(SO2[1:7,:],axis=0),'-',color="k",linewidth=0.3)
ax3.plot(year,np.nansum(SO2[1:8,:],axis=0),'-',color="k",linewidth=0.3)
ax3.fill_between(year,SO2[1,:],0, where=SO2[1,:]>=0,facecolor='m',label='India',alpha=0.7)
ax3.fill_between(year,np.nansum(SO2[1:3,:],axis=0),np.nansum(SO2[1:2,:],axis=0),facecolor="r",label='China',alpha=0.7)
ax3.fill_between(year,np.nansum(SO2[1:4,:],axis=0),np.nansum(SO2[1:3,:],axis=0),facecolor="y",label='Europe',alpha=0.7)
ax3.fill_between(year,np.nansum(SO2[1:5,:],axis=0),np.nansum(SO2[1:4,:],axis=0),facecolor="g",label='USA',alpha=0.7)
ax3.fill_between(year,np.nansum(SO2[1:6,:],axis=0),np.nansum(SO2[1:5,:],axis=0),facecolor="c",label='Brazil',alpha=0.7)
ax3.fill_between(year,np.nansum(SO2[1:7,:],axis=0),np.nansum(SO2[1:6,:],axis=0),facecolor="b",label='Australia',alpha=0.7)
ax3.fill_between(year,np.nansum(SO2[1:8,:],axis=0),np.nansum(SO2[1:7,:],axis=0),facecolor="orange",label='Southern Africa',alpha=0.7)
ax3.fill_between(year,SO2[0,:],np.nansum(SO2[1:8,:],axis=0),facecolor="wheat",label='Others',alpha=0.7)

ax4 =plt.subplot2grid((2, 4), (0, 3), colspan=1);ax4.axis('off')

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
r_patch = mpatches.Patch(color='r', alpha=0.7,label='India')
m_patch = mpatches.Patch(color='m', alpha=0.7,label='China')
y_patch = mpatches.Patch(color='y', alpha=0.7,label='Europe')
g_patch = mpatches.Patch(color='g', alpha=0.7,label='USA')
c_patch = mpatches.Patch(color='c', alpha=0.7,label='Brazil')
b_patch = mpatches.Patch(color='b', alpha=0.7,label='Australia')
o_patch = mpatches.Patch(color='orange', alpha=0.7,label='Southern Africa')
w_patch = mpatches.Patch(color='wheat', alpha=0.7,label='Others')
lines = [w_patch,o_patch,b_patch,c_patch,g_patch,y_patch,r_patch,m_patch]
labels = ['Others','Southern Africa','Australia','Brazil','USA','Europe','China','India']
legend = plt.legend(lines,labels,ncol=1,loc=10,labelspacing=1.5,markerscale =10)
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)


lon,lat,AOD=get_AOD();
plot_setup();pad= 5;
ax =plt.subplot2grid((2, 4), (1, 0), colspan=2)
colormap='RdYlBu';colormap = reverse_colourmap(colormap);colorbar_min=-0.1;colorbar_max=0.1
p_value = np.zeros((np.shape(AOD[0,:,:])))
ax.annotate(r'(d) $\mathrm{\mathsf{\Delta}}$AOD (2006-2025)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
spatial_figure(ax,stats.nanmean(AOD[0:20,:,:],axis=0),lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)

ax =plt.subplot2grid((2, 4), (1, 2), colspan=2)
ax.annotate(r'(e) $\mathrm{\mathsf{\Delta}}$AOD (2081-2100)',xy=(0.02,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)	
colormesh1=spatial_figure(ax,stats.nanmean(AOD[75:95,:,:],axis=0),lon,lat,colormap,colorbar_min,colorbar_max,p_value,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.2, 0.05, 0.60, 0.025])
char = fig.colorbar(colormesh1,orientation='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.01,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.04, bottom=0.08, right=0.97, top=0.95, wspace=0.15, hspace=0.20);
plt.savefig('Fig1_BC_OC_SO2_countries.png', format='png', dpi=600)
