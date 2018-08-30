"""
This is to plot the ozone concentration difference between 2010 and 1970 
the 1970 data is from ozone_1.9x2.5_L26_1850-2005_c090803.nc
the 2010 data is from ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202.nc
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
from scipy.stats import mannwhitneyu as man_test
lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


def get_LENS_OZONE():
	input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/emissions_PP/ozone/'
	ozone1970 = input_path+'ozone_1.9x2.5_L26_1850-2005_c090803.nc'
	nc_fid = nc4.Dataset(ozone1970,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	lev = nc_fid.variables['lev'][:]
	o1970 = nc_fid.variables['O3'][156:168,:,:,:]
	nc_fid.close()
	ozone2010 = input_path+'ozone_rcp85_v1_1.9x2.5_L26_1995-2105_c100202.nc'
	nc_fid = nc4.Dataset(ozone2010,mode='r')
	o2010 = nc_fid.variables['O3'][48:60,:,:,:]
	nc_fid.close()
	ozone_lat_alt_cross = stats.nanmean(stats.nanmean(o2010-o1970,axis=3),axis=0)
	ozone_lat_alt_cross = np.ma.masked_invalid(ozone_lat_alt_cross)
	ozone_lat_alt_cross = np.ma.masked_equal(ozone_lat_alt_cross,0)
	return lat,lev,ozone_lat_alt_cross
	
		
def data_readin(variable,exp):	
	def data_netcdf(exp,scenario,variable,region_key='All'):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'+exp+'/'
		var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lev = nc_fid.variables['lev'][:]
		days = nc_fid.variables['time'][:]; 
		data = nc_fid.variables['O3'][:]
		nc_fid.close()
		return lev,lat,data
		
	if exp=='full_chem':
		lev,lat,EdgRef = data_netcdf(exp,'EdgRefP',variable);
		_,_,EdgRef70AP = data_netcdf(exp,'EdgRef70AP',variable);
		_,_,EdgEne = data_netcdf(exp,'EdgEneP',variable);
		_,_,EdgTech = data_netcdf(exp,'EdgTechP',variable);		
	return lev,lat,EdgRef70AP,EdgEne,EdgRef,EdgTech
	
def spa_pat_reg_mean(variable):	
	def mannwhitneyu_test(vairable1,variable2):
		p_threshold=0.05
		size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
		p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
		# from scipy.stats import ttest_ind as test
		for x in range(size[0]):
			for y in range(size[1]):
				cache1 = vairable1[:,x,y]
				cache2 = variable2[:,x,y]
				if (cache2==cache1).all(): p_value[x,y] = np.nan
				else: _,p_value[x,y] = man_test(cache1,cache2);
		p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
		return p_value
	####calculate the difference and their significance
	def diff_sig(vairable1,variable2):
		dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
		sig = mannwhitneyu_test(vairable1, variable2)
		# sig_dif = np.multiply(dif,sig)
		return dif,sig	
	lev,lat,EdgRef70AP,EdgEne,EdgRef,EdgTech = data_readin(variable,exp='full_chem');    # fully coupled total respons	
	BEoA,BEoA_s = diff_sig(EdgRef,EdgRef70AP)
	Ene,Ene_s = diff_sig(EdgEne, EdgRef70AP)
	Tech,Tech_s = diff_sig(EdgTech, EdgRef70AP)	
	CAMchem = np.stack((BEoA,Ene,Tech),axis=0)
	return lev,lat,CAMchem

lat,lev,LENS = get_LENS_OZONE();
lev1,lat1,CAMchem =spa_pat_reg_mean('O3.YZmean');
x,y= np.meshgrid(lat,lev)
x1,y1= np.meshgrid(lat1,lev1);


fig = plt.figure(facecolor='White',figsize=[15,6]);plot_setup();pad= 5;cmap = 'RdBu_r';

###CAM5
ax =plt.subplot2grid((3, 4), (0, 0), rowspan=1);colorbar_min=-0.2;colorbar_max=0.2;
cmap = discrete_cmap(16,cmap)
colormesh1 = plt.pcolormesh(x[0:13,:],np.flipud(y[0:13,:]),np.flipud(LENS[0:13,:])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([0,193]);ax.set_xticks([]);
plt.gca().invert_yaxis()

ax.annotate('(a) CAM5',xy=(0.05,1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax =plt.subplot2grid((3, 4), (1, 0), rowspan=2);colorbar_min=-0.02;colorbar_max=0.02;
colormesh1 = plt.pcolormesh(x[13:-1,:],np.flipud(y[13:-1,:]),np.flipud(LENS[13:-1,:])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([266,1000]);ax.set_xticks(range(-90,91,30));ax.set_ylabel('');
plt.gca().invert_yaxis()


### CAM5-CHEM BEST ESTIMATE
ax =plt.subplot2grid((3, 4), (0, 1), rowspan=1);colorbar_min=-0.2;colorbar_max=0.2;
colormesh1 = plt.pcolormesh(x1[0:13,:],np.flipud(y1[0:13,:]),np.flipud(CAMchem[0,0:13,:,0])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([0,193]);ax.set_xticks([]);
plt.gca().invert_yaxis()

ax.annotate('(b) Best Estimate',xy=(0.05,1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax =plt.subplot2grid((3, 4), (1, 1), rowspan=2);colorbar_min=-0.02;colorbar_max=0.02;
colormesh1 = plt.pcolormesh(x1[13:-1,:],np.flipud(y1[13:-1,:]),np.flipud(CAMchem[0,13:-1,:,0])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([266,1000]);ax.set_xticks(range(-90,91,30));ax.set_ylabel('');
plt.gca().invert_yaxis()

### CAM5-CHEM Ene
ax =plt.subplot2grid((3, 4), (0, 2), rowspan=1);colorbar_min=-0.2;colorbar_max=0.2;
colormesh1 = plt.pcolormesh(x1[0:13,:],np.flipud(y1[0:13,:]),np.flipud(CAMchem[1,0:13,:,0])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([0,193]);ax.set_xticks([]);
plt.gca().invert_yaxis()
ax.annotate('(c) Energy Consumption',xy=(0.05,1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax =plt.subplot2grid((3, 4), (1, 2), rowspan=2);colorbar_min=-0.02;colorbar_max=0.02;
colormesh1 = plt.pcolormesh(x1[13:-1,:],np.flipud(y1[13:-1,:]),np.flipud(CAMchem[1,13:-1,:,0])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([266,1000]);ax.set_xticks(range(-90,91,30));ax.set_ylabel('');
plt.gca().invert_yaxis()

### CAM5-CHEM Tech
ax =plt.subplot2grid((3, 4), (0, 3), rowspan=1);colorbar_min=-0.2;colorbar_max=0.2;
colormesh1 = plt.pcolormesh(x1[0:13,:],np.flipud(y1[0:13,:]),np.flipud(CAMchem[2,0:13,:,0])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([0,193]);ax.set_xticks([]);
plt.gca().invert_yaxis()
cbar_ax = fig.add_axes([0.95, 0.67, 0.015, 0.27])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.25)*(colorbar_max-colorbar_min)+colorbar_min,2))


ax.annotate('(d) Technology Advancements',xy=(0.05,1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
ax =plt.subplot2grid((3, 4), (1, 3), rowspan=2);colorbar_min=-0.02;colorbar_max=0.02;
colormesh1 = plt.pcolormesh(x1[13:-1,:],np.flipud(y1[13:-1,:]),np.flipud(CAMchem[2,13:-1,:,0])*10**(6),cmap=cmap,vmin=colorbar_min, vmax=colorbar_max);
ax.set_xlim([-90,90]);ax.set_ylim([266,1000]);ax.set_xticks(range(-90,91,30));ax.set_ylabel('');
# ax.set_xlabel('Lattitude (degree)');
plt.gca().invert_yaxis()
cbar_ax = fig.add_axes([0.95, 0.07, 0.015, 0.57])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.25)*(colorbar_max-colorbar_min)+colorbar_min,2))


plt.subplots_adjust(left=0.05, bottom=0.06, right=0.93, top=0.95, wspace=0.2, hspace=0.02); 
plt.savefig('Ozone_COMPARISON.png', format='png', dpi=600)







