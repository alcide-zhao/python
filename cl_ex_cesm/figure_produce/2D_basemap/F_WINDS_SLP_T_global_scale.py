
"""
This is to show the variation in  
		1. The surface temperature
		2. the sea level pressure with 850pa winds overlapped
		3. The 200hpa winds 
"""
import numpy as np
from scipy import stats
import math
import os; import site
from scipy.stats import ttest_ind as test
import netCDF4 as nc4
from scipy.interpolate import interp2d  as interp2d
from scipy import stats
lib_path = os.path.join(
	os.path.realpath(
		os.path.dirname(__file__)
	), 
	os.path.pardir, 
	os.path.pardir,
    # 'lib'
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio
oceanmask=sio.loadmat('//home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan


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

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55.1],'EA':[100,145,20,50],'SA':[65,100,5,30]}
region = rergion_dic['GLOBE'][:]

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_P.mat')
PS_rcp85 = data['PS_rcp85'];PS_fixa = data['PS_fixa']; 
del data

data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_TS_PSL_Q_UV(Q)_LHFLX_FSNS(C).mat')
lon = data['lon'][0,:];lat = data['lat'][0,:];#lev = data['lev'][0,:];
TS_rcp85 = data['TS_rcp85'];TS_fixa = data['TS_fixa'];TS = TS_rcp85-TS_fixa;
PSL_rcp85 = data['PSL_rcp85'];PSL_fixa = data['PSL_fixa']; 
U_rcp85 = data['U850_rcp85'];U_fixa = data['U850_fixa'];
V_rcp85 = data['V850_rcp85'];V_fixa = data['V850_fixa'];
del data

PSL = np.multiply(PSL_rcp85-PSL_fixa,1);
U = np.multiply(U_rcp85-U_fixa,1);
V = np.multiply(V_rcp85-V_fixa,1);


p_mask = stats.nanmean(PSL_rcp85-PS_rcp85,axis=0)/100;  p_mask[:] = 1
# # p_mask[p_mask<=300] =1 ;p_mask[p_mask>300] =0 
_,_,p_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,p_mask)

def spatial_diff_sig(region,lon,lat,variable_D):
	lons,lats,variable = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable_D)
	lon_interp_AM = np.arange(0,180.1,4);lat_interp_AM = np.arange(-10,80.1,3);
	f = interp2d(lons,lats,p_mask,kind='linear')
	p_ma = f(lon_interp_AM, lat_interp_AM); p_ma[p_ma==0] =np.nan
	variable_m =stats.nanmean(variable[0:10,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_0615_D = np.multiply(f(lon_interp_AM, lat_interp_AM),p_ma)
	variable_m = stats.nanmean(variable[25:45,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_3150_D = np.multiply(f(lon_interp_AM, lat_interp_AM),p_ma)
	variable_m = stats.nanmean(variable[75:95,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_8100_D = np.multiply(f(lon_interp_AM, lat_interp_AM),p_ma)
	return lon_interp_AM,lat_interp_AM,att_0615_D,att_3150_D,att_8100_D

lon_interp_AM,lat_interp_AM,PSL_0615_D,PSL_3150_D,PSL_8100_D = spatial_diff_sig(region,lon,lat,PSL)
_,_,U_0615_D,U_3150_D,U_8100_D = spatial_diff_sig(region,lon,lat,U)
_,_,V_0615_D,V_3150_D,V_8100_D = spatial_diff_sig(region,lon,lat,V)

def spatial_diff_sig(region,lon,lat,variable_D):
	lons,lats,variable = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable_D)
	lon_interp_AM = np.arange(0,180.1,4);lat_interp_AM = np.arange(-10,80.1,3);
	variable_m =stats.nanmean(variable[0:10,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_0615_D = np.multiply(f(lon_interp_AM, lat_interp_AM),1)
	variable_m = stats.nanmean(variable[25:45,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_3150_D = np.multiply(f(lon_interp_AM, lat_interp_AM),1)
	variable_m = stats.nanmean(variable[75:95,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_8100_D = np.multiply(f(lon_interp_AM, lat_interp_AM),1)
	return lon_interp_AM,lat_interp_AM,att_0615_D,att_3150_D,att_8100_D


lon,lat,AOD=get_AOD();
_,_,TS_0615_D,TS_3150_D,TS_8100_D = spatial_diff_sig(region,lon,lat,TS)
_,_,AOD_0615_D,AOD_3150_D,AOD_8100_D = spatial_diff_sig(region,lon,lat,AOD)


fig = plt.figure(facecolor='White',figsize=(9, 12));pad= 5 

colorbar_min=-0.1;colorbar_max=0.1;colormap='RdYlBu_r';p_value=np.zeros((np.shape(TS_0615_D)))
ax = plt.subplot(3,1,1)
ax.annotate(r'(a) $\mathrm{\mathsf{\Delta}}$AOD',xy=(0.02,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
colormesh1 = spatial_figure_norm(ax,AOD_8100_D,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=False )
cbar_ax = fig.add_axes([0.86, 0.67, 0.02, 0.26])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2)); # 
char.set_label('#')


ax = plt.subplot(3,1,2)
colorbar_min=-1.6;  colorbar_max =1.6; colormap ='RdYlBu_r';
colormesh1 = spatial_figure_norm(ax,TS_8100_D,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=False )
ax.annotate('(b) SAT',xy=(0.02, 1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
cbar_ax = fig.add_axes([0.86, 0.39, 0.02, 0.24])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2)); # 
char.set_label('K')
				
ax = plt.subplot(3,1,3)
colorbar_min=-50.1;  colorbar_max =50; 	
ax.annotate('(c) Winds850+SLP',xy=(0.02, 1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
qk_caption = r'$\mathrm{\mathsf{1\/m\/s^{-1}}}$'; qk_scale =1
colormesh1 = spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,PSL_8100_D, U_8100_D, V_8100_D,qk_scale,qk_caption,qk_is = True, tb_lef=True,tb_bot=True )
cbar_ax = fig.add_axes([0.86, 0.09, 0.02, 0.24])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(-50,50.1,20))); # 
char.set_label('Pa')


plt.subplots_adjust(left=0.05, bottom=0.05, right=0.81, top=0.95, wspace=0.000000000001, hspace=0.15);
plt.savefig('winds_T.png', format='png', dpi=1200)

