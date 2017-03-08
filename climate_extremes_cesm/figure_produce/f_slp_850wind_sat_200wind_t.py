
"""
This is to show the variation in  
		1. The surface temperature
		2. the sea level pressure with 850pa winds overlapped
		3. The 200hpa winds 
"""
import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import os
from scipy import stats
import site
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
# from scipy.interpolate import interp2d

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)

site.addsitedir(lib_path)

from lib import *



# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}

region = rergion_dic['ASIA'][:]



def Aerosol_TS_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name):
	'''
	This fUnction is to process thw diagnostic field into JJA mean from the CESM
	monthly data
	The inpUt is the file which contains both rcp85 and fixA simUlations
	The oUtpUt is shoUld be normaly lons, lats, time_series and range clipped data
	the lev information willl also be returned for 4D field
	'''
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	rcp85 = nc_fid.variables['rcp85'][:] 
	RCP85_MV = nc_fid.variables['rcp85'].missing_value
	rcp85[rcp85 == RCP85_MV] = np.nan
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:] 
	RCP85_fixA_MV = nc_fid.variables['rcp85_fixA'].missing_value
	rcp85_fixA[rcp85_fixA == RCP85_fixA_MV] = np.nan
	units =  nc_fid.variables['rcp85_fixA'].units
	# long_name =  nc_fid.variables['rcp85_fixA'].long_name
	long_name= 'Surface Temperature (radiative)'
	
	lons,lats,rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85)
	lons,lats,fixA_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,rcp85_fixA)	
	size = np.shape(rcp85_clipped)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	###JJA mean
	if (np.rank(rcp85_clipped) == 4):
		levs = nc_fid.variables['lev'][:]
		rcp85_JJAmean = np.empty((95,30,size[2],size[3])); rcp85_JJAmean[:] =np.nan
		fixa_JJAmean = np.empty((95,30,size[2],size[3])); fixa_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = rcp85_clipped[layer_b:layer_e+1,:,:,:]
			rcp85_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
			cache = fixA_clipped[layer_b:layer_e+1,:,:,:]
			fixa_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
	else:
		levs=np.nan
		rcp85_JJAmean = np.empty((95,size[1],size[2])); rcp85_JJAmean[:] =np.nan
		fixa_JJAmean = np.empty((95,size[1],size[2])); fixa_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = rcp85_clipped[layer_b:layer_e+1,:,:]
			rcp85_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
			cache = fixA_clipped[layer_b:layer_e+1,:,:]
			fixa_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	nc_fid.close()
	return lons, lats, levs, year_series,rcp85_JJAmean,fixa_JJAmean,units,long_name

#TS
variable = 'TS'
imput_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = imput_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, TS_rcp85_JJAmean,TS_fixa_JJAmean,units,long_name= Aerosol_TS_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

TS_0615_D = stats.nanmean((TS_rcp85_JJAmean-TS_fixa_JJAmean)[0:10,:,:],axis=0)
_, TS_0615_pvalue =scipy.stats.ttest_1samp(TS_rcp85_JJAmean[0:10,:,:], stats.nanmean(TS_fixa_JJAmean[0:10,:,:],axis=0), axis=0)
TS_0615_pvalue[TS_0615_pvalue>0.05] =np.nan; TS_0615_pvalue[TS_0615_pvalue<=0.05] =1

TS_3645_D = stats.nanmean((TS_rcp85_JJAmean-TS_fixa_JJAmean)[30:40,:,:],axis=0)
_, TS_3645_pvalue =scipy.stats.ttest_1samp(TS_rcp85_JJAmean[30:40,:,:], stats.nanmean(TS_fixa_JJAmean[30:40,:,:],axis=0), axis=0)
TS_3645_pvalue[TS_3645_pvalue>0.05] =np.nan; TS_3645_pvalue[TS_3645_pvalue<=0.05] =1

TS_9100_D = stats.nanmean((TS_rcp85_JJAmean-TS_fixa_JJAmean)[85:95,:,:],axis=0)
_, TS_9100_pvalue =scipy.stats.ttest_1samp(TS_rcp85_JJAmean[85:95,:,:], stats.nanmean(TS_fixa_JJAmean[85:95,:,:],axis=0), axis=0)
TS_9100_pvalue[TS_9100_pvalue>0.05] =np.nan; TS_9100_pvalue[TS_9100_pvalue<=0.05] =1

#PS
variable = 'PS'
imput_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = imput_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, PS_rcp85_JJAmean,PS_fixa_JJAmean,units,long_name= Aerosol_TS_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
PS_0615_D = stats.nanmean((PS_rcp85_JJAmean-PS_fixa_JJAmean)[0:10,:,:],axis=0)
_, PS_0615_pvalue =scipy.stats.ttest_1samp(PS_rcp85_JJAmean[0:10,:,:], stats.nanmean(PS_fixa_JJAmean[0:10,:,:],axis=0), axis=0)
PS_0615_pvalue[PS_0615_pvalue>0.05] =np.nan; PS_0615_pvalue[PS_0615_pvalue<=0.05] =1

PS_3645_D = stats.nanmean((PS_rcp85_JJAmean-PS_fixa_JJAmean)[30:40,:,:],axis=0)
_, PS_3645_pvalue =scipy.stats.ttest_1samp(PS_rcp85_JJAmean[30:40,:,:], stats.nanmean(PS_fixa_JJAmean[30:40,:,:],axis=0), axis=0)
PS_3645_pvalue[PS_3645_pvalue>0.05] =np.nan; PS_3645_pvalue[PS_3645_pvalue<=0.05] =1

PS_9100_D = stats.nanmean((PS_rcp85_JJAmean-PS_fixa_JJAmean)[85:95,:,:],axis=0)
_, PS_9100_pvalue =scipy.stats.ttest_1samp(PS_rcp85_JJAmean[85:95,:,:], stats.nanmean(PS_fixa_JJAmean[85:95,:,:],axis=0), axis=0)
PS_9100_pvalue[PS_9100_pvalue>0.05] =np.nan; PS_9100_pvalue[PS_9100_pvalue<=0.05] =1


#U
variable = 'U'
imput_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = imput_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, U_rcp85,U_fixa,units,long_name= Aerosol_TS_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
U_fixa_850 =U_fixa[:,23,:,:]; U_rcp85_850 =U_rcp85[:,23,:,:]
U_fixa_200 =U_fixa[:,12,:,:]; U_rcp85_200 =U_rcp85[:,12,:,:]

U_0615_D_850 = stats.nanmean((U_rcp85_850-U_fixa_850)[0:10,:,:],axis=0)
U_3645_D_850 = stats.nanmean((U_rcp85_850-U_fixa_850)[30:40,:,:],axis=0)
U_9100_D_850 = stats.nanmean((U_rcp85_850-U_fixa_850)[85:95,:,:],axis=0)

U_0615_D_200 = stats.nanmean((U_rcp85_200-U_fixa_200)[0:10,:,:],axis=0)
U_3645_D_200 = stats.nanmean((U_rcp85_200-U_fixa_200)[30:40,:,:],axis=0)
U_9100_D_200 = stats.nanmean((U_rcp85_200-U_fixa_200)[85:95,:,:],axis=0)

#V
variable = 'V'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, V_rcp85,V_fixa,units,long_name= Aerosol_TS_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
V_fixa_850 =U_fixa[:,23,:,:]; V_rcp85_850 =U_rcp85[:,23,:,:]
V_fixa_200 =U_fixa[:,12,:,:]; V_rcp85_200 =U_rcp85[:,12,:,:]

V_0615_D_850 = stats.nanmean((V_rcp85_850-V_fixa_850)[0:10,:,:],axis=0)
V_3645_D_850 = stats.nanmean((V_rcp85_850-V_fixa_850)[30:40,:,:],axis=0)
V_9100_D_850 = stats.nanmean((V_rcp85_850-V_fixa_850)[85:95,:,:],axis=0)

V_0615_D_200 = stats.nanmean((V_rcp85_200-V_fixa_200)[0:10,:,:],axis=0)
V_3645_D_200 = stats.nanmean((V_rcp85_200-V_fixa_200)[30:40,:,:],axis=0)
V_9100_D_200 = stats.nanmean((V_rcp85_200-V_fixa_200)[85:95,:,:],axis=0)

# now plotting
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(6, 10), facecolor='White')
fig.tight_layout()
# tight_layouts doesn't take these labels into account. We'll need 
# to make some room. These numbers are are manually tweaked. 
# You could automatically calculate them, but it's a pain.
fig.subplots_adjust(left=0.15, top=0.9)


colormap ='jet'; 
pad= 5 
#TS
colorbar_min=-0.4; colorbar_max = 1.6;
ax = plt.subplot(3,3,1)
spatial_figure(TS_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,TS_0615_pvalue)
ax.annotate('2006_2015',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax.annotate('TS',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(3,3,2)
spatial_figure(TS_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,TS_3645_pvalue)
ax.annotate('2036_2045',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(3,3,3)
colormesh1=spatial_figure(TS_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,TS_9100_pvalue)
ax.annotate('2090_2010',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
cbar_ax = fig.add_axes([0.93, 0.65, 0.01, 0.25])
char = fig.colorbar(colormesh1, cax=cbar_ax,extend='both')
char.set_label('K',fontsize=15)

#PS and 850hpa winds
colormap ='BrBG'; 	
colorbar_min=-60; colorbar_max = 60;		
ax=plt.subplot(3,3,4)
spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,PS_0615_D,U_0615_D_850,V_0615_D_850,projection='cyl')
ax.annotate('SLP\n850hPa Wind',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(3,3,5)
spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,PS_3645_D,U_3645_D_850,V_3645_D_850,projection='cyl')
plt.subplot(3,3,6)
colormesh2 = spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,PS_9100_D,U_9100_D_850,V_9100_D_850,projection='cyl')
cbar_ax = fig.add_axes([0.93, 0.35, 0.01, 0.25])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both')
char.set_label('Pa',fontsize=15)

# U
colormap ='gist_rainbow'; 
colorbar_min=-0.2; colorbar_max = 2	
ax=plt.subplot(3,3,7)
W_0615_D_200 = np.sqrt(np.power(U_0615_D_200,2)+np.power(V_0615_D_200,2))
spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,W_0615_D_200,U_0615_D_200,V_0615_D_200,projection='cyl')
ax.annotate('200hPa\nWind',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(3,3,8)
W_3645_D_200 = np.sqrt(np.power(U_3645_D_200,2)+np.power(V_3645_D_200,2))
spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,W_3645_D_200,U_3645_D_200,V_3645_D_200,projection='cyl')
plt.subplot(3,3,9)
W_9100_D_200 = np.sqrt(np.power(U_9100_D_200,2)+np.power(V_9100_D_200,2))
colormesh3 = spatial_scaler_vector(lons,lats,colormap,colorbar_min,colorbar_max,W_9100_D_200,U_9100_D_200,V_9100_D_200,projection='cyl')
cbar_ax = fig.add_axes([0.93, 0.05, 0.01, 0.25])
char = fig.colorbar(colormesh3, cax=cbar_ax,extend='both')
char.set_label('m/s',fontsize=15)

plt.show()
