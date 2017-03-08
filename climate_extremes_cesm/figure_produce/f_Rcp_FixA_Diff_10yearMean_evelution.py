# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 2017
This is to produce the time evelution of 10-years averaged spatial precipitation extremes 
Here using the rx5day
@author: Alcide.Zhao
"""
import site
import os
import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from scipy import stats
from matplotlib import rcParams
from scipy.interpolate import interp2d


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)
site.addsitedir(lib_path)

def time_seeries_of_spatial_mean(time_s, time_e, time, data):
	"""
	for the input array, values of the interesed geographical domain is to be averaged
	the time seris of the averaged value is then produced 
	"""
	time_series = np.empty((time_e-time_s+1))
	
	layer_s = [layer for layer in range(len(time)) if time[layer] == time_s]
	layer_e = [layer for layer in range(len(time)) if time[layer] == time_e]
	layer_s = layer_s[0]
	layer_e = layer_e[0]
	for layer in range(layer_s , layer_e+1):
		# print data[layer,:,:]
		time_series[layer-layer_s] = stats.nanmean(stats.nanmean(data[layer,:,:]))
	return time_series



oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
time_s = 2006
time_e = 2100

# precipitation extremes
att_ref = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
att_list=['rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot']
unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm','mean_precip':'Mm'}

# temperature extremes
# att_ref stores the reference value which is used to remove climatology
# att_ref = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_dic = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_list = ['txx', 'txn', 'tnx', 'tnn', 'dtr', 'fd0', 'su25','id0','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
# unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}


rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}

region = rergion_dic['ASIA'][:]

# #######################################################
# # 1. spatial features                                 #
# #######################################################

# att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}

att_name='rx5day'
colorbar_unit = unit_dic[att_name]
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
# series if the starting points of the time slices

file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

time = nc_fid.variables['time'][:]
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]
att_value = nc_fid.variables[att_name][:]
att_value = np.multiply(att_value,oceanmask)
lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
nc_fid.close()
_,_,ocean_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,oceanmask)


file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
# time = nc_fid.variables['time'][:]
# lat = nc_fid.variables['lat'][:]
# lon = nc_fid.variables['lon'][:]
att_value = nc_fid.variables[att_name][:]
att_value = np.multiply(att_value,oceanmask)
lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
nc_fid.close()



"""
colorbar_min = -5
colorbar_max = 20
year_b = 2010 #beginning point of the time slice 
year_e = 2100
year_step = 10 # moving bins every tenyears
year_span = 20 #length of the time slice]
spst = [2,int(figure_no/2)] # subplot_style define the rows and the columns of the s being plotted
year_s = range(year_b,year_e-year_span+year_step,year_step)
figure_no = np.shape(year_s)[0]

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	plt.figure(figure_n, facecolor='White',figsize=[10,4])
	plt.subplot(spst[0],spst[1],subplot_no)
	# the PD_F refers to the beggining 10 years of the simulation under pcp_fixA scenario
	PD_F =stats.nanmean(att_clipped_r85F[0:9,:,:],axis=0)
	statistic,p_value=scipy.stats.ttest_1samp(att_clipped_r85F[year-2006:year-2006+year_span,:,:], PD_F, axis=0)
	att_annual_F = stats.nanmean(att_clipped_r85F[year-2006:year-2006+year_span,:,:],axis=0)-PD_F
	figure_title = str(year)+'_'+str(year+year_span)+'_FixA'
	p_value[p_value>=0.1] =np.nan
	p_value[p_value<0.1] =1
	p_value[:] =np.nan
	
	colormap ='jet'
	att_annual_F = np.multiply(att_annual_F,ocean_mask)
	# spatial_figure(np.divide(att_annual_F,PD_F)*100,lons,lats,figure_title,colormap,-12,50,'%',p_value)
	spatial_figure(att_annual_F,lons,lats,figure_title,colormap,-1,5,unit_dic[att_name],p_value)
	
	# spatial evelution under the RCP85 scenario
	figure_n = 2;
	# subplot_no =subplot_no+figure_no
	plt.figure(figure_n, facecolor='White',figsize=[10,4])
	plt.subplot(spst[0],spst[1],subplot_no)
	statistic,p_value=scipy.stats.ttest_1samp(att_clipped_r85[year-2006:year-2006+year_span,:,:],stats.nanmean(att_clipped_r85[0:9,:,:]) , axis=0)
	p_value[p_value>=0.1] =np.nan
	p_value[p_value<0.1] =1
	p_value[:] =np.nan
	att_annual_R = stats.nanmean(att_clipped_r85[year-2006:year-2006+year_span,:,:],axis=0)-PD_F
	figure_title = str(year)+'_'+str(year+year_span)+'_RCP8.5'
	
	colormap ='jet'
	att_annual_R = np.multiply(att_annual_R,ocean_mask)
	# spatial_figure(np.divide(att_annual_R,PD_F)*100,lons,lats,figure_title,colormap,-12,50,'%',p_value)	
	spatial_figure(att_annual_R,lons,lats,figure_title,colormap,-1,5,unit_dic[att_name],p_value)
	figure_n = 3;
	# subplot_no =subplot_no+figure_no
	plt.figure(figure_n, facecolor='White',figsize=[10,4])
	# spatial evelution of RCP85-FixA 
	ax=plt.subplot(spst[0],spst[1],subplot_no)
	att_annual_D = att_annual_R-att_annual_F
	figure_title = 'Rcp-FixA_'+str(year)+'_'+str(year+year_span)
	statistic,p_value=scipy.stats.ttest_1samp(att_clipped_r85[year-2006:year-2006+year_span,:,:]-PD_F,att_annual_F,axis=0)
	p_value[p_value>=0.1]=  np.nan
	p_value[p_value<0.1] =1
	p_value[:] =np.nan
	colormap ='seismic'
	att_annual_D = np.multiply(att_annual_D,ocean_mask)
	# spatial_figure(np.divide(att_annual_D,PD_F)*100,lons,lats,figure_title,colormap,-30,30,'%',p_value)
	spatial_figure(att_annual_D,lons,lats,figure_title,colormap,-0.4,1.2,unit_dic[att_name],p_value)
plt.show()
	# plt.savefig(input_path+att_name+'.tif')
	# figure_n+ figure_n+1	
"""


# follw yang's figures 


att_0615_R =stats.nanmean(att_clipped_r85[0:10,:,:],axis=0)
att_0615_F =stats.nanmean(att_clipped_r85F[0:10,:,:],axis=0)
att_0615_D = att_0615_R -att_0615_F
att_3645_pd_R = stats.nanmean(att_clipped_r85[30:40,:,:],axis=0)
att_3645_pd_F = stats.nanmean(att_clipped_r85F[30:40,:,:],axis=0)
att_3645_D = att_3645_pd_R -att_3645_pd_F
att_9100_pd_R = stats.nanmean(att_clipped_r85[85:95,:,:],axis=0)
att_9100_pd_F = stats.nanmean(att_clipped_r85F[85:95,:,:],axis=0)
att_9100_D = att_9100_pd_R -att_9100_pd_F


# now plotting
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 15), facecolor='White')
fig.tight_layout()
# tight_layout doesn't take these labels into account. We'll need 
# to make some room. These numbers are are manually tweaked. 
# You could automatically calculate them, but it's a pain.
fig.subplots_adjust(left=0.15, top=0.9)

		
# plt.setp(axes.flat, xlabel='X-label', ylabel='Y-label')
# rcp85 and fixa
colormap ='GnBu'; p_value = np.zeros((np.shape(att_0615_R))); colorbar_min=0;  colorbar_max = 250;
pad= 5 
ax = plt.subplot(3,3,1)
spatial_figure(att_0615_R,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('2006_2015',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax.annotate('RCP8.5',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(3,3,2)
spatial_figure(att_3645_pd_R,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('2036_2045',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
ax = plt.subplot(3,3,3)
ax.annotate('2090_2010',xy=(0.5, 1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

spatial_figure(att_9100_pd_R,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax=plt.subplot(3,3,4)
spatial_figure(att_0615_F,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('FixA',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
plt.subplot(3,3,5)
spatial_figure(att_3645_pd_F,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
plt.subplot(3,3,6)
colormesh1 =spatial_figure(att_9100_pd_F,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
#   Rcp and Fixa share the same colorbar while the diff plots use its own colorbar
cbar_ax = fig.add_axes([0.93, 0.35, 0.01, 0.55])
char = fig.colorbar(colormesh1, cax=cbar_ax,extend='both')
char.set_label('Day',fontsize=20)

# rcp85 minus fixa
colormap = 'BrBG'; colorbar_min=-10;  colorbar_max = 30; #BrBG seismic


ax=plt.subplot(3,3,7)
statistic,p_value=scipy.stats.ttest_1samp(att_clipped_r85[0:9,:,:], att_0615_F,axis=0)
p_value[p_value>=0.1] =np.nan; p_value[p_value<0.1] =1
# Zm = np.ma.masked_where(Z > 1.2, Z)
spatial_figure_norm(att_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
ax.annotate('Diff',xy=(-0.3, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')
				
plt.subplot(3,3,8)
statistic,p_value=scipy.stats.ttest_1samp(att_clipped_r85[30:40,:,:], att_3645_pd_F,axis=0)
p_value[p_value>=0.1] =np.nan; p_value[p_value<0.1] =1
spatial_figure_norm(att_3645_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)
plt.subplot(3,3,9)
statistic,p_value=scipy.stats.ttest_1samp(att_clipped_r85[85:95,:,:], att_9100_pd_F,axis=0)
p_value[p_value>=0.1] =np.nan; p_value[p_value<0.1] =1
# print p_value
colormesh2 =spatial_figure_norm(att_9100_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value)

cbar_ax = fig.add_axes([0.93, 0.05, 0.01, 0.25])
char = fig.colorbar(colormesh2, cax=cbar_ax,extend='both')
char.set_label('Day',fontsize=20)


plt.show()