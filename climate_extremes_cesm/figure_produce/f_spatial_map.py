
import site
import os
import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt

from scipy import stats
from matplotlib import rcParams
from scipy.interpolate import interp2d
from lib import *


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
	
# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/'

oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
time_s = 2006
time_e = 2100

precipitation extremes
att_ref = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
att_list=['rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot']
unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

# temperature extremes
# att_ref stores the reference value which is used to remove climatology
# att_ref = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_dic = {'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
# att_list = ['txx', 'txn', 'tnx', 'tnn', 'dtr', 'fd0', 'su25','id0','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
# unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}


rergion_dic={'GLOBE':[0,360,-90,90],'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,6.,17.]\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.],\
'ASIA':[0,360,0,90]}


region = 'SI'
# #######################################################
# # 1. spatial features                                 #
# #######################################################

# att_dic = {'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}

figure_no = 1
for att_name in ['rx5day']:
	colorbar_unit = unit_dic[att_name]
	
	
	colorbar_min = 0
	colorbar_max = 40
	colormap ='jet'
	
	
	year_b = 2006 #beginning point of the time slice 
	year_step = 10 # moving bins every tenyears
	year_span = 20 #length of the time slice]
	year_e = ((2100-year_span-2006)/year_step)*year_step+2006+year_step # the last begining point of the time sclice
	
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	
	year_s = range(year_b,year_e,year_step)
	
	spst = [1,3] # subplot_style define the rows and the columns of the s being plotted
	subplot_no = 1
	figure_no = 1
	plt.figure(figure_no, facecolor='White',figsize=[12,8])
	for year in year_s:
	
		subplot_no =subplot_no+1
		file_name = 'Spatial_ensumble_mean_Prep_extremes_global_2006_2100_fixA.nc'
		nc_fid = nc4.Dataset(input_path+file_name,mode='r')
		ax=plt.subplot(spst[0],spst[1],subplot_no)
		time = nc_fid.variables['time'][:]
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		att_value = nc_fid.variables[att_name][:]
		att_value = np.multiply(att_value,oceanmask)
		lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
		# the PD_F refers to the beggining 10 years of the simulation under pcp_fixA scenario
		PD_F =stats.nanmean(att_clipped_r85F[0:9,:,:],axis=0)
		att_annual_F = stats.nanmean(att_clipped_r85F[year-2006:year-2006+20,:,:],axis=0)-PD_F	
		figure_title = str(year)+'_'+str(year+20)+'_'+'_FixA'
		statistic,p_value=scipy.stats.ttest_1samp(att_annual_F, PD_F, axis=0)
		p_value[p_value>=0.05] =np.nan
		p_value[p_value<0.05] =1
		# print p_value[p_value==1]
		# print p_vale(p_value==1)
		spatial_figure(att_annual_F,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name],p_value)
		nc_fid.close()
		
		file_name = 'Spatial_ensumble_mean_Prep_extremes_global_2006_2100_Rcp8.5.nc'
		nc_fid = nc4.Dataset(input_path+file_name,mode='r')
		ax=plt.subplot(spst[0],spst[1],subplot_no)
		time = nc_fid.variables['time'][:]
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		att_value = nc_fid.variables[att_name][:]
		att_value = np.multiply(att_value,oceanmask)
		lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
		att_annual_R = stats.nanmean(att_clipped_r85[year-2006:year-2006+20,:,:],axis=0)-PD_F
		figure_title = str(year)+'_'+str(year+20)+'_'+'_RCP8.5'
		statistic,p_value=scipy.stats.ttest_1samp(att_annual_R, PD_F, axis=0)
		p_value[p_value>=0.05] =np.nan
		p_value[p_value<0.05] =1
		spatial_figure(att_annual_R,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name],p_value)	
		nc_fid.close()
		
		
		subplot_no =subplot_no+1
		ax=plt.subplot(spst[0],spst[1],subplot_no)
		att_annual_D = att_annual_R-att_annual_F
		figure_title = 'Rcp-FixA_'str(year)+'_'+str(year+20)+'_'+'_Diff'
		statistic,p_value=scipy.stats.ttest_1samp(att_annual_F, att_annual_R, axis=0)
		p_value[p_value>=0.05] =np.nan
		p_value[p_value<0.05] =1
		spatial_figure(att_annual_D,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name],p_value)
		subplot_no =subplot_no+1
	plt.show()
	# plt.savefig(input_path+att_name+'.tif')
	figure_no+ figure_no+1	
	'''	
	# follw yang's figures 
	ax=plt.subplot(2,3,1)
	file_name = 'spatial_ensumble_mean_Temp_extremes_global_2006_2080.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lats = nc_fid.variables['lat'][:]
	lons = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[att_name][:]
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	PD_R =stats.nanmean(att_value[0:9,:,:],axis=0)
	att_3050_pd_R = stats.nanmean(att_value[25:44,:,:],axis=0)-PD_R
	att_3050_pd_R = np.multiply(att_3050_pd_R,oceanmask)
	figure_title = att_name+'_2031-2050-PD'+'_RCP8.5'
	spatial_figure(att_3050_pd_R,lons,lats,figure_title,'jet',colorbar_min,colorbar_max,unit_dic[att_name])
	nc_fid.close()

	ax1=plt.subplot(2,3,4)
	file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2081_2100.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lats = nc_fid.variables['lat'][:]
	lons = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[att_name][:]
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	att_8000_pd_R = stats.nanmean(att_value[0:19,:,:],axis=0)-PD_R
	att_8000_pd_R = np.multiply(att_8000_pd_R,oceanmask)
	figure_title = att_name+'_2081-2100-PD'+'_RCP8.5'
	spatial_figure(att_8000_pd_R,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name])
	nc_fid.close()
	

	file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100.nc'
	ax1=plt.subplot(2,3,2)
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lats = nc_fid.variables['lat'][:]
	lons = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[att_name][:]
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	PD_F = stats.nanmean(att_value[0:9,:,:],axis=0)
	att_3050_pd_F = stats.nanmean(att_value[25:44,:,:],axis=0)-PD_F
	att_3050_pd_F = np.multiply(att_3050_pd_F,oceanmask)
	figure_title = att_name+'_2031-2050-PD'+'_RCP8.5_FxixA'
	spatial_figure(att_3050_pd_F,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name])
	nc_fid.close()
	
	ax1=plt.subplot(2,3,5)

	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	att_8000_pd_F = stats.nanmean(att_value[-20:-1,:,:],axis=0)-PD_F
	att_8000_pd_F = np.multiply(att_8000_pd_F,oceanmask)
	figure_title = att_name+'_2081-2100-PD'+'_RCP8.5_FixA'
	spatial_figure(att_8000_pd_F,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name])

	
	ax1=plt.subplot(2,3,3)
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	diff_3050_R_F = att_3050_pd_R-att_3050_pd_F
	figure_title = att_name+'_Rcp-FixA_2031_2050'
	spatial_figure(diff_3050_R_F,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name])
		
	
	ax1=plt.subplot(2,3,6)
	# lons,lats,att_value = range_clip(60,150,0,55,lons,lats,att_value)
	diff_8000_R_F = att_8000_pd_R-att_8000_pd_F
	figure_title = att_name+'_Rcp-FixA_2031_2050'
	spatial_figure(diff_3050_R_F,lons,lats,figure_title,colormap,colorbar_min,colorbar_max,unit_dic[att_name])
	
	
	
	
	plt.show()
	# plt.savefig(input_path+att_name+'.tif')
	figure_no+ figure_no+1
	'''