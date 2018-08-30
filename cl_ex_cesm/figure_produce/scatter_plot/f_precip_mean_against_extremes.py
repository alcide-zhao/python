# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to show the correlations between mean and extreme precipitation evelutions
@author: Alcide.Zhao
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as spio
import glob  
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale as scale


lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)

site.addsitedir(lib_path)
# print lib_path
# from lib import ncdump

from lib import range_clip



rergion_dic={'GLOBE':[0,360,-90,90],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.]}
indices_list =('cdd','cwd','r10','r20','r95p','r99p','rx1day','rx5day','sdii','precptot','total_precip')

RCP85_globe_TimeSeries_dic ={'cdd':np.zeros((95,30)),'cwd':np.zeros((95,30)),'r10':np.zeros((95,30)),\
'r20':np.zeros((95,30)),'r95p':np.zeros((95,30)),'r99p':np.zeros((95,30)),'rx1day':np.zeros((95,30)),\
'rx5day':np.zeros((95,30)),'sdii':np.zeros((95,30)),'precptot':np.zeros((95,30)),'total_precip':np.zeros((95,30))}
RCP85_SI_TimeSeries_dic ={'cdd':np.zeros((95,30)),'cwd':np.zeros((95,30)),'r10':np.zeros((95,30)),\
'r20':np.zeros((95,30)),'r95p':np.zeros((95,30)),'r99p':np.zeros((95,30)),'rx1day':np.zeros((95,30)),\
'rx5day':np.zeros((95,30)),'sdii':np.zeros((95,30)),'precptot':np.zeros((95,30)),'total_precip':np.zeros((95,30))}
RCP85_CI_TimeSeries_dic ={'cdd':np.zeros((95,30)),'cwd':np.zeros((95,30)),'r10':np.zeros((95,30)),\
'r20':np.zeros((95,30)),'r95p':np.zeros((95,30)),'r99p':np.zeros((95,30)),'rx1day':np.zeros((95,30)),\
'rx5day':np.zeros((95,30)),'sdii':np.zeros((95,30)),'precptot':np.zeros((95,30)),'total_precip':np.zeros((95,30))}

oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/rcp85/*.nc'

files=sorted(glob.glob(input_path))

# globe_mean_TimeSeries=np.zeros((95,30));SI_mean_TimeSeries=np.zeros((95,30));CI_mean_TimeSeries=np.zeros((95,30));
layer =0
for file in files:
	nc_fid = nc4.Dataset(file,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	# time = nc_fid.variables['time'][:]
	# value = nc_fid.variables['total_precip'][:]/91
	# value = np.multiply(value,oceanmask)
	# _,_,globe_mean_3d = range_clip(rergion_dic['GLOBE'][0],rergion_dic['GLOBE'][1],rergion_dic['GLOBE'][2],rergion_dic['GLOBE'][3],lon,lat,value)
	# _,_,SI_mean_3d = range_clip(rergion_dic['SI'][0],rergion_dic['SI'][1],rergion_dic['SI'][2],rergion_dic['SI'][3],lon,lat,value)
	# _,_,CI_mean_3d = range_clip(rergion_dic['CI'][0],rergion_dic['CI'][1],rergion_dic['CI'][2],rergion_dic['CI'][3],lon,lat,value)	
	# globe_mean_TimeSeries[:,layer] = stats.nanmean(stats.nanmean(globe_mean_3d,axis=2),axis=1)
	# SI_mean_TimeSeries[:,layer] = stats.nanmean(stats.nanmean(SI_mean_3d,axis=2),axis=1)
	# CI_mean_TimeSeries[:,layer] = stats.nanmean(stats.nanmean(CI_mean_3d,axis=2),axis=1)
	# print CI_mean_TimeSeries
	for variable in indices_list:
		value = nc_fid.variables[variable][:]
		value = np.multiply(value,oceanmask)
		_,_,globe_3d = range_clip(rergion_dic['GLOBE'][0],rergion_dic['GLOBE'][1],rergion_dic['GLOBE'][2],rergion_dic['GLOBE'][3],lon,lat,value)
		# _,_,SI_3d = range_clip(rergion_dic['SI'][0],rergion_dic['SI'][1],rergion_dic['SI'][2],rergion_dic['SI'][3],lon,lat,value)
		_,_,CI_3d = range_clip(rergion_dic['CI'][0],rergion_dic['CI'][1],rergion_dic['CI'][2],rergion_dic['CI'][3],lon,lat,value)	
		# RCP85_globe_TimeSeries_dic[variable][:,layer] = stats.nanmean(stats.nanmean(globe_3d,axis=2),axis=1)
		# RCP85_SI_TimeSeries_dic[variable][:,layer] = stats.nanmean(stats.nanmean(SI_3d,axis=2),axis=1)
		RCP85_CI_TimeSeries_dic[variable][:,layer] = stats.nanmean(stats.nanmean(CI_3d,axis=2),axis=1)
	layer+=1
# plt.figure(1, facecolor='White',figsize=[20,5])
# spst = [2,11]   # subplot style
# sub_plot = 1
# for variable in indices_list:

	# ax=plt.subplot(spst[0],spst[1],sub_plot)
	# # print np.shape(movingaverage(time_series,5))
	# print SI_mean_TimeSeries
	# print SI_TimeSeries_dic[variable]
	# x = scale(stats.nanmean(globe_mean_TimeSeries,axis=1),feature_range=(0,1))
	# y = scale(stats.nanmean(globe_TimeSeries_dic[variable],axis=1),feature_range=(0,1))
	# ax.scatter(x,y,marker='s',color="red",label='Global')
	# x = scale(stats.nanmean(SI_mean_TimeSeries,axis=1),feature_range=(0,1))
	# y = scale(stats.nanmean(SI_TimeSeries_dic[variable],axis=1),feature_range=(0,1))	
	# ax.scatter(x,y,marker="o",color="blue",label='SOuth India')
	# x = scale(stats.nanmean(CI_mean_TimeSeries,axis=1),feature_range=(0,1))
	# y = scale(stats.nanmean(SI_TimeSeries_dic[variable],axis=1),feature_range=(0,1))		
	# ax.scatter(x,y,marker='v',color="green",label='Central India')
	# ax.set_title(variable.title(),fontsize=15)
	# sub_plot = sub_plot+1
		
		
# input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/fixa/*.nc'

# files=sorted(glob.glob(input_path))
# rergion_dic={'GLOBE':[0,360,-90,90],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.]}
# indices_list =('cdd','cwd','r10','r20','r95p','r99p','rx1day','rx5day','sdii','precptot')
# FIXA_globe_TimeSeries_dic ={'cdd':np.zeros((95,30)),'cwd':np.zeros((95,30)),'r10':np.zeros((95,30)),\
# 'r20':np.zeros((95,30)),'r95p':np.zeros((95,30)),'r99p':np.zeros((95,30)),'rx1day':np.zeros((95,30)),\
# 'rx5day':np.zeros((95,30)),'sdii':np.zeros((95,30)),'precptot':np.zeros((95,30)),'total_precip':np.zeros((95,30))}
# FIXA_SI_TimeSeries_dic ={'cdd':np.zeros((95,30)),'cwd':np.zeros((95,30)),'r10':np.zeros((95,30)),\
# 'r20':np.zeros((95,30)),'r95p':np.zeros((95,30)),'r99p':np.zeros((95,30)),'rx1day':np.zeros((95,30)),\
# 'rx5day':np.zeros((95,30)),'sdii':np.zeros((95,30)),'precptot':np.zeros((95,30)),'total_precip':np.zeros((95,30))}
# FIXA_CI_TimeSeries_dic ={'cdd':np.zeros((95,30)),'cwd':np.zeros((95,30)),'r10':np.zeros((95,30)),\
# 'r20':np.zeros((95,30)),'r95p':np.zeros((95,30)),'r99p':np.zeros((95,30)),'rx1day':np.zeros((95,30)),\
# 'rx5day':np.zeros((95,30)),'sdii':np.zeros((95,30)),'precptot':np.zeros((95,30)),'total_precip':np.zeros((95,30))}

# # globe_mean_TimeSeries=np.zeros((95,30));SI_mean_TimeSeries=np.zeros((95,30));CI_mean_TimeSeries=np.zeros((95,30));
# layer =0
# for file in files:
	# nc_fid = nc4.Dataset(file,mode='r')
	# lat = nc_fid.variables['lat'][:]
	# lon = nc_fid.variables['lon'][:]
	# # time = nc_fid.variables['time'][:]
	# # value = nc_fid.variables['total_precip'][:]/91
	# # mean_precp = np.multiply(value,oceanmask)
	# # _,_,globe_mean_3d = range_clip(rergion_dic['GLOBE'][0],rergion_dic['GLOBE'][1],rergion_dic['GLOBE'][2],rergion_dic['GLOBE'][3],lon,lat,value)
	# # _,_,SI_mean_3d = range_clip(rergion_dic['SI'][0],rergion_dic['SI'][1],rergion_dic['SI'][2],rergion_dic['SI'][3],lon,lat,value)
	# # _,_,CI_mean_3d = range_clip(rergion_dic['CI'][0],rergion_dic['CI'][1],rergion_dic['CI'][2],rergion_dic['CI'][3],lon,lat,value)	
	# # FIXA_globe_mean_TimeSeries[:,layer] = stats.nanmean(stats.nanmean(globe_mean_3d,axis=2),axis=1)
	# # FIXA_SI_mean_TimeSeries[:,layer] = stats.nanmean(stats.nanmean(SI_mean_3d,axis=2),axis=1)
	# # FIXA_CI_mean_TimeSeries[:,layer] = stats.nanmean(stats.nanmean(CI_mean_3d,axis=2),axis=1)
	# # print CI_mean_TimeSeries
	# for variable in indices_list:
		# value = nc_fid.variables[variable][:]
		# value = np.multiply(value,oceanmask)
		# _,_,globe_3d = range_clip(rergion_dic['GLOBE'][0],rergion_dic['GLOBE'][1],rergion_dic['GLOBE'][2],rergion_dic['GLOBE'][3],lon,lat,value)
		# _,_,SI_3d = range_clip(rergion_dic['SI'][0],rergion_dic['SI'][1],rergion_dic['SI'][2],rergion_dic['SI'][3],lon,lat,value)
		# _,_,CI_3d = range_clip(rergion_dic['CI'][0],rergion_dic['CI'][1],rergion_dic['CI'][2],rergion_dic['CI'][3],lon,lat,value)	
		# FIXA_globe_TimeSeries_dic[variable][:,layer] = stats.nanmean(stats.nanmean(globe_3d,axis=2),axis=1)
		# FIXA_SI_TimeSeries_dic[variable][:,layer] = stats.nanmean(stats.nanmean(SI_3d,axis=2),axis=1)
		# FIXA_CI_TimeSeries_dic[variable][:,layer] = stats.nanmean(stats.nanmean(CI_3d,axis=2),axis=1)
	# layer+=1
'''
plt.figure(1, facecolor='White',figsize=[20,5])
spst = [2,5]   # subplot style
sub_plot = 1
indices_list = ('cdd','cwd','r10','r20','r95p','r99p','rx1day','rx5day','sdii','precptot')
for variable in indices_list:
	ax=plt.subplot(spst[0],spst[1],sub_plot)
	# print np.shape(movingaverage(time_series,5))
	# print SI_mean_TimeSeries
	# print SI_TimeSeries_dic[variable]
	x = scale(stats.nanmean(RCP85_globe_TimeSeries_dic['total_precip'],axis=1),feature_range=(0,1))
	y = scale(stats.nanmean(RCP85_globe_TimeSeries_dic[variable],axis=1),feature_range=(0,1))
	cc = np.corrcoef(x,y)
	ax.scatter(x,y,marker='s',color="red",label='Gl '+str(round(cc[0,1],2)))
	x = scale(stats.nanmean(RCP85_SI_TimeSeries_dic['total_precip'],axis=1),feature_range=(0,1))
	y = scale(stats.nanmean(RCP85_SI_TimeSeries_dic[variable],axis=1),feature_range=(0,1))	
	cc = np.corrcoef(x,y)
	ax.scatter(x,y,marker="o",color="blue",label='SI '+str(round(cc[0,1],2)))
	x = scale(stats.nanmean(RCP85_CI_TimeSeries_dic['total_precip'],axis=1),feature_range=(0,1))
	y = scale(stats.nanmean(RCP85_CI_TimeSeries_dic[variable],axis=1),feature_range=(0,1))		
	cc = np.corrcoef(x,y)
	ax.scatter(x,y,marker='v',color="green",label='CI '+str(round(cc[0,1],2)))
	ax.set_title(variable.title(),fontsize=10)
	legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(0.25, 0.85),fontsize=15)	
	sub_plot = sub_plot+1
'''
# spst = [2,11]   # subplot style
# sub_plot = 12
# indices_list = ('cdd','cwd','r10','r20','r95p','r99p','rx1day','rx5day','sdii','precptot')
# for variable in indices_list:
	# ax=plt.subplot(spst[0],spst[1],sub_plot)
	# # print np.shape(movingaverage(time_series,5))
	# # print SI_mean_TimeSeries
	# # print SI_TimeSeries_dic[variable]
	# x = scale(stats.nanmean(RCP85_globe_TimeSeries_dic['total_precip']-FIXA_globe_TimeSeries_dic['total_precip'],axis=1),feature_range=(0,1))
	# y = scale(stats.nanmean(RCP85_globe_TimeSeries_dic[variable]-FIXA_globe_TimeSeries_dic[variable],axis=1),feature_range=(0,1))
	# cc = np.corrcoef(x,y)
	# ax.scatter(x,y,marker='s',color="red",label='Gl '+str(round(cc[0,1],2)))
	# x = scale(stats.nanmean(RCP85_SI_TimeSeries_dic['total_precip']-FIXA_SI_TimeSeries_dic['total_precip'],axis=1),feature_range=(0,1))
	# y = scale(stats.nanmean(RCP85_SI_TimeSeries_dic[variable]-FIXA_SI_TimeSeries_dic[variable],axis=1),feature_range=(0,1))	
	# cc = np.corrcoef(x,y)
	# ax.scatter(x,y,marker="o",color="blue",label='SI '+str(round(cc[0,1],2)))
	# x = scale(stats.nanmean(RCP85_CI_TimeSeries_dic['total_precip']-FIXA_CI_TimeSeries_dic['total_precip'],axis=1),feature_range=(0,1))
	# y = scale(stats.nanmean(RCP85_CI_TimeSeries_dic[variable]-FIXA_CI_TimeSeries_dic[variable],axis=1),feature_range=(0,1))	
	# cc = np.corrcoef(x,y)
	# ax.scatter(x,y,marker='v',color="green",label='CI '+str(round(cc[0,1],2)))
	# ax.set_title(variable.title(),fontsize=15)
	# legend = ax.legend(loc=10, shadow=True,bbox_to_anchor=(1, 0.9),fontsize=10)	
	# sub_plot = sub_plot+1	
	
# plt.show()	

"""
#########################################
#PCA
#########################################
indices_list = ('cdd','cwd','r10','r20','r95p','r99p','rx1day','rx5day','sdii','precptot')
array =np.zeros((95*30,10))
column=0
for variable in indices_list:
	cache= np.reshape(RCP85_CI_TimeSeries_dic[variable][:],(95*30))
	print cache
	array[:,column] =cache
	print column
	column=column+1

from matplotlib.mlab import PCA as mlabPCA
mlab_pca = mlabPCA(array)
print 'weight'
print mlab_pca.Wt
print 'the proportion of variance of each of the principal components'
print mlab_pca.fracs 
"""

