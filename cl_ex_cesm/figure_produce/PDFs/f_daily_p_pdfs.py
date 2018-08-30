# -*- coding: utf-8 -*-
"""
This is to plot the daily precipitation PDFs of the following three time period
	1986-2005
	rcp8.5 2081-2100
	rcp8.5_FixA 2081-2100
	note all are JJA: 92 days each year
"""
import netCDF4 as nc4
import numpy as np
import scipy.io as spio
from scipy import stats
import matplotlib.pyplot as plt
from numpy import linspace
from scipy.stats.kde import gaussian_kde
import os; import site
import scipy.io as sio
from scipy import integrate

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}
oceanmask=spio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

"""

def data_netcdf(scenario,period):
	no_member =6
	if period =='6190': lent=30
	else: lent =20
	def data_reshape (scenario,time,data,lon,lat,period):

		if period=='6190':
			year_series = range(1961,1990)
		elif period=='8605':
			year_series = range(1986,2005)
		else :
			year_series = range(2081,2100)
		SA_CACHE = np.empty((lent,92,27,29));
		EA_CACHE = np.empty((lent,92,33,37));
		for iyear in year_series:
			layer_b = [layer for layer in range(len(time)) if time[layer] == iyear*10000+601][0]  
			layer_e = [layer for layer in range(len(time)) if time[layer] == iyear*10000+831][0]
			region = rergion_dic['SA']
			_,_,SA_CACHE[iyear-year_series[0],:,:,:] = range_clip(region[0],region[1],region[2],region[3],lon,lat,data[layer_b:layer_e+1,:,:])
			region = rergion_dic['EA']
			_,_,EA_CACHE[iyear-year_series[0],:,:,:] = range_clip(region[0],region[1],region[2],region[3],lon,lat,data[layer_b:layer_e+1,:,:]); 
		return  np.reshape(SA_CACHE,(lent*92,27,29)),np.reshape(EA_CACHE,(lent*92,33,37))
	
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/'+scenario+'/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + '/file_list.txt', "r")
	text_content = text_file.readlines()
	SA = np.empty((no_member,lent*92,27,29));EA = np.empty((no_member,lent*92,33,37));
	
	for en_no in range(0,no_member): 
		nc_f = text_content[en_no][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		time = nc_fid.variables['date'][:]
		att_value = np.multiply(nc_fid.variables['PRECT'],oceanmask)
		att_value [att_value<=0] = np.nan;
		att_value=att_value*24*60*60*1000  # from m/s to mm/day
		nc_fid.close()
		SA[en_no,:,:,:],EA[en_no,:,:,:] = data_reshape(scenario,time,att_value,lon,lat,period)
		print 'en_no '+str(en_no) + '  finished'
	return EA, SA
	
	
scenario = 'his'; period = '8605';EA_8605,SA_8605 = data_netcdf(scenario,period)
scenario = 'his'; period = '6190';EA_6190,SA_6190 = data_netcdf(scenario,period)
scenario = 'rcp85';period = '8100'; EA_rcp8100,SA_rcp8100 = data_netcdf(scenario,period)
scenario = 'fixa';period = '8100'; EA_fixa8100,SA_fixa8100 = data_netcdf(scenario,period)

file_out ='/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/'

sio.savemat(file_out+'EA_SA_DailyP_reshaped_10en_20yrs.mat', {'EA_8605':EA_8605,'EA_6190':EA_6190,'SA_6190':SA_6190,\
'SA_8605':SA_8605,'EA_rcp8100':EA_rcp8100,'SA_rcp8100':SA_rcp8100,'EA_fixa8100':EA_fixa8100,\
'SA_fixa8100':SA_fixa8100})

"""
def get_r95p_mean(SA_6190,EA_6190):
	SA_6190=np.reshape(np.nanmean(np.nanmean(SA_6190,axis=3),axis=2),(6,30*92))
	EA_6190=np.reshape(np.nanmean(np.nanmean(EA_6190,axis=3),axis=2),(6,30*92))
	EA = np.nanpercentile(EA_6190,95)
	SA = np.nanpercentile(SA_6190,95)
	return EA,SA
		
def kde_bins(data,bin_min=0,bin_max=30):
	data[data>200] = np.nan
	data = np.nanmean(np.nanmean(data,axis=3),axis=2)
	size_m =np.shape(data);size = size_m[0]*size_m[1]  #*size_m[2]*size_m[3]
	data = np.reshape(data,(size));
	data_m = data[~np.isinf(data)];
	data_mm = data_m[~np.isnan(data_m)];  # print np.nanmin(data_mm),np.nanmax(data_mm),
	bins = np.arange(bin_min,bin_max,.5)
	pdfs,_ = np.histogram(data_mm,bins=bins,density=True) #range=(bin_min,bin_max),
	pdfs[pdfs==0]=np.nan
	# kde = gaussian_kde( np.ma.masked_invalid(data_mm) );
	# pdfs = kde(bins); print np.sum(pdfs)
	mean = np.nanmean(data_mm); #print np.shape(bins[:-1]),np.shape(pdfs)
	return bins[:-1],pdfs,mean  #integrate.cumtrapz(hist, initial=0.)

def plot_pdfs(ax,x,y1,y2,y3,v1,v2,v3,rc9p,remove_mean =False):
	if remove_mean:
		ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
		ax.plot(x+0.25,y1,'-',color="k",linewidth=2,label = 'Baseline')
		ax.axvline(x=v1,color='k',linewidth = 2,ls='--');
		ax.plot(x+0.25-(v2-v1),y2,'-',color="r",linewidth=2,label = 'RCP8.5')
		# ax.axvline(x=v2-(v2-v1),color='r',linewidth = 2,ls='--');
		ax.plot(x+0.25-(v3-v1),y3,'-',color="b",linewidth=2,label = 'RCP8.5_FixA')
		# ax.axvline(x=v3-(v3-v1),color='b',linewidth = 2,ls='--');
		ax.axvline(x=10,color='g',linewidth = 2,ls='-');
		ax.axvline(x=rc9p,color='k',linewidth = 2,ls='-');
	else:
		ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
		ax.plot(x+0.25,y1,'-',color="k",linewidth=2,label = 'Baseline')
		ax.axvline(x=v1,color='k',linewidth = 2,ls='--');
		ax.plot(x+0.25,y2,'-',color="r",linewidth=2,label = 'RCP8.5')
		ax.axvline(x=v2,color='r',linewidth = 2,ls='--');
		ax.plot(x+0.25,y3,'-',color="b",linewidth=2,label = 'RCP8.5_FixA')
		ax.axvline(x=v3,color='b',linewidth = 2,ls='--');
		ax.axvline(x=10,color='g',linewidth = 2,ls='-');
		ax.axvline(x=rc9p,color='k',linewidth = 2,ls='-');		
	return ax
	
mat = sio.loadmat('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/EA_SA_DailyP_reshaped_10en_20yrs.mat')
SA_6190 = mat['SA_6190'][:];EA_6190 = mat['EA_6190'][:];	
EA95,SA95 = get_r95p_mean(SA_6190,EA_6190)
fig = plt.figure(figsize=(12, 6), facecolor='White');plot_setup();pad=5
bin_min =0;bin_max=24
######SA
region = 'SA'
ax = plt.subplot(2,2,1);
ax.annotate('(a) South Asia',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
SA_8605 = mat['SA_8605'][:];SA_rcp8100 = mat['SA_rcp8100'][:];SA_fixa8100 = mat['SA_fixa8100'][:];	
SA_bins,SA_8605_PDFs,SA_8605_mean = kde_bins(SA_8605);print 'bug', SA_8605_mean
_,SA_rcp8100_PDFs,SA_rcp8100_mean = kde_bins(SA_rcp8100)
_,SA_fixa8100_PDFs,SA_fixa8100_mean = kde_bins(SA_fixa8100); #print SA_bins,SA_fixa8100_PDFs
ax = plot_pdfs(ax,SA_bins,SA_8605_PDFs,SA_rcp8100_PDFs,SA_fixa8100_PDFs,SA_8605_mean,SA_rcp8100_mean,SA_fixa8100_mean,SA95)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/6));
ax.set_ylim([0,0.2]);ax.set_yticks(np.arange(0,.21,.05));
legend1 = ax.legend(shadow=False,ncol=1,loc='upper left')	 
legend1.get_frame().set_facecolor('None');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.3)

axins =  inset_axes(ax, width="30%",height="70%", loc=1) # zoom = 6
x1, x2, y1, y2 = 16, 24, 0, 0.015
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# plt.xticks(visible=False)
# plt.yticks(visible=False)
axins = plot_pdfs(axins,SA_bins,SA_8605_PDFs,SA_rcp8100_PDFs,SA_fixa8100_PDFs,SA_8605_mean,SA_rcp8100_mean,SA_fixa8100_mean,SA95)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.4")

ax = plt.subplot(2,2,2);
ax.annotate('(b) South Asia (mean rmoved)',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
SA_8605 = mat['SA_8605'][:];SA_rcp8100 = mat['SA_rcp8100'][:];SA_fixa8100 = mat['SA_fixa8100'][:];	
SA_bins,SA_8605_PDFs,SA_8605_mean = kde_bins(SA_8605)
_,SA_rcp8100_PDFs,SA_rcp8100_mean = kde_bins(SA_rcp8100)
_,SA_fixa8100_PDFs,SA_fixa8100_mean = kde_bins(SA_fixa8100); #print SA_bins,SA_fixa8100_PDFs
# print SA_rcp8100_PDFs
# print SA_fixa8100_PDFs
ax = plot_pdfs(ax,SA_bins,SA_8605_PDFs,SA_rcp8100_PDFs,SA_fixa8100_PDFs,SA_8605_mean,SA_rcp8100_mean,SA_fixa8100_mean,SA95,remove_mean=True)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/6));
ax.set_ylim([0,0.2]);ax.set_yticks(np.arange(0,.21,.05));
legend1 = ax.legend(shadow=False,ncol=1,loc='upper left')	 
legend1.get_frame().set_facecolor('None');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.3)

axins =  inset_axes(ax, width="30%",height="70%", loc=1) # zoom = 6
x1, x2, y1, y2 = 16, 24, 0, 0.015
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# plt.xticks(visible=False)
# plt.yticks(visible=False)
axins = plot_pdfs(axins,SA_bins,SA_8605_PDFs,SA_rcp8100_PDFs,SA_fixa8100_PDFs,SA_8605_mean,SA_rcp8100_mean,SA_fixa8100_mean,SA95,remove_mean=True)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.4")


######EA
ax = plt.subplot(2,2,3);bin_min =0;bin_max=18
ax.annotate('(c) East Asia',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
EA_8605 = mat['EA_8605'][:];EA_rcp8100 = mat['EA_rcp8100'][:];EA_fixa8100 = mat['EA_fixa8100'][:];
EA_bins,EA_8605_PDFs,EA_8605_mean = kde_bins(EA_8605)
_,EA_rcp8100_PDFs,EA_rcp8100_mean = kde_bins(EA_rcp8100)
_,EA_fixa8100_PDFs,EA_fixa8100_mean = kde_bins(EA_fixa8100)

EA_rcp = np.concatenate((EA_rcp8100_PDFs[0:12],EA_rcp8100_PDFs[12:30]),axis=0)
EA_fixa = np.concatenate((EA_fixa8100_PDFs[0:12],EA_fixa8100_PDFs[12:30]),axis=0)

ax = plot_pdfs(ax,EA_bins,EA_8605_PDFs,EA_rcp8100_PDFs,EA_fixa8100_PDFs,EA_8605_mean,EA_rcp8100_mean,EA_fixa8100_mean,EA95)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/6));
ax.set_ylim([0,0.30]);ax.set_yticks(np.arange(0,.31,.10));

axins =  inset_axes(ax, width="30%",height="70%", loc=1) # zoom = 6
x1, x2, y1, y2 = 12, 18, 0, 0.006
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins = plot_pdfs(axins,EA_bins,EA_8605_PDFs,EA_rcp8100_PDFs,EA_fixa8100_PDFs,EA_8605_mean,EA_rcp8100_mean,EA_fixa8100_mean,EA95)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.4")



ax = plt.subplot(2,2,4);bin_min =0;bin_max=18
ax.annotate('(d) East Asia (mean rmoved)',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
EA_8605 = mat['EA_8605'][:];EA_rcp8100 = mat['EA_rcp8100'][:];EA_fixa8100 = mat['EA_fixa8100'][:];
EA_bins,EA_8605_PDFs,EA_8605_mean = kde_bins(EA_8605)
_,EA_rcp8100_PDFs,EA_rcp8100_mean = kde_bins(EA_rcp8100)
_,EA_fixa8100_PDFs,EA_fixa8100_mean = kde_bins(EA_fixa8100)

EA_rcp = np.concatenate((EA_rcp8100_PDFs[0:12],EA_rcp8100_PDFs[12:30]),axis=0)
EA_fixa = np.concatenate((EA_fixa8100_PDFs[0:12],EA_fixa8100_PDFs[12:30]),axis=0)

ax = plot_pdfs(ax,EA_bins,EA_8605_PDFs,EA_rcp8100_PDFs,EA_fixa8100_PDFs,EA_8605_mean,EA_rcp8100_mean,EA_fixa8100_mean,EA95,remove_mean=True)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/6));
ax.set_ylim([0,0.30]);ax.set_yticks(np.arange(0,.31,.10));

axins =  inset_axes(ax, width="30%",height="70%", loc=1) # zoom = 6
x1, x2, y1, y2 = 12, 18, 0, 0.006
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins = plot_pdfs(axins,EA_bins,EA_8605_PDFs,EA_rcp8100_PDFs,EA_fixa8100_PDFs,EA_8605_mean,EA_rcp8100_mean,EA_fixa8100_mean,EA95,remove_mean=True)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.4")



			
plt.subplots_adjust(left=0.07, bottom=0.06, right=0.97, top=0.92, wspace=0.1, hspace=0.30);
plt.savefig('PDFs_daily_p_his_and_81_10.png', format='png', dpi=1000)

