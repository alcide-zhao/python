# -*- coding: utf-8 -*-
"""
Created on Feb 01 2017
This is to plot the RCP8.5 and RCP8.5_FIXA means and model variabilities of following fields
(1)TS
(2)TOtal precipitation
(3) RX5DAY
(4) r95p
(5) CWD
(6) SDII
"""
import scipy
import numpy as np
from scipy import stats
import os; import site; 
import netCDF4 as nc4
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


region_dic={'ASIA':[60,155,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}
region_key = 'EA';region = region_dic[region_key]

# ocean land mask
import scipy.io as sio
ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan

def time_series_range(variable, ref):
	from numpy import linspace
	from scipy.stats.kde import gaussian_kde
	size = np.shape(variable)[1]
	variable_median = stats.nanmean(variable,axis = 0)-ref
	# variable [variable <=0]= np.nan; mask = ~np.isnan(variable); variable = variable[mask]
	variable_upper = np.empty((size));variable_upper[:]= np.nan
	variable_loweer = np.empty((size));variable_loweer[:]= np.nan
	for itime in range(size):	
		data = variable[:,itime];
		bins = linspace(min(data), max(data), 50); 
		gkde=gaussian_kde(data)
		kdepdf = gkde.evaluate(bins)
		kde = 100*kdepdf/np.sum(kdepdf)
		variable_upper[itime] = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]-ref
		variable_loweer[itime] = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == 10 or (sum(kde[0:i])<=10 and sum(kde[0:i+1])>=10))][0]-ref
	return variable_median,variable_upper,variable_loweer
#############################
# climatology references
#############################
file_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
file = file_path+'TS_P_850UV_model_eval_1971_2000.mat'
data = sio.loadmat(file)
TS_JJA_mean_CESM = np.multiply(data['TS_JJA_mean_CESM'][:],ocean_mask_CESM);   
lon_CESM = data['lon_CESM'][0,:];lat_CESM = data['lat_CESM'][0,:];

# mean surface temperature and total precipitation climatology
_,_,TS_JJA_mean = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TS_JJA_mean_CESM)
TS_ref = stats.nanmean(stats.nanmean(stats.nanmean(TS_JJA_mean,axis=2),axis=1),axis=0)

# precipitation indices climatology : reference time period is 1971-2000 
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
TP = np.multiply(nc_fid.variables['total_precip'][66:86,:,:],ocean_mask_CESM)
rx5day = np.multiply(nc_fid.variables['rx5day'][66:86,:,:],ocean_mask_CESM)
r95p = np.multiply(nc_fid.variables['r95p'][66:86,:,:],ocean_mask_CESM)
sdii =np.multiply(nc_fid.variables['r10'][66:86,:,:],ocean_mask_CESM)
cwd = np.multiply(nc_fid.variables['cwd'][66:86,:,:],ocean_mask_CESM)
nc_fid.close()
_,_,TP_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TP)
TP_ref = stats.nanmean(stats.nanmean(stats.nanmean(TP_range,axis=2),axis=1),axis=0)/92
_,_,rx5day_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,rx5day)
rx5day_ref = stats.nanmean(stats.nanmean(stats.nanmean(rx5day_range,axis=2),axis=1),axis=0)
_,_,r95p_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,r95p)
r95p_ref = stats.nanmean(stats.nanmean(stats.nanmean(r95p_range,axis=2),axis=1),axis=0)
_,_,sdii_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,sdii)
sdii_ref = stats.nanmean(stats.nanmean(stats.nanmean(sdii_range,axis=2),axis=1),axis=0)
_,_,cwd_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,cwd)
cwd_ref = stats.nanmean(stats.nanmean(stats.nanmean(cwd_range,axis=2),axis=1),axis=0)
# _,_,precptot_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TP)
# precptot_ref = stats.nanmean(stats.nanmean(stats.nanmean(precptot_range,axis=2),axis=1),axis=0)
# r95pr_ref = r95p_ref/(TP_ref*92)

# TS
file = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/ensenmble_TS_200602_210101.nc'
nc_fid = nc4.Dataset(file,mode='r')
lon = nc_fid.variables['lon'][:]
lat = nc_fid.variables['lat'][:]
TS_fixa = nc_fid.variables['TS_fixa'][:]
TS_rcp85 = nc_fid.variables['TS_rcp85'][:]
nc_fid.close()
# mask off the TS data as this has been done when preprocessing
TS_fixa = np.multiply(TS_fixa, ocean_mask_CESM);
TS_rcp85 = np.multiply(TS_rcp85, ocean_mask_CESM);

# region clip
_,_,TS_rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,TS_rcp85)
_,_,TS_fixa_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,TS_fixa)

# plt.imshow(TS_rcp85_clipped[0,0,:,:],origin = 'lower')
# plt.show()
TS_rcp85_spatial_mean = stats.nanmean(stats.nanmean(TS_rcp85_clipped,axis = 3),axis=2)
TS_fixa_spatial_mean = stats.nanmean(stats.nanmean(TS_fixa_clipped,axis = 3),axis=2)

# precipitation 
TS_mean_rcp85_TiSe,TS_upper_rcp85_TiSe,TS_lower_rcp85_TiSe = time_series_range(TS_rcp85_spatial_mean,TS_ref  )
TS_mean_fixa_TiSe,TS_upper_fixa_TiSe,TS_lower_fixa_TiSe  = time_series_range(TS_fixa_spatial_mean,TS_ref )


print 'TS',np.nanmean(TS_mean_rcp85_TiSe[-21:-1]), np.nanmean(TS_mean_fixa_TiSe[-21:-1])
# RCP85
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/TiSe_PEI_rcp85_'+region_key+'_2006_2100.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
year_series = nc_fid.variables['time'][:]
rx5day_rcp85 = nc_fid.variables['rx5day'][:]
r95p_rcp85= nc_fid.variables['r95p'][:]
sdii_rcp85 = nc_fid.variables['r10'][:]
cwd_rcp85 = nc_fid.variables['cwd'][:]
TP_rcp85 = nc_fid.variables['total_precip'][:]/92
# precptot_rcp85 = nc_fid.variables['total_precip'][:]
# r95pr_rcp85= np.divide(r95p_rcp85,TP_rcp85*92);r95pr_rcp85[np.isinf(r95pr_rcp85)]=np.nan
# nc_fid.close()

TP_mean_rcp85_TiSe,TP_upper_rcp85_TiSe,TP_lower_rcp85_TiSe = time_series_range(TP_rcp85 ,TP_ref)
rx5day_mean_rcp85_TiSe,rx5day_upper_rcp85_TiSe,rx5day_lower_rcp85_TiSe = time_series_range(rx5day_rcp85,rx5day_ref)
cwd_mean_rcp85_TiSe,cwd_upper_rcp85_TiSe,cwd_lower_rcp85_TiSe = time_series_range(cwd_rcp85,cwd_ref )
sdii_mean_rcp85_TiSe,sdii_upper_rcp85_TiSe,sdii_lower_rcp85_TiSe = time_series_range(sdii_rcp85,sdii_ref )
r95pr_mean_rcp85_TiSe,r95pr_upper_rcp85_TiSe,r95pr_lower_rcp85_TiSe = time_series_range(r95p_rcp85,r95p_ref )

# FIXA
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/TiSe_PEI_fixa_'+region_key+'_2006_2100.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
year_series = nc_fid.variables['time'][:]
rx5day_fixa = nc_fid.variables['rx5day'][:]
r95p_fixa  = nc_fid.variables['r95p'][:]
sdii_fixa  = nc_fid.variables['r10'][:]
cwd_fixa  = nc_fid.variables['cwd'][:]
TP_fixa  = nc_fid.variables['total_precip'][:]/92
# precptot_fixa = nc_fid.variables['total_precip'][:]
# r95pr_fixa=np.divide(r95p_fixa,TP_fixa*92);r95pr_fixa[np.isinf(r95pr_fixa)]=np.nan
nc_fid.close()

TP_mean_fixa_TiSe,TP_upper_fixa_TiSe,TP_lower_fixa_TiSe = time_series_range(TP_fixa ,TP_ref )
rx5day_mean_fixa_TiSe,rx5day_upper_fixa_TiSe,rx5day_lower_fixa_TiSe = time_series_range(rx5day_fixa ,rx5day_ref )
cwd_mean_fixa_TiSe,cwd_upper_fixa_TiSe,cwd_lower_fixa_TiSe = time_series_range(cwd_fixa,cwd_ref  )
sdii_mean_fixa_TiSe,sdii_upper_fixa_TiSe,sdii_lower_fixa_TiSe = time_series_range(sdii_fixa,sdii_ref)
# r95p_mean_fixa_TiSe,r95p_upper_fixa_TiSe,r95p_lower_fixa_TiSe = time_series_range(r95p_fixa,r95p_ref)) #
r95pr_mean_fixa_TiSe,r95pr_upper_fixa_TiSe,r95pr_lower_fixa_TiSe = time_series_range(r95p_fixa,r95p_ref) #

print 'TP',np.nanmean(TP_mean_rcp85_TiSe[-21:-1]), np.nanmean(TP_mean_fixa_TiSe[-21:-1])

# T_test
p_threshold = 0.1
# from scipy.stats import ttest_ind as test

##The Mann-Whistney U-test
def mannwhitneyu_test(data1,data2):
	times = np.shape(data2)[1]
	p_value = np.empty((times));p_value[:]=np.nan
	from scipy.stats import mannwhitneyu as test
	for itime in range(times):
		cache1 = data1[:,itime]
		cache2 = data2[:,itime]
		_,p_value[itime] = test(cache1,cache2);
	p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
	return p_value

p_TS = mannwhitneyu_test(TS_rcp85_spatial_mean,TS_fixa_spatial_mean );
p_TP = mannwhitneyu_test(TP_rcp85,TP_fixa );
p_rx5day = mannwhitneyu_test(rx5day_rcp85,rx5day_fixa); 
p_r95pr = mannwhitneyu_test(r95p_rcp85,r95p_fixa); 
p_cwd = mannwhitneyu_test(cwd_rcp85,cwd_fixa); 
p_sdii = mannwhitneyu_test(sdii_rcp85,sdii_fixa ); 


####plotting
import matplotlib.pyplot as plt
fig = plt.figure(facecolor='White',figsize=[6,6]);plot_setup();

ax = plt.subplot(3,2,1);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))#EA
# ax.set_xlim([2000,2100]);ax.set_ylim([-1,5.5]);  #SA
ax.set_xticklabels([]);#ax.locator_params(axis='y',nbins=6);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(a) SAT (K)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
lins1 = ax.plot(year_series,TS_mean_rcp85_TiSe,'-',color="r",label='RCP8.5',linewidth=1)
lins2 = ax.fill_between(year_series,TS_upper_rcp85_TiSe,TS_lower_rcp85_TiSe, where=TS_upper_rcp85_TiSe>=TS_lower_rcp85_TiSe, color="r", alpha = 0.1) 
lins3 = ax.plot(year_series,TS_mean_fixa_TiSe,'-',color="b",label='RCP8.5_FixA',linewidth=1)
lins4 = ax.fill_between(year_series,TS_upper_fixa_TiSe,TS_lower_fixa_TiSe, where=TS_upper_fixa_TiSe>=TS_lower_fixa_TiSe, color="b", alpha = 0.1) 
lins5 = ax.plot(year_series,p_TS*-0.0,'*',color="g",markersize=3, markeredgecolor='none')
# lins5 = ax.plot(year_series,p_TS*(-0.5),'*',color="g",markersize=3, markeredgecolor='none')


legend = ax.legend(shadow=False,ncol=1,loc ='upper left')	 
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)

ax = plt.subplot(3,2,2);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_ylim([-1,2]);ax.set_yticks(np.arange(-1,2.1,0.75))#EA
# ax.set_xlim([2000,2100]);ax.set_ylim([-100,200]);
ax.set_xticklabels([]);#ax.locator_params(axis='y',nbins=6);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate(r'(b) PMEAN $\mathrm{\mathsf{(mm\/day^{-1})}}$',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
lins1 = ax.plot(year_series,TP_mean_rcp85_TiSe,'-',color="r",label='RCP8.5',linewidth=1)
lins2 = ax.fill_between(year_series,TP_upper_rcp85_TiSe,TP_lower_rcp85_TiSe, where=TP_upper_rcp85_TiSe>=TP_lower_rcp85_TiSe, color="r", alpha = 0.1) 
lins3 = ax.plot(year_series,TP_mean_fixa_TiSe,'-',color="b",label='RCP8.5_FixA',linewidth=1)
lins4 = ax.fill_between(year_series,TP_upper_fixa_TiSe,TP_lower_fixa_TiSe, where=TP_upper_fixa_TiSe>=TP_lower_fixa_TiSe, color="b", alpha = 0.1) 
lins5 = ax.plot(year_series,p_TP*-0.85,'*',color="g",markersize=3, markeredgecolor='none')
# lins5 = ax.plot(year_series,p_TP*-80,'*',color="g",markersize=3, markeredgecolor='none')


ax = plt.subplot(3,2,3);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_ylim([-20,80]);ax.set_yticks(np.arange(-20,80.1,25))#EA
# ax.set_xlim([2000,2100]);ax.set_ylim([-30,100]);
ax.set_xticklabels([]);#ax.locator_params(axis='y',nbins=6);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(c) RX5DAY (mm)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
lins1 = ax.plot(year_series,rx5day_mean_rcp85_TiSe,'-',color="r",label='RCP8.5',linewidth=1)
lins2 = ax.fill_between(year_series,rx5day_upper_rcp85_TiSe,rx5day_lower_rcp85_TiSe, where=rx5day_upper_rcp85_TiSe>=rx5day_lower_rcp85_TiSe, color="r", alpha = 0.1) 
lins3 = ax.plot(year_series,rx5day_mean_fixa_TiSe,'-',color="b",label='RCP8.5_FixA',linewidth=1)
lins4 = ax.fill_between(year_series,rx5day_upper_fixa_TiSe,rx5day_lower_fixa_TiSe, where=rx5day_upper_fixa_TiSe>=rx5day_lower_fixa_TiSe, color="b", alpha = 0.1) 
lins5 = ax.plot(year_series,p_rx5day*-12,'*',color="g",markersize=3, markeredgecolor='none')
# lins5 = ax.plot(year_series,p_rx5day*-20,'*',color="g",markersize=3, markeredgecolor='none')

ax = plt.subplot(3,2,4);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
# ax.set_xlim([2000,2100]);ax.set_ylim([-1.5,3]);
ax.set_xlim([2000,2100]);ax.set_ylim([-2,5.2]);ax.set_yticks(np.arange(-2.0,5.3,1.8))#EA
ax.set_xticklabels([]);#ax.locator_params(axis='y',nbins=6);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(d) R10 (days)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
lins1 = ax.plot(year_series,sdii_mean_rcp85_TiSe,'-',color="r",label='RCP8.5',linewidth=1)
lins2 = ax.fill_between(year_series,sdii_upper_rcp85_TiSe,sdii_lower_rcp85_TiSe, where=sdii_upper_rcp85_TiSe>=sdii_lower_rcp85_TiSe, color="r", alpha = 0.1) 
lins3 = ax.plot(year_series,sdii_mean_fixa_TiSe,'-',color="b",label='RCP8.5_FixA',linewidth=1)
lins4 = ax.fill_between(year_series,sdii_upper_fixa_TiSe,sdii_lower_fixa_TiSe, where=sdii_upper_fixa_TiSe>=sdii_lower_fixa_TiSe, color="b", alpha = 0.1) 
lins5 = ax.plot(year_series,p_sdii*-1.55,'*',color="g",markersize=3, markeredgecolor='none')
# lins5 = ax.plot(year_series,p_sdii*-0.7,'*',color="g",markersize=3, markeredgecolor='none')

ax = plt.subplot(3,2,5);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_ylim([-8,4]);ax.set_yticks(np.arange(-8,4.1,3))#EA
# ax.set_xlim([2000,2100]);ax.set_ylim([-12,5]);ax.locator_params(axis='y',nbins=6);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(e) CWD (days)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
lins1 = ax.plot(year_series,cwd_mean_rcp85_TiSe,'-',color="r",label='RCP8.5',linewidth=1)
lins2 = ax.fill_between(year_series,cwd_upper_rcp85_TiSe,cwd_lower_rcp85_TiSe, where=cwd_upper_rcp85_TiSe>=cwd_lower_rcp85_TiSe, color="r", alpha = 0.1) 
lins3 = ax.plot(year_series,cwd_mean_fixa_TiSe,'-',color="b",label='RCP8.5_FixA',linewidth=1)
lins4 = ax.fill_between(year_series,cwd_upper_fixa_TiSe,cwd_lower_fixa_TiSe, where=cwd_upper_fixa_TiSe>=cwd_lower_fixa_TiSe, color="b", alpha = 0.1) 
lins5 = ax.plot(year_series,p_cwd*-7,'*',color="g",markersize=3, markeredgecolor='none')
# lins5 = ax.plot(year_series,p_cwd*-11,'*',color="g",markersize=3, markeredgecolor='none')

ax = plt.subplot(3,2,6);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_ylim([-20,100]);ax.set_yticks(np.arange(-20,100.1,30))#EA
# ax.set_xlim([2000,2100]);ax.set_ylim([-0.05,0.15]); #ax.locator_params(axis='y',nbins=6);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(f) R95P (mm)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
lins1 = ax.plot(year_series,r95pr_mean_rcp85_TiSe,'-',color="r",label='RCP8.5',linewidth=1)
lins2 = ax.fill_between(year_series,r95pr_upper_rcp85_TiSe,r95pr_lower_rcp85_TiSe, where=r95pr_upper_rcp85_TiSe>=r95pr_lower_rcp85_TiSe, color="r", alpha = 0.1) 
lins3 = ax.plot(year_series,r95pr_mean_fixa_TiSe,'-',color="b",label='RCP8.5_FixA',linewidth=1)
lins4 = ax.fill_between(year_series,r95pr_upper_fixa_TiSe,r95pr_lower_fixa_TiSe, where=r95pr_upper_fixa_TiSe>=r95pr_lower_fixa_TiSe, color="b", alpha = 0.1) 
lins5 = ax.plot(year_series,p_r95pr*-12.5,'*',color="g",markersize=3, markeredgecolor='none') #-9
# lins5 = ax.plot(year_series,p_r95pr*-0.04,'*',color="g",markersize=3, markeredgecolor='none')

plt.subplots_adjust(left=0.06, bottom=0.05, right=0.92, top=0.95, wspace=0.20, hspace=0.20); 
plt.savefig(region_key+'.png', format='png', dpi=1000)
# plt.show()
