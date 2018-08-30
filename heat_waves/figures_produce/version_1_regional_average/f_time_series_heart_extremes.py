# -*- coding: utf-8 -*-
"""
This is to plot the time ecvelution of the intensity, duration and frequeency of 


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



###########READOUT DATA
CATEGORY = 1  # 0 FOR 90TH AND 1 FOR 99TH 
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/TD_TN_HW/Heat_TD_TN_HW_EastChina_India_rcp85.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
year_series = nc_fid.variables['year'][:]
HW_CEC_r =  stats.nanmean(nc_fid.variables['HW_CEC'][:],axis=0)[:,CATEGORY,:]
HW_SEC_r =  stats.nanmean(nc_fid.variables['HW_SEC'][:],axis=0)[:,CATEGORY,:]
HW_NID_r =  stats.nanmean(nc_fid.variables['HW_NID'][:],axis=0)[:,CATEGORY,:]
HW_SID_r =  stats.nanmean(nc_fid.variables['HW_SID'][:],axis=0)[:,CATEGORY,:]
HE_CEC_r =  stats.nanmean(nc_fid.variables['TD_CEC'][:],axis=0)[:,CATEGORY,:]
HE_SEC_r =  stats.nanmean(nc_fid.variables['TD_SEC'][:],axis=0)[:,CATEGORY,:]
HE_NID_r =  stats.nanmean(nc_fid.variables['TD_NID'][:],axis=0)[:,CATEGORY,:]
HE_SID_r =  stats.nanmean(nc_fid.variables['TD_SID'][:],axis=0)[:,CATEGORY,:]
TN_CEC_r =  stats.nanmean(nc_fid.variables['TN_CEC'][:],axis=0)[:,CATEGORY,:]
TN_SEC_r =  stats.nanmean(nc_fid.variables['TN_SEC'][:],axis=0)[:,CATEGORY,:]
TN_NID_r =  stats.nanmean(nc_fid.variables['TN_NID'][:],axis=0)[:,CATEGORY,:]
TN_SID_r =  stats.nanmean(nc_fid.variables['TN_SID'][:],axis=0)[:,CATEGORY,:]
nc_fid.close()
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/TD_TN_HW/Heat_TD_TN_HW_EastChina_India_fixa.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
year_series = nc_fid.variables['year'][:]
HW_CEC_f =  stats.nanmean(nc_fid.variables['HW_CEC'][:],axis=0)[:,CATEGORY,:]
HW_SEC_f =  stats.nanmean(nc_fid.variables['HW_SEC'][:],axis=0)[:,CATEGORY,:]
HW_NID_f =  stats.nanmean(nc_fid.variables['HW_NID'][:],axis=0)[:,CATEGORY,:]
HW_SID_f =  stats.nanmean(nc_fid.variables['HW_SID'][:],axis=0)[:,CATEGORY,:]
HE_CEC_f =  stats.nanmean(nc_fid.variables['TD_CEC'][:],axis=0)[:,CATEGORY,:]
HE_SEC_f =  stats.nanmean(nc_fid.variables['TD_SEC'][:],axis=0)[:,CATEGORY,:]
HE_NID_f =  stats.nanmean(nc_fid.variables['TD_NID'][:],axis=0)[:,CATEGORY,:]
HE_SID_f =  stats.nanmean(nc_fid.variables['TD_SID'][:],axis=0)[:,CATEGORY,:]
TN_CEC_f =  stats.nanmean(nc_fid.variables['TN_CEC'][:],axis=0)[:,CATEGORY,:]
TN_SEC_f =  stats.nanmean(nc_fid.variables['TN_SEC'][:],axis=0)[:,CATEGORY,:]
TN_NID_f =  stats.nanmean(nc_fid.variables['TN_NID'][:],axis=0)[:,CATEGORY,:]
TN_SID_f =  stats.nanmean(nc_fid.variables['TN_SID'][:],axis=0)[:,CATEGORY,:]
nc_fid.close()

####plotting
import matplotlib.pyplot as plt
fig = plt.figure(facecolor='White',figsize=[15,8]);plot_setup();

ax = plt.subplot(3,4,1);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(a)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('SEC',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax.annotate('ratio',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)

lins1 = ax.plot(year_series,HE_CEC_r[:,0],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_CEC_f[:,0],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,2);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(b)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('CEC',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
lins1 = ax.plot(year_series,HE_SEC_r[:,0],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_SEC_f[:,0],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,3);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(c)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('SID',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
lins1 = ax.plot(year_series,HE_SID_r[:,0],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_SID_f[:,0],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,4);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(d)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('NID',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
lins1 = ax.plot(year_series,HE_NID_r[:,0],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_NID_f[:,0],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,5);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(e)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('iNTENSITY',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)

lins1 = ax.plot(year_series,HE_SEC_r[:,2],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_SEC_f[:,2],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,6);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(f)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins1 = ax.plot(year_series,HE_CEC_r[:,2],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_CEC_f[:,2],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,7);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(g)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins1 = ax.plot(year_series,HE_SID_r[:,2],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_SID_f[:,2],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,8);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(h)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

lins1 = ax.plot(year_series,HE_NID_r[:,2],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_NID_f[:,2],'-',color="b",label='RCP8.5_FixA',linewidth=1)


ax = plt.subplot(3,4,9);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(i)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('Intensity',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)

lins1 = ax.plot(year_series,HE_SEC_r[:,4],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_SEC_f[:,4],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,10);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(j)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins1 = ax.plot(year_series,HE_CEC_r[:,4],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_CEC_f[:,4],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,11);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(k)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins1 = ax.plot(year_series,HE_SID_r[:,4],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_SID_f[:,4],'-',color="b",label='RCP8.5_FixA',linewidth=1)

ax = plt.subplot(3,4,12);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([2000,2100]);ax.set_xticklabels([]);
# ax.set_ylim([-0.4,6.0]); ax.set_yticks(np.arange(-0.4,6.1,1.6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axvspan(2031, 2050, alpha=0.2, color='gray');ax.axvspan(2081, 2100, alpha=0.2, color='gray');
ax.annotate('(l)',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

lins1 = ax.plot(year_series,HE_NID_r[:,4],'-',color="r",label='RCP8.5',linewidth=1)
lins1 = ax.plot(year_series,HE_NID_f[:,4],'-',color="b",label='RCP8.5_FixA',linewidth=1)

plt.show()




# precptot_fixa = nc_fid.variables['total_precip'][:]
# r95pr_fixa=np.divide(r95p_fixa,TP_fixa*92);r95pr_fixa[np.isinf(r95pr_fixa)]=np.nan
