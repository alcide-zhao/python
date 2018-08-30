# -*- coding: utf-8 -*-
"""
Created on Feb 9 2017
This is to show the relationships betwen India surface Temperature and precipitation evelutions
@author: Alcide.Zhao
"""
import numpy as np
import netCDF4 as nc4
from scipy import stats

import os; import site
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


#######################################
# 0.1 varibale definition
#######################################
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}
def time_series_range(variable, ref,region):
	_,_,variable_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,variable)
	variable_mean = 100*(stats.nanmean(stats.nanmean(variable_clipped,axis = 2),axis=1)-ref)/ref
	return variable_mean

def time_series_absolute(variable, ref,region):
	_,_,variable_clipped = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,variable)
	variable_mean = stats.nanmean(stats.nanmean(variable_clipped,axis = 2),axis=1)-ref
	return variable_mean
	
# ocean land mask
import scipy.io as sio
ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask_CESM.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan

#############################
# climatology references
#############################
file_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
file = file_path+'TS_P_850UV_model_eval_1971_2000.mat'
data = sio.loadmat(file)
TS_JJA_mean_CESM =np.multiply(data['TS_JJA_mean_CESM'][:],ocean_mask_CESM);   
lon_CESM = data['lon_CESM'][0,:];lat_CESM = data['lat_CESM'][0,:];

# precipitation indices climatology : reference time period is 1971-2000 
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/historical_ensumble_mean_PEI_global_1920_2005.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
PMEAN = np.multiply(nc_fid.variables['total_precip'][66:86,:,:],ocean_mask_CESM)/92
rx5day = np.multiply(nc_fid.variables['rx5day'][66:86,:,:],ocean_mask_CESM)
r95p = np.multiply(nc_fid.variables['r95p'][66:86,:,:],ocean_mask_CESM)
r30 =np.multiply(nc_fid.variables['rnm'][66:86,:,:],ocean_mask_CESM)
cwd = np.multiply(nc_fid.variables['cwd'][66:86,:,:],ocean_mask_CESM)
nc_fid.close()


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

# RCP85
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/ensumble_mean_PEI_global_2006_2100_rcp85.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
year_series = nc_fid.variables['time'][:]
rx5day_rcp85 = np.multiply(nc_fid.variables['rx5day'][:],ocean_mask_CESM)
r95p_rcp85= np.multiply(nc_fid.variables['r95p'][:],ocean_mask_CESM)
r30_rcp85 = np.multiply(nc_fid.variables['rnm'][:],ocean_mask_CESM)
cwd_rcp85 = np.multiply(nc_fid.variables['cwd'][:],ocean_mask_CESM)
PMEAN_rcp85 = np.multiply(nc_fid.variables['total_precip'][:],ocean_mask_CESM)/92

nc_fid.close()

# FIXA
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
year_series = nc_fid.variables['time'][:]
rx5day_fixa = np.multiply(nc_fid.variables['rx5day'][:],ocean_mask_CESM)
r95p_fixa  = np.multiply(nc_fid.variables['r95p'][:],ocean_mask_CESM)
r30_fixa  =np.multiply( nc_fid.variables['rnm'][:],ocean_mask_CESM)
cwd_fixa  =np.multiply( nc_fid.variables['cwd'][:],ocean_mask_CESM)
PMEAN_fixa  = np.multiply(nc_fid.variables['total_precip'][:],ocean_mask_CESM)/92

nc_fid.close()

#############################
# EA
#############################
region = rergion_dic['EA'][:]
# mean surface temperature and total precipitation climatology
_,_,TS_JJA_mean = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TS_JJA_mean_CESM)
TS_ref = stats.nanmean(stats.nanmean(stats.nanmean(TS_JJA_mean,axis=2),axis=1),axis=0)
_,_,PMEAN_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,PMEAN)
# plt.imshow(stats.nanmean(PMEAN_range,axis=0))
# plt.show()
PMEAN_ref = stats.nanmean(stats.nanmean(stats.nanmean(PMEAN_range,axis=2),axis=1),axis=0)
_,_,rx5day_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,rx5day)
rx5day_ref = stats.nanmean(stats.nanmean(stats.nanmean(rx5day_range,axis=2),axis=1),axis=0)
_,_,r95p_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,r95p)
r95p_ref = stats.nanmean(stats.nanmean(stats.nanmean(r95p_range,axis=2),axis=1),axis=0)
_,_,r30_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,r30)
r30_ref = stats.nanmean(stats.nanmean(stats.nanmean(r30_range,axis=2),axis=1),axis=0)
_,_,cwd_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,cwd)
cwd_ref = stats.nanmean(stats.nanmean(stats.nanmean(cwd_range,axis=2),axis=1),axis=0)
print PMEAN_ref

# region clip
_,_,TS_rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,TS_rcp85)-TS_ref
_,_,TS_fixa_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,TS_fixa)-TS_ref
TS_rcp85_TiSe = stats.nanmean(stats.nanmean(stats.nanmean(TS_rcp85_clipped,axis = 3),axis=2),axis=0)
TS_fixa_TiSe = stats.nanmean(stats.nanmean(stats.nanmean(TS_fixa_clipped,axis = 3),axis=2),axis=0)

PMEAN_rcp85_TiSe= time_series_range(PMEAN_rcp85 ,PMEAN_ref,region)
rx5day_rcp85_TiSe= time_series_range(rx5day_rcp85,rx5day_ref,region)
cwd_rcp85_TiSe= time_series_range(cwd_rcp85,cwd_ref,region )
r30_rcp85_TiSe= time_series_range(r30_rcp85,r30_ref,region )
r95p_rcp85_TiSe= time_series_range(r95p_rcp85,r95p_ref,region )

PMEAN_fixa_TiSe = time_series_range(PMEAN_fixa ,PMEAN_ref,region)
rx5day_fixa_TiSe= time_series_range(rx5day_fixa ,rx5day_ref,region)
cwd_fixa_TiSe= time_series_range(cwd_fixa,cwd_ref,region)
r30_fixa_TiSe= time_series_range(r30_fixa,r30_ref,region)
r95p_fixa_TiSe= time_series_range(r95p_fixa,r95p_ref,region) #

PMEAN_fixa_abso = time_series_absolute(PMEAN_fixa ,PMEAN_ref,region)
PMEAN_rcp85_abso= time_series_absolute(PMEAN_rcp85 ,PMEAN_ref,region)
#############################
# plotting
#############################
import matplotlib.pyplot as plt
fig = plt.figure(facecolor='White',figsize=[6,9]);plot_setup();

stats_para = {}
# P against T 
ax = plt.subplot(5,2,1);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([0,6]);ax.locator_params(axis='x',nbins=6);ax.set_xticklabels([]);
ax.set_ylim([-5,25]); ax.set_yticks(np.arange(-5,25.1,6))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(a)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate(r'$\mathrm{\mathsf{\Delta}}$PMEAN (%)',xy=(-0.20, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=10)					
x_rcp85 = TS_rcp85_TiSe; y_rcp = PMEAN_rcp85_TiSe; 
slope, intercept, r_value, p_value, std_err = stats.linregress(x_rcp85,y_rcp)
stats_para.update(p_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x_rcp85*slope+intercept;
x_fixa = TS_fixa_TiSe;y_fixa = PMEAN_fixa_TiSe;
slope, intercept, r_value, p_value, std_err = stats.linregress(x_fixa,y_fixa)
stats_para.update(p_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x_fixa*slope+intercept;
ax.scatter(x_rcp85,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5 (EA)')
ax.plot(x_rcp85,y_rcp_fit,'-',color="r",linewidth=1,label='Rcp8.5 fitting (EA)')
ax.scatter(x_fixa,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA (EA)')
ax.plot(x_fixa,y_fix_fit,'-',color="b",linewidth=1,label='Rcp8.5_FixA fitting (EA)')		

legend = ax.legend(shadow=False,ncol=1,bbox_to_anchor=(2.1, 1.0)) 
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.1)

# rx5day against T 
ax = plt.subplot(5,2,3);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([0,6]);ax.locator_params(axis='x',nbins=6);ax.set_xticklabels([]);
ax.set_ylim([-5,45]); ax.set_yticks(np.arange(-5,45.1,10))

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(b)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate(r'$\mathrm{\mathsf{\Delta}}$Rx5DAY (%)',xy=(-0.20, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=10)	
				
x = TS_rcp85_TiSe; y_rcp = rx5day_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(rx5day_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = rx5day_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(rx5day_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)	


# r30 against T 
ax = plt.subplot(5,2,5);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([0,6]);ax.locator_params(axis='x',nbins=6);ax.set_xticklabels([]);
ax.set_ylim([-10,60]); ax.set_yticks(np.arange(-10,60.1,14))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(c)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate(r'$\mathrm{\mathsf{\Delta}}$R30 (%)',xy=(-0.20, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=10)	
			
x = TS_rcp85_TiSe; y_rcp = r30_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(r30_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = r30_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(r30_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)	

# CWD against T 
ax = plt.subplot(5,2,7);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([0,6]);ax.locator_params(axis='x',nbins=6);ax.set_xticklabels([]);
ax.set_ylim([-15,10]); ax.set_yticks(np.arange(-15,10.1,5))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(d)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate(r'$\mathrm{\mathsf{\Delta}}$CWD (%)',xy=(-0.20, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=10)	

x = TS_rcp85_TiSe; y_rcp = cwd_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(cwd_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = cwd_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(cwd_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)	

# R99PR against T 
ax = plt.subplot(5,2,9);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([0,6]);ax.locator_params(axis='x',nbins=6);
ax.set_ylim([-10,40]); ax.set_yticks(np.arange(-10,40.1,10))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(e)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate(r'$\mathrm{\mathsf{\Delta}}$R95P (%)',xy=(-0.20, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical',fontsize=10)	
ax.annotate(r'$\mathrm{\mathsf{\Delta}}$SAT (degree)',xy=(0.5, -0.3), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal',fontsize=10)					
x = TS_rcp85_TiSe;y_rcp = r95p_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(r95p_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = r95p_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(r95p_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)
stats_para_P = {}

# rx5day against T 
ax = plt.subplot(5,2,4);pad= 5;	
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([-0.5,1.5]);ax.locator_params(axis='x',nbins=4);ax.set_xticklabels([]);
ax.set_ylim([-5,45]);ax.set_yticks(np.arange(-5,45.1,10))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(f)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.set_xticklabels([]);ax.locator_params(axis='x',nbins=6);		
x = PMEAN_rcp85_abso; y_rcp = rx5day_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(rx5day_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = rx5day_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(rx5day_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)	

# r30 against T 
ax = plt.subplot(5,2,6);pad= 5;	
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([-0.5,1.5]);ax.locator_params(axis='x',nbins=4);ax.set_xticklabels([]);
ax.set_ylim([-10,60]); ax.set_yticks(np.arange(-10,60.1,14))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(g)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.set_xticklabels([]);ax.locator_params(axis='x',nbins=6);		
x = PMEAN_rcp85_abso; y_rcp = r30_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(r30_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = r30_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(r30_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)	

# CWD against T 	
ax = plt.subplot(5,2,8);pad= 5;	
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_xlim([-0.5,1.5]);ax.locator_params(axis='x',nbins=4);ax.set_xticklabels([]);
ax.set_ylim([-15,10]); ax.set_yticks(np.arange(-15,10.1,5))
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(h)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.set_xticklabels([]);ax.locator_params(axis='x',nbins=6);
x = PMEAN_rcp85_abso; y_rcp = cwd_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(cwd_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = cwd_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(cwd_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)	

# R99PR against T 
ax = plt.subplot(5,2,10);pad= 5;
ax.set_xlim([-0.5,1.5]);ax.locator_params(axis='x',nbins=4);
ax.set_ylim([-10,40]); ax.set_yticks(np.arange(-10,40.1,10))
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=True, pad=5)
ax.annotate('(i)',xy=(0.02, 0.85), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.locator_params(axis='x',nbins=6);
ax.annotate(r'$\mathrm{\mathsf{\Delta}PMEAN (mm\/day^{-1})}$',xy=(0.5, -0.3), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal',fontsize=10)
x = PMEAN_rcp85_abso;y_rcp = r95p_rcp85_TiSe;
ax.scatter(x,y_rcp, marker=',',color="r",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(r95p_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="r",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = r95p_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="b",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(r95p_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="b",label='RCP8.5_FixA fitting',linewidth=1)


print '#######################################'
print 'EA'
print '#######################################'
print 'temp'
for key in stats_para.keys():
	print key
	print np.round(stats_para[key][:],3)

	
	
	
#############################
# SA
#############################
region = rergion_dic['SA'][:]
# mean surface temperature and total precipitation climatology
_,_,TS_JJA_mean = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,TS_JJA_mean_CESM)
TS_ref = stats.nanmean(stats.nanmean(stats.nanmean(TS_JJA_mean,axis=2),axis=1),axis=0)
_,_,PMEAN_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,PMEAN)
# plt.imshow(stats.nanmean(PMEAN_range,axis=0))
# plt.show()
PMEAN_ref = stats.nanmean(stats.nanmean(stats.nanmean(PMEAN_range,axis=2),axis=1),axis=0)
_,_,rx5day_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,rx5day)
rx5day_ref = stats.nanmean(stats.nanmean(stats.nanmean(rx5day_range,axis=2),axis=1),axis=0)
_,_,r95p_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,r95p)
r95p_ref = stats.nanmean(stats.nanmean(stats.nanmean(r95p_range,axis=2),axis=1),axis=0)
_,_,r30_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,r30)
r30_ref = stats.nanmean(stats.nanmean(stats.nanmean(r30_range,axis=2),axis=1),axis=0)
_,_,cwd_range = range_clip(region[0],region[1],region[2],region[3],lon_CESM,lat_CESM,cwd)
cwd_ref = stats.nanmean(stats.nanmean(stats.nanmean(cwd_range,axis=2),axis=1),axis=0)


# region clip
_,_,TS_rcp85_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,TS_rcp85)-TS_ref
_,_,TS_fixa_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,TS_fixa)-TS_ref
TS_rcp85_TiSe = stats.nanmean(stats.nanmean(stats.nanmean(TS_rcp85_clipped,axis = 3),axis=2),axis=0)
TS_fixa_TiSe = stats.nanmean(stats.nanmean(stats.nanmean(TS_fixa_clipped,axis = 3),axis=2),axis=0)

PMEAN_rcp85_TiSe= time_series_range(PMEAN_rcp85 ,PMEAN_ref,region)
rx5day_rcp85_TiSe= time_series_range(rx5day_rcp85,rx5day_ref,region)
cwd_rcp85_TiSe= time_series_range(cwd_rcp85,cwd_ref,region )
r30_rcp85_TiSe= time_series_range(r30_rcp85,r30_ref,region )
r95p_rcp85_TiSe= time_series_range(r95p_rcp85,r95p_ref,region )

PMEAN_fixa_TiSe = time_series_range(PMEAN_fixa ,PMEAN_ref,region)
rx5day_fixa_TiSe= time_series_range(rx5day_fixa ,rx5day_ref,region)
cwd_fixa_TiSe= time_series_range(cwd_fixa,cwd_ref,region)
r30_fixa_TiSe= time_series_range(r30_fixa,r30_ref,region)
r95p_fixa_TiSe= time_series_range(r95p_fixa,r95p_ref,region) #

PMEAN_rcp85_abso= time_series_absolute(PMEAN_rcp85 ,PMEAN_ref,region)
PMEAN_fixa_abso = time_series_absolute(PMEAN_fixa ,PMEAN_ref,region)


#############################
# plotting
#############################
stats_para = {}
# P against T 
ax = plt.subplot(5,2,1);				
x_rcp85 = TS_rcp85_TiSe; y_rcp = PMEAN_rcp85_TiSe; 
slope, intercept, r_value, p_value, std_err = stats.linregress(x_rcp85,y_rcp)
stats_para.update(p_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x_rcp85*slope+intercept;
x_fixa = TS_fixa_TiSe;y_fixa = PMEAN_fixa_TiSe;
slope, intercept, r_value, p_value, std_err = stats.linregress(x_fixa,y_fixa)
stats_para.update(p_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x_fixa*slope+intercept;
ax.scatter(x_rcp85,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5 (SA)')
ax.plot(x_rcp85,y_rcp_fit,'-',color="m",linewidth=1,label='Rcp8.5 fitting (SA)')
ax.scatter(x_fixa,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA (SA)')
ax.plot(x_fixa,y_fix_fit,'-',color="g",linewidth=1,label='Rcp8.5_FixA fitting (SA)')		
legend = ax.legend(shadow=False,ncol=1,bbox_to_anchor=(2.1, 1.0))	 
legend.get_frame().set_facecolor('gray');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.1)

# rx5day against T 	
ax = plt.subplot(5,2,3);			
x = TS_rcp85_TiSe; y_rcp = rx5day_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(rx5day_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = rx5day_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(rx5day_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)	


# r30 against T
ax = plt.subplot(5,2,5); 		
x = TS_rcp85_TiSe; y_rcp = r30_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(r30_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = r30_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(r30_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)	

# CWD against T 
ax = plt.subplot(5,2,7);
x = TS_rcp85_TiSe; y_rcp = cwd_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(cwd_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = cwd_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(cwd_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)	

# R99PR against T 
ax = plt.subplot(5,2,9);			
x = TS_rcp85_TiSe;y_rcp = r95p_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para.update(r95p_t_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = TS_fixa_TiSe;y_fixa = r95p_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para.update(r95p_t_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)
stats_para_P = {}

# rx5day against T 
ax = plt.subplot(5,2,4);
x = PMEAN_rcp85_abso; y_rcp = rx5day_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(rx5day_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = rx5day_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(rx5day_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)	

# r30 against T 
ax = plt.subplot(5,2,6);
x = PMEAN_rcp85_abso; y_rcp = r30_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(r30_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = r30_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(r30_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)	

# CWD against T 	
ax = plt.subplot(5,2,8);
x = PMEAN_rcp85_abso; y_rcp = cwd_rcp85_TiSe; 
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(cwd_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = cwd_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(cwd_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)	

# R99PR against T 
ax = plt.subplot(5,2,10);
x = PMEAN_rcp85_abso;y_rcp = r95p_rcp85_TiSe;
ax.scatter(x,y_rcp, marker=',',color="m",s=3, edgecolor='none',label='Rcp8.5');
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_rcp)
stats_para_P.update(r95p_p_rcp=[slope, intercept, r_value, p_value, std_err])
y_rcp_fit = x*slope+intercept;
ax.plot(x,y_rcp_fit,'-',color="m",label='RCP8.5 fitting',linewidth=1)
x = PMEAN_fixa_abso;y_fixa = r95p_fixa_TiSe;
ax.scatter(x,y_fixa, marker=',',color="g",s=3, edgecolor='none',label='Rcp8.5_FixA')
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y_fixa)
stats_para_P.update(r95p_p_fix=[slope, intercept, r_value, p_value, std_err])
y_fix_fit = x*slope+intercept;
ax.plot(x,y_fix_fit,'-',color="g",label='RCP8.5_FixA fitting',linewidth=1)

print '#######################################'
print 'SA'
print '#######################################'
print 'temp'
for key in stats_para.keys():
	print key
	print np.round(stats_para[key][:],3)

print 'prep'
for key in stats_para_P.keys():
	print key
	print np.round(stats_para_P[key][:],3)		

plt.subplots_adjust(left=0.12, bottom=0.07, right=0.95, top=0.96, wspace=0.2, hspace=0.1);
plt.savefig('fig11.png', format='png', dpi=1000)
plt.show()

