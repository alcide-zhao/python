# -*- coding: utf-8 -*-
"""
This is to produce the PDFs of precipitatio mean and extrmes over 2081-2100
The baseline and the net chanegs are plotted. 
"""
import netCDF4 as nc4
import numpy as np
import scipy.io as spio
from scipy import stats
import matplotlib.pyplot as plt
from numpy import linspace
from scipy.stats.kde import gaussian_kde
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

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}
oceanmask=spio.loadmat('//home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

#############################################
### aerosol differebces
#############################################
	
def time_series_range(variable):
	size = np.shape(variable)[1]
	variable_median = stats.nanmean(variable,axis = 0)
	variable_upper = np.empty((size));variable_upper[:]= np.nan
	variable_lower = np.empty((size));variable_lower[:]= np.nan
	for itime in range(size):	
		data = variable[:,itime];
		bins = linspace(min(data), max(data), 50); 
		gkde=gaussian_kde(data)
		kdepdf = gkde.evaluate(bins)
		kde = 100*kdepdf/np.sum(kdepdf)
		variable_upper[itime] = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
		variable_lower[itime] = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == 10 or (sum(kde[0:i])<=10 and sum(kde[0:i+1])>=10))][0]
	return variable_median,variable_upper,variable_lower



##The Mann-Whistney U-test
def mannwhitneyu_test(data1,data2):
	p_threshold = 0.025
	times = np.shape(data2)[1]
	p_value = np.empty((times));p_value[:]=np.nan
	from scipy.stats import mannwhitneyu as test
	for itime in range(times):
		cache1 = data1[:,itime]
		cache2 = data2[:,itime]
		_,p_value[itime] = test(cache1,cache2);
	p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
	return p_value

def kde_bins(region, variable):
	from scipy.stats import ttest_ind as test
	######## historical references
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	file_name = 'historical_ensumble_mean_PEI_global_1920_2005.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	att_ref = att_clipped_r85F[66:86,:,:];size_m =np.shape(att_ref);size = size_m[0]*size_m[1]*size_m[2]
	dis_min =np.nanmin(att_ref[:]);dis_max = np.nanmax(att_ref[:])
	hist_CESM = np.reshape(att_ref,(size,1))
	mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
	ds_ref = linspace( dis_min,dis_max, 100);kde = gaussian_kde( hist_CESM );
	kde_ref = 100*kde(ds_ref)
	nc_fid.close()
	
	######## rcp8.5
	input_path = '//exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/rcp85/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + '/file_list.txt', "r")
	text_content = text_file.readlines()
	kde_rcp85 = np.empty((30,100)); kde_rcp85[:]=np.nan
	for ensumble_member in range(0,30): 
		nc_f = text_content[ensumble_member][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		att_value = nc_fid.variables[variable][75:95,:,:]
		att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
		lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
		hist_CESM=np.reshape(att_clipped_r85,(1,size))
		mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
		kde = gaussian_kde( hist_CESM );
		kde_rcp85[ensumble_member,:] = 100*kde(ds_ref)
		nc_fid.close()
	######## fixa
	input_path = '//exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/fixa/'
	os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
	text_file = open(input_path + '/file_list.txt', "r")
	text_content = text_file.readlines()
	kde_fixa = np.empty((12,100)); kde_fixa[:]=np.nan
	for ensumble_member in range(0,12): 
		nc_f = text_content[ensumble_member][:-1]
		nc_fid = nc4.Dataset(nc_f,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		att_value = nc_fid.variables[variable][75:95,:,:]
		att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
		lons,lats,att_clipped_fixa = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
		hist_CESM=np.reshape(att_clipped_fixa,(1,size))
		mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
		kde = gaussian_kde( hist_CESM );
		kde_fixa[ensumble_member,:] = 100*kde(ds_ref)
		nc_fid.close()
	
	rcp85_mean,rcp85_upper,rcp85_lower= time_series_range(kde_rcp85)
	rcp85_mean=rcp85_mean-kde_ref;rcp85_upper=rcp85_upper-kde_ref;rcp85_lower=rcp85_lower-kde_ref;
	AA_diff = rcp85_mean - 	stats.nanmean(kde_fixa,axis = 0)+kde_ref
	significance = mannwhitneyu_test(kde_rcp85,kde_fixa )	
	return ds_ref,kde_ref,rcp85_mean,rcp85_upper,rcp85_lower,AA_diff,significance


	
######SA
region = rergion_dic['SA'][:]
ptot_ds,ptot_ref,ptot_rcp85_mean,ptot_rcp85_upper,ptot_rcp85_lower,ptot_AA,ptot_sig= kde_bins(region, variable = 'total_precip')
rx5day_ds,rx5day_ref,rx5day_rcp85_mean,rx5day_rcp85_upper,rx5day_rcp85_lower,rx5day_AA,rx5day_sig= kde_bins(region, variable = 'rx5day')
sdii_ds,sdii_ref,sdii_rcp85_mean,sdii_rcp85_upper,sdii_rcp85_lower,sdii_AA,sdii_sig = kde_bins(region, variable = 'r10')
cwd_ds,cwd_ref,cwd_rcp85_mean,cwd_rcp85_upper,cwd_rcp85_lower,cwd_AA,cwd_sig = kde_bins(region, variable = 'cwd')
r95p_ds,r95p_ref,r95p_rcp85_mean,r95p_rcp85_upper,r95p_rcp85_lower,r95p_AA,r95p_sig= kde_bins(region,variable = 'r95p')	


fig = plt.figure(figsize=(8, 8), facecolor='White');plot_setup()

ax1 = plt.subplot(5,2,1);ax1.axhline(0, color='k',alpha=1);pad= 5;
ax1.annotate(r'(a) PMEAN ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.annotate('South Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
# ax1.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',rotation = 90)
ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax1.set_xlim([0,18]); ax1.set_xticks(np.arange(0,18.1,6));
ax1.set_ylim([-0.08,0.08]); ax1.set_yticks(np.arange(-0.08,0.081,0.04))
ax11 = ax1.twinx();
ax11.set_xlim([0,18]);ax11.set_xticks(np.arange(0,18.1,6));
ax11.set_ylim([0,0.2]);ax11.set_yticks(np.arange(0,0.21,0.05))
align_yaxis(ax1,-0.08,ax11,0)

ax11.plot(ptot_ds/92, ptot_ref,'k',alpha=1,label='Baseline')
ax1.plot(ptot_ds/92, ptot_rcp85_mean,'r-',alpha=1,label='RCP8.5')
ax1.fill_between(ptot_ds/92,ptot_rcp85_upper,ptot_rcp85_lower, where=ptot_rcp85_upper>=ptot_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax1.plot(ptot_ds/92, ptot_AA,'r--',alpha=1,label='AAs')
ax1.plot(ptot_ds/92,ptot_sig*-0.06,"g-",  markersize=3, markeredgecolor='none')
	
ax3 = plt.subplot(5,2,3);ax3.axhline(0, color='k',alpha=1);
ax3.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax3.set_xlim([0,450]); ax3.set_xticks(np.arange(0,450.1,150));
ax3.set_ylim([-0.6,0.60]);ax3.set_yticks(np.arange(-0.6,0.61,0.3))
ax3.annotate('(b) RX5DAY (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
# ax3.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',rotation = 90)	
ax33 = ax3.twinx();
ax33.set_xlim([0,450]);ax33.set_xticks(np.arange(0,450.1,150));
ax33.set_ylim([0,1.00]);ax33.set_yticks(np.arange(0,1.01,0.25))
align_yaxis(ax3,-0.6,ax33,0)

ax33.plot(rx5day_ds, rx5day_ref,'k',alpha=1,label='Baseline')
ax3.plot(rx5day_ds, rx5day_rcp85_mean,'r-',alpha=1,label='Baseline')
ax3.fill_between(rx5day_ds,rx5day_rcp85_upper,rx5day_rcp85_lower, where=rx5day_rcp85_upper>=rx5day_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax3.plot(rx5day_ds, rx5day_AA,'r--',alpha=1,label='Baseline')
ax3.plot(rx5day_ds,rx5day_sig*-0.45,"g-",  markersize=3, markeredgecolor='none')
			
ax5 = plt.subplot(5,2,5);ax5.axhline(0, color='k',alpha=1);
ax5.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax5.set_xlim([0,42]); ax5.set_xticks(np.arange(0,42.1,14));
ax5.set_ylim([-1.8,1.8]);ax5.set_yticks(np.arange(-1.8,1.81,0.9))
ax5.annotate('(c) R10 (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)							
ax5.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)	
ax55 = ax5.twinx();
ax55.set_xlim([0,42]);ax55.set_xticks(np.arange(0,42.1,14));
ax55.set_ylim([0,4.8]);ax55.set_yticks(np.arange(0,4.81,1.2));
align_yaxis(ax5,-1.8,ax55,0)

ax55.plot(sdii_ds, sdii_ref,'k',alpha=1,label='Baseline')
ax5.plot(sdii_ds, sdii_rcp85_mean,'r-',alpha=1,label='Baseline')
ax5.fill_between(sdii_ds,sdii_rcp85_upper,sdii_rcp85_lower, where=sdii_rcp85_upper>=sdii_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax5.plot(sdii_ds, sdii_AA,'r--',alpha=1,label='Baseline')
ax5.plot(sdii_ds,sdii_sig*-1.35,"g-",  markersize=3, markeredgecolor='none')
			
ax7 = plt.subplot(5,2,7);ax7.axhline(0, color='k',alpha=1);
ax7.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax7.set_xlim([0,90]); ax7.set_xticks(np.arange(0,90.1,30));
ax7.set_ylim([-1.2,1.2]);ax7.set_yticks(np.arange(-1.2,1.21,0.6))
ax7.annotate('(d) CWD (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
# ax7.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',rotation = 90)	
ax77 = ax7.twinx();
ax77.set_xlim([0,90]);ax77.set_xticks(np.arange(0,90.1,30));
ax77.set_ylim([0,2.8]);ax77.set_yticks(np.arange(0,2.81,0.7));
align_yaxis(ax7,-1.2,ax77,0)

ax77.plot(cwd_ds, cwd_ref,'k',alpha=1,label='Baseline')
ax7.plot(cwd_ds, cwd_rcp85_mean,'r-',alpha=1,label='Baseline')
ax7.fill_between(cwd_ds,cwd_rcp85_upper,cwd_rcp85_lower, where=cwd_rcp85_upper>=cwd_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax7.plot(cwd_ds, cwd_AA,'r--',alpha=1,label='Baseline')
ax7.plot(cwd_ds,cwd_sig*-0.9,"g-",  markersize=3, markeredgecolor='none')
				
ax9 = plt.subplot(5,2,9);ax9.axhline(0, color='k',alpha=1);
ax9.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax9.set_xlim([0,600]); ax9.set_xticks(np.arange(0,600.1,200));
ax9.set_ylim([-0.4,0.4]);ax9.set_yticks(np.arange(-0.4,0.41,0.2))
ax9.annotate('(e) R95P (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)		
# ax9.annotate(r'$\mathrm{\mathsf{\Delta}}$PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='center',rotation = 90)	
ax99 = ax9.twinx();
ax99.set_xlim([0,600]);ax99.set_xticks(np.arange(0,600.1,200));
ax99.set_ylim([0,0.72]);ax99.set_yticks(np.arange(0,0.721,0.18));
align_yaxis(ax9,-0.4,ax99,0)	

ax99.plot(r95p_ds, r95p_ref,'k',alpha=1,label='Baseline')
ax9.plot(r95p_ds, r95p_rcp85_mean,'r-',alpha=1,label='Baseline')
ax9.fill_between(r95p_ds,r95p_rcp85_upper,r95p_rcp85_lower, where=r95p_rcp85_upper>=r95p_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax9.plot(r95p_ds, r95p_AA,'r--',alpha=1,label='Baseline')
ax9.plot(r95p_ds,r95p_sig*-0.3,"g-",  markersize=3, markeredgecolor='none')
	

region = rergion_dic['EA'][:]	
ptot_ds,ptot_ref,ptot_rcp85_mean,ptot_rcp85_upper,ptot_rcp85_lower,ptot_AA,ptot_sig= kde_bins(region, variable = 'total_precip')
rx5day_ds,rx5day_ref,rx5day_rcp85_mean,rx5day_rcp85_upper,rx5day_rcp85_lower,rx5day_AA,rx5day_sig= kde_bins(region, variable = 'rx5day')
sdii_ds,sdii_ref,sdii_rcp85_mean,sdii_rcp85_upper,sdii_rcp85_lower,sdii_AA,sdii_sig = kde_bins(region, variable = 'r10')
cwd_ds,cwd_ref,cwd_rcp85_mean,cwd_rcp85_upper,cwd_rcp85_lower,cwd_AA,cwd_sig = kde_bins(region, variable = 'cwd')
r95p_ds,r95p_ref,r95p_rcp85_mean,r95p_rcp85_upper,r95p_rcp85_lower,r95p_AA,r95p_sig= kde_bins(region,variable = 'r95p')	
		
ax2 = plt.subplot(5,2,2);ax2.axhline(0, color='k',alpha=1);
ax2.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax2.set_xlim([0,15]); ax2.set_xticks(np.arange(0,15.1,5));
ax2.set_ylim([-0.08,0.08]); ax2.set_yticks(np.arange(-0.08,0.081,0.04))
ax2.annotate(r'(f) PMEAN ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax2.annotate('East Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax22 = ax2.twinx();
ax22.set_xlim([0,15]);ax22.set_xticks(np.arange(0,15.1,5));
ax22.set_ylim([0,0.28]);ax22.set_yticks(np.arange(0,0.281,0.07))
align_yaxis(ax2,-0.08,ax22,0)

ax22.plot(ptot_ds/92, ptot_ref,'k',alpha=1,label='Baseline')
ax2.plot(ptot_ds/92, ptot_rcp85_mean,'r-',alpha=1,label='Rcp8.5')
ax2.fill_between(ptot_ds/92,ptot_rcp85_upper,ptot_rcp85_lower, where=ptot_rcp85_upper>=ptot_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax2.plot(ptot_ds/92, ptot_AA,'r--',alpha=1,label='AAs')
ax2.plot(ptot_ds/92,ptot_sig*-0.06,"g-",  markersize=3, markeredgecolor='none')
							
ax4 = plt.subplot(5,2,4);ax4.axhline(0, color='k',alpha=1);
ax4.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax4.set_xlim([0,240]);ax4.set_xticks(np.arange(0.0,250.1,80))
ax4.set_ylim([-0.4,0.4]);ax4.set_yticks(np.arange(-0.4,0.41,0.2));
ax4.annotate('(g) RX5DAY (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax44 = ax4.twinx();
ax44.set_xlim([0,270]);ax44.set_xticks(np.arange(0,270.1,90));
ax44.set_ylim([0,1.2]);ax44.set_yticks(np.arange(0,1.21,0.3));
align_yaxis(ax4,-0.4,ax44,0)

ax44.plot(rx5day_ds, rx5day_ref,'k',alpha=1,label='Baseline')
ax4.plot(rx5day_ds, rx5day_rcp85_mean,'r-',alpha=1,label='Baseline')
ax4.fill_between(rx5day_ds,rx5day_rcp85_upper,rx5day_rcp85_lower, where=rx5day_rcp85_upper>=rx5day_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax4.plot(rx5day_ds, rx5day_AA,'r--',alpha=1,label='Baseline')
ax4.plot(rx5day_ds,rx5day_sig*-0.3,"g-",  markersize=3, markeredgecolor='none')
					
ax6 = plt.subplot(5,2,6);ax6.axhline(0, color='k',alpha=1);
ax6.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax6.set_xlim([0,30]);ax6.set_xticks(np.arange(0,30.1,10));
ax6.set_ylim([-3,3]);ax6.set_yticks(np.arange(-3,3.1,1.5))
ax6.annotate('(h) R10 (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
ax6.annotate('PDF (%)',xy=(1.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)					
ax66 = ax6.twinx();
ax66.set_xlim([0,30]);ax66.set_xticks(np.arange(0,30.1,10));
ax66.set_ylim([0,8]);ax66.set_yticks(np.arange(0,8.1,2));
align_yaxis(ax6,-3,ax66,0)

ax66.plot(sdii_ds, sdii_ref,'k',alpha=1,label='Baseline')
ax6.plot(sdii_ds, sdii_rcp85_mean,'r-',alpha=1,label='Baseline')
ax6.fill_between(sdii_ds,sdii_rcp85_upper,sdii_rcp85_lower, where=sdii_rcp85_upper>=sdii_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax6.plot(sdii_ds, sdii_AA,'r--',alpha=1,label='Baseline')
ax6.plot(sdii_ds,sdii_sig*-2.25,"g-",  markersize=3, markeredgecolor='none')
					
ax8 = plt.subplot(5,2,8);ax8.axhline(0, color='k',alpha=1)	;
ax8.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax8.set_xlim([0,24]);ax8.set_xticks(np.arange(0,24.1,8));
ax8.set_ylim([-4,4]);ax8.set_yticks(np.arange(-4,4.1,2))
ax8.annotate('(i) CWD (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)					
ax88 = ax8.twinx();
ax88.set_xlim([0,30]);ax88.set_xticks(np.arange(0,30.1,10));
ax88.set_ylim([0,10]);ax88.set_yticks(np.arange(0,10.1,2.5));
align_yaxis(ax8,-4,ax88,0)

ax88.plot(cwd_ds, cwd_ref,'k',alpha=1,label='Baseline')
ax8.plot(cwd_ds, cwd_rcp85_mean,'r-',alpha=1,label='Baseline')
ax8.fill_between(cwd_ds,cwd_rcp85_upper,cwd_rcp85_lower, where=cwd_rcp85_upper>=cwd_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax8.plot(cwd_ds, cwd_AA,'r--',alpha=1,label='Baseline')
ax8.plot(cwd_ds,cwd_sig*-3,"g-",  markersize=3, markeredgecolor='none')
			
ax10 = plt.subplot(5,2,10);ax10.axhline(0, color='k',alpha=1);	
ax10.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='off',pad =3,which= 'major')
ax10.set_xlim([0,600]);ax10.set_xticks(np.arange(0,600.1,200));
ax10.set_ylim([-0.3,0.3]);ax10.set_yticks(np.arange(-0.3,0.31,0.15))
ax10.annotate('(j) R95P (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)			
ax00 = ax10.twinx();
ax00.set_xlim([0,600]);ax00.set_xticks(np.arange(0,600.1,200));
ax00.set_ylim([0,0.72]);ax00.set_yticks(np.arange(0,0.721,0.18));
align_yaxis(ax10,-0.3,ax00,0)

ax00.plot(r95p_ds, r95p_ref,'k',alpha=1,label='Baseline')
ax10.plot(r95p_ds, r95p_rcp85_mean,'r-',alpha=1,label='Baseline')
ax10.fill_between(r95p_ds,r95p_rcp85_upper,r95p_rcp85_lower, where=r95p_rcp85_upper>=r95p_rcp85_lower, color="r", alpha = 0.3,edgecolor='none') 
ax10.plot(r95p_ds, r95p_AA,'r--',alpha=1,label='Baseline')
ax10.plot(r95p_ds,r95p_sig*-0.215,"g-",  markersize=3, markeredgecolor='none')
	

legend1 = ax1.legend(shadow=False,ncol=1,bbox_to_anchor=(0.32, 0.93))	 
legend1.get_frame().set_facecolor('None');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.3)
legend11 = ax11.legend(shadow=False,ncol=1,bbox_to_anchor=(0.35, 1.05))	 
legend11.get_frame().set_facecolor('None');legend11.get_frame().set_edgecolor('None');legend11.get_frame().set_alpha(0.1)


							
plt.subplots_adjust(left=0.1, bottom=0.03, right=0.90, top=0.95, wspace=0.25, hspace=0.2);
plt.savefig('Fig10.png', format='png', dpi=1000)
plt.show()