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
oceanmask=spio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
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

def kde_bins(region_key,variable,bin_min,bin_max):
	if region_key =='SA': length =15660
	else: length =24420
	region =rergion_dic[region_key][:];
	def readin_reshape_pdf(scenario,region,variable,no_member,bin_min,bin_max):
		if scenario =='his':
			input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
			file_name = 'historical_ensumble_mean_PEI_global_1920_2005.nc'
			nc_fid = nc4.Dataset(input_path+file_name,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			att_value = nc_fid.variables[variable][:]
			att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
			lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
			att_ref = att_clipped_r85F[66:86,:,:];
			nc_fid.close()
			size_m =np.shape(att_ref);size = size_m[0]*size_m[1]*size_m[2]
			data = np.reshape(att_ref,(size*no_member,1))
			mask = ~np.isnan(data); data = data[mask]					
		else:
			input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'+scenario+'/'
			os.system('find ' + os.path.abspath(input_path) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path + '/file_list.txt')
			text_file = open(input_path + '/file_list.txt', "r")
			text_content = text_file.readlines()
			data = np.empty((no_member,length))
			for ensumble_member in range(0,no_member): 
				nc_f = text_content[ensumble_member][:-1]
				nc_fid = nc4.Dataset(nc_f,mode='r')
				lat = nc_fid.variables['lat'][:]
				lon = nc_fid.variables['lon'][:]
				att_value = nc_fid.variables[variable][75:95,:,:]
				att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
				lons,lats,att_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
				size_m =np.shape(att_clipped);size = size_m[0]*size_m[1]*size_m[2]
				data[ensumble_member,:]=np.reshape(att_clipped,(1,size))
				nc_fid.close()
			data = np.reshape(data,(1,no_member*length))
			mask = ~np.isnan(data); data = data[mask]					
		
		bins = linspace( bin_min,bin_max, 100);kde = gaussian_kde( data );
		pdfs = 100*kde(bins)
		mean = np.nanmean(data)
		return bins,pdfs,mean
	scenario='his';no_member=1;bins,his_pdf,his_mean = readin_reshape_pdf(scenario,region,variable,no_member,bin_min,bin_max)		
	scenario='rcp85';no_member=30;_,rcp85_pdf,rcp85_mean = readin_reshape_pdf(scenario,region,variable,no_member,bin_min,bin_max)		
	scenario='fixa';no_member=12;_,fixa_pdf,fixa_mean = readin_reshape_pdf(scenario,region,variable,no_member,bin_min,bin_max)		
	return bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean

def plot_pdfs(ax,x,y1,y2,y3,v1,v2,v3):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.plot(x,y1,'-',color="k",linewidth=1,label = 'Baseline')
	ax.axvline(x=v1,color='k',linewidth = 1,ls='--');
	index =[index for index in range(len(y1)) if y1[index]==np.nanmedian(y1)]; print index
	# ax.axvline(x=x[index],color='k',linewidth = 1,ls='-.');
	ax.plot(x,y2,'-',color="r",linewidth=1,label = 'RCP8.5')
	ax.axvline(x=v2,color='r',linewidth = 1,ls='--');
	index=[index for index in range(len(y2)) if y2[index]==np.nanmedian(y2)]; print index
	# ax.axvline(x=x[index],color='r',linewidth = 1,ls='-.');
	ax.plot(x,y3,'-',color="b",linewidth=1,label = 'RCP8.5_FixA')
	ax.axvline(x=v3,color='b',linewidth = 1,ls='--');
	index =[index for index in range(len(y3)) if y3[index]==np.nanmedian(y3)]; print index
	# ax.axvline(x=x[index],color='b',linewidth = 1,ls='-.');
	return ax
	

fig = plt.figure(figsize=(8, 8), facecolor='White');plot_setup();pad=5
######SA
region = 'SA'

ax = plt.subplot(5,2,1);
ax.annotate(r'(a) PMEAN ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('South Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
variable ='total_precip';bin_min=0;bin_max=20*92
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));ax.set_xticklabels(np.arange(bin_min,bin_max/92+.1,(bin_max-bin_min)/5/92));
ax.set_ylim([0,.25]);ax.set_yticks(np.arange(0,.26,.05));
legend1 = ax.legend(shadow=False,ncol=1,loc='center right')	 
legend1.get_frame().set_facecolor('None');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.3)
	
ax = plt.subplot(5,2,3);
			
ax.annotate('(b) RX5DAY (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
variable ='rx5day';bin_min=0;bin_max=450
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
ax.set_ylim([0,1.2]);ax.set_yticks(np.arange(0,1.21,.2));

		
ax = plt.subplot(5,2,5);
ax.annotate('(c) R10 (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)							
ax.annotate('PDF (%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)	
variable ='r10';bin_min=0;bin_max=42
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/6));
ax.set_ylim([0,8]);ax.set_yticks(np.arange(0,8.1,2));


ax = plt.subplot(5,2,7);
ax.annotate('(d) CWD (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
variable ='cwd';bin_min=0;bin_max=90
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
# ax.set_ylim([0,4.8]);ax.set_yticks(np.arange(0,4.81,1.2));


ax = plt.subplot(5,2,9);
ax.annotate('(e) R95P (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)		
variable ='r95p';bin_min=0;bin_max=600
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
ax.set_ylim([0,.8]);ax.set_yticks(np.arange(0,.81,.2));

region = 'EA'

ax = plt.subplot(5,2,2);
ax.annotate(r'(f) PMEAN ($\mathrm{\mathsf{mm\/day^{-1}}}$)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('East Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
variable ='total_precip';bin_min=0;bin_max=15*92
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));ax.set_xticklabels(np.arange(bin_min,bin_max/92+.1,(bin_max-bin_min)/5/92));

ax = plt.subplot(5,2,4);
			
ax.annotate('(g) RX5DAY (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
variable ='rx5day';bin_min=0;bin_max=279
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
# ax.set_ylim([0,4.8]);ax.set_yticks(np.arange(0,4.81,1.2));

		
ax = plt.subplot(5,2,6);
ax.annotate('(h) R10 (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)								
variable ='r10';bin_min=0;bin_max=30
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
ax.set_ylim([0,8]);ax.set_yticks(np.arange(0,8.1,2));


ax = plt.subplot(5,2,8);
ax.annotate('(i) CWD (days)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)	
variable ='cwd';bin_min=0;bin_max=30
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
# ax.set_ylim([0,4.8]);ax.set_yticks(np.arange(0,4.81,1.2));


ax = plt.subplot(5,2,10);
ax.annotate('(j) R95P (mm)',(0.98, 0.82), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)		
variable ='r95p';bin_min=0;bin_max=600
bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean = kde_bins(region,variable,bin_min,bin_max)
ax = plot_pdfs(ax,bins,his_pdf,rcp85_pdf,fixa_pdf,his_mean,rcp85_mean,fixa_mean)
ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/5));
ax.set_ylim([0,.8]);ax.set_yticks(np.arange(0,.81,.2));
			
plt.subplots_adjust(left=0.1, bottom=0.03, right=0.97, top=0.95, wspace=0.25, hspace=0.2);
plt.savefig('Fig10_a.png', format='png', dpi=1000)
plt.show()