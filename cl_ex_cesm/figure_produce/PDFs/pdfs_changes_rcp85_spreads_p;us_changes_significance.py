from numpy import linspace
from scipy.stats.kde import gaussian_kde
def time_series_range(variable):
	size = np.shape(variable)[1]
	variable_median = stats.nanmean(variable,axis = 0)
	variable_upper = np.empty((size));variable_upper[:]= np.nan
	variable_loweer = np.empty((size));variable_loweer[:]= np.nan
	for itime in range(size):	
		data = variable[:,itime];
		bins = linspace(min(data), max(data), 50); 
		gkde=gaussian_kde(data)
		kdepdf = gkde.evaluate(bins)
		kde = 100*kdepdf/np.sum(kdepdf)
		variable_upper[itime] = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
		variable_loweer[itime] = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == 10 or (sum(kde[0:i])<=10 and sum(kde[0:i+1])>=10))][0]
	return variable_median,variable_upper,variable_loweer



##The Mann-Whistney U-test
def mannwhitneyu_test(data1,data2):
	p_threshold = 0.1
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
	dis_min =0.75*np.nanmin(att_ref[:]);dis_max = 1.25*np.nanmax(att_ref[:])
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
	
	rcp85_mean,rcp85_upper,rcp85_loweer= time_series_range(kde_rcp85)
	AA_diff = rcp85_mean - 	stats.nanmean(kde_fixa,axis = 0)
	significance = mannwhitneyu_test(kde_rcp85,kde_fixa )	
	return rcp85_mean,rcp85_upper,rcp85_loweer,AA_diff,significance
	
	
	
	
	
	
	