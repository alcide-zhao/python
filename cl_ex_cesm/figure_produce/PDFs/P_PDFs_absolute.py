	
def spatial_diff_sig(region, variable):
	from scipy.stats import ttest_ind as test	
	input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'
	#fixa
	file_name = 'historical_ensumble_mean_PEI_global_1920_2005.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask); att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85F = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	#RCP8.5
	file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	
	att_3150_R = att_clipped_r85[25:45,:,:]
	att_3150_F = att_clipped_r85F[66:86,:,:]
	att_3150_D = att_3150_R -att_3150_F
	att_8100_R = att_clipped_r85[75:95,:,:]
	att_8100_F = att_clipped_r85F[66:86,:,:]
	att_8100_D = att_8100_R -att_8100_F
	return lons,lats,att_3150_D,att_8100_D


# INDICES REPRODUCING
region = rergion_dic['EA'][:]	
lons,lats,r95pr_3150_D,r95pr_8100_D= spatial_diff_sig(region, variable = 'r95pr')
_,_,cwd_3150_D,cwd_8100_D = spatial_diff_sig(region, variable = 'cwd')
_,_,sdii_3150_D,sdii_8100_D = spatial_diff_sig(region, variable = 'sdii')
_,_,rx5day_3150_D,rx5day_8100_D= spatial_diff_sig(region, variable = 'rx5day')
_,_,PTOT_3150_D,PTOT_8100_D= spatial_diff_sig(region, variable = 'total_precip')	

		
ax = plt.subplot(5,2,2);
size =np.shape(PTOT_3150_D)
hist_CESM = np.reshape(PTOT_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
dist_space = linspace( -4, 4, 100)
kde = gaussian_kde( hist_CESM/92 )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1)
size =np.shape(PTOT_8100_D); 
hist_CESM = np.reshape(PTOT_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM/92 )	
plt.plot(dist_space, kde(dist_space),'r-',alpha=1 )

							
ax = plt.subplot(5,2,4);
dist_space = linspace( -80, 80, 100)
size =np.shape(rx5day_3150_D)
hist_CESM = np.reshape(rx5day_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1)
size =np.shape(rx5day_8100_D)
hist_CESM = np.reshape(rx5day_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)

					
ax = plt.subplot(5,2,6);
dist_space = linspace( -4, 4, 100)
size =np.shape(sdii_3150_D)
hist_CESM = np.reshape(sdii_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1 )
size =np.shape(sdii_8100_D)
hist_CESM = np.reshape(sdii_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)

							
ax = plt.subplot(5,2,8);
dist_space = linspace( -10, 10, 100)
size =np.shape(cwd_3150_D)
hist_CESM = np.reshape(cwd_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1 )
size =np.shape(cwd_8100_D)
hist_CESM = np.reshape(cwd_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1 )

			
				
ax = plt.subplot(5,2,10);
dist_space = linspace( -0.6,0.6,100)
size =np.shape(r95pr_3150_D)
hist_CESM = np.reshape(r95pr_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1 )
size =np.shape(r95pr_8100_D)
hist_CESM = np.reshape(r95pr_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)


region = rergion_dic['SA'][:]
# INDICES REPRODUCING

lons,lats,r95pr_3150_D,r95pr_8100_D= spatial_diff_sig(region, variable = 'r95pr')
_,_,cwd_3150_D,cwd_8100_D = spatial_diff_sig(region, variable = 'cwd')
_,_,sdii_3150_D,sdii_8100_D = spatial_diff_sig(region, variable = 'sdii')
_,_,rx5day_3150_D,rx5day_8100_D= spatial_diff_sig(region, variable = 'rx5day')
_,_,PTOT_3150_D,PTOT_8100_D= spatial_diff_sig(region, variable = 'total_precip')


ax = plt.subplot(5,2,1)
dist_space = linspace( -4, 4, 100)
size =np.shape(PTOT_3150_D)
hist_CESM = np.reshape(PTOT_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM/92 )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1,label ='ALL_3150')

size =np.shape(PTOT_8100_D); 
hist_CESM = np.reshape(PTOT_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM/92 )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1,label ='ALL_8100')

legend=plt.legend(loc=2,shadow=False);
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.3)
		
ax = plt.subplot(5,2,3);
dist_space = linspace( -150,150,100)
size =np.shape(rx5day_3150_D)
hist_CESM = np.reshape(rx5day_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1)

size =np.shape(rx5day_8100_D)
hist_CESM = np.reshape(rx5day_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)
				
ax = plt.subplot(5,2,5)
dist_space = linspace( -6, 6, 100)
size =np.shape(sdii_3150_D)
hist_CESM = np.reshape(sdii_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1)

size =np.shape(sdii_8100_D)
hist_CESM = np.reshape(sdii_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)
				
ax = plt.subplot(5,2,7);	
dist_space = linspace( -20,20,100)
size =np.shape(cwd_3150_D)
hist_CESM = np.reshape(cwd_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1)

size =np.shape(cwd_8100_D)
hist_CESM = np.reshape(cwd_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)

				
ax = plt.subplot(5,2,9);
dist_space = linspace( -0.6,0.6,100)
size =np.shape(r95pr_3150_D)
hist_CESM = np.reshape(r95pr_3150_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'k-',alpha=1)

size =np.shape(r95pr_8100_D)
hist_CESM = np.reshape(r95pr_8100_D,(size[0]*size[1]*size[2],1))
mask = ~np.isnan(hist_CESM); hist_CESM = hist_CESM[mask]					
kde = gaussian_kde( hist_CESM )
plt.plot(dist_space, kde(dist_space),'r-',alpha=1)
							
plt.subplots_adjust(left=0.2, bottom=0.03, right=0.98, top=0.95, wspace=0.15, hspace=0.15);
plt.savefig('Fig14.png', format='png', dpi=1000)
plt.show()