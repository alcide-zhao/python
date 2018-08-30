
'''
# precipitation extremes
prep_att_ref = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_dic = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_list=['rx5day','r99p','sdii','cwd']
prep_unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

# temperature extremes
# temp_att_ref stores the reference value which is used to remove climatology
temp_att_ref = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_dic = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_list = ['txn','su25','dtr','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
temp_unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}

# precipitation extremes
prep_att_ref = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_dic = {'rx5day':[],'r99p':[],'sdii':[],'cwd':[]}
prep_att_list=['rx5day','r99p','sdii','cwd']
prep_unit_dic = {'rx1day':'Mm', 'rx5day':'Mm', 'sdii':'Mm/day', 'r10':'Days', 'r20':'Days', 'rnm':'Days', 'cdd':'Days','cwd':'Days', 'r95p':'Mm', 'r99p':'Mm', 'precptot':'Mm'}

# temperature extremes
# temp_att_ref stores the reference value which is used to remove climatology
temp_att_ref = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_dic = {'txn':[],'su25':[],'dtr':[],'tr20':[]}
temp_att_list = ['txn','su25','dtr','tr20'] #,'tn10p', 'tn90p','tx10p','tx90p'
temp_unit_dic = {'txx':'Degree', 'txn':'Degree', 'tnx':'Degree', 'tnn':'Degree', 'dtr':'Degree', 'fd0':'Days', 'su25':'Days','id0':'Days', 'tr20':'Days','tn10p':'%','tn90p':'%','tx10p':'%','tx90p':'%'}

#######################################################
# 1. time evolution of extrmes indices                #
#######################################################

spst = [2,2]   # subplot style
"""
Time series of temperature extreme indecies 
"""

input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/'
plt.figure(1, facecolor='White',figsize=[120,30])

file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100_Rcp8.5.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

sub_plot=1
for att_name in temp_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	temp_att_ref[att_name] = time_series[0]
	time_series = time_series -temp_att_ref[att_name]
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	# print np.shape(movingaverage(time_series,5))
	ax1.plot(time,movingaverage(time_series),'-',color="red",linewidth=2,label='RCP8.5')
	sub_plot = sub_plot+1
nc_fid.close()


file_name = 'Spatial_ensumble_mean_Temp_extremes_global_2006_2100_fixA.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	
sub_plot=1	
for att_name in temp_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	# temp_att_ref[att_name] = time_series[0]
	time_series = time_series -temp_att_ref[att_name]
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	ax1.plot(time,movingaverage(time_series),'-',color="blue",linewidth=2,label='RCP8.5_fixA')
	ax1.set_xlabel('Year',fontsize=10)
	ax1.set_xlim([2000,2100])
	# ax1.set_ylim([4,7.5])
	ax1.set_ylabel(temp_unit_dic[att_name].title(),fontsize=15)
	# ax1.set(aspect=10)
	ax1.set_title(att_name.title(),fontsize=15)
	ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
	sub_plot = sub_plot+1
nc_fid.close()

"""
Time series of precipitation ectreme indecies 
"""
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'


file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')
spst = [2,4]
sub_plot=1
for att_name in prep_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	prep_att_ref[att_name] = time_series[0]
	time_series = time_series -prep_att_ref[att_name]
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	ax1.plot(time,movingaverage(time_series),'-',color="red",linewidth=2,label='RCP8.5')
	sub_plot = sub_plot+1
nc_fid.close()


file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(input_path+file_name,mode='r')

sub_plot=1	
for att_name in prep_att_list:
	time = nc_fid.variables['time'][:]
	lons = nc_fid.variables['lon'][:]
	lats = nc_fid.variables['lat'][:]
	att_value = nc_fid.variables[att_name][:]
	att_value = np.multiply(att_value,oceanmask)
	lons,lats,att_value_clipped = range_clip(region[0],region[1],region[2],region[3],lons,lats,att_value)
	time_series = time_seeries_of_spatial_mean(time_s, time_e, time, att_value_clipped)
	# prep_att_ref[att_name] = time_series[0]
	time_series = time_series -prep_att_ref[att_name]
	
	ax1=plt.subplot(spst[0],spst[1],sub_plot)
	ax1.plot(time,movingaverage(time_series),'-',color="blue",linewidth=2,label='RCP8.5_fixA')
	ax1.set_xlabel('Year',fontsize=10)
	ax1.set_xlim([2000,2100])
	# ax1.set_ylim([4,7.5])
	ax1.set_ylabel(prep_unit_dic[att_name].title(),fontsize=15)
	# ax1.set(aspect=10)
	ax1.set_title(att_name.title(),fontsize=15)
	ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
	ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
	sub_plot = sub_plot+1
nc_fid.close()

# legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(0.95, 0.8),fontsize=20)	
# plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/time_series_global.tif')

'''

