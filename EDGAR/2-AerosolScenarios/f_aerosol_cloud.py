"""
This is to plot the TS  from the EDGAR six experiments
data inputs are forty years of monthly TS
	first step is to process the monthly TS into monthly mean
	second step is to process annual mean of monthly mean
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d  as interp2d
from scipy.stats import mannwhitneyu as man_test
from scipy.stats import ttest_ind as student_test
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

def spatial_figure(axs,data,p_value,lons,lats,colormap,colorbar_min,colorbar_max,tb_lef=True,tb_bot=True): #c_bad,c_under,c_over,c_number=20,
	"""
	input : all parameters and data rel;ated to the figure you want to plot_title
	output : a spatial map of the data
	"""
	# plt.imshow(p_value);plt.show()
	lons[lons>180]-=360; 
	# calculate the origin of the map
	lon_0 = lons.mean(); 
	lat_0 = lats.mean(); 
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 60; lat_bin = 30
	map = Basemap(lat_0=0, lon_0=0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs,projection='cyl')
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# s = map.pcolor(xi, yi, data)
	
	if tb_lef:
		map.drawparallels(np.arange(round(lat_b,0)-lat_bin, round(lat_e,0)+lat_bin, lat_bin), labels=[1,0,0,0],linewidth=0.0,fontsize=8)
	if tb_bot:
		map.drawmeridians(np.arange(round(lon_b,0), round(lon_e,0)+lon_bin, lon_bin), labels=[0,0,0,1],linewidth=0.0,fontsize=8)
	# Add Coastlines, States, and Country Boundaries
	map.drawcoastlines(); #map.drawcountries() #map.drawstates(); #
	masked_obj = np.ma.masked_where(np.isnan(data), data)
	# masked_obj = maskoceans(lon,lat,masked_obj)
	cmap = discrete_cmap(10,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); cmap.set_over('k'); cmap.set_under('darkblue');
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)
	pm = np.ma.masked_not_equal(p_value, 1)
	map.pcolor(xi, yi, np.multiply(pm,masked_obj), hatch='+', alpha=0.,lw=0.9,latlon=True)
	return colormesh
		

def monthly_mean_unceer(variable,zonal_mean,zonal_mean_annual):
	def day2datetime(scenario,days):
		"""
		# convert days from a reference into int datetime 
		# do not take leap years into account
		"""
		date_int = np.empty((len(days)));date_int[:]=np.nan
		if scenario =='T1970C': start_year =1970
		else: start_year =2010
		start =(start_year*365)
		ith=0	
		for iday in days:
			month_days =np.array([31,28,31,30,31,30,31,31,30,31,30,31])
			calendar_days = np.array([0,31,59,90,120,151,181,212,243,273,304,334,365])
			total_days = int(iday) + start; 
			year = total_days//365; 
			remainder =  total_days%365
			if remainder ==0: year=year-1;month=12;day=31
			else: 
				month = 1+[layer for layer in range(len(calendar_days)) if calendar_days[layer]< remainder and calendar_days[layer+1]>=remainder][0]
				day = int(remainder - calendar_days[month-1])
				if day == 0: day = month_days[month-1]
			date_int[ith] = year*10000+month*100+day
			ith=ith+1
		return date_int.astype(int)
			
	def data_readin(variable):	
		def annual_and_monthly (scenario,time,data):
			# plt.imshow(data[3,:,:]);plt.show()		
			annual_month=np.empty((30,13,np.shape(data)[1],np.shape(data)[2]));annual_month[:]=np.nan
			# plt.imshow(annual_month[0,3,:,:]);plt.show()
			if scenario=='T1970RCP':
				year_series = range(2020,2050)
			elif scenario=='EdgEne':
				year_series = range(2200,2230)
			elif scenario=='Edg70GO':
				year_series = range(2070,2100)
			else:
				year_series = range(2130,2160)
			for iyear in year_series:
				# print iyear*100+12,time//100 
				if (iyear == year_series[0] and time[0]//100 > year_series[0] *100+1):
					layer_b=0
				else:
					layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  #January
				# print layer_b
				if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
					layer_e=-2
				else:
					layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+12][0]   #december
				data_cache = data[layer_b:layer_e+1,:,:]
				annual_month[iyear-year_series[0],0,:,:] = np.nanmean(data_cache,axis=0)
				# plt.imshow(annual_month[iyear-year_series[0],0,:,:]);plt.show()
				for imonth in range(1,13):
					# print time
					layer = [layer for layer in range(len(time)) if (time[layer]//100  == iyear*100+imonth)]
					# print time[layer]
					annual_month[iyear-year_series[0],imonth,:,:] = data[layer,:,:]
			return annual_month	

		def data_netcdf(scenario,variable):
			input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
			# print var_path
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			days = nc_fid.variables['time'][:];print days
			
			time = day2datetime(scenario,days);print time
			data = nc_fid.variables[variable][:]
			if variable =='CDNUMC' :   #CDNUMC ACTREL TGCLDLWP CLDTOT
				data =data*10**(-9)  # scaled by 10**(-10)
			elif variable =='TGCLDLWP':
				data =data*10**(3)  # kg/m2 to g/m2
			elif variable =='CLDTOT':  
				data =data*100   # convert to percent  
			elif variable.startswith('BURDEN'):
				data =  10**(6)*data
			elif variable =='AODVIS': 
				data[data>1] =np.nan;data[data<0] =np.nan;
			nc_fid.close()
			monthly_map = annual_and_monthly(scenario,time,data)
			return lon,lat,monthly_map
		
		lon,lat,Edg70GO = data_netcdf('Edg70GO',variable)
		lon,lat,T1970 = data_netcdf('T1970RCP',variable)
		_,_,EdgRef = data_netcdf('EdgRef',variable)
		_,_,Edg70Oz = data_netcdf('Edg70Oz',variable)
		_,_,EdgEne = data_netcdf('EdgEne',variable)
		_,_,EdgTech = data_netcdf('EdgTech',variable)
		return lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech
	
	"""
	This function block is to produce a weighted_mask for specific regions (either administrative or box)
	and then produce the spatial mean (weighted)
	"""
	def AreaWeight(lon1,lon2,lat1,lat2):
		'''
		calculate the earth radius in m2
		'''
		radius = 6371000;
		area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
		(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
		return area
		
	def box_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,mask):
		"""
		fill the range outside the box with 0
		"""
		lon = np.array(lon)
		lat = np.array(lat)
		colum_s = [index for index in range(len(lon)) if np.abs(lon-lon_s)[index] == np.min(np.abs(lon-lon_s))][0]
		colum_e = [index for index in range(len(lon)) if np.abs(lon-lon_e)[index] == np.min(np.abs(lon-lon_e))][0]
		row_s = [index for index in range(len(lat)) if np.abs(lat-lat_s)[index] == np.min(np.abs(lat-lat_s))][0]
		row_e = [index for index in range(len(lat)) if np.abs(lat-lat_e)[index] == np.min(np.abs(lat-lat_e))][0]
		if (colum_s> colum_e):
			cache = colum_e; colum_e = colum_s; colum_s = cache;
		if (row_s> row_e):
			cache = row_e; row_e = row_s; row_s = cache;
		mask[:,0:colum_s] =0; mask[:,colum_e:-1] =0
		mask[0:row_s,:] =0; mask[row_e:-1,:] =0
		return mask

	def mask_weight(region_key,lon,lat):
		"""
		Read in the country mask
		interpolate it to the required resolution grids with lon_interp,lat_interp 
		
		"""
		lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
		lons,lats = np.meshgrid (lon,lat)
		area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
		if region_key == 'All':
			mask=area
			mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
		else:
			##OCEAN_MASKS FOR COUNTRIES
			ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
			lon_mask = ocean_mask['lon'][0,:];
			lat_mask = ocean_mask['lat'][0,:];
			box_region_dic={'Land':[0,360,-90,90],'Arctic':[0,360,60,90],'ASIA':[65,145,5,45],'EUS':[265,280,30,50],'EA':[100,145,20,50],'SA':[65,100,5,30],'SESA':[295,315,-40,-25]}
			if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'Globe'):
				mask= ocean_mask[region_key][:]
			elif (region_key == 'Land' or region_key == 'ASIA' or region_key == 'EA' or region_key == 'SA' or region_key == 'SESA' or region_key == 'EUS'):
				mask= ocean_mask['Globe'][:]
				box = box_region_dic[region_key]
				mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
			elif region_key == 'Arctic':
				mask= ocean_mask['Globe'][:]; mask[:] = 1;
				box = box_region_dic[region_key]
				mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
				# plt.imshow(mask);plt.show()
			else:
				print "error region name"
			# interpolate from 360*720 to 192*288
			mask[np.isnan(mask)]=0;	mask[mask>0]=1;
			f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
			mask = f(lon, lat);
			mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
			# weight each grid cell by its area weight against the total area
			mask=np.multiply(mask,area);  
			mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
			# print np.nansum(np.nansum(mask_weighted,axis=1),axis=0)
		return mask_weighted	
	
	def mannwhitneyu_test(vairable1,variable2):
		p_threshold=0.05
		size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
		p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
		# from scipy.stats import ttest_ind as test
		for x in range(size[0]):
			for y in range(size[1]):
				cache1 = vairable1[:,x,y]
				cache2 = variable2[:,x,y]
				if (cache2==cache1).all(): p_value[x,y] = np.nan
				else: _,p_value[x,y] = man_test(cache1,cache2);
		p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
		return p_value
	####calculate the difference and their significance
	def diff_sig(vairable1,variable2):
		dif = np.nanmean(vairable1,axis=0)-np.nanmean(variable2,axis=0)
		sig = mannwhitneyu_test(vairable1, variable2)
		return dif,sig
		
	if variable == "CF":	  # cloud forcing
		lon,lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_readin('SWCF');
		lon,lat,T1970_L,Edg70GO_L,Edg70Oz_L,EdgRef_L,EdgEne_L,EdgTech_L = data_readin('LWCF');
		T1970=(T1970_S+T1970_L)
		Edg70GO=(Edg70GO_S+Edg70GO_L) 
		Edg70Oz=(Edg70Oz_S+Edg70Oz_L)
		EdgRef=(EdgRef_S+EdgRef_L)
		EdgEne=(EdgEne_S+EdgEne_L)
		EdgTech=(EdgTech_S+EdgTech_L)
	elif variable == "SFCSN":
		index = 'FSNSC'; lon,lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_readin(index);
		index = 'FLNSC'; _,lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_readin(index);
		T1970 = T1970_S-T1970_l
		Edg70GO = Edg70GO_S-Edg70GO_l
		EdgRef = EdgRef_S-EdgRef_l
		EdgEne = EdgEne_S-EdgEne_l
		EdgTech = EdgTech_S-EdgTech_l	
		
	elif variable == "SFCF":     # SURFACE CLOUD FORCING
		index = 'FSNS'; lon,lat,T1970_S,Edg70GO_S,Edg70Oz_S,EdgRef_S,EdgEne_S,EdgTech_S = data_readin(index);
		index = 'FLNS'; _,lat,T1970_l,Edg70GO_l,Edg70Oz_l,EdgRef_l,EdgEne_l,EdgTech_l = data_readin(index);
		index = 'FSNSC'; _,lat,T1970_SC,Edg70GO_SC,Edg70Oz_SC,EdgRef_SC,EdgEne_SC,EdgTech_SC = data_readin(index);
		index = 'FLNSC'; _,lat,T1970_lC,Edg70GO_lC,Edg70Oz_lC,EdgRef_lC,EdgEne_lC,EdgTech_lC = data_readin(index);
		T1970 = (T1970_S-T1970_l) -(T1970_SC-T1970_lC)
		Edg70GO = (Edg70GO_S-Edg70GO_l)-(Edg70GO_SC-Edg70GO_lC)
		EdgRef =( EdgRef_S-EdgRef_l)-(EdgRef_SC-EdgRef_lC)
		EdgEne =( EdgEne_S-EdgEne_l)-(EdgEne_SC-EdgEne_lC)
		EdgTech =( EdgTech_S-EdgTech_l)-(EdgTech_SC-EdgTech_lC)
	else:
		lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_readin(variable);
	if zonal_mean:
		if zonal_mean_annual:
			def decompose_mean_std(var1,var2):
				diff = np.nanmean(var1 - var2,axis=2)
				# print "diff", np.shape(diff)
				mean = np.nanmean(diff,axis=0)
				P25 = mean - np.abs(mean-np.nanpercentile(diff,25,axis=0))/1.8
				P75 = mean + np.abs(np.nanpercentile(diff,75,axis=0)-mean)/1.8
				return mean,P25,P75
			lon,lat,T1970,Edg70GO,Edg70Oz,EdgRef,EdgEne,EdgTech = data_readin(variable); #,EdgEne,EdgTech
			mean_TAs,P25_TAs,P75_TAs = 	decompose_mean_std(EdgRef[:,0,:,:],T1970[:,0,:,:])
			mean_AAs,P25_AAs,P75_AAs = 	decompose_mean_std(Edg70GO[:,0,:,:],T1970[:,0,:,:])
			mean_Ene,P25_Ene,P75_Ene = 	decompose_mean_std(EdgRef[:,0,:,:],EdgEne[:,0,:,:])
			mean_Tech,P25_Tech,P75_Tech = decompose_mean_std(EdgRef[:,0,:,:],EdgTech[:,0,:,:])
			
			ref1970,_,_ = decompose_mean_std(T1970,0)
			
			TAs = np.concatenate((mean_TAs,P25_TAs,P75_TAs),axis=0)
			AAs = np.concatenate((mean_AAs,P25_AAs,P75_AAs),axis=0)
			Ene = np.concatenate((mean_Ene,P25_Ene,P75_Ene),axis=0)
			Tech = np.concatenate((mean_Tech,P25_Tech,P75_Tech),axis=0)
			return lat, AAs, Ene,Tech
		else:
			AAs = np.nanmean(Edg70GO-T1970,axis=3)
			ENE= np.nanmean(EdgRef-EdgEne,axis=3)
			TECH= np.nanmean(EdgRef-EdgTech,axis=3)
			return lat,AAs,ENE,TECH
	else:
		def spatial_diff_sig(vairable1,variable2):
			"""
			####calculate spatial the difference and their significance
			"""
			def mannwhitneyu_test(vairable1,variable2):
				p_threshold=0.05
				size = np.array([np.shape(variable2)[1],np.shape(variable2)[2]]); 
				p_value = np.empty((size[0],size[1]));p_value[:]=np.nan
				# from scipy.stats import ttest_ind as test
				for x in range(size[0]):
					for y in range(size[1]):
						cache1 = vairable1[:,x,y]
						cache2 = variable2[:,x,y]
						if (cache2==cache1).all(): p_value[x,y] = np.nan
						else: _,p_value[x,y] = man_test(cache1,cache2);
				p_value[p_value>p_threshold]=np.nan;p_value[p_value<=p_threshold]=1
				return p_value
			dif = np.nanmean(vairable1-variable2,axis=0)
			sig = mannwhitneyu_test(vairable1, variable2)
			return dif,sig
		AAs,AAs_s = spatial_diff_sig(Edg70GO[:,0,:,:],T1970[:,0,:,:])
		ENE,ENE_s= spatial_diff_sig(EdgRef[:,0,:,:],EdgEne[:,0,:,:])
		TECH,TECH_s= spatial_diff_sig(EdgRef[:,0,:,:],EdgTech[:,0,:,:])
		return lon,lat,AAs,ENE,TECH,AAs_s,ENE_s,TECH_s


	
def plot_zonal_mean_uncertainty(ax,x,y1,y2,y3,ylim):
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
	ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
	ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
	ax.set_xlim([-90,90]);ax.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
	ax.set_xticklabels(('90S','70S','50S','30S','10S','EQ','10N','30N','50N','70N','90N'));
	ax.set_ylim([-ylim,ylim]);ax.set_yticks(np.arange(-ylim,ylim+.00001,ylim/4.0))
	# ax.plot(x,y1[0:192],'-',color="k",linewidth=3,label = 'Historical (2010-1970)')
	ax.plot(x,y1[0:192],'-',color="b",linewidth=3,label = 'BEoA (2010-1970)')
	ax.plot(x,y2[0:192],'-',color="r",linewidth=3,label = 'Energy Consumption')
	ax.plot(x,y3[0:192],'-',color="g",linewidth=3,label = 'Technology Advancements')
	
	# ax.fill_between(x,y1[192:384],y1[384:576],facecolor="k",alpha=0.15)
	# ax.fill_between(x,y1[192:384],y1[384:576],facecolor="b",alpha=0.15)
	# ax.fill_between(x,y2[192:384],y2[384:576],facecolor="r",alpha=0.15)
	# ax.fill_between(x,y3[192:384],y3[384:576],facecolor="g",alpha=0.15)
	
	ax2 = ax.twinx();
	# # ax2.plot(x,y4,'-',color="k",linewidth=3, label = "Historical (2010-1970")
	# # ax2.plot(x,ref,'-',color="m",linewidth=4,label = 'E1970')
	# # ax2.fill_between(x,ref[192:384],ref[384:576],facecolor="k",alpha=0.15)
	ax2.set_xlim([-90,90]);ax2.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
	ax2.set_ylim([-ylim,ylim]);ax2.set_yticks(np.arange(-ylim,ylim+.000001,ylim/4.0))
	return ax
	

	
def calculation_R2(a,b):
	iteral = np.shape(a)[-1]
	r2 = np.empty((iteral));r2[:] =np.nan;
	p_value= np.empty((iteral));p_value[:] =np.nan;
	for itimes in range(iteral):
		a_ = np.reshape(a[:,:,itimes],12*40); 
		b_ = np.reshape(b[:,:,itimes],12*40)
		mask = ~np.ma.masked_invalid(a_).mask
		# r2[itimes] = np.corrcoef(a_[mask],b_[mask])[0,1]
		r2[itimes],p_value[itimes] = stats.pearsonr(a_[mask],b_[mask])
	p_value[p_value>=0.0001]=np.nan;
	p_value[p_value<0.0001] =1
	return r2,p_value


	
##################zonal mean cloud features ##################
"""
fig = plt.figure(facecolor='White',figsize=[18.5,12.5]);pad= 5  #plot_setup();

index ='CDNUMC'; lat, AAs, Ene,Tech = monthly_mean_unceer(index,zonal_mean=True,zonal_mean_annual=True);
print np.shape(AAs)
ax = plt.subplot(2,2,2);
ax.annotate('(b) Vertically-integrated droplet concentration'+r'($\mathrm{\mathsf{10^{9}}}$)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
ax = plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,ylim=4)	
legend = ax.legend(shadow=False,ncol=1,loc ='lower left',fontsize=20)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)
		
index ='TGCLDLWP'; lat, AAs, Ene,Tech = monthly_mean_unceer(index,zonal_mean=True,zonal_mean_annual=True);
ax = plt.subplot(2,2,3);
ax.annotate('(c) Total grid-box cloud liquid water path '+r'($\mathrm{\mathsf{g\/m^{-2}}}$)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,ylim=4)

index ='ACTREL'; lat, AAs, Ene,Tech = monthly_mean_unceer(index,zonal_mean=True,zonal_mean_annual=True);
ax = plt.subplot(2,2,4);
ax.annotate('(d) Cloud top effective droplet radius (micron)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,ylim=0.1)

index ='CLDTOT'; lat, AAs, Ene,Tech = monthly_mean_unceer(index,zonal_mean=True,zonal_mean_annual=True);
ax = plt.subplot(2,2,1);
ax.annotate('(a) Cloud fraction (%)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
plot_zonal_mean_uncertainty(ax,lat, AAs, Ene,Tech,ylim=3)
plt.subplots_adjust(left=0.05, bottom=0.03, right=0.95, top=0.95, wspace=0.2, hspace=0.15); 
plt.savefig('zonal_mean_aerosol_cloud.png', format='png', dpi=1000)
"""
##################################### CC  feedbacks ###############################
"""
ax = plt.subplot(2,2,2);
index ='SFCSN';     AAs_R,Ene_R,Tech_R = monthly_mean_unceer(index,zonal_mean=True); 
index ='AODVIS'; AAs_w,Ene_w,Tech_w = monthly_mean_unceer(index,zonal_mean=True); 


ax.annotate('CC of Clear-sky surface net radiation vs. AOD',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
AAs,AAs_p = calculation_R2(AAs_w,AAs_R);
Ene,Ene_p = calculation_R2(Ene_w,Ene_R);
Tech,Tech_p= calculation_R2(Tech_w,Tech_R);

sig = np.multiply(np.multiply(AAs_p,Ene_p),Tech_p); #print sig
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
ax.set_xticklabels(('90S','70S','50S','30S','10S','EQ','10N','30N','50N','70N','90N'));
ax.set_ylim([-1,1]);ax.set_yticks(np.arange(-1,1.1,0.25))
# ax.axvspan(lat[20],lat[-14], alpha=0.2, color='grey',edgecolor='none')

ax.plot(lat,AAs,'-',color="b",linewidth=3,label = 'BEoA (2010-1970)')
ax.plot(lat,Ene,'-',color="r",linewidth=3,label = 'Energy Consumption')
ax.plot(lat,Tech,'-',color="g",linewidth=3,label = 'Technology Advancements')
ax.plot(np.multiply(lat,sig),sig*0.8,marker="*",markeredgecolor='g',color="g",ms=8)
ax2 = ax.twinx();
ax2.set_ylim([-1,1]);ax2.set_yticks(np.arange(-1,1.1,0.25))
ax = plt.subplot(2,2,4);
index ='SFCF';     AAs_R,Ene_R,Tech_R = monthly_mean_unceer(index,zonal_mean=True); 
index ='TGCLDLWP'; AAs_w,Ene_w,Tech_w = monthly_mean_unceer(index,zonal_mean=True); 

ax.annotate('CC of surface cloud forcing vs. cloud liquid water path',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)				
AAs,AAs_p = calculation_R2(AAs_w,AAs_R);
Ene,Ene_p = calculation_R2(Ene_w,Ene_R);
Tech,Tech_p= calculation_R2(Tech_w,Tech_R);

sig = np.multiply(np.multiply(AAs_p,Ene_p),Tech_p)
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=0,color='k',linewidth = 1);
ax.axvline(x=-30,ls='--',color='k',linewidth = 1);ax.axvline(x=30,ls='--',color='k',linewidth = 1);
ax.axvline(x=-60,ls='--',color='k',linewidth = 1);ax.axvline(x=60,ls='--',color='k',linewidth = 1);
ax.set_xlim([-90,90]);ax.set_xticks((-90,-70,-50,-30,-10,0,10,30,50,70,90));
ax.set_xticklabels(('90S','70S','50S','30S','10S','EQ','10N','30N','50N','70N','90N'));
ax.set_ylim([-1,1]);ax.set_yticks(np.arange(-1,1.1,0.25))
# ax.axvspan(lat[20],lat[-14], alpha=0.2, color='grey',edgecolor='none')

ax.plot(lat,AAs,'-',color="b",linewidth=3,label = 'BEoA (2010-1970)')
ax.plot(lat,Ene,'-',color="r",linewidth=3,label = 'Energy Consumption')
ax.plot(lat,Tech,'-',color="g",linewidth=3,label = 'Technology Advancements')
ax.plot(np.multiply(lat,sig),sig*0.8,marker="*",markeredgecolor='g',color="g",ms=8)
ax2 = ax.twinx();
ax2.set_ylim([-1,1]);ax2.set_yticks(np.arange(-1,1.1,0.25))


plt.subplots_adjust(left=0.05, bottom=0.03, right=0.95, top=0.95, wspace=0.2, hspace=0.15); 
plt.savefig('zonal_mean_aerosol_response_vs_feedbacks.png', format='png', dpi=1000)	
"""
"""
#####################
##  Scatter plots   #
#####################
def select_bands(a,b,band_s,band_e):
	a_ = np.reshape(a[:,:,band_s:band_e],12*40*22)
	b_ = np.reshape(b[:,:,band_s:band_e],12*40*22)
	mask = ~np.multiply(np.ma.masked_invalid(a_).mask,np.ma.masked_invalid(a_).mask)
	return a_[mask],b_[mask],round(np.corrcoef(a_[mask],b_[mask])[0,1],1)
fig = plt.figure(facecolor='White',figsize=[16.5,12.5]);pad= 5  #plot_setup();
index ='SFCSN';  AAs_R,Ene_R,Tech_R = monthly_mean_unceer(index,zonal_mean=True); 
index ='AODVIS'; AAs_w,Ene_w,Tech_w = monthly_mean_unceer(index,zonal_mean=True); 

ax = plt.subplot(2,2,1);
			
AAs_x,AAs_y,AA_cc = select_bands(AAs_w,AAs_R,85,107); 
Ene_x,Ene_y,Ene_cc = select_bands(Ene_w,Ene_R,85,107);
Tech_x,Tech_y,Tech_cc= select_bands(Tech_w,Tech_R,85,107);
ax.annotate('(a) FSNS vs. AOD over Tropics (10S - 10N)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
cc = np.array([AA_cc,Ene_cc,Tech_cc]); colors=['b','r','g']			
ax.annotate('CC: ',xy=(0.60,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color='k',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)			
for i in range(3):
	ax.annotate(' '+str(cc[i]),xy=(0.68+i*0.10,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color=colors[i],
					ha='left', va='baseline',rotation='horizontal',fontsize=15)	

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.plot(AAs_x,AAs_y, ls='',marker='s',color="b",ms=5,markeredgecolor='b',label = 'BEoA (2010-1970)')
ax.plot(Ene_x,Ene_y, ls='',marker='s',color="r",ms=5,markeredgecolor='r',label = 'Energy Consumption')
ax.plot(Tech_x,Tech_y,ls='',marker='s',color="g",ms=5,markeredgecolor='g',label = 'Technology Advancements')
ax.axhline(y=0,color='k',linewidth = 3);ax.axvline(x=0,color='k',linewidth = 3);
plt.xlabel('Aerosol optical depth',fontsize = 15);plt.ylabel('Clear-sky net radiative flux at surface '+r'($\mathrm{\mathsf{W\/m^{-2}}}$)',fontsize = 15)
ax.set_xlim([-0.15,0.15]);ax.set_ylim([-20,20]);

legend = ax.legend(shadow=False,ncol=1,loc ='upper left',fontsize=15)	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0)


ax = plt.subplot(2,2,2);
			
AAs_x,AAs_y,AA_cc = select_bands(AAs_w,AAs_R,170,192); 
Ene_x,Ene_y,Ene_cc = select_bands(Ene_w,Ene_R,170,192);
Tech_x,Tech_y,Tech_cc= select_bands(Tech_w,Tech_R,170,192);
ax.annotate('(b) FSNS vs. AOD Arctic (70 - 90N)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)	
cc = np.array([AA_cc,Ene_cc,Tech_cc]); colors=['b','r','g']			
ax.annotate('CC: ',xy=(0.60,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color='k',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)			
for i in range(3):
	ax.annotate(' '+str(cc[i]),xy=(0.68+i*0.10,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color=colors[i],
					ha='left', va='baseline',rotation='horizontal',fontsize=15)	

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.plot(AAs_x,AAs_y, ls='',marker='s',color="b",ms=5,markeredgecolor='b',label = 'BEoA (2010-1970)')
ax.plot(Ene_x,Ene_y, ls='',marker='s',color="r",ms=5,markeredgecolor='r',label = 'Energy Consumption')
ax.plot(Tech_x,Tech_y,ls='',marker='s',color="g",ms=5,markeredgecolor='g',label = 'Technology Advancements')
ax.axhline(y=0,color='k',linewidth = 3);ax.axvline(x=0,color='k',linewidth = 3);
plt.xlabel('Aerosol optical depth',fontsize = 15);plt.ylabel('Clear-sky net radiative flux at surface '+r'($\mathrm{\mathsf{W\/m^{-2}}}$)',fontsize = 15)
ax.set_xlim([-0.03,0.03]);ax.set_ylim([-100,100]);


index ='SFCF';     AAs_R,Ene_R,Tech_R = monthly_mean_unceer(index,zonal_mean=True); 
index ='TGCLDLWP'; AAs_w,Ene_w,Tech_w = monthly_mean_unceer(index,zonal_mean=True);
 
ax = plt.subplot(2,2,3);
				
AAs_x,AAs_y,AA_cc = select_bands(AAs_w,AAs_R,85,107); 
Ene_x,Ene_y,Ene_cc = select_bands(Ene_w,Ene_R,85,107);
Tech_x,Tech_y,Tech_cc= select_bands(Tech_w,Tech_R,85,107);
ax.annotate('(c) CRFNS vs. LWP over Tropics (10S - 10N)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)
cc = np.array([AA_cc,Ene_cc,Tech_cc]); colors=['b','r','g']			
ax.annotate('CC: ',xy=(0.60,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color='k',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)			
for i in range(3):
	ax.annotate(' '+str(cc[i]),xy=(0.68+i*0.10,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color=colors[i],
					ha='left', va='baseline',rotation='horizontal',fontsize=15)	
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.plot(AAs_x,AAs_y, ls='',marker='s',color="b",ms=5,markeredgecolor='b',label = 'BEoA (2010-1970)')
ax.plot(Ene_x,Ene_y, ls='',marker='s',color="r",ms=5,markeredgecolor='r',label = 'Energy Consumption')
ax.plot(Tech_x,Tech_y,ls='',marker='s',color="g",ms=5,markeredgecolor='g',label = 'Technology Advancements')
ax.axhline(y=0,color='k',linewidth = 3);ax.axvline(x=0,color='k',linewidth = 3);
plt.xlabel('Cloud liquid water parth '+r'($\mathrm{\mathsf{g\/m^{-2}}}$)',fontsize = 15);plt.ylabel('Cloud radiative forcing at surface '+r'($\mathrm{\mathsf{W\/m^{-2}}}$)',fontsize = 15)
ax.set_xlim([-40,40]);ax.set_xticks(np.arange(-40,41,20));
ax.set_ylim([-40,40]);ax.set_yticks(np.arange(-40,41,20));
			
AAs_x,AAs_y,AA_cc = select_bands(AAs_w,AAs_R,170,192); 
Ene_x,Ene_y,Ene_cc = select_bands(Ene_w,Ene_R,170,192);
Tech_x,Tech_y,Tech_cc= select_bands(Tech_w,Tech_R,170,192);

ax = plt.subplot(2,2,4);
ax.annotate('(d) CRFNS vs. LWP over Arctic (70 - 90N)',xy=(0.02,1.01), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=15)

cc = np.array([AA_cc,Ene_cc,Tech_cc]); colors=['b','r','g']			
ax.annotate('CC: ',xy=(0.60,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color='k',
					ha='left', va='baseline',rotation='horizontal',fontsize=15)			
for i in range(3):
	ax.annotate(' '+str(cc[i]),xy=(0.68+i*0.10,0.90), xytext=(0, pad),   #\mathrm{\mathsf{\Delta} VT\/( K\/m\/s^{-1})}
					xycoords='axes fraction', textcoords='offset points',color=colors[i],
					ha='left', va='baseline',rotation='horizontal',fontsize=15)	

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.plot(AAs_x,AAs_y, ls='',marker='s',color="b",ms=5,markeredgecolor='b',label = 'BEoA (2010-1970)')
ax.plot(Ene_x,Ene_y, ls='',marker='s',color="r",ms=5,markeredgecolor='r',label = 'Energy Consumption')
ax.plot(Tech_x,Tech_y,ls='',marker='s',color="g",ms=5,markeredgecolor='g',label = 'Technology Advancements')
ax.axhline(y=0,color='k',linewidth = 3);ax.axvline(x=0,color='k',linewidth = 3);
plt.xlabel('Cloud liquid water parth '+r'($\mathrm{\mathsf{g\/m^{-2}}}$)',fontsize = 15);plt.ylabel('Cloud radiative forcing at surface '+r'($\mathrm{\mathsf{W\/m^{-2}}}$)',fontsize = 15)
ax.set_xlim([-40,40]);ax.set_xticks(np.arange(-40,41,20));
ax.set_ylim([-40,40]);ax.set_yticks(np.arange(-40,41,20));
# x=np.arange(1,366); y =10*x/366
# ax.plot(x,y,'-',color="w",linewidth=0)
# ax.set_xlim([1,366]);ax.set_ylim([0,10]);
# # line_k = mlines.Line2D(x,y,ls='-',color='k',lw=2)
# line_b = mlines.Line2D(x,y,ls='-',color='b',lw=2)
# line_r = mlines.Line2D(x,y,ls='-',color='r',lw=2)
# line_g = mlines.Line2D(x,y,ls='-',color='g',lw=2)

# lines = [line_b,line_r,line_g]
# labels = ['BEoA (2010 - 1970)','Energy consumption','Technology advancements']
# legend = plt.legend(lines,labels,ncol=1,loc ='center',labelspacing=1.5,markerscale =10)
# legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)


plt.subplots_adjust(left=0.10, bottom=0.07, right=0.95, top=0.95, wspace=0.25, hspace=0.25); 
plt.savefig('scatterplots_raiation_aerosol.png', format='png', dpi=1000)	
"""

######################spatial maps ###########################################

fig = plt.figure(facecolor='White',figsize=[22.5,14]);pad= 5  #plot_setup();
colormap='RdBu_r';  #colormap = reverse_colourmap(colormap);


index ='CLDTOT'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT
colorbar_min=-3;colorbar_max=3;
lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s = monthly_mean_unceer(index,zonal_mean=False,zonal_mean_annual=False);

ax = plt.subplot(4,3,1);
ax.annotate('(a)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('BEoA (2010-1970)',xy=(0.5,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)	
# ax.annotate(r'CDNC ($\mathrm{\mathsf{10^{10}\/m^{-2}}}$)',xy=(-0.1,0.5), xytext=(0, pad),
ax.annotate('Cloud (%)',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)	
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,2);
ax.annotate('(e)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Energy Consumption',xy=(0.5,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)		
colormesh1 = spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,3);
ax.annotate('(i)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('Technology advancements',xy=(0.5,1.03), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=20)		
colormesh1 = spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.92, 0.75, 0.015, 0.20])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))


index ='CDNUMC'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT
colorbar_min= -5 ;colorbar_max=5;
lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s = monthly_mean_unceer(index,zonal_mean=False,zonal_mean_annual=False);
ax = plt.subplot(4,3,4);
ax.annotate('(b)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('CDNC' +r'($10^{9}$)',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,5);
ax.annotate('(f)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,6);
ax.annotate('(j)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.92, 0.51, 0.015, 0.20])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

index ='TGCLDLWP'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT
colorbar_min=-10;colorbar_max=10;
lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s= monthly_mean_unceer(index,zonal_mean=False,zonal_mean_annual=False);
ax = plt.subplot(4,3,7);
ax.annotate('(c)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate(r'CLWP ($\mathrm{\mathsf{g\/m^{-2}}}$)',xy=(-0.1, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)	
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,8);
ax.annotate('(g)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,9);
ax.annotate('(k)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

cbar_ax = fig.add_axes([0.92, 0.27, 0.015, 0.20])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))

index ='ACTREL'	  #F CDNUMC ACTREL TGCLDLWP CLDTOT
colorbar_min=-0.15;colorbar_max=0.15;
lon,lat,AEROSOL,ENE,TECH,AEROSOL_s,ENE_s,TECH_s = monthly_mean_unceer(index,zonal_mean=False,zonal_mean_annual=False);
ax = plt.subplot(4,3,10);
ax.annotate('(d)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)	
ax.annotate('CTER (micron)',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90,fontsize=20)
spatial_figure(ax,AEROSOL,AEROSOL_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,11);
ax.annotate('(h)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = spatial_figure(ax,ENE,ENE_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)

ax = plt.subplot(4,3,12);
ax.annotate('(l)',xy=(0.02,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=20)			
colormesh1 = spatial_figure(ax,TECH,TECH_s,lon,lat,colormap,colorbar_min,colorbar_max,tb_lef=False,tb_bot=False)
cbar_ax = fig.add_axes([0.92, 0.04, 0.015, 0.20])
char = fig.colorbar(colormesh1,orientation='vertical',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.1,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))



plt.subplots_adjust(left=0.07, bottom=0.03, right=0.90, top=0.95, wspace=0.04, hspace=0.15); 
plt.savefig('cloud_interactions_Decompose.png', format='png', dpi=1000)

