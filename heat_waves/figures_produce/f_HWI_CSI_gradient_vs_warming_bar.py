# -*- coding: utf-8 -*-
'''
This is to plot the gradients of HW and CS magnitude against one degree of warming from 
	hiS (1920-2005) GHG(RCp8.5_FixA) GHG+AEROSOL (RCP8.5) and AEROSOL (RCP8.5-FixA)
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import math
import scipy.io as sio
from scipy.interpolate import interp2d  as interp2d

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_720_360.mat')
lon_mask = ocean_mask['lon'][0,:];


def get_baseline_data():

	directory ='CalDayThr'
	baseline_data =data_path+'/Temp_pp/'+directory+'/HWI_CSI_1961_1990_0595_'+directory+'.nc' #_InterOsi
	nc_fid = nc4.Dataset(baseline_data,mode='r')
	HWIM= stats.nanmean(nc_fid.variables['HWIM'],axis=0);
	HWIS= stats.nanmean(nc_fid.variables['HWIS'],axis=0);
	return HWIM,HWIS
	

def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# plt.imshow(area);plt.show()
	return area

def mask_match(country_key,lon,lat):
	"""
	Read in the country mask
	interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	"""
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_StAf_USA_AUS_BRICS_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	mask= ocean_mask[country_key][:]
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	# print np.unique(mask)
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); 
	mask = f(lon, lat);
	mask[mask >= 1] = 1;mask[mask < 1] = np.nan;mask[0:27,:]=np.nan
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid (lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)
	mask=np.multiply(mask,area);  #crop the interested region
	mask=np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
	return mask
	

def HWI_CSI_gradient_pp(country_key,Regional_T=False):
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	def FilterAnomoly(data):
		for stats in range(3):
			for ilayer in range(np.shape(data)[1]):
				cache = np.reshape(data[stats,ilayer,:,:],192*288)
				P99= np.nanpercentile(cache,90)
				cache = data[stats,ilayer,:,:];cache[cache>P99]=np.nan
				data[stats,ilayer,:,:] =cache
		return data
	#rcp85
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_rcp85.nc',mode='r')
	lon = nc_fid.variables['lon'][:]
	lat = nc_fid.variables['lat'][:]
	mask = mask_match(country_key,lon,lat)
	# plt.imshow(mask,origin='lower');plt.show()
	mask_globe = mask_match('Globe',lon,lat)
	HWI_M_R = np.nansum(np.nansum(np.multiply(FilterAnomoly(nc_fid.variables['HWI_M'][[0,4,5],:,:,:]),mask),axis=3),axis=2)
	CSI_M_R = np.nansum(np.nansum(np.multiply(FilterAnomoly(nc_fid.variables['CSI_M'][[0,4,5],:,:,:]),mask),axis=3),axis=2)
	# plt.imshow(nc_fid.variables['HWI_M'][0,94,:,:],origin='lower');plt.show()
	# print country_key,HWI_M_R[0]
	nc_fid.close()
	#fixa
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_fixa.nc',mode='r')
	HWI_M_F = np.nansum(np.nansum(np.multiply(FilterAnomoly(nc_fid.variables['HWI_M'][[0,4,5],:,:,:]),mask),axis=3),axis=2)
	CSI_M_F = np.nansum(np.nansum(np.multiply(FilterAnomoly(nc_fid.variables['CSI_M'][[0,4,5],:,:,:]),mask),axis=3),axis=2)
	nc_fid.close()
	
	## aerosol
	HWI_M_A= HWI_M_R-HWI_M_F; #max=np.nanmax(temp,axis=0);min=np.nanmin(temp,axis=0);median=np.nanmedian(temp,axis=0);HWI_M_A=np.array([median,max,min])
	CSI_M_A= CSI_M_R-CSI_M_F; #max=np.nanmax(temp,axis=0);min=np.nanmin(temp,axis=0);median=np.nanmedian(temp,axis=0);CSI_M_A=np.array([median,max,min])

	#TS
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_TS_200602_210101.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	R85=nc_fid.variables['rcp85'];
	FiX=nc_fid.variables['rcp85_fixA'];
	rcp85_T=np.empty((95,192,288)); rcp85_T[:]=np.nan
	FixA_T=np.empty((95,192,288)); FixA_T[:]=np.nan
	for iyear in range(0,95):
		rcp85_T[iyear,:,:] = stats.nanmean(R85[iyear*12:(iyear+1)*12,:,:],axis=0)
		FixA_T[iyear,:,:] = stats.nanmean(FiX[iyear*12:(iyear+1)*12,:,:],axis=0)
	Aerosol_T = rcp85_T-FixA_T
	
	# use global sat or regional sat to regress
	if Regional_T: mask_TS= mask
	else:  mask_TS = mask_globe
		
	TS_R = np.nansum(np.nansum(np.multiply(rcp85_T,mask_TS),axis=2),axis=1)      #RCP85
	TS_F = np.nansum(np.nansum(np.multiply(FixA_T,mask_TS),axis=2),axis=1)       #FIXA
	TS_A = np.nansum(np.nansum(np.multiply(Aerosol_T,mask_TS),axis=2),axis=1)  # AEROSOL 
	nc_fid.close();

	### historical (1920-2005)	
	file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/ensumble_mean_TS_192002_200512.nc'
	nc_fid = nc4.Dataset(file_name,mode='r');
	his=nc_fid.variables['TS'];
	##TS
	his_T=np.empty((86,192,288)); his_T[:]=np.nan
	for iyear in range(0,86):
		his_T[iyear,:,:] = stats.nanmean(his[iyear*12:(iyear+1)*12,:,:],axis=0)
	TS_H = np.nansum(np.nansum(np.multiply(his_T,mask_TS),axis=2),axis=1)[51:86]      #R
	nc_fid.close();	
	file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/CalDayThr/'
	nc_fid = nc4.Dataset(file_path+'HWI_CSI_EnPP_his.nc',mode='r')
	HWI_M_H = np.nansum(np.nansum(np.multiply(nc_fid.variables['HWI_M'][[0,4,5],-36:-1,:,:],mask),axis=3),axis=2)
	CSI_M_H = np.nansum(np.nansum(np.multiply(nc_fid.variables['CSI_M'][[0,4,5],-36:-1,:,:],mask),axis=3),axis=2)
	nc_fid.close()

	TS =np.array([TS_R,TS_F,TS_A])
	HWI_M =np.array([HWI_M_R,HWI_M_F,HWI_M_A])
	CSI_M =np.array([CSI_M_R,CSI_M_F,CSI_M_A])
	return TS_H,HWI_M_H,CSI_M_H,TS,HWI_M,CSI_M
	
def movingaverage (values, window=5):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.nanmean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.nanmean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.nanmean(values[index-boudary:index+boudary+1])
	return result

def statistics_print(country_key,Regional_T=False):
	def li_fit(TS,variable):
			x = movingaverage (TS, window=5); y =movingaverage(variable[0,:], window=5); mask = ~np.isnan(x) & ~np.isnan(y); #ensemble mean
			slope_m, intercept_m, _, _, _ = stats.linregress(x[mask],y[mask]);
			x = movingaverage (TS, window=5); y =movingaverage (variable[1,:], window=5); mask = ~np.isnan(x) & ~np.isnan(y); #75p
			slope_u, intercept_u, _, _, _ = stats.linregress(x[mask],y[mask]);
			x = movingaverage (TS, window=5); y = movingaverage(variable[2,:], window=01); mask = ~np.isnan(x) & ~np.isnan(y);  #25p
			slope_l, intercept_l, _, _, _ = stats.linregress(x[mask],y[mask]);	
			slope=np.array([slope_m,slope_u,slope_l]);
			# intercept=np.array([intercept_m,intercept_u,intercept_l]);
			# print np.shape(stat)
			return x,slope
	def fit_print(TS,variable):
		_,slope_rcp=li_fit(TS[0,:],variable[0,:,:])	
		# print 'GHG';print slope_rcp
		_,slope_fix=li_fit(TS[1,:],variable[1,:,:])	
		# print 'GHG';print slope_fix
		_,slope_aer=li_fit(TS[2,:],variable[2,:,:])
		# print 'AER';print slope_aer
		return slope_rcp,slope_fix,slope_aer
	his_T,HWI_M_H,CSI_M_H,TS,HWI_M,CSI_M = HWI_CSI_gradient_pp(country_key,Regional_T) 
	_,slope_his=li_fit(his_T,HWI_M_H)
	# print '-----HW---------';
	slope_rcp,slope_fix,slope_aer=fit_print(TS,HWI_M);
	HW_SLOPE=np.array([slope_his,slope_rcp,slope_fix,slope_aer])
	_,slope_his=li_fit(his_T,CSI_M_H)
	# print '-----CS---------';
	slope_rcp,slope_fix,slope_aer=fit_print(TS,CSI_M);
	CS_SLOPE=np.array([slope_his,slope_rcp,slope_fix,slope_aer])
	# print slope_rcp
	return HW_SLOPE,CS_SLOPE

def autolabel(X_pos,values,height_lift):
	"""
	Attach a text label above each bar displaying its height
	"""
	height= np.round(np.nan_to_num(values),2);y_pos = height_lift*height
	for i in range(len(height)):
		ax.text(X_pos[i],y_pos[i],'%4.2f' % height[i], ha='center', va='bottom',size=4)

fig = plt.figure(facecolor='White',figsize=[6,6]);plot_setup();pad= 5;
ax = plt.subplot(2,1,1);

	
# print '*************globe*************'	
HW_SLOPE_GLO,CS_SLOPE_GLO = statistics_print('Globe',Regional_T=False)
# print '*************Australia*************'
HW_SLOPE_AUS,CS_SLOPE_AUS = statistics_print('Australia',Regional_T=False)
# print '*************China*************'
HW_SLOPE_CHA,CS_SLOPE_CHA = statistics_print('China',Regional_T=False)
# print '*************Europe*************
HW_SLOPE_EUR,CS_SLOPE_EUR = statistics_print('Europe',Regional_T=False)
# print '*************India*************'
HW_SLOPE_IND,CS_SLOPE_IND = statistics_print('India',Regional_T=False)
# print '*************RUASSIA*************'
HW_SLOPE_BRA,CS_SLOPE_BRA = statistics_print('Brazil',Regional_T=False)
# print '*************USA*************'
HW_SLOPE_USA,CS_SLOPE_USA = statistics_print('USA',Regional_T=False)
# print '*************USA*************'
HW_SLOPE_SAF,CS_SLOPE_SAF = statistics_print('Southern African',Regional_T=False)

HW_SLOPE = np.array([HW_SLOPE_GLO,HW_SLOPE_AUS,HW_SLOPE_BRA,HW_SLOPE_CHA,HW_SLOPE_EUR,HW_SLOPE_IND,HW_SLOPE_SAF,HW_SLOPE_USA])
CS_SLOPE = -1*np.array([CS_SLOPE_GLO,CS_SLOPE_AUS,CS_SLOPE_BRA,CS_SLOPE_CHA,CS_SLOPE_EUR,CS_SLOPE_IND,CS_SLOPE_SAF,CS_SLOPE_USA])

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
# ax.spines['right'].set_visible(False);ax.spines['top'].set_visible(False);ax.spines['top'].set_visible(False);
ax.set_xlim([-0.5,20]);ax.set_xticks(np.arange(0.8,20.5,2.5));ax.set_xticklabels(('GLO','AUS','BRA','CHA','EUR','IND','SAF','USA'));
ax.set_ylim([0,4]);ax.set_yticks(np.arange(0,4.1,1))
X_pos=np.arange(0,20,2.5);dis=0.5
ax.annotate(r'(a) $\mathrm{\mathsf{\Delta}\/Heatwave\/magnitude\/against\/global\/mean\/{\Delta}SAT\/(K^{-1})}$',xy=(0.05,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
				
y_value = HW_SLOPE[:,0,0];yerr=abs(np.array([HW_SLOPE[:,0,0]-HW_SLOPE[:,0,2],HW_SLOPE[:,0,1]- HW_SLOPE[:,0,0]])); 
rects1=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='k', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos,y_value,1.1)

y_value = HW_SLOPE[:,1,0];yerr=abs(np.array([HW_SLOPE[:,1,0]-HW_SLOPE[:,1,2],HW_SLOPE[:,1,1]- HW_SLOPE[:,1,0]])); 
rects2=ax.bar(X_pos+dis, y_value, yerr=yerr, align='center',color='r', ecolor='r',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis,y_value,1.2)

y_value = HW_SLOPE[:,2,0];yerr=abs(np.array([HW_SLOPE[:,2,0]-HW_SLOPE[:,2,2],HW_SLOPE[:,2,1]- HW_SLOPE[:,2,0]])); 
rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr,align='center',color='b', ecolor='b',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis*2,y_value,1.1)

y_value = HW_SLOPE[:,3,0];yerr=abs(np.array([HW_SLOPE[:,3,0]-HW_SLOPE[:,3,2],HW_SLOPE[:,3,1]- HW_SLOPE[:,3,0]])); 
rects4=ax.bar(X_pos+dis*3, y_value, yerr=yerr, align='center',color='m', ecolor='m',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis*3,y_value,1.05)

print np.divide(HW_SLOPE[:,3,0],HW_SLOPE[:,2,0])
print np.divide(HW_SLOPE[:,3,1],HW_SLOPE[:,2,1])
print np.divide(HW_SLOPE[:,3,2],HW_SLOPE[:,2,2])

legend1=ax.legend((rects1[0], rects2[0],rects3[3], rects4[4]),('His (1970-2005)', 'GHGs+AAs (2006-2100)','GHGs (2006-2100)','AAs (2006-2100)'), shadow=False,ncol=1,bbox_to_anchor=(0.4, 1.01))
legend1.get_frame().set_facecolor('white');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)


ax = plt.subplot(2,1,2);

# print '*************globe*************'	
HW_SLOPE_GLO,CS_SLOPE_GLO = statistics_print('Globe',Regional_T=True)
# print '*************Australia*************'
HW_SLOPE_AUS,CS_SLOPE_AUS = statistics_print('Australia',Regional_T=True)
# print '*************China*************'
HW_SLOPE_CHA,CS_SLOPE_CHA = statistics_print('China',Regional_T=True)
# print '*************Europe*************'
HW_SLOPE_EUR,CS_SLOPE_EUR = statistics_print('Europe',Regional_T=True)
# print '*************India*************'
HW_SLOPE_IND,CS_SLOPE_IND = statistics_print('India',Regional_T=True)
# print '*************RUASSIA*************'
HW_SLOPE_BRA,CS_SLOPE_BRA = statistics_print('Brazil',Regional_T=True)	
# print '*************USA*************'
HW_SLOPE_USA,CS_SLOPE_USA = statistics_print('USA',Regional_T=True)
# print '*************USA*************'
HW_SLOPE_SAF,CS_SLOPE_SAF = statistics_print('Southern African',Regional_T=True)

HW_SLOPE = np.array([HW_SLOPE_GLO,HW_SLOPE_AUS,HW_SLOPE_BRA,HW_SLOPE_CHA,HW_SLOPE_EUR,HW_SLOPE_IND,HW_SLOPE_SAF,HW_SLOPE_USA])
CS_SLOPE = -1*np.array([CS_SLOPE_GLO,CS_SLOPE_AUS,CS_SLOPE_BRA,CS_SLOPE_CHA,CS_SLOPE_EUR,CS_SLOPE_IND,CS_SLOPE_SAF,CS_SLOPE_USA])

ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=False,top=False, pad=5)
# ax.spines['right'].set_visible(False);ax.spines['top'].set_visible(False);ax.spines['top'].set_visible(False);
ax.set_xlim([-0.5,20]);ax.set_xticks(np.arange(0.8,20.5,2.5));ax.set_xticklabels(('GLO','AUS','BRA','CHA','EUR','IND','SAF','USA'));
ax.set_ylim([0,4]);ax.set_yticks(np.arange(0,4.1,1))
ax.annotate(r'(b) $\mathrm{\mathsf{\Delta}\/Heatwave\/magnitude\/against\/regional\/mean\/{\Delta}SAT\/(K^{-1})}$',xy=(0.05,1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
X_pos=np.arange(0,20,2.5);dis=0.5

y_value = HW_SLOPE[:,0,0];yerr=abs(np.array([HW_SLOPE[:,0,0]-HW_SLOPE[:,0,2],HW_SLOPE[:,0,1]- HW_SLOPE[:,0,0]])); 
rects1=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='k', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos,y_value,1.1)

y_value = HW_SLOPE[:,1,0];yerr=abs(np.array([HW_SLOPE[:,1,0]-HW_SLOPE[:,1,2],HW_SLOPE[:,1,1]- HW_SLOPE[:,1,0]])); 
rects2=ax.bar(X_pos+dis, y_value, yerr=yerr, align='center',color='r', ecolor='r',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis,y_value,1.2)

y_value = HW_SLOPE[:,2,0];yerr=abs(np.array([HW_SLOPE[:,2,0]-HW_SLOPE[:,2,2],HW_SLOPE[:,2,1]- HW_SLOPE[:,2,0]])); 
rects3=ax.bar(X_pos+dis*2, y_value, yerr=yerr,align='center',color='b', ecolor='b',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis*2,y_value,1.1)

y_value = HW_SLOPE[:,3,0];yerr=abs(np.array([HW_SLOPE[:,3,0]-HW_SLOPE[:,3,2],HW_SLOPE[:,3,1]- HW_SLOPE[:,3,0]])); 
rects4=ax.bar(X_pos+dis*3, y_value, yerr=yerr, align='center',color='m', ecolor='m',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis*3,y_value,1.05)
print np.divide(HW_SLOPE[:,3,0],HW_SLOPE[:,2,0])
print np.divide(HW_SLOPE[:,3,1],HW_SLOPE[:,2,1])
print np.divide(HW_SLOPE[:,3,2],HW_SLOPE[:,2,2])

plt.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95, wspace=0.2, hspace=0.2); 
# plt.show()
plt.savefig('Fig6_HW_gradient_bar_local.png', format='png', dpi=1000)
















"""
ax = plt.subplot(2,1,2);
ax.axhline(y=0,color='k',lw=1)
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=False,bottom=False,top=False, pad=5)
ax.spines['right'].set_visible(False);ax.spines['top'].set_visible(False);ax.spines['top'].set_visible(False);
ax.set_xlim([-0.5,20]);ax.set_xticks(np.arange(0.8,20.5,2.5));ax.set_xticklabels(('GLO','AUS','BRA','CHA','EUR','IND','SAF','USA'));
ax.set_ylim([-1.0,0.5]);ax.set_yticks(np.arange(-1.0,0.51,0.5));
X_pos=np.arange(0,20,2.5);dis=0.5

y_value = CS_SLOPE[:,0,0];yerr=abs(np.array([CS_SLOPE[:,0,0]-CS_SLOPE[:,0,2],CS_SLOPE[:,0,1]- CS_SLOPE[:,0,0]])); 
rects=ax.bar(X_pos, y_value, yerr=yerr, align='center',color='k', ecolor='k',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos,y_value,2.0)

y_value = CS_SLOPE[:,1,0];yerr=abs(np.array([CS_SLOPE[:,1,0]-CS_SLOPE[:,1,2],CS_SLOPE[:,1,1]- CS_SLOPE[:,1,0]])); 
rects=ax.bar(X_pos+dis, y_value, yerr=yerr, align='center',color='r', ecolor='r',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis,y_value,1.5)

y_value = CS_SLOPE[:,2,0];yerr=abs(np.array([CS_SLOPE[:,2,0]-CS_SLOPE[:,2,2],CS_SLOPE[:,2,1]- CS_SLOPE[:,2,0]])); 
rects=ax.bar(X_pos+dis*2, y_value, yerr=yerr,align='center',color='b', ecolor='b',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis*2,y_value,2.0)

y_value = CS_SLOPE[:,3,0];yerr=abs(np.array([CS_SLOPE[:,3,0]-CS_SLOPE[:,3,2],CS_SLOPE[:,3,1]- CS_SLOPE[:,3,0]])); 
rects=ax.bar(X_pos+dis*3, y_value, yerr=yerr,align='center',color='m', ecolor='m',capsize=0,alpha=0.7,width=0.5,lw=0)
autolabel(X_pos+dis*3,y_value,1.0)
"""



"""
def scatterplot_ensemble(TS, variable):
	ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)

	# x_rcp,stats_rcp=li_fit(TS[0,:],variable[0,:,:])
	# height = variable[0,0,:]; yerr=np.array([variable[0,0,:]- variable[0,2,:],variable[0,1,:]-variable[0,0,:]]); 
	# ax.errorbar(x_rcp, height,yerr=yerr,fmt='-',color='r',ecolor='r')
	# ax.scatter(x_rcp,variable[0,0,:], marker='*',color="r",s=3, edgecolor='r')
	# ax.scatter(x_rcp,variable[0,1,:], marker='o',color="r",s=3, edgecolor='r')
	# ax.scatter(x_rcp,variable[0,2,:], marker='s',color="r",s=3, edgecolor='r')
	# ax.plot(x_rcp,x_rcp*stats_rcp[0,0]+stats_rcp[1,0], color="r")


	x_fix,stats_fix=li_fit(TS[1,:],variable[1,:,:])
	height = variable[1,0,:]; yerr=np.array([variable[1,0,:]- variable[1,2,:],variable[1,1,:]-variable[1,0,:]]); 
	ax.errorbar(x_fix, height,yerr=yerr,fmt='-',color='b',ecolor='b')
	
	# ax.scatter(x_fix,variable[1,0,:], marker='*',color="b",s=3, edgecolor='b')
	# ax.scatter(x_fix,variable[1,1,:], marker='o',color="b",s=3, edgecolor='b')
	# ax.scatter(x_fix,variable[1,2,:], marker='s',color="b",s=3, edgecolor='b')
	# ax.plot(x_fix,x_fix*stats_fix[0,0]+stats_fix[1,0], color="b")
	
	x_aer,stats_aer=li_fit(TS[2,:],variable[2,:,:])
	height = variable[2,0,:]; yerr=np.array([variable[2,0,:]- variable[2,2,:],variable[2,1,:]-variable[2,0,:]]); 
	ax.errorbar(x_aer, height,yerr=yerr,fmt='-',color='k',ecolor='k')
	# ax.scatter(x_aer,variable[2,0,:], marker='*',color="k",s=3, edgecolor='k')
	# ax.scatter(x_aer,variable[2,1,:], marker='o',color="k",s=3, edgecolor='k')
	# ax.scatter(x_aer,variable[2,2,:], marker='s',color="k",s=3, edgecolor='k')
	# ax.plot(x_aer,x_aer*stats_aer[0,0]+stats_aer[1,0], color="k")
	# ax.fill_between(x_aer,variable[2,1,:],variable[2,2,:],facecolor="k",alpha=0.7)	


fig = plt.figure(facecolor='White',figsize=[8,13]);plot_setup();pad= 5;
ax = plt.subplot(4,2,1);
ax.annotate('Intensity\n($^\circ$C)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('Heat Waves',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(a)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
scatterplot_ensemble(TS_CHA, HWInM_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));
# ax.set_ylim([0,110]);ax.set_yticks(np.arange(0,110,20));		


ax = plt.subplot(4,2,3);
ax.annotate('Duration\n(days/event)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
scatterplot_ensemble(TS_CHA, HWDuM_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));
# ax.set_ylim([0,110]);ax.set_yticks(np.arange(0,110,20));		
ax = plt.subplot(4,2,5);
ax.annotate('Frequency\n(events/year)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(e)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
scatterplot_ensemble(TS_CHA, HWI_NO_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));	
			
ax = plt.subplot(4,2,7);
ax.annotate('Total days\n(days/year)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(g)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
scatterplot_ensemble(TS_CHA, HWI_M_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));	

ax = plt.subplot(4,2,2);
ax.annotate('Cold Spells',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(b)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
scatterplot_ensemble(TS_CHA, CSInM_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));

ax = plt.subplot(4,2,4);
ax.annotate('(d)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)			
scatterplot_ensemble(TS_CHA, CSDuM_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));
		
ax = plt.subplot(4,2,6);
ax.annotate('(f)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
scatterplot_ensemble(TS_CHA, CSI_NO_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));	

ax = plt.subplot(4,2,8);
ax.annotate('(h)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)			
scatterplot_ensemble(TS_CHA, CSI_M_CHA)
ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));	

ax.set_xlim([0,6]);ax.set_xticks(np.arange(0,6.1,1));
"""

