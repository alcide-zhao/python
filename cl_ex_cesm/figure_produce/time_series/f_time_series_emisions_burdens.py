# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to show the time evelution of selected ectreme indices 
				the evelution of aerosol emissions, burden and AOD 
@author: Alcide.Zhao
"""
import scipy

import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import math
import site
import scipy.io as sio
# from scipy.interpolate import interp2d

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
# 0.0 data input
#######################################

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}

#######################################################
# 0.2 functions and subroutines                       #
#######################################################
def movingaverage (values, window=3):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result
def getArea(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	return area
	
def time_series_region_total(variable,region,lon,lat):
	lons,lats,value = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable)
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid (lons,lats)
	area = getArea(lons,lons+lon_res,lats,lats+lat_res)
	TSRT = np.nansum(np.nansum(np.multiply(value,area),axis=2),axis=1) # time series of region total
	return TSRT

################################################
# Loading the JJA data from the .mat files 
# process the data into region total 
################################################

# emissions
file_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
emission_file = file_path+'RCP_JJA_emissions_2005_2100.mat'
data = sio.loadmat(emission_file)
time_E = data['time'][0,:];lon_E = data['lon'][0,:];lat_E = data['lat'][0,:];
BC_E = data['BC'][:];OC_E = data['OC'][:];SO2_E = data['SO2'][:]; 


# burdens
burden_file = file_path+'Aerosol_Burden_AOD_JJA_2006_2100.mat'
data = sio.loadmat(burden_file)
time_B = data['time'][0,:];lon_B = data['lon'][0,:];lat_B = data['lat'][0,:];
BC_RCP85 = data['BC_RCP85'];BC_FIXA = data['BC_FIXA'];
OC_RCP85 = data['OC_RCP85'];OC_FIXA = data['OC_FIXA'];
# AOD_RCP85 = data['AOD_RCP85'];AOD_FIXA = data['AOD_FIXA'];
SO4_RCP85 = data['SO4_RCP85'];SO4_FIXA = data['SO4_FIXA'];


#################################
# plotting the figures
################################
# fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(6, 6), facecolor='White')
fig1 = plt.figure(facecolor='White',figsize=[6,6]);plot_setup();

# subplot_1 the emissions over Asia domain
region = rergion_dic['ASIA'][:]
BC_E_TSRT = time_series_region_total(BC_E,region,lon_E,lat_E)
OC_E_TSRT = time_series_region_total(OC_E,region,lon_E,lat_E)
SO2_E_TSRT = time_series_region_total(SO2_E,region,lon_E,lat_E)

BC_B_TSRT = time_series_region_total(BC_RCP85,region,lon_B,lat_B)-time_series_region_total(BC_FIXA,region,lon_B,lat_B)
OC_B_TSRT = time_series_region_total(OC_RCP85,region,lon_B,lat_B)-time_series_region_total(OC_FIXA,region,lon_B,lat_B)
SO4_B_TSRT = time_series_region_total(SO4_RCP85,region,lon_B,lat_B)-time_series_region_total(SO4_FIXA,region,lon_B,lat_B)
"""
region = rergion_dic['ASIA'][:]
BC_B_TSRT = time_series_region_total(SO4_RCP85,region,lon_B,lat_B)-time_series_region_total(SO4_FIXA,region,lon_B,lat_B)
print BC_B_TSRT[94]*10**(-8)
BC_B_TSRT =np.divide(BC_B_TSRT,time_series_region_total(SO4_FIXA,region,lon_B,lat_B))*100
print BC_B_TSRT[94]
region = rergion_dic['EA'][:]
BC_B_TSRT = time_series_region_total(SO4_RCP85,region,lon_B,lat_B)-time_series_region_total(SO4_FIXA,region,lon_B,lat_B)
print BC_B_TSRT[94]*10**(-8)
BC_B_TSRT =np.divide(BC_B_TSRT,time_series_region_total(SO4_FIXA,region,lon_B,lat_B))*100
print BC_B_TSRT[94]
region = rergion_dic['SA'][:]
BC_B_TSRT = time_series_region_total(SO4_RCP85,region,lon_B,lat_B)-time_series_region_total(SO4_FIXA,region,lon_B,lat_B)
print BC_B_TSRT[44]*10**(-8)
BC_B_TSRT =np.divide(BC_B_TSRT,time_series_region_total(SO4_FIXA,region,lon_B,lat_B))*100
print BC_B_TSRT[44]


"""
## left y-asix the emissions 
ax1 = plt.subplot(3,1,1);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_xlim([2000,2100]);ax1.set_xticklabels([])
ax1.set_ylim([-1.5,0.3]);ax1.set_yticks(np.arange(-1.5,0.31,0.45))
ax1.axvspan(2031, 2050, alpha=0.3, color='pink');ax1.axvspan(2081, 2100, alpha=0.3, color='pink');
ax1.annotate('(a) AMR land',xy=(-0.17, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
## right y-axis the burrdens
ax2 = ax1.twinx();ax2.set_ylim([-3.0,0.6]);ax2.set_yticks(np.arange(-3.0,0.61,0.9))
lins1 = ax1.plot(time_B,10**(-2)*movingaverage(BC_E_TSRT[1:]-BC_E_TSRT[0]),'-.',color="b",label='BC*10',linewidth=1.5)
lins2 = ax1.plot(time_B,10**(-3)*movingaverage(OC_E_TSRT[1:]-OC_E_TSRT[0]),'-.',color="g",label='OC',linewidth=1.5)
lins3 = ax1.plot(time_B,10**(-3)*movingaverage(SO2_E_TSRT[1:]-SO2_E_TSRT[0]),'-.',color="r",label='SO2',linewidth=1.5)
lins4 = ax2.plot(time_B,10**(-7)*movingaverage(BC_B_TSRT),'-',color="b",label='BC*10',linewidth=1.5)
lins5 = ax2.plot(time_B,10**(-8)*movingaverage(OC_B_TSRT),'-',color="g",label='POA+SOA',linewidth=1.5)
lins6 = ax2.plot(time_B,10**(-8)*movingaverage(SO4_B_TSRT),'-',color="r",label='SO4',linewidth=1.5)
lins7 = ax2.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
align_yaxis(ax1,0,ax2,0)
# subplot_2 the emissions over East Aian Summer Monsoon Region
region = rergion_dic['EA'][:]
BC_E_TSRT = time_series_region_total(BC_E,region,lon_E,lat_E)
OC_E_TSRT = time_series_region_total(OC_E,region,lon_E,lat_E)
SO2_E_TSRT = time_series_region_total(SO2_E,region,lon_E,lat_E)

BC_B_TSRT = time_series_region_total(BC_RCP85,region,lon_B,lat_B)-time_series_region_total(BC_FIXA,region,lon_B,lat_B)
OC_B_TSRT = time_series_region_total(OC_RCP85,region,lon_B,lat_B)-time_series_region_total(OC_FIXA,region,lon_B,lat_B)
SO4_B_TSRT = time_series_region_total(SO4_RCP85,region,lon_B,lat_B)-time_series_region_total(SO4_FIXA,region,lon_B,lat_B)

## left y-asix the emissions 
ax1 = plt.subplot(3,1,2);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([-1.3,0.3]);ax1.set_yticks(np.arange(-1.3,0.31,0.4))
ax1.set_xlim([2000,2100]);ax1.set_xticklabels([])
ax1.axvspan(2031, 2050, alpha=0.3, color='pink');ax1.axvspan(2081, 2100, alpha=0.3, color='pink');
ax1.annotate('(b) EA land',xy=(-0.17, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{Mg\/s^{-1}}}$'+')', color='k')

## the right y-axis the burrdens
ax2 = ax1.twinx();ax2.set_ylim([-1.60,0.37]);ax2.set_yticks(np.arange(-1.60,0.38,0.49))
ax2.set_ylabel('Aerosol burdens (Tg)', color='k')

lins1 = ax1.plot(time_B,10**(-2)*movingaverage(BC_E_TSRT[1:]-BC_E_TSRT[0]),'-.',color="b",label='BC*10',linewidth=1.5)
lins2 = ax1.plot(time_B,10**(-3)*movingaverage(OC_E_TSRT[1:]-OC_E_TSRT[0]),'-.',color="g",label='OC',linewidth=1.5)
lins3 = ax1.plot(time_B,10**(-3)*movingaverage(SO2_E_TSRT[1:]-SO2_E_TSRT[0]),'-.',color="r",label='SO2',linewidth=1.5)
lins4 = ax2.plot(time_B,10**(-7)*movingaverage(BC_B_TSRT),'-',color="b",label='BC*10',linewidth=1.5)
lins5 = ax2.plot(time_B,10**(-8)*movingaverage(OC_B_TSRT),'-',color="g",label='POA+SOA',linewidth=1.5)
lins6 = ax2.plot(time_B,10**(-8)*movingaverage(SO4_B_TSRT),'-',color="r",label='SO4',linewidth=1.5)
lins7 = ax2.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
align_yaxis(ax1,0,ax2,0)


# subplot_3 the emissions over South Aian Summer Monsoon Region
region = rergion_dic['SA'][:]
BC_E_TSRT = time_series_region_total(BC_E,region,lon_E,lat_E)
OC_E_TSRT = time_series_region_total(OC_E,region,lon_E,lat_E)
SO2_E_TSRT = time_series_region_total(SO2_E,region,lon_E,lat_E)
BC_B_TSRT = time_series_region_total(BC_RCP85,region,lon_B,lat_B)-time_series_region_total(BC_FIXA,region,lon_B,lat_B)
OC_B_TSRT = time_series_region_total(OC_RCP85,region,lon_B,lat_B)-time_series_region_total(OC_FIXA,region,lon_B,lat_B)
SO4_B_TSRT = time_series_region_total(SO4_RCP85,region,lon_B,lat_B)-time_series_region_total(SO4_FIXA,region,lon_B,lat_B)
## left y-asix the emissions 
ax1 = plt.subplot(3,1,3);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_xlim([2000,2100]);ax1.set_xlabel('Year', color='k');
ax1.axvspan(2031, 2050, alpha=0.3, color='pink');ax1.axvspan(2081, 2100, alpha=0.3, color='pink');
ax1.set_ylim([-0.2,0.2]);ax1.set_yticks(np.arange(-0.2,0.21,0.1))
ax1.annotate('(c) SA land',xy=(-0.17, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
## the right y-axis the burrdens
ax2 = ax1.twinx();ax2.set_ylim([-0.5,0.5]);ax2.set_yticks(np.arange(-0.5,0.51,0.25))
lins1 = ax1.plot(time_B,10**(-2)*movingaverage(BC_E_TSRT[1:]-BC_E_TSRT[0]),'-.',color="b",label='BC'+r'$\times$'+'10',linewidth=1.5)
lins2 = ax1.plot(time_B,10**(-3)*movingaverage(OC_E_TSRT[1:]-OC_E_TSRT[0]),'-.',color="g",label='           OC',linewidth=1.5)
lins3 = ax1.plot(time_B,10**(-3)*movingaverage(SO2_E_TSRT[1:]-SO2_E_TSRT[0]),'-.',color="r",label='SO2',linewidth=1.5)
lins4 = ax2.plot(time_B,10**(-7)*movingaverage(BC_B_TSRT),'-',color="b",label='BC'+r'$\times$'+'10',linewidth=1.5)
lins5 = ax2.plot(time_B,10**(-8)*movingaverage(OC_B_TSRT),'-',color="g",label='POA+SOA',linewidth=1.5)
lins6 = ax2.plot(time_B,10**(-8)*movingaverage(SO4_B_TSRT),'-',color="r",label='SO4',linewidth=1.5)
lins7 = ax2.plot(range(2000,2101),0*np.ones((101)),'-',color="k")
align_yaxis(ax1, 0, ax2, 0)

##legends
legend1 = ax1.legend(shadow=False,ncol=3,bbox_to_anchor=(0.75, 0.288))	 
legend1.get_frame().set_facecolor('gray');legend1.get_frame().set_edgecolor('None');legend1.get_frame().set_alpha(0.1)
legend2 = ax2.legend(shadow=False,ncol=3,bbox_to_anchor=(0.75, 0.17))	 
legend2.get_frame().set_facecolor('gray');legend2.get_frame().set_edgecolor('None');legend2.get_frame().set_alpha(0.1)


plt.subplots_adjust(left=0.15, bottom=0.08, right=0.85, top=0.96, wspace=0, hspace=0.1);
plt.savefig('fig4_time_eccelution_emissions_burdens.png', format='png', dpi=800)
plt.show()
"""
"""
# emissions_sum
sub_plot=5
spst = [2,4]
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_sum[1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_sum[0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('Emissions',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
# aerosol_burden sum
sub_plot=6

ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(burden_sum_fixA)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)

ax1.plot(year_series,movingaverage(burden_sum_rcp85)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

#aerosol optical depth
sub_plot=7
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(AODVIS_fixA_jja_time_series),'-',color="blue",linewidth=2,label='RCP8.5_fixA')
ax1.plot(year_series,movingaverage(AODVIS_rcp85_jja_time_series),'-',color="red",linewidth=2,label='RCP8.5')
# ax1.hold(True)
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
# ax1.set_ylabel('AOD550'.title(),fontsize=30)
# ax1.set(aspect=10)
ax1.set_title('AOD550',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(1.5, 0.8),fontsize=15)	

spst =[2,4]
plt.figure(2, facecolor='White',figsize=[10,10])

# BC emissions
sub_plot=1
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_TiSe_dic['BC'][1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_TiSe_dic['BC'][0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('BC',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# OC emissions
sub_plot=2
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_TiSe_dic['OC'][1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_TiSe_dic['OC'][0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('OC',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# BC emissions
sub_plot=3
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(emission_TiSe_dic['SO2'][1:]),'-',color="red",linewidth=2,label='RCP8.5') #*1000 considering the value scale
ax1.plot(year_series,movingaverage(emission_TiSe_dic['SO2'][0]*np.ones((95))),'-',color="blue",linewidth=2,label='RCP8.5_FixA')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2/s'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('SO2',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(1.5, 0.8),fontsize=15)	

# BC BURDEN
sub_plot=5
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENBC_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENBC_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('BC Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# POM+SOA BURDEN
sub_plot = 6
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENPOM_fixA_jja_time_series+BURDENSOA_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENPOM_rcp85_jja_time_series+BURDENSOA_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('POA+SOA Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)
legend = ax1.legend(loc=10, shadow=True,bbox_to_anchor=(1.5, 0.9),fontsize=15)	


# SO4 BURDEN
sub_plot = 7
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENSO4_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
ax1.hold(True)
ax1=plt.subplot(spst[0],spst[1],sub_plot)
ax1.plot(year_series,movingaverage(BURDENSO4_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
ax1.set_xlabel('Year',fontsize=10)
ax1.set_xlim([2000,2100])
# ax1.set_ylim([4,7.5])
ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# ax1.set(aspect=10)
ax1.set_title('SO4 Aerosol Burden',fontsize=15)
ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# SEASALT  BURDEN
# sub_plot = 8
# ax1=plt.subplot(spst[0],spst[1],sub_plot)
# ax1.plot(year_series,movingaverage(BURDENSEASALT_fixA_jja_time_series+BURDENDUST_fixA_jja_time_series)*10000,'-',color="blue",linewidth=2,label='RCP8.5_fixA') #*1000 considering the value scale
# ax1.hold(True)
# ax1=plt.subplot(spst[0],spst[1],sub_plot)
# ax1.plot(year_series,movingaverage(BURDENSEASALT_rcp85_jja_time_series+BURDENDUST_rcp85_jja_time_series)*10000,'-',color="red",linewidth=2,label='RCP8.5')
# ax1.set_xlabel('Year',fontsize=10)
# ax1.set_xlim([2000,2100])
# # ax1.set_ylim([4,7.5])
# ax1.set_ylabel('Kg/m^2*10^(-4)'.title(),fontsize=15)
# # ax1.set(aspect=10)
# ax1.set_title('SEASALT Aerosol Burden',fontsize=15)
# ax1.grid(True,color='k', linestyle='-', linewidth=1,alpha=0.1) 
# ax1.tick_params(labelsize=15,axis='both', which='major', pad=5)

# plt.savefig('/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/temp/time_series_global.tif')
