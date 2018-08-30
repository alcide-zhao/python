# -*- coding: utf-8 -*-
'''
This is to produce a diagram to show the prod=cedures to calculate the HW indices
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import time as clock
import math as math

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio

# functions
def zeros_lookups(data):
    # Create an array that is 1 where data is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(data, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where abs diff is 1.
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges

def get_all_baseline_data(scenario,directory,lat_point,lon_point):
	"""
	reading all the necessary date including the thresholds, 
	the 1961-1990 baselines and the inter annual/seasonal osculations upon requet
	"""
	data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp'
	# TmTn 95percentiles for 1961-1995
	threshold_data =data_path+'/Temp_pp/TmTn_percentile_calender_enMean.nc'
	nc_fid = nc4.Dataset(threshold_data,mode='r')
	lat = nc_fid.variables['lat'][lat_point]
	lon = nc_fid.variables['lon'][lon_point]
	TX95P = nc_fid.variables['TX95P'][:,lat_point,lon_point]+273.15
	TN95P = nc_fid.variables['TN95P'][:,lat_point,lon_point]+273.15
	nc_fid.close()
	
	# baseline data
	baseline_data =data_path+'/Temp_pp/'+directory+'/HWI_CSI_1961_1990_0595_'+directory+'.nc' #_InterOsi
	nc_fid = nc4.Dataset(baseline_data,mode='r')
	HWIM= stats.nanmean(nc_fid.variables['HWIM'],axis=0)[lat_point,lon_point];
	HWIS= stats.nanmean(nc_fid.variables['HWIS'],axis=0)[lat_point,lon_point];
	nc_fid.close()
	return lat,lon,TX95P,TN95P,HWIM,HWIS

################
#     Main     #
################
scenario = 'his';directory='CalDayThr' # _InterOsi
lat_point =112;lon_point = 118   #50N 10E
HW_index = np.ones(365);
scenario_dic ={'his':[1,84,85]}   ##2003
ensemble = scenario_dic[scenario][0]; layer_s = scenario_dic[scenario][1]; layer_e = scenario_dic[scenario][2]
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/'+scenario+'/'
##get all data input
lat,lon,TX95P,TN95P,HWIM,HWIS = get_all_baseline_data(scenario,directory,lat_point,lon_point) 
print 'lon:',lon,'lat:',lat, 'HWIM:',HWIM,'HWIS:',HWIS
# print TX95P
# print TX95P
# print TN95P

nc_f = input_path +'TmTn_his_01.nc'
print nc_f
nc_fid = nc4.Dataset(nc_f,mode='r')
year = nc_fid.variables['year'][layer_s:layer_e]
TX = nc_fid.variables['TX'][layer_s:layer_e,:,lat_point,lon_point]+273.15;
TN = nc_fid.variables['TN'][layer_s:layer_e,:,lat_point,lon_point]+273.15; 

TX_HW = (TX-TX95P)[0,:]; TN_HW = (TN-TN95P)[0,:]; 
TXTN = ((TX_HW+TN_HW)/2)

##### TXTN
index =  np.ones(365);
tag = [item for item in range(len(TXTN)) if (TX_HW[item]>0 and TN_HW[item]>0)];HW_index[tag] = 0;
ranges_all = zeros_lookups(HW_index);


HWDC = 365-np.sum(HW_index)  #Total days of heatwaves
print ranges_all
print 'HWDC:' ,HWDC

# print np.shape(ranges)
for ino in range(np.shape(ranges_all)[0]): 
	# print ino
	if (ranges_all[ino,1]-ranges_all[ino,0] < 3): 
		# print ranges[ino,1],ranges[ino,0]
		HW_index[ranges_all[ino,0]:ranges_all[ino,1]] = 1
if (len(TXTN)- np.sum(HW_index)==0):  # no cays meet the conditions
	HWI_M=np.nan;HWI_N=np.nan;HWI_X=np.nan; HWI_NO=np.nan;
	HWCM=np.nan;HWCN=np.nan;HWCD=np.nan;HWCS=np.nan; HWCE=np.nan;
	HWCV=np.nan;HWCP=np.nan;HWCU=np.nan;RF_M=np.nan; RF_N=np.nan; RF_X=np.nan;
	Du_X=np.nan; Du_M=np.nan;In_X=np.nan; In_M=np.nan;
else:
	ranges = zeros_lookups(HW_index);
	print ranges
	duraton= ranges[:,1]-ranges[:,0]; Du_X=np.nanmax(duraton); 
	event_no = 0
	HW_events = np.array([]);HW_intensity = np.array([]);
	for ino in range(np.shape(ranges)[0]):
		event_no=event_no+1
		HW_events=np.append(HW_events,np.sum(TXTN[ranges[ino,0]:ranges[ino,1]]))
		HW_intensity=np.append(HW_intensity,stats.nanmean(TXTN[ranges[ino,0]:ranges[ino,1]]))
	In_X=np.nanmax(HW_intensity);
	HWI_NO = event_no
	print 'duration  :',duraton
	print 'DU_X :' ,Du_X
	print 'IN_X :' ,In_X
	print 'intensity :', HW_intensity
	print 'HW_event :', HW_events
	## HW magnitude
	HW_scaled=np.divide(HW_events-HWIM,HWIS)

fig = plt.figure(facecolor='White',figsize=[10,6]);plot_setup();pad= 5;

ax =plt.subplot2grid((4, 1), (0, 0), rowspan=3);
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
# ax.axhline(y=0,color='k',linewidth = 1);ax.axvline(x=12,color='k',ls = '--',linewidth = 1);

ax.set_xlim([1,366]);ax.set_xticks(np.array([31,59,90,120,151,181,212,243,273,304,334,365]));
ax.set_xticklabels(('31 Jan','28 Feb','31 Mar','30 April','31 May','30 June','31 July','31 Aug',' 30 Sep','31 Oct','30 Nov','31 Dec'));
ax.set_ylabel('Temperature (K)')
x=np.arange(1,366)
##line
ax.plot(x,TX[0,:],'-',color="m",linewidth=2,label = 'TX')
ax.plot(x,TX95P,'-.',color="m",linewidth=2,label = 'TX95P (1961-1990)')
ax.plot(x,TN[0,:],'-',color="g",linewidth=2,label = 'TN')
ax.plot(x,TN95P,'-.',color="g",linewidth=2,label = 'TN95P (1961-1990)')
##### TX excess#
index =  np.ones(365);
tag = [item for item in range(len(TXTN)) if (TX_HW[item]>0 and TN_HW[item]<=0)];index[tag] = 0;
TX_ranges = zeros_lookups(index);

for ino in range(np.shape(TX_ranges)[0]):
	ax.axvspan(TX_ranges[ino,0],TX_ranges[ino,1], alpha=0.5, color='m',edgecolor='none',lw=0);

##### TN excess#
index =  np.ones(365);
tag=[item for item in range(len(TXTN)) if (TX_HW[item]<=0 and TN_HW[item]>0)];index[tag] = 0;
TN_ranges = zeros_lookups(index);
for ino in range(np.shape(TN_ranges)[0]):
	ax.axvspan(TN_ranges[ino,0],TN_ranges[ino,1], alpha=0.5, color='g',edgecolor='none',lw=0);
	
#### Tx>95p and TN>95p but foe less than 3 days	
tag = [item for item in range(len(TXTN)) if (TX_HW[item]>0 and TN_HW[item]>0)];HW_index[tag] = 0;
for ino in range(np.shape(ranges_all)[0]): 
	if (ranges_all[ino,1]-ranges_all[ino,0] >= 3): 
		HW_index[ranges_all[ino,0]:ranges_all[ino,1]] = 1	
ranges_all = zeros_lookups(HW_index);
for ino in range(np.shape(ranges_all)[0]):
	ax.axvspan(ranges_all[ino,0],ranges_all[ino,1], alpha=0.5, color='y',edgecolor='none',lw=0);
###### heatwaves detected
for ino in range(np.shape(ranges)[0]):
	ax.axvspan(ranges[ino,0],ranges[ino,1], alpha=0.7, color='r',edgecolor='none',lw=0);
for ino in range(np.shape(ranges)[0]):	
	x = range(ranges[ino,0]+1,ranges[ino,1]+1)
	y_lower = TX95P[ranges[ino,0]:ranges[ino,1]]
	y_upper = TX[0,ranges[ino,0]:ranges[ino,1]]
	ax.fill_between(x,y_lower,y_upper,facecolor='k',alpha=1)	
	y_lower = TN95P[ranges[ino,0]:ranges[ino,1]]
	y_upper = TN[0,ranges[ino,0]:ranges[ino,1]]
	ax.fill_between(x,y_lower,y_upper,facecolor="k",alpha=1)	

ax.annotate('TX excess', xy=(ranges[2,0]+4, TX95P[ranges[2,0]+1]+0.05), xytext=(ranges[2,0]-50, 29.5+273.15),
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"),fontsize=10
            )
ax.annotate('TN excess', xy=(ranges[2,0]+4, TN95P[ranges[2,0]+1]+0.1), xytext=(ranges[2,0]-50, 28.5+273.15),
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"),fontsize=10
            )
			
ax.annotate(r'$ Intensity =\frac{1}{4}\sum_ {i=1}^{4} \frac{[TX(i)-TX95P(i)]+[TN(i)-TN95P(i)]}{2} = 0.23 K$', xy=(ranges[2,0]+4, 26+273.15), xytext=(ranges[2,0]-100, 25.1+273.15),
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"),fontsize=13,bbox=dict(boxstyle="round4,pad=.5", fc="0.9"),
            )
ax.annotate(r'$ Magnitude =\sum_ {i=1}^{4} \frac{[TX(i)-TX95P(i)]+[TN(i)-TN95P(i)]}{2} = 0.93 K days$', xy=(ranges[2,0]+2, 26+273.15), xytext=(ranges[2,0]-30, 24.3+273.15),
             arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"),fontsize=13,bbox=dict(boxstyle="round4,pad=.5", fc="0.9"),
            )					
# r'$\sum_{i=0}^\infty x_i$'			
ann = ax.annotate('Duration = 4 days', (ranges[2,0]-2, 29.5),(ranges[2,1]+3, 29.5+273.15),
                  arrowprops={'arrowstyle':'<->'},fontsize=10)
				  

ax2 =plt.subplot2grid((4, 1), (3, 0), rowspan=1);ax2.axis('off')

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

x=np.arange(1,366); y =10*x/366
ax2.plot(x,y,'-',color="w",linewidth=0)
ax2.set_xlim([1,366]);ax2.set_ylim([0,10]);

TX = mlines.Line2D(x,y,ls='-',color='m',lw=2)
TX95 = mlines.Line2D(x,y,ls='-.',color='m',lw=2)
TN = mlines.Line2D(x,y,ls='-',color='g',lw=2)
TN95 = mlines.Line2D(x,y,ls='-.',color='g',lw=2)
m_patch = mpatches.Patch(color='m', alpha=0.5)
g_patch = mpatches.Patch(color='g', alpha=0.5)
y_patch = mpatches.Patch(color='y', alpha=0.5)
r_patch = mpatches.Patch(color='r', alpha=0.7)

lines = [TX,TN,TX95,TN95,m_patch,g_patch,y_patch,r_patch]
labels = ['TX','TN','TX95P (1961-1990)','TN95P (1961-1990)','TX > TX95P but TN <= TN95P','TN > TN95P but TX <= TX95P','TX > TX95P and TN > TN95P, but for less than 3 dyas','Heatwave (TX > TX95P and TN > TN95P, for at least 3 days)']
legend = plt.legend(lines,labels,ncol=2,loc='lower left',labelspacing=1.5,markerscale =10)
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)

ax2.annotate('Annual Frequency = 3 events/year', xy=(225, 8),fontsize=8)
ax2.annotate('Annual Peak Intensity = MAX{0.11,0.12,0.23} = 0.23 K', xy=(225, 6),fontsize=8)
ax2.annotate('Annual Peak Duration = MAX{3,6,4} = 6 days/event', xy=(225, 4),fontsize=8)
ax2.annotate('Annual Total Hot Days =\n        SUM{1,3,2,2,1,1,1,1,1,1,1,6,1,4,1,1,2,1} = 31 days/year', xy=(225, 1),fontsize=8)

plt.subplots_adjust(left=0.06, bottom=0.03, right=0.98, top=0.98, wspace=0.15, hspace=0.2); 
plt.savefig('HW_DF_diagram.png', format='png', dpi=1000)