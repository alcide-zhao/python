# -*- coding: utf-8 -*-
"""
This is to show aerosol induced cpercentage changes in different categories of rainfall
	light raindall: <1mm/day
	moderate rainfall: 1 - 95 th percentiles
	heavy : > 95th percentile
"""
import scipy
import numpy as np
from scipy import stats
import os; import site; 
import netCDF4 as nc4
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
import math

region_dic={'ASIA':[60,155,-5,55],'EA':[100,145,20,50],'SA':[65,100,5,30]}
# ocean land masr
import scipy.io as sio
oceanmask = sio.loadmat('/home/s1667168/coding/python/climate_extremes_cesm/external_data/landoceanmask_CESM.mat')['landoceanmask']
oceanmask[oceanmask==0]=np.nan


#############################
# climatology references
#############################
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/'

def movingaverage (values, window=5):
	boudary = int(math.floor(window/2))
	result = np.empty((len(values))); result[:]=np.nan
	for index in range(0,boudary):
		result[index] = np.mean(values[index:index+window])
	for index in range(-1*boudary-1,0):
		result[index] = np.mean(values[index-window:index])
	for index in range(boudary,len(values)-boudary):
		result[index] = np.mean(values[index-boudary:index+boudary+1])
	return result
	
def RCP_TiSe(region, variable):
	file_name = 'ensumble_mean_PEI_global_2006_2100_rcp85.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	attt_TiSe = stats.nanmean(stats.nanmean(att_clipped_r85,axis=2),axis=1)
	return time,attt_TiSe

def fix_TiSe(region, variable):
	file_name = 'ensumble_mean_PEI_global_2006_2100_fixa.nc'
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	att_value = nc_fid.variables[variable][:]
	att_value = np.multiply(att_value,oceanmask);att_value[att_value<= 0] = np.nan
	lons,lats,att_clipped_r85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,att_value)
	nc_fid.close()
	attt_TiSe = stats.nanmean(stats.nanmean(att_clipped_r85,axis=2),axis=1)
	return time,attt_TiSe
	

fig = plt.figure(figsize=(8, 6), facecolor='White');plot_setup(); pad=5

region_rey = 'SA';region = region_dic[region_rey]

time,r95p_RCP = RCP_TiSe(region, 'r95p')
time,ptot_RCP = RCP_TiSe(region, 'total_precip')
time,precptot_RCP = RCP_TiSe(region, 'precptot')

time,r95p_fix = fix_TiSe(region, 'r95p')
time,ptot_fix = fix_TiSe(region, 'total_precip')
time,precptot_fix = fix_TiSe(region, 'precptot')

h_f = r95p_fix; h_r = r95p_RCP;
m_f = precptot_fix-r95p_fix; m_r = precptot_RCP - r95p_RCP;
l_f = ptot_fix-precptot_fix; l_r = ptot_RCP - precptot_RCP;

ax1 = plt.subplot(3,2,1);ax1.axhline(0, color='k',alpha=1);

y1 =100*(movingaverage (l_f/ptot_fix, window=5)- movingaverage (l_f/ptot_fix, window=5)[0])
y2 =100*(movingaverage (l_r/ptot_RCP, window=5)- movingaverage (l_r/ptot_RCP, window=5)[0])
baseline = y2[0];base_label = '  baseline = '+str(round(100.000*(l_r/ptot_RCP)[0],2))+'%'

#ax1.grid(False,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.axvspan(2031, 2050, alpha=0.2, color='gray');ax1.axvspan(2081, 2100, alpha=0.2, color='gray');
#ax1.plot(time, y1,color='lightgreen', marker = '+')
ax1.plot(time, y2 ,'k',alpha=1,label ='Rcp8.5',marker = 's',markersize = 1.5)
ax1.fill_between(time, baseline, y1, where=y1 <= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, baseline, y1, where=y1 >= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, y1, y2, where= y1 >= y2, facecolor='green',edgecolor='None', interpolate=True,alpha=0.5,label ='Negative AAs offset')
ax1.fill_between(time, y1, y2, where=y1 <= y2, facecolor='red',edgecolor='None', interpolate=True,alpha=0.5,label ='Positive AAs offset')

ax1.annotate('(a)'+base_label,(0.46,0.42), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.annotate('South Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax1.annotate(r'LPP ($\mathrm{\mathsf{\Delta}}$%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)
ax1.set_xlim([2000,2100]); ax1.set_xticks([]);
ax1.set_ylim([-0.68,0.04]); ax1.set_yticks(np.arange(-0.68,0.041,0.18))
print 'light precipitaiton over SA'
print stats.nanmean(y2[25:45]-y1[25:45]),stats.nanmean(y1[25:45])
print stats.nanmean(y2[75:95]-y1[75:95]),stats.nanmean(y1[75:95])

#  legend for fill_between and line
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
black_line = mlines.Line2D([], [], color='k',label='Rcp8.5',marker = 's',markersize = 1.5)
green_cross = mlines.Line2D([], [],marker='+', color='lightgreen')
# lightgreen_patch = mpatches.Patch(color = 'lightgreen',alpha = 0.1, hatch = 'k+',linewidth =0.0)
red_patch = mpatches.Patch(color='red', alpha=0.5,label='AAs negative')
green_patch = mpatches.Patch(color='green', alpha=0.5,label='AAs positive')
lines = [black_line,green_cross,green_patch,  red_patch]
labels = ['RCP8.5','RCP8.5_FixA','Negative AAs influence','Positive AAs influence']
legend2 = plt.legend(lines,labels,ncol=1,loc=3,bbox_to_anchor=(0.43, 0.001))
legend2.get_frame().set_facecolor('white');legend2.get_frame().set_edgecolor('k');legend2.get_frame().set_alpha(1)


ax1 = plt.subplot(3,2,3);ax1.axhline(0, color='k',alpha=1);
y1 =100*(movingaverage (m_f/ptot_fix, window=5)- movingaverage (m_f/ptot_fix, window=5)[0])
y2 =100*(movingaverage (m_r/ptot_RCP, window=5)- movingaverage (m_r/ptot_RCP, window=5)[0])
baseline = y2[0];base_label = '  baseline = '+str(round(100.000*(m_r/ptot_RCP)[0],2)+15)+'%'
# #ax1.grid(False,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.axvspan(2031, 2050, alpha=0.2, color='gray');ax1.axvspan(2081, 2100, alpha=0.2, color='gray');
#ax1.plot(time, y1,color='lightgreen', marker = '+')
ax1.plot(time, y2 ,'k',alpha=1,label ='Rcp8.5',marker = 's',markersize = 1.5)
ax1.fill_between(time, baseline, y1, where=y1 <= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, baseline, y1, where=y1 >= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, y1, y2, where= y1 >= y2, facecolor='green',edgecolor='None', interpolate=True,alpha=0.5,label ='Negative AAs offset')
ax1.fill_between(time, y1, y2, where=y1 <= y2, facecolor='red',edgecolor='None', interpolate=True,alpha=0.5,label ='Positive AAs offset')

ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='on',pad =3,which= 'major')
ax1.annotate('(b)'+base_label,(0.49,0.42), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.annotate(r'MHP ($\mathrm{\mathsf{\Delta}}$%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)
ax1.set_xlim([2000,2100]); ax1.set_xticks([]);
ax1.set_ylim([-4.6,0.4]); ax1.set_yticks(np.arange(-4.6,0.41,1.25))
print 'Moderate precipitaiton over SA'
print stats.nanmean(y2[25:45]-y1[25:45]),stats.nanmean(y1[25:45])
print stats.nanmean(y2[75:95]-y1[75:95]),stats.nanmean(y1[75:95])		
ax1 = plt.subplot(3,2,5);ax1.axhline(0, color='k',alpha=1);
y1 =100*(movingaverage (h_f/ptot_fix, window=5)- movingaverage (h_f/ptot_fix, window=5)[0])
y2 =100*(movingaverage (h_r/ptot_RCP, window=5)- movingaverage (h_r/ptot_RCP, window=5)[0])
baseline = y2[0];base_label = '  baseline = '+str(round(100.000*(h_r/ptot_RCP)[0],2)-15)+'%'

#ax1.grid(False,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.axvspan(2031, 2050, alpha=0.2, color='gray');ax1.axvspan(2081, 2100, alpha=0.2, color='gray');
#ax1.plot(time, y1,color='lightgreen', marker = '+')
ax1.plot(time, y2 ,'k',alpha=1,label ='Rcp8.5',marker = 's',markersize = 1.5)
ax1.fill_between(time, baseline, y1, where=y1 <= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, baseline, y1, where=y1 >= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, y1, y2, where= y1 >= y2, facecolor='green',edgecolor='None', interpolate=True,alpha=0.5,label ='Negative AAs offset')
ax1.fill_between(time, y1, y2, where=y1 <= y2, facecolor='red',edgecolor='None', interpolate=True,alpha=0.5,label ='Positive AAs offset')
print 'extreme precipitaiton over SA'
print stats.nanmean(y2[25:45]-y1[25:45]),stats.nanmean(y1[25:45])
print stats.nanmean(y2[75:95]-y1[75:95]),stats.nanmean(y1[75:95])
ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='on',pad =3,which= 'major')
ax1.annotate('(c)'+base_label,(0.46,0.42), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.annotate(r'EPP ($\mathrm{\mathsf{\Delta}}$%)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation = 90)
ax1.set_xlim([2000,2100]); ax1.set_xticks(np.arange(2000,2100.1,20));
ax1.set_ylim([-0.3,4.5]); ax1.set_yticks(np.arange(-0.3,4.6,1.2))
				
region_rey = 'EA';region = region_dic[region_rey]

time,r95p_RCP = RCP_TiSe(region, 'r95p')
time,ptot_RCP = RCP_TiSe(region, 'total_precip')
time,precptot_RCP = RCP_TiSe(region, 'precptot')

time,r95p_fix = fix_TiSe(region, 'r95p')
time,ptot_fix = fix_TiSe(region, 'total_precip')
time,precptot_fix = fix_TiSe(region, 'precptot')

h_f = r95p_fix; h_r = r95p_RCP;
m_f = precptot_fix-r95p_fix; m_r = precptot_RCP - r95p_RCP;
l_f = ptot_fix-precptot_fix; l_r = ptot_RCP - precptot_RCP;

ax1 = plt.subplot(3,2,2);ax1.axhline(0, color='k',alpha=1);
y1 =100*(movingaverage (l_f/ptot_fix, window=5)- movingaverage (l_f/ptot_fix, window=5)[0])
y2 =100*(movingaverage (l_r/ptot_RCP, window=5)- movingaverage (l_r/ptot_RCP, window=5)[0])
baseline = y2[0];base_label = '  baseline = '+str(round(100.000*(l_r/ptot_RCP)[0],2))+'%'
#ax1.grid(False,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.axvspan(2031, 2050, alpha=0.2, color='gray');ax1.axvspan(2081, 2100, alpha=0.2, color='gray');
#ax1.plot(time, y1,color='lightgreen', marker = '+')
ax1.plot(time, y2 ,'k',alpha=1,label ='Rcp8.5',marker = 's',markersize = 1.5)
ax1.fill_between(time, baseline, y1, where=y1 <= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, baseline, y1, where=y1 >= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, y1, y2, where= y1 >= y2, facecolor='green',edgecolor='None', interpolate=True,alpha=0.5,label ='Negative AAs offset')
ax1.fill_between(time, y1, y2, where=y1 <= y2, facecolor='red',edgecolor='None', interpolate=True,alpha=0.5,label ='Positive AAs offset')

ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='on',pad =3,which= 'major')
ax1.annotate('(d)'+base_label,(0.46,0.42), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.annotate('East Asia',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline')
ax1.set_xlim([2000,2100]); ax1.set_xticks([]);
ax1.set_ylim([-0.68,0.04]); ax1.set_yticks(np.arange(-0.68,0.041,0.18))
print 'light precipitaiton over EA'
print stats.nanmean(y2[25:45]-y1[25:45]),stats.nanmean(y1[25:45])
print stats.nanmean(y2[75:95]-y1[75:95]),stats.nanmean(y1[75:95])

ax1 = plt.subplot(3,2,4);ax1.axhline(0, color='k',alpha=1);
y1 =100*(movingaverage (m_f/ptot_fix, window=5)- movingaverage (m_f/ptot_fix, window=5)[0])
y2 =100*(movingaverage (m_r/ptot_RCP, window=5)- movingaverage (m_r/ptot_RCP, window=5)[0])
baseline = y2[0];base_label = '  baseline = '+str(24+round(100.000*(m_r/ptot_RCP)[0],2))+'%'
#ax1.grid(False,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.axvspan(2031, 2050, alpha=0.2, color='gray');ax1.axvspan(2081, 2100, alpha=0.2, color='gray');
#ax1.plot(time, y1,color='lightgreen', marker = '+')
ax1.plot(time, y2 ,'k',alpha=1,label ='Rcp8.5',marker = 's',markersize = 1.5)
ax1.fill_between(time, baseline, y1, where=y1 <= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, baseline, y1, where=y1 >= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, y1, y2, where= y1 >= y2, facecolor='green',edgecolor='None', interpolate=True,alpha=0.5,label ='Negative AAs offset')
ax1.fill_between(time, y1, y2, where=y1 <= y2, facecolor='red',edgecolor='None', interpolate=True,alpha=0.5,label ='Positive AAs offset')

ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='on',pad =3,which= 'major')
ax1.annotate('(e)'+base_label,(0.46,0.42), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.set_xlim([2000,2100]); ax1.set_xticks([]);
ax1.set_ylim([-4.6,0.4]); ax1.set_yticks(np.arange(-4.6,0.41,1.25))
print 'moderate precipitaiton over EA'
print stats.nanmean(y2[25:45]-y1[25:45]),stats.nanmean(y1[25:45])
print stats.nanmean(y2[75:95]-y1[75:95]),stats.nanmean(y1[75:95])			
ax1 = plt.subplot(3,2,6);ax1.axhline(0, color='k',alpha=1);
y1 =100*(movingaverage (h_f/ptot_fix, window=5)- movingaverage (h_f/ptot_fix, window=5)[0])
y2 =100*(movingaverage (h_r/ptot_RCP, window=5)- movingaverage (h_r/ptot_RCP, window=5)[0])
baseline = y2[0];base_label = '  baseline = '+str(round(100.000*(h_r/ptot_RCP)[0],2)-24)+'%'
#ax1.grid(False,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.axvspan(2031, 2050, alpha=0.2, color='gray');ax1.axvspan(2081, 2100, alpha=0.2, color='gray');
#ax1.plot(time, y1,color='lightgreen', marker = '+')
ax1.plot(time, y2 ,'k',alpha=1,label ='Rcp8.5',marker = 's',markersize = 1.5)
ax1.fill_between(time, baseline, y1, where=y1 <= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, baseline, y1, where=y1 >= baseline, hatch = '++',edgecolor='green',color ='None',linewidth =0.0)
ax1.fill_between(time, y1, y2, where= y1 >= y2, facecolor='green',edgecolor='None', interpolate=True,alpha=0.5,label ='Negative AAs offset')
ax1.fill_between(time, y1, y2, where=y1 <= y2, facecolor='red',edgecolor='None', interpolate=True,alpha=0.5,label ='Positive AAs offset')
print 'extreme precipitaiton over EA'
print stats.nanmean(y2[25:45]-y1[25:45]),stats.nanmean(y1[25:45])
print stats.nanmean(y2[75:95]-y1[75:95]),stats.nanmean(y1[75:95])
ax1.tick_params(axis='both',direction ='in',bottom = 'on', top='off', left='on', right='on',pad =3,which= 'major')
ax1.annotate('(f)'+base_label,(0.46,0.42), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='right', va='baseline',rotation='horizontal',fontsize=10)
ax1.set_xlim([2000,2100]); ax1.set_xticks(np.arange(2000,2100.1,20));
ax1.set_ylim([-0.3,4.5]); ax1.set_yticks(np.arange(-0.3,4.6,1.2))
plt.subplots_adjust(left=0.10, bottom=0.05, right=0.98, top=0.95, wspace=0.15, hspace=0.10);
plt.savefig('Fig15.pdf', format='pdf', dpi=1000)
plt.show()
