# -*- coding: utf-8 -*-
'''
This is to compare the HWI and CSI from NCEP and CESM ensemble mean ()
'''
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy import stats
import scipy.io as io
io
import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

file_path ='/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/ano/'

NCEP_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_NCEP.nc',mode='r')
NCEP_lon = NCEP_fid.variables['lon'][:]
NCEP_lat = NCEP_fid.variables['lat'][:]
NCEP_year = NCEP_fid.variables['year'][:]
NCEP_HWDC = stats.nanmean(stats.nanmean(NCEP_fid.variables['HWDC'][:],axis=2),axis=1)
NCEP_HWDuM =stats.nanmean(stats.nanmean(NCEP_fid.variables['HWDuM'][:],axis=2),axis=1)
NCEP_HWI_NO = stats.nanmean(stats.nanmean(NCEP_fid.variables['HWI_NO'][:],axis=2),axis=1)
NCEP_HWInM = stats.nanmean(stats.nanmean(NCEP_fid.variables['HWInM'][:],axis=2),axis=1)
NCEP_CSDC = stats.nanmean(stats.nanmean(NCEP_fid.variables['CSDC'][:],axis=2),axis=1)
NCEP_CSDuM =stats.nanmean(stats.nanmean(NCEP_fid.variables['CSDuM'][:],axis=2),axis=1)
NCEP_CSI_NO = stats.nanmean(stats.nanmean(NCEP_fid.variables['CSI_NO'][:],axis=2),axis=1)
NCEP_CSInM = stats.nanmean(stats.nanmean(NCEP_fid.variables['CSInM'][:],axis=2),axis=1)
NCEP_fid.close()


his_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_his.nc',mode='r')
his_lon = his_fid.variables['lon'][:]
his_lat = his_fid.variables['lat'][:]
his_year = his_fid.variables['year'][:]
his_HWDC = stats.nanmean(stats.nanmean(his_fid.variables['HWDC'][:],axis=3),axis=2)
his_HWDuM = stats.nanmean(stats.nanmean(his_fid.variables['HWDuM'][:],axis=3),axis=2)
his_HWI_NO = stats.nanmean(stats.nanmean( his_fid.variables['HWI_NO'][:],axis=3),axis=2)
his_HWInM = stats.nanmean(stats.nanmean(his_fid.variables['HWInM'][:],axis=3),axis=2)
his_CSDC = stats.nanmean(stats.nanmean(his_fid.variables['CSDC'][:],axis=3),axis=2)
his_CSDuM = stats.nanmean(stats.nanmean(his_fid.variables['CSDuM'][:],axis=3),axis=2)
his_CSI_NO = stats.nanmean(stats.nanmean( his_fid.variables['CSI_NO'][:],axis=3),axis=2)
his_CSInM = stats.nanmean(stats.nanmean(his_fid.variables['CSInM'][:],axis=3),axis=2)
his_fid.close()


rcp85_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_rcp85.nc',mode='r')
rcp85_lon = rcp85_fid.variables['lon'][:]
rcp85_lat = rcp85_fid.variables['lat'][:]
rcp85_year = rcp85_fid.variables['year'][:]
rcp85_HWDC = stats.nanmean(stats.nanmean(rcp85_fid.variables['HWDC'][:],axis=3),axis=2)
rcp85_HWDuM = stats.nanmean(stats.nanmean(rcp85_fid.variables['HWDuM'][:],axis=3),axis=2)
rcp85_HWI_NO = stats.nanmean(stats.nanmean( rcp85_fid.variables['HWI_NO'][:],axis=3),axis=2)
rcp85_HWInM = stats.nanmean(stats.nanmean(rcp85_fid.variables['HWInM'][:],axis=3),axis=2)
rcp85_CSDC = stats.nanmean(stats.nanmean(rcp85_fid.variables['CSDC'][:],axis=3),axis=2)
rcp85_CSDuM = stats.nanmean(stats.nanmean(rcp85_fid.variables['CSDuM'][:],axis=3),axis=2)
rcp85_CSI_NO = stats.nanmean(stats.nanmean( rcp85_fid.variables['CSI_NO'][:],axis=3),axis=2)
rcp85_CSInM = stats.nanmean(stats.nanmean(rcp85_fid.variables['CSInM'][:],axis=3),axis=2)
rcp85_fid.close()

fixa_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_fixa.nc',mode='r')
fixa_lon = fixa_fid.variables['lon'][:]
fixa_lat = fixa_fid.variables['lat'][:]
fixa_year = fixa_fid.variables['year'][:]
fixa_HWDC = stats.nanmean(stats.nanmean(fixa_fid.variables['HWDC'][:],axis=3),axis=2)
fixa_HWDuM = stats.nanmean(stats.nanmean(fixa_fid.variables['HWDuM'][:],axis=3),axis=2)
fixa_HWI_NO = stats.nanmean(stats.nanmean( fixa_fid.variables['HWI_NO'][:],axis=3),axis=2)
fixa_HWInM = stats.nanmean(stats.nanmean(fixa_fid.variables['HWInM'][:],axis=3),axis=2)
fixa_CSDC = stats.nanmean(stats.nanmean(fixa_fid.variables['CSDC'][:],axis=3),axis=2)
fixa_CSDuM = stats.nanmean(stats.nanmean(fixa_fid.variables['CSDuM'][:],axis=3),axis=2)
fixa_CSI_NO = stats.nanmean(stats.nanmean( fixa_fid.variables['CSI_NO'][:],axis=3),axis=2)
fixa_CSInM = stats.nanmean(stats.nanmean(fixa_fid.variables['CSInM'][:],axis=3),axis=2)
fixa_fid.close()

# rcp45_fid= nc4.Dataset(file_path+'HWI_CSI_EnPP_rcp45.nc',mode='r')
# rcp45_lon = rcp45_fid.variables['lon'][:]
# rcp45_lat = rcp45_fid.variables['lat'][:]
# rcp45_year = rcp45_fid.variables['year'][:]
# rcp45_HWDC = stats.nanmean(stats.nanmean(rcp45_fid.variables['HWDC'][:],axis=3),axis=2)
# rcp45_HWDuM = stats.nanmean(stats.nanmean(rcp45_fid.variables['HWDuM'][:],axis=3),axis=2)
# rcp45_HWI_NO = stats.nanmean(stats.nanmean( rcp45_fid.variables['HWI_NO'][:],axis=3),axis=2)
# rcp45_HWInM = stats.nanmean(stats.nanmean(rcp45_fid.variables['HWInM'][:],axis=3),axis=2)
# rcp45_CSDC = stats.nanmean(stats.nanmean(rcp45_fid.variables['CSDC'][:],axis=3),axis=2)
# rcp45_CSDuM = stats.nanmean(stats.nanmean(rcp45_fid.variables['CSDuM'][:],axis=3),axis=2)
# rcp45_CSI_NO = stats.nanmean(stats.nanmean( rcp45_fid.variables['CSI_NO'][:],axis=3),axis=2)
# rcp45_CSInM = stats.nanmean(stats.nanmean(rcp45_fid.variables['CSInM'][:],axis=3),axis=2)
# rcp45_fid.close()

fig = plt.figure(facecolor='White',figsize=[8.5,10]);plot_setup();pad= 5;
ax = plt.subplot(4,2,1);ax.set_xlim([1920,2100]);ax.set_xticklabels([]);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.annotate('Intensity\n($^\circ$C)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('HeatWaves',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(a)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
lins0 = ax.plot(NCEP_year,NCEP_HWInM[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_HWInM[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_HWInM[2,:],his_HWInM[3,:], where=his_HWInM[2,:]>=his_HWInM[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_HWInM[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_HWInM[2,:],rcp85_HWInM[3,:], where=rcp85_HWInM[2,:]>=rcp85_HWInM[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_HWInM[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_HWInM[2,:],rcp45_HWInM[3,:], where=rcp45_HWInM[2,:]>=rcp45_HWInM[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_HWInM[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_HWInM[2,:],fixa_HWInM[3,:], where=fixa_HWInM[2,:]>=fixa_HWInM[3,:], color="b", alpha = 0.1) 

ax = plt.subplot(4,2,3);ax.set_xlim([1920,2100]);ax.set_xticklabels([]);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.annotate('Duration\n(days/event)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(c)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
lins0 = ax.plot(NCEP_year,NCEP_HWDuM[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_HWDuM[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_HWDuM[2,:],his_HWDuM[3,:], where=his_HWDuM[2,:]>=his_HWDuM[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_HWDuM[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_HWDuM[2,:],rcp85_HWDuM[3,:], where=rcp85_HWDuM[2,:]>=rcp85_HWDuM[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_HWDuM[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_HWDuM[2,:],rcp45_HWDuM[3,:], where=rcp45_HWDuM[2,:]>=rcp45_HWDuM[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_HWDuM[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_HWDuM[2,:],fixa_HWDuM[3,:], where=fixa_HWDuM[2,:]>=fixa_HWDuM[3,:], color="b", alpha = 0.1) 
		
ax = plt.subplot(4,2,5);ax.set_xlim([1920,2100]);ax.set_xticklabels([]);
ax.annotate('Frequency\n(events/year)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(e)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
lins0 = ax.plot(NCEP_year,NCEP_HWI_NO[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_HWI_NO[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_HWI_NO[2,:],his_HWI_NO[3,:], where=his_HWI_NO[2,:]>=his_HWI_NO[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_HWI_NO[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_HWI_NO[2,:],rcp85_HWI_NO[3,:], where=rcp85_HWI_NO[2,:]>=rcp85_HWI_NO[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_HWI_NO[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_HWI_NO[2,:],rcp45_HWI_NO[3,:], where=rcp45_HWI_NO[2,:]>=rcp45_HWI_NO[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_HWI_NO[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_HWI_NO[2,:],fixa_HWI_NO[3,:], where=fixa_HWI_NO[2,:]>=fixa_HWI_NO[3,:], color="b", alpha = 0.1) 

ax = plt.subplot(4,2,7);ax.set_xlim([1920,2100]);
ax.annotate('Total days\n(days/year)',xy=(-0.15, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='vertical')
ax.annotate('(g)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
lins0 = ax.plot(NCEP_year,NCEP_HWDC[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_HWDC[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_HWDC[4,:],his_HWDC[5,:], where=his_HWDC[4,:]>=his_HWDC[5,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_HWDC[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_HWDC[2,:],rcp85_HWDC[3,:], where=rcp85_HWDC[4,:]>=rcp85_HWDC[5,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_HWDC[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_HWDC[4,:],rcp45_HWDC[5,:], where=rcp45_HWDC[4,:]>=rcp45_HWDC[5,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_HWDC[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_HWDC[4,:],fixa_HWDC[5,:], where=fixa_HWDC[4,:]>=fixa_HWDC[5,:], color="b", alpha = 0.1) 
	
ax = plt.subplot(4,2,2);ax.set_xlim([1920,2100]);ax.set_xticklabels([]);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.annotate('Cold Spells',xy=(0.5, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='center',rotation='horizontal')
ax.annotate('(b)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)				
lins0 = ax.plot(NCEP_year,NCEP_CSInM[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_CSInM[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_CSInM[2,:],his_CSInM[3,:], where=his_CSInM[2,:]>=his_CSInM[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_CSInM[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_CSInM[2,:],rcp85_CSInM[3,:], where=rcp85_CSInM[2,:]>=rcp85_CSInM[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_CSInM[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_CSInM[2,:],rcp45_CSInM[3,:], where=rcp45_CSInM[2,:]<=rcp45_CSInM[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_CSInM[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_CSInM[2,:],fixa_CSInM[3,:], where=fixa_CSInM[2,:]>=fixa_CSInM[3,:], color="b", alpha = 0.1) 

ax = plt.subplot(4,2,4);ax.set_xlim([1920,2100]);ax.set_xticklabels([]);
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
ax.annotate('(d)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)			
lins0 = ax.plot(NCEP_year,NCEP_CSDuM[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_CSDuM[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_CSDuM[2,:],his_CSDuM[3,:], where=his_CSDuM[2,:]>=his_CSDuM[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_CSDuM[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_CSDuM[2,:],rcp85_CSDuM[3,:], where=rcp85_CSDuM[2,:]>=rcp85_CSDuM[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_CSDuM[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_CSDuM[2,:],rcp45_CSDuM[3,:], where=rcp45_CSDuM[2,:]>=rcp45_CSDuM[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_CSDuM[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_CSDuM[2,:],fixa_CSDuM[3,:], where=fixa_CSDuM[2,:]>=fixa_CSDuM[3,:], color="b", alpha = 0.1) 
		
ax = plt.subplot(4,2,6);ax.set_xlim([1920,2100]);ax.set_xticklabels([]);
ax.annotate('(f)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)						
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
lins0 = ax.plot(NCEP_year,NCEP_CSI_NO[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_CSI_NO[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_CSI_NO[2,:],his_CSI_NO[3,:], where=his_CSI_NO[2,:]>=his_CSI_NO[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_CSI_NO[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_CSI_NO[2,:],rcp85_CSI_NO[3,:], where=rcp85_CSI_NO[2,:]>=rcp85_CSI_NO[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_CSI_NO[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_CSI_NO[2,:],rcp45_CSI_NO[3,:], where=rcp45_CSI_NO[2,:]>=rcp45_CSI_NO[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_CSI_NO[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_CSI_NO[2,:],fixa_CSI_NO[3,:], where=fixa_CSI_NO[2,:]>=fixa_CSI_NO[3,:], color="b", alpha = 0.1) 

ax = plt.subplot(4,2,8);ax.set_xlim([1920,2100]);
ax.annotate('(h)',xy=(0.02,0.89), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)		
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5)
lins0 = ax.plot(NCEP_year,NCEP_CSDC[:],color='k',linewidth=1)
lins1 = ax.plot(his_year,his_CSDC[1,:],'-',color='r',linewidth=1,label='his')
lins2 = ax.fill_between(his_year,his_CSDC[2,:],his_CSDC[3,:], where=his_CSDC[2,:]>=his_CSDC[3,:], color='r', alpha = 0.1) 
lins3 = ax.plot(rcp85_year,rcp85_CSDC[1,:],'-',color="r",linewidth=1,label='RCP85')
lins4 = ax.fill_between(rcp85_year,rcp85_CSDC[2,:],rcp85_CSDC[3,:], where=rcp85_CSDC[2,:]>=rcp85_CSDC[3,:], color="r", alpha = 0.1) 
# lins5 = ax.plot(rcp45_year,rcp45_CSDC[1,:],'-',color="g",linewidth=1,label='RCP45')
# lins6 = ax.fill_between(rcp45_year,rcp45_CSDC[2,:],rcp45_CSDC[3,:], where=rcp45_CSDC[2,:]>=rcp45_CSDC[3,:], color="g", alpha = 0.1) 
lins7 = ax.plot(fixa_year,fixa_CSDC[1,:],'-', color="b",linewidth=1,label='RCP85_FixA')
lins8 = ax.fill_between(fixa_year,fixa_CSDC[2,:],fixa_CSDC[3,:], where=fixa_CSDC[2,:]>=fixa_CSDC[3,:], color="b", alpha = 0.1) 


plt.subplots_adjust(left=0.1, bottom=0.03, right=0.98, top=0.95, wspace=0.1, hspace=0.1); 
plt.savefig('Temp.png', format='png', dpi=1000)

