# -*- coding: utf-8 -*-
"""
This is to plot the PDFs of daily P and T to reveal the roles of GHGs, Aerosols, Trop. O3 and Strat O3 in changing these varaibles
interested regions include Global (ocean+land) global land, E.Asia, Europe and America
v2 includes GPCC data but failed to produce becaue GPCC has monthly data
"""
import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
from numpy import linspace
from scipy.interpolate import interp2d  as interp2d
import scipy.io as sio
import math as math
import cPickle as pkl

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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
oceanmask=spio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
"""
def AreaWeight(lon1,lon2,lat1,lat2):
	'''
	calculate the earth radius in m2
	'''
	radius = 6371000;
	area = (math.pi/180)*np.power(radius,2)*np.abs(lon1-lon2)*\
	(np.abs(np.sin(np.radians(lat1))-np.sin(np.radians(lat2))))
	# print np.nansum(np.nansum(area,axis=1),axis=0)
	return area
	
def box_clip(lon_s,lon_e,lat_s,lat_e,lon,lat,mask):
	# fill the range outside the box with 0
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
	mask[:,:colum_s] =0; mask[:,colum_e:] =0
	mask[:row_s,:] =0; mask[row_e:,:] =0
	return mask

def mask_weight(region_key,lon,lat,reverse=False):
	# Read in the country mask
	# interpolate it to the required resolution grids with lon_interp,lat_interp 
	
	lon_res = lon[1] - lon[0];lat_res = lat[1] - lat[0];
	lons,lats = np.meshgrid(lon,lat)
	area = AreaWeight(lons,lons+lon_res,lats,lats+lat_res)

	##OCEAN_MASKS FOR COUNTRIES
	ocean_mask = sio.loadmat('/home/s1667168/coding/python/external_data/Euro_USA_AUS_BRICS_STA_720_360.mat')
	lon_mask = ocean_mask['lon'][0,:];
	lat_mask = ocean_mask['lat'][0,:];
	box_region_dic={'All':[0,360,-90,90],'ASIA':[65,145,5,45],'US':[240,290,30,50],'ARCTIC':[0,360,60,90],'TROPICS':[0,360,-28,28],'EUROPE':[0,40,30,70],}
	if (region_key == 'USA' or region_key == 'Europe' or region_key == 'India' or region_key == 'China' or region_key == 'GloLand'):
		mask= ocean_mask[region_key][:]
	elif  region_key in box_region_dic:
		mask= ocean_mask['All'][:]
		box = box_region_dic[region_key]
		mask = box_clip(box[0],box[1],box[2],box[3],lon_mask,lat_mask,mask)
	else:
		print "error region name"
	
	# interpolate from 360*720 to 192*288
	mask[np.isnan(mask)]=0;	mask[mask>0]=1;
	f = interp2d(lon_mask, lat_mask, mask,kind='linear'); mask = f(lon, lat);
	# plt.imshow(mask,origin='lower');plt.show()
	mask[mask >= 1] = 1;mask[mask < 1] = 0;
	# weight each grid cell by its area weight against the total area
	if reverse:
		mask=1-mask
	mask[mask==0] = np.nan
	mask=np.multiply(mask,area); 
	mask_weighted = np.divide(mask,np.nansum(np.nansum(mask,axis=1),axis=0))
	return mask_weighted

def data_readin(variable):	
	def day2datetime(scenario,days):
		###
		# convert days from a reference into int datetime 
		# do not take leap years into account
		###
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
		
		
	def RegionMean(scenario,time,lon,lat,data):
		regions ={'All':np.empty((30*365)),'GloLand':np.empty((30*365)),'ASIA':np.empty((30*365)),'Europe':np.empty((30*365)),'USA':np.empty((30*365))}  #,'USA','SESA','EA','India'
		for region_key in regions:
			mask = mask_weight(region_key,lon,lat); #print mask
			regions[region_key][0:np.shape(data)[0]] = np.nansum(np.nansum(np.multiply(mask,data),axis=2),axis=1)
		return regions
	
	def data_netcdf(scenario,variable):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/FullCp/'
		if variable == "TS":
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.'+variable+'.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			data = nc_fid.variables[variable][:]-273.15
			data[data<-100] = np.nan;data[data>100] = np.nan;
			nc_fid.close()
		elif variable == "Pr":
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.PRECL.nc'
			nc_fid = nc4.Dataset(var_path,mode='r')
			PRECL = nc_fid.variables['PRECL'][:]
			var_path = input_path+scenario+'/mon/atm/'+scenario+'.atm.mon.PRECC.nc'
			nc_fid.close()
			nc_fid = nc4.Dataset(var_path,mode='r')
			lat = nc_fid.variables['lat'][:]
			lon = nc_fid.variables['lon'][:]
			days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
			PRECC = nc_fid.variables['PRECC'][:]		
			nc_fid.close()
			data = (PRECC+PRECL)*24*60*60*1000
			data[data<0] = np.nan;data[data>200] = np.nan;
		RegionM = RegionMean(scenario,time,lon,lat,data)
		# print RegionM
		return RegionM
	
	Edg70GO = data_netcdf('Edg70GO',variable)
	T1970 = data_netcdf('T1970RCP',variable)
	EdgRef = data_netcdf('EdgRef',variable)
	Edg70Oz = data_netcdf('Edg70Oz',variable)
	Edg70T10SOZ = data_netcdf('Edg70T10SOZ',variable)

	return T1970,Edg70GO,Edg70Oz,Edg70T10SOZ,EdgRef

	
	
def GPCC_PP():
	def GPCC_pr(period):
		def RegionMeanGPCC(lon,lat,data):
			regions ={'All':np.empty((15*12)),'GloLand':np.empty((15*12)),'ASIA':np.empty((15*12)),'Europe':np.empty((15*12)),'USA':np.empty((30*12))}  #,'USA','SESA','EA','India'
			for region_key in regions:
				mask = mask_weight(region_key,lon,lat); #print mask
				regions[region_key][0:np.shape(data)[0]] = np.nansum(np.nansum(np.multiply(mask,data),axis=2),axis=1)
			return regions
		def day2datetime(days):
			from datetime import date, timedelta
			date_int = np.empty((len(days)));date_int[:]=np.nan;
			start = date(1800,1,1)
			ith = 0
			for day in days:
				delta = timedelta(int(day))
				offset = start + delta
				date_int[ith] = offset.year*10000+offset.month*100+offset.day
				ith=ith+1
			return date_int			
		
		input_path ='/exports/csce/datastore/geos/users/s1667168/obs/Pr/'
		var_path = input_path+'GPCCprecip.mon.total.1x1.v7.nc'
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat =nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		time1 = nc_fid.variables['time'][:];
		data1 = nc_fid.variables['precip'][:]
		
		input_path ='/exports/csce/datastore/geos/users/s1667168/obs/Pr/'
		var_path = input_path+'GPCCprecip.monitor.mon.total.1x1.v4.nc'
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat =nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		time2 = nc_fid.variables['time'][:];
		data2 = nc_fid.variables['precip'][84:-2,:,:]
		data = np.concatenate((data1,data2),axis=0)
		time = np.concatenate((time1,time2),axis=0)
		data[data<0]= np.nan; data[data>3000]= np.nan; 
		data=data/30
		nc_fid.close()

		time=day2datetime(time);print time
		if period == '1970':
			layer_b = [layer for layer in range(len(time)) if time[layer]//100 == 1963*100+1][0]
			layer_e = [layer for layer in range(len(time)) if time[layer]//100 == 1977*100+12][0]
			Pr15yr = data[layer_b:layer_e+1,:,:]
		elif period == '2010':
			layer_b = [layer for layer in range(len(time)) if time[layer]//100 == 2003*100+1][0]
			layer_e = [layer for layer in range(len(time)) if time[layer]//100 == 2017*100+12][0]
			Pr15yr = data[layer_b:layer_e+1,:,:]
		RegionM = RegionMeanGPCC(lon,lat,Pr15yr)
		return RegionM
	Pr1670 = GPCC_pr('1970');
	Pr2010 = GPCC_pr('2010');
	return Pr1670,Pr2010
	
	
	
def TotDecomposePDFs(variable):
	####
	# 1.. input regional mean of variables in a dic from all the differenct simulations
	# 2.. decompose into GHGs, AAs and ozones and recombine one-by-one-by-one
	# 3.. produce PDFs
	####
	def kde_bins(data,bins,variable):
		if variable == 'Pr':
			data[data>500] = np.nan;
		elif variable == 'TS':
			data[data>100] = np.nan;data[data<-100] = np.nan;
		data_m = data[~np.isinf(data)];
		data_mm = data_m[~np.isnan(data_m)];  # print np.nanmin(data_mm),np.nanmax(data_mm),
		pdfs,_ = np.histogram(data_mm,bins=bins,density=True) #range=(bin_min,bin_max),
		pdfs[pdfs==0]=np.nan			
		mean = np.nanmean(data_mm); #print np.shape(bins[:-1]),np.shape(pdfs)
		P25 =  np.nanpercentile(data_mm,25);
		P75 =  np.nanpercentile(data_mm,75);
		return pdfs,mean,P25,P75  #integrate.cumtrapz(hist, initial=0.)		
	BASE,B_AA,B_AA_GHG,B_AA_GHG_StrO3,All = data_readin(variable)
	GPCC70,GPCC10 = GPCC_PP()
	if  variable == 'TS':
		bin_min = -50;bin_max = 50;
		bins = np.arange(bin_min,bin_max,.5)
	elif variable == 'Pr':
		bin_min = 0;bin_max = 30;
		bins = np.arange(bin_min,bin_max,.1)
	
	PDFs = {'All':np.empty((7,len(bins)-1)),'GloLand':np.empty((7,len(bins)-1)),'ASIA':np.empty((7,len(bins)-1)),'Europe':np.empty((7,len(bins)-1)),'USA':np.empty((7,len(bins)-1))}
	means = {'All':np.empty((7)),'GloLand':np.empty((7)),'ASIA':np.empty((7)),'Europe':np.empty((7)),'USA':np.empty((7))}
	P25s = {'All':np.empty((7)),'GloLand':np.empty((7)),'ASIA':np.empty((7)),'Europe':np.empty((7)),'USA':np.empty((7))}
	P75s = {'All':np.empty((7)),'GloLand':np.empty((7)),'ASIA':np.empty((7)),'Europe':np.empty((7)),'USA':np.empty((7))}
	
	for region_key in PDFs:
		PDFs[region_key][0,:],means[region_key][0],P25s[region_key][0],P75s[region_key][0] = kde_bins(BASE[region_key],bins,variable)
		PDFs[region_key][1,:],means[region_key][1],P25s[region_key][1],P75s[region_key][1] = kde_bins(B_AA[region_key],bins,variable)
		PDFs[region_key][2,:],means[region_key][2],P25s[region_key][2],P75s[region_key][2] = kde_bins(B_AA_GHG[region_key],bins,variable)
		PDFs[region_key][3,:],means[region_key][3],P25s[region_key][3],P75s[region_key][3] = kde_bins(B_AA_GHG_StrO3[region_key],bins,variable)
		PDFs[region_key][4,:],means[region_key][4],P25s[region_key][4],P75s[region_key][4] = kde_bins(All[region_key],bins,variable)
		PDFs[region_key][5,:],means[region_key][5],P25s[region_key][5],P75s[region_key][5] = kde_bins(GPCC70[region_key],bins,variable)
		PDFs[region_key][6,:],means[region_key][6],P25s[region_key][6],P75s[region_key][6] = kde_bins(GPCC10[region_key],bins,variable)
				
	return bins[:-1],PDFs,means,P25s,P75s	
					

input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/analysis/PDFs/'
# variable = 'TS'; bins_TS,PDFs_TS,means_TS,P25_TS,P75_TS = TotDecomposePDFs(variable)
# TS_DICS = [bins_TS, PDFs_TS,means_TS,P25_TS,P75_TS]
# pkl.dump( TS_DICS, open(input_path+"TS_PDF_mean_var_DICS.p", "wb" ) )

variable = 'Pr'; bins_Pr,PDFs_Pr,means_Pr,P25_Pr,P75_Pr = TotDecomposePDFs(variable)
Pr_DICS = [bins_Pr, PDFs_Pr,means_Pr,P25_Pr,P75_Pr]
pkl.dump( Pr_DICS, open(input_path+"Pr_monthly_PDF_mean_var_Dics.p", "wb" ) )
"""

def plot_pdfs(ax,x,PDFs,means,region_key,smth,abs_change =False):
	if abs_change:
		ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
		ax.plot(x+0.25,PDFs[region_key][0,:],'-',color="b",linewidth=1,label = 'Baseline (E1970)')
		ax.plot(x+0.25,PDFs[region_key][1,:],'-',color="r",linewidth=1,label = 'AAs')
		ax.plot(x+0.25,PDFs[region_key][2,:],'-',color="g",linewidth=1,label = 'AAs+GHGs')
		ax.plot(x+0.25,PDFs[region_key][3,:],'-',color="m",linewidth=1,label = 'AAs+GHGs+Str.O3')
		ax.plot(x+0.25,PDFs[region_key][4,:],'-',color="k",linewidth=1,label = 'AAs+GHGs+Str.O3+Tro.O3 (E2010)')
		
	else:
		ax.axhline(y=0,color='k',linewidth = 1,ls='-');
		ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
		# ax.plot(x+0.25,smooth(PDFs[region_key][6,:]-PDFs[region_key][5,:],smth),'--',color="k",linewidth=1.5,label = 'GPCC')
		ax.plot(x+0.25,smooth(PDFs[region_key][4,:]-PDFs[region_key][0,:],smth),'-',color="k",linewidth=1,label = 'All')
		ax.plot(x+0.25,smooth(PDFs[region_key][2,:]-PDFs[region_key][1,:],smth),'-',color="g",linewidth=1,label = 'GHGs')
		ax.plot(x+0.25,smooth(PDFs[region_key][1,:]-PDFs[region_key][0,:],smth),'-',color="r",linewidth=1,label = 'AAs')
		ax.plot(x+0.25,smooth(PDFs[region_key][4,:]-PDFs[region_key][3,:],smth),'-',color="m",linewidth=1,label = 'Trop. O3')
		ax.plot(x+0.25,smooth(PDFs[region_key][3,:]-PDFs[region_key][2,:],smth),'-',color="b",linewidth=1,label = 'Strat. O3')	
		plt.xlabel('Pr '+r'($\mathrm{\mathsf{\Delta}\/mm\/day^{-1}}$)');plt.ylabel(r'$\mathrm{\mathsf{\Delta}PDF}$'+' (%)')
		ax2 = ax.twinx();
		ax2.plot(x+0.25,PDFs[region_key][0,:],'--',color="k",linewidth=2,label = 'Baseline (E1970)')
	return ax,ax2

def plot_mean_var(ax,means,P25s,P75s,region_key,percent=False):
	ax.axvline(x=0,color='k',linewidth = 1,ls='-');
	means_all = means[region_key][:];
	# P25s = P25s[region_key][:]-means;print P25s
	# P75s = P75s[region_key][:]-P25s;print P75s
	means = means_all-means_all[0]; 
	# var = P75s - P25s
	ax.tick_params(axis='both', which='major', direction = 'out',left=True,right=True,bottom=True,top=False, pad=5)
	mean_decomposed = np.array([means[2]-means[1],means[1]-means[0],means[4]-means[3],means[3]-means[2]])
	if percent: mean_decomposed =mean_decomposed*100/means_all[0]
	# var_decomposed = np.array([means[2]-means[1],means[1]-means[0],means[4]-means[3],means[3]-means[2]])	
	ax.barh(0.5,mean_decomposed[0],color='g',align='center',height=0.3)
	ax.barh(1.5,mean_decomposed[1],color='r',align='center',height=0.3)
	ax.barh(2.5,mean_decomposed[2],color='m',align='center',height=0.3)
	ax.barh(3.5,mean_decomposed[3],color='b',align='center',height=0.3)
	
	ax.set_ylim([0,4]);ax.set_yticks(np.arange(0.5,4,1));ax.set_yticklabels(['GHGs','AAs','Trop. O3','Strat. O3']);

input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/analysis/PDFs/'
# TS_DICS = pkl.load(open(input_path+"TS_PDF_mean_var_DICS.p", "rb") )
# bins_TS = TS_DICS[0];PDFs_TS = TS_DICS[1];means_TS = TS_DICS[2];P25_TS = TS_DICS[3];P75_TS = TS_DICS[4];

Pr_DICS = pkl.load(open(input_path+"Pr_PDF_mean_var_Dics_v1.p", "rb") )
bins_Pr = Pr_DICS[0];PDFs_Pr = Pr_DICS[1];means_Pr = Pr_DICS[2];P25_Pr = Pr_DICS[3];P75_Pr = Pr_DICS[4];

# ax = plt.subplot2grid((3, 4), (0, 0), colspan=3)		
fig = plt.figure(figsize=(16, 8), facecolor='White');pad=5
bin_min =-50;bin_max=50

print np.shape(PDFs_Pr)

######SA
region_key = 'All'
ax = plt.subplot2grid((2, 10), (0, 0), colspan=3)			
ax.annotate('(a) Global area-weighted mean daily precipitation',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=18)			
ax,ax2 = plot_pdfs(ax,bins_Pr,PDFs_Pr,means_Pr,region_key,smth=1,abs_change =False)

ax.set_xlim([2.5,4.0]);ax.set_xticks(np.arange(2.5,4.0+.1,0.5));
ax.set_ylim([-1,1]);ax.set_yticks(np.arange(-1,1.01,.5));

ax = plt.subplot2grid((2, 10), (0, 3), colspan=2)
plot_mean_var(ax,means_Pr,P25_Pr,P75_Pr,region_key,percent=True)
ax.set_xlim([-2,3]);ax.set_xticks(np.arange(-2,3.1,1));
plt.xlabel(r'$\mathrm{\mathsf{\Delta}}$'+' Pr (%)');

region_key = 'ASIA'
ax = plt.subplot2grid((2, 10), (0, 5), colspan=3)		
ax.annotate('(c) Asia area-weighted mean daily precipitation',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=18)			
ax,ax2 = plot_pdfs(ax,bins_Pr,PDFs_Pr,means_Pr,region_key,smth=20,abs_change =False)
ax.set_xlim([1,8]);ax.set_xticks(np.arange(1,8+.1,1));
ax.set_ylim([-0.02,0.02]);ax.set_yticks(np.arange(-0.02,0.021,.01));
ax = plt.subplot2grid((2, 10), (0, 8), colspan=2)
plot_mean_var(ax,means_Pr,P25_Pr,P75_Pr,region_key,percent=True)
ax.set_xlim([-5,5]);ax.set_xticks(np.arange(-5,5.1,2.5));
plt.xlabel(r'$\mathrm{\mathsf{\Delta}}$'+' Pr (%)');

region_key = 'Globe'
ax = plt.subplot2grid((2, 10), (1, 0), colspan=3)		
ax.annotate('(b) Global land area-weighted mean daily precipitation',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=18)			
ax,ax2 = plot_pdfs(ax,bins_Pr,PDFs_Pr,means_Pr,region_key,smth=1,abs_change =False)
ax.set_xlim([2.0,3.5]);ax.set_xticks(np.arange(2.0,3.5+.1,0.5));
ax.set_ylim([-0.3,0.3]);ax.set_yticks(np.arange(-0.3,0.31,.1));
ax = plt.subplot2grid((2, 10), (1, 3), colspan=2)
plot_mean_var(ax,means_Pr,P25_Pr,P75_Pr,region_key,percent=True)
ax.set_xlim([-2,3]);ax.set_xticks(np.arange(-2,3.1,1));
plt.xlabel(r'$\mathrm{\mathsf{\Delta}}$'+' Pr (%)');


region_key = 'Europe'
ax = plt.subplot2grid((2, 10), (1, 5), colspan=3)		
ax.annotate('(d) Europe area-weighted mean daily precipitation',(0.02, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=18)			
ax,ax2 = plot_pdfs(ax,bins_Pr,PDFs_Pr,means_Pr,region_key,smth=6,abs_change =False)
ax.set_xlim([0.5,5.5]);ax.set_xticks(np.arange(0.5,5.5+.1,1));
ax.set_ylim([-0.09,0.09]);ax.set_yticks(np.arange(-0.09,0.091,.03));
ax = plt.subplot2grid((2, 10), (1, 8), colspan=2)
plot_mean_var(ax,means_Pr,P25_Pr,P75_Pr,region_key,percent=True)
ax.set_xlim([-5,5]);ax.set_xticks(np.arange(-5,5.1,2.5));
plt.xlabel(r'$\mathrm{\mathsf{\Delta}}$'+' Pr (%)');


ax = fig.add_axes([0.20, -0.00, 0.70, 0.10])
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
kd_patch = mlines.Line2D([],[],color='k',ls='--',lw=2,label='Clim')
k_patch = mlines.Line2D([],[],color='k', lw=1.5,label='All')
g_patch = mlines.Line2D([],[],color='g',lw=1.5,label='GHGs')
r_patch = mlines.Line2D([],[],color='r', lw=1.5,label='AAs')
m_patch = mlines.Line2D([],[],color='m', lw=1.5,label='Trop. O3')
b_patch = mlines.Line2D([],[],color='b', lw=1.5,label='Strat. O3')
lines = [kd_patch,k_patch,g_patch,r_patch,m_patch,b_patch]
labels = ['Climatology','All','GHGs','AAs','Trop. O3','Strat. O3']
legend = plt.legend(lines,labels,ncol=6,loc='lower left',labelspacing=0.2)
legend.get_frame().set_facecolor('none');legend.get_frame().set_edgecolor('none');legend.get_frame().set_alpha(1)
ax.axis('off')

# region_key = 'Globe'
# ax = plt.subplot2grid((2, 10), (1, 0), colspan=3)			
# ax.annotate('(c) Global land area-mean daily temperature',(0.9, 1.05), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='center', va='baseline',rotation='horizontal',fontsize=15)			
# ax = plot_pdfs(ax,bins_TS,PDFs_TS,means_TS,region_key,smth=1,abs_change =False)
# plt.xlabel('SAT (K)');plt.ylabel(r'$\mathrm{\mathsf{\Delta}PDF}$'+' (%)')
# # ax.set_xlim([bin_min,bin_max]);ax.set_xticks(np.arange(bin_min,bin_max+.1,(bin_max-bin_min)/6));
# # ax.set_ylim([0,0.2]);ax.set_yticks(np.arange(0,.21,.05));

# ax = plt.subplot2grid((2, 10), (1, 3), colspan=2)
# plot_mean_var(ax,means_TS,P25_TS,P75_TS,region_key)
# ax.set_xlim([-0.5,1.5]);ax.set_xticks(np.arange(-0.5,1.51,0.5));
# plt.xlabel(r'$\mathrm{\mathsf{\Delta}SAT\/(K)}$');

# plt.show()
plt.subplots_adjust(left=0.07, bottom=0.13, right=0.98, top=0.94, wspace=4.0, hspace=0.50);
plt.savefig('Global_PDFs_daily_Pr.png', format='png', dpi=1000)
