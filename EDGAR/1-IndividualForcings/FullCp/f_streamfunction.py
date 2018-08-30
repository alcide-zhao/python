"""
This is to process the monthly fields into zonal mean and plot 
	currently, SAT ans precipitation
	the model internal variability (25-t5 th percentile) and mean are both plotted
"""
import site
import os
import numpy as np
import netCDF4 as nc4
from scipy import stats
import scipy.io as sio
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d  as interp1d
from scipy import integrate

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir,os.path.pardir,os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *


	
def data_read_annual_mean(variable,levs):	
	"""
	calculate the earth radius in m2
	"""	
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
		
		
	def mon_mean2annual_mean(scenario,time,data,levs):
		
		calendar_day = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
		if scenario=='T1970RCP':
			year_series = range(2020,2050)
		elif scenario=='EdgEne':
			year_series = range(2200,2230)
		elif scenario=='Edg70GO':
			year_series = range(2070,2100)
		else:
			year_series = range(2130,2160)
		for iyear in year_series:
			if (iyear == year_series[0] and time[0]//100 >= year_series[0] *100+1):
				layer_b=0
			else:
				layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+1][0]  #June01
			if (iyear == year_series[-1] and time[-1]//100 <= year_series[-1] *100+12):
				layer_e=-2
			else:
				layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+12][0]  #August 31
			# print iyear,layer_b
			if levs:  # 4d DATA
				annual_mean=np.empty((40,30,192,288));annual_mean[:]=np.nan
				data_cache = data[layer_b:layer_e+1,:,:,:]
				annual_mean[iyear-year_series[0],:,:,:] = stats.nanmean(data_cache,axis=0)
			else:
				annual_mean=np.empty((40,192,288));annual_mean[:]=np.nan
				data_cache = data[layer_b:layer_e+1,:,:]
				annual_mean[iyear-year_series[0],:,:] = stats.nanmean(data_cache,axis=0)
		annual_mean = np.nanmean(annual_mean,axis=0)
		return annual_mean

	def data_netcdf(scenario,FREQ,variable,levs):
		input_path ='/exports/csce/datastore/geos/users/s1667168/CESM_EDGAR/ModelOutput/'
		var_path = input_path+scenario+'/'+FREQ+'/atm/'+scenario+'.atm.'+FREQ+'.'+variable+'.nc'
		# print var_path
		nc_fid = nc4.Dataset(var_path,mode='r')
		lat = nc_fid.variables['lat'][:]
		lon = nc_fid.variables['lon'][:]
		days = nc_fid.variables['time'][:]; time = day2datetime(scenario,days);#print time
		data = nc_fid.variables[variable][:]# 850hpa
		VarAnMean = mon_mean2annual_mean(scenario,time,data,levs)
		if levs: lev = nc_fid.variables['lev'][:]
		else : lev =np.nan	
		nc_fid.close()
		return lev,lat,lon,VarAnMean
		
	FREQ = 'mon'
	lev,lat,lon,Edg70GO = data_netcdf('Edg70GO',FREQ,variable,levs)
	lev,lat,lon,T1970 = data_netcdf('T1970RCP',FREQ,variable,levs)
	lev,lat,lon,EdgRef = data_netcdf('EdgRef',FREQ,variable,levs)
	lev,lat,lon,Edg70Oz = data_netcdf('Edg70Oz',FREQ,variable,levs)
	#lev,lat,lon,EdgEne = data_netcdf('EdgEne',FREQ,variable,levs)
	#lev,lat,lon,EdgTech = data_netcdf('EdgTech',FREQ,variable,levs)
	return lev,lat,lon,T1970,Edg70GO,Edg70Oz,EdgRef  #,EdgEne,EdgTech

	
def hybrid2p_interp(PS,midpoint=True,interface=False):
	"""
	This is to convert the CESM hybrid vertical coordinate to the pressure system
	datain must be an array of 3, 4, or 5 dimensions. Needs to contain a level dimension in hybrid coordinates. 
	The order of the dimensions is specific. The three rightmost dimensions must be 
	level x lat x lon [e.g. T(time,lev,lat,lon)]. 
	The order of the level dimension must be top-to-bottom
	"""
	
	PO = 1000   #surface reference pressure
	if midpoint:
		hya = np.array([ 0.00364346569404006, 0.00759481964632869, 0.0143566322512925,\
			0.0246122200042009, 0.0382682997733355, 0.0545954797416925,\
			0.0720124505460262, 0.0878212302923203, 0.103317126631737,\
			0.121547240763903, 0.142994038760662, 0.168225079774857,\
			0.178230673074722, 0.170324325561523, 0.161022908985615,\
			0.150080285966396, 0.137206859886646, 0.122061938047409,\
			0.104244712740183, 0.0849791541695595, 0.0665016956627369,\
			0.0501967892050743, 0.037188658490777, 0.028431948274374,\
			0.0222089774906635, 0.016407382208854, 0.0110745579004288,\
			0.00625495356507599, 0.00198940909467638, 0])
		hyb = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0196774136275053,\
			0.062504293397069, 0.112887907773256, 0.172161616384983,\
			0.241894043982029, 0.323930636048317, 0.420442461967468,\
			0.524799540638924, 0.624887734651566, 0.713207691907883,\
			0.783669710159302, 0.831102818250656, 0.864811271429062,\
			0.896237164735794, 0.92512384057045, 0.951230525970459,\
			0.974335998296738, 0.992556095123291])
	elif interface:
		hya = np.array([ 0.00225523952394724, 0.00503169186413288, 0.0101579474285245,\
			0.0185553170740604, 0.0306691229343414, 0.0458674766123295,\
			0.0633234828710556, 0.0807014182209969, 0.0949410423636436,\
			0.11169321089983, 0.131401270627975, 0.154586806893349,\
			0.181863352656364, 0.17459799349308, 0.166050657629967,\
			0.155995160341263, 0.14416541159153, 0.130248308181763,\
			0.113875567913055, 0.0946138575673103, 0.0753444507718086,\
			0.0576589405536652, 0.0427346378564835, 0.0316426791250706,\
			0.0252212174236774, 0.0191967375576496, 0.0136180268600583,\
			0.00853108894079924, 0.00397881818935275, 0, 0])
		hyb = np.array([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0393548272550106,\
			0.0856537595391273, 0.140122056007385, 0.204201176762581,\
			0.279586911201477, 0.368274360895157, 0.47261056303978,\
			0.576988518238068, 0.672786951065063, 0.753628432750702,\
			0.813710987567902, 0.848494648933411, 0.881127893924713,\
			0.911346435546875, 0.938901245594025, 0.963559806346893,\
			0.985112190246582, 1])
	Pnew = np.empty((len(hya),np.shape(PS)[0],np.shape(PS)[1])); Pnew[:] = np.nan
	
	for ilev in range(np.shape(hya)[0]):
		Pnew[ilev,:,:] = PO * hya[ilev]+ hyb[ilev] *PS/100
	
	# p_interp=np.array([10, 20, 30, 50, 70,100,150,200,250,\
                     # 300,400,500,600,700,850,925,1000]
	# data_new = np.array([np.shape(Pnew)[0],np.shape(PS)[0],np.shape(PS)[1]]); P[:] = np.nan
	# ## interpolate to the desired pressure levels
	# for ilat in range(np.shape(PS)[0]):
		# for ilon in range(np.shape(PS)[1]):
			# data_cache = datain[:,ilat,ilon]
			# pressure = Pnew[:,ilat,ilon]
			# f = interp1d(pressure, data_cache); 
			# data_interp[:,ilat,ilon] = f(p_interp);	 
	return Pnew       #, p_interp, data_interp


 
def streamfunction_zonal_mean(V,lev,lat,lon):
	'''Return the zonal mean Reconstruction with input V'''
	a = 6.373*10**6
	g = 10
	lat_rad = np.deg2rad(lat)
	factor = 2 * np.math.pi *a*np.cos(lat_rad)/g
	lev_new = np.array([10, 20, 30, 50,70,100,150,200,250,\
                     300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])
	stf = np.empty((len(lev_new),len(lat),len(lon)))
	for ilat in range(len(lat)):
		lat_rad = np.deg2rad(lat[ilat])
		factor = 2 * np.math.pi *a*np.cos(lat_rad)/g
		for ilon in range(len(lon)):
			integration = integrate.cumtrapz(np.flip(V[:,ilat,ilon],axis=0),x=np.flip(lev[:,ilat,ilon],axis=0)*100, initial=0.)  ## integrate from surface to a pressure level
			cache  = np.flip(integration,axis=0) * factor
			f = interp1d(lev[:,ilat,ilon], cache,bounds_error=False,fill_value = 0); 
			stf[:,ilat,ilon] = f(lev_new);
	stf =stf*10**(-9)     # 10**(9) kg/s
	# stf_zm = np.nanmean(stf,axis=2)
	return lev_new,stf
	
def nterpolate2pressure(DataIn,lev,lat,lon):
	lev_new = np.array([10, 20, 30, 50,70,100,150,200,250,
		300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000])
	interp = np.empty((len(lev_new),len(lat),len(lon)));interp[:] = np.nan
	for ilat in range(len(lat)):
		for ilon in range(len(lon)):
			f = interp1d(lev[:,ilat,ilon], DataIn[:,ilat,ilon],bounds_error=False,fill_value = 0); 
			interp[:,ilat,ilon] = f(lev_new); 
	return  lev_new,interp 

def decompose_stf(filed):
	variable='PS'; lev,lat,lon,PS_T1970,PS_Edg70GO,PS_Edg70Oz,PS_EdgRef = data_read_annual_mean(variable,levs=False)
	P_T1970 = hybrid2p_interp(PS_T1970,midpoint=True,interface=False)
	P_Edg70GO = hybrid2p_interp(PS_Edg70GO,midpoint=True,interface=False)
	P_Edg70Oz = hybrid2p_interp(PS_Edg70Oz,midpoint=True,interface=False)
	P_EdgRef = hybrid2p_interp(PS_EdgRef,midpoint=True,interface=False)
	
	if filed == "stf":
		variable='V'; lev,lat,lon,V_T1970,V_Edg70GO,V_Edg70Oz,V_EdgRef = data_read_annual_mean(variable,levs=True)
		lev_new,T1970 = streamfunction_zonal_mean(V_T1970,P_T1970,lat,lon)
		lev_new,Edg70GO = streamfunction_zonal_mean(V_Edg70GO,P_Edg70GO,lat,lon)
		lev_new,Edg70Oz = streamfunction_zonal_mean(V_Edg70Oz,P_Edg70Oz,lat,lon)
		lev_new,EdgRef = streamfunction_zonal_mean(V_EdgRef,P_EdgRef,lat,lon)
	else:
		variable=filed; lev,lat,lon,T1970,Edg70GO,Edg70Oz,EdgRef = data_read_annual_mean(variable,levs=True)
		lev_new,T1970 = nterpolate2pressure(T1970,P_T1970,lat,lon)
		lev_new,Edg70GO = nterpolate2pressure(Edg70GO,P_Edg70GO,lat,lon)
		lev_new,Edg70Oz = nterpolate2pressure(Edg70Oz,P_Edg70Oz,lat,lon)
		lev_new,EdgRef = nterpolate2pressure(EdgRef,P_EdgRef,lat,lon)
	
	total = np.nanmean(EdgRef  - T1970,axis=2)
	GHG =   np.nanmean( Edg70Oz - Edg70GO,axis=2)
	AAs =   np.nanmean( Edg70GO - T1970,axis=2)
	O3 =    np.nanmean( EdgRef  - Edg70Oz,axis=2)	
	return lev_new,lat,total,GHG,AAs,O3,np.nanmean(T1970,axis=2)



# filed = 'T';lev_new,lat,T_total,T_GHG,T_AAs,T_O3,T_1970 = decompose_stf(filed)
filed = 'VQ';lev_new,lat,VQ_total,VQ_GHG,VQ_AAs,VQ_O3,VQ_1970 = decompose_stf(filed)
VQ_total=VQ_total*10**(3);VQ_GHG=VQ_GHG*10**(3);VQ_AAs=VQ_AAs*10**(3);VQ_O3=VQ_O3*10**(3);VQ_1970 =VQ_1970*10**(3)   # from kg/kg/s to g/kg m/s

filed = 'stf';lev_new,lat,stf_total,stf_GHG,stf_AAs,stf_O3,stf_1970 = decompose_stf(filed)

def plot_shade_plot(ax,x,y,shade,ContourData,colormap,colorbar_min,colorbar_max,lev_sel):
	
	matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
	xi, yi = np.meshgrid(x, y)
	# y = np.flipud(y); shade = np.flipud(shade);ContourData = np.flipud(ContourData);
	masked_obj = np.ma.masked_where(np.isnan(shade), shade)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over('r'); #cmap.set_under('darkmagenta');
	CS1 = plt.pcolor(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)#
	# clevs = np.round(np.arange(np.nanmin(ContourData[:]),np.nanmax(ContourData[:]),(np.nanmax(ContourData[:])-np.nanmin(ContourData[:]))/8),decimals=0)
	if lev_sel == "STF":
		clevs=np.array([-100,-80,-60,-40,-20,20,40,60,80,100])
	elif lev_sel =="VQ":
		clevs=np.array([-8,-6,-4,-2,-1,1,2,4,6,-8])
	CS2 = plt.contour(xi,yi,ContourData,levels = clevs,linewidths=1,colors='k')
	plt.clabel(CS2,clevs,inline=False,fontsize=10,fmt= '%1.0f',inline_spacing=10) #
	ax.set_xlim([-70,70]);ax.set_xticks(np.arange(-60,61,30));ax.set_xticklabels(('60S','30S','EQ','30N','60N'));
	ax.set_ylim([200,1000]);ax.set_yticks([200,300,400,500,600,700,800,850,925,1000]);
	ax.axvline(x=0,c='k',ls=':',lw=3)
	
	plt.gca().invert_yaxis()
	return CS1

	
fig = plt.figure(facecolor='White',figsize=[12,14]);plot_setup();pad= 5
colormap='coolwarm'


####STREAM FUNCITION
colorbar_min=-20;colorbar_max=20;
ax = plt.subplot(4,2,1);
ax.annotate('(a)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
ax.annotate('2010 - 1970',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)						
plot_shade_plot(ax,lat,lev_new,stf_total,stf_1970,colormap,colorbar_min,colorbar_max,lev_sel ="STF")

ax = plt.subplot(4,2,3);
ax.annotate('(b)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
ax.annotate('GHGs',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)						
plot_shade_plot(ax,lat,lev_new,stf_GHG,stf_1970,colormap,colorbar_min,colorbar_max,lev_sel ="STF")

ax = plt.subplot(4,2,5);
ax.annotate('(c)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
ax.annotate('AAs',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)						
plot_shade_plot(ax,lat,lev_new,stf_AAs,stf_1970,colormap,colorbar_min,colorbar_max,lev_sel ="STF")

ax = plt.subplot(4,2,7);
ax.annotate('(d)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
ax.annotate('O3',xy=(-0.1,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)					
colormesh1 = plot_shade_plot(ax,lat,lev_new,stf_O3,stf_1970,colormap,colorbar_min,colorbar_max,lev_sel ="STF")
cbar_ax = fig.add_axes([0.06, 0.03, 0.44, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/Tg\/s^{-1}}$',xy=(0.5,1.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)		


####VQ
colorbar_min=-1;colorbar_max=1;
ax = plt.subplot(4,2,2);
ax.annotate('(e)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
plot_shade_plot(ax,lat,lev_new,VQ_total,VQ_1970,colormap,colorbar_min,colorbar_max,lev_sel ="VQ")

ax = plt.subplot(4,2,4);
ax.annotate('(f)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
plot_shade_plot(ax,lat,lev_new,VQ_GHG,VQ_1970,colormap,colorbar_min,colorbar_max,lev_sel ="VQ")

ax = plt.subplot(4,2,6);
ax.annotate('(g)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
plot_shade_plot(ax,lat,lev_new,VQ_AAs,VQ_1970,colormap,colorbar_min,colorbar_max,lev_sel ="VQ")

ax = plt.subplot(4,2,8);
ax.annotate('(h)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
colormesh1 = plot_shade_plot(ax,lat,lev_new,VQ_O3,VQ_1970,colormap,colorbar_min,colorbar_max,lev_sel ="VQ")
cbar_ax = fig.add_axes([0.54, 0.03, 0.44, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.1)*(colorbar_max-colorbar_min)+colorbar_min,2))
cbar_ax.annotate(r'$\mathrm{\mathsf{\Delta}\/g\/kg^{-1}\/m\/s^{-1}}$',xy=(0.5,1.2), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='center', va='baseline',rotation='horizontal',fontsize=15)	
				
plt.subplots_adjust(left=0.06, bottom=0.08, right=0.98, top=0.96, wspace=0.1, hspace=0.15); 
plt.savefig('./GHG_AEROSOL_OZONE/latitude_altitude_STF_VQ.png', format='png', dpi=1000)



"""
####T
colorbar_min=-5;colorbar_max=5;
ax = plt.subplot(4,3,1);
ax.annotate('2010-1970',xy=(-0.15,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)
ax.annotate(r'$\mathrm{\mathsf{\Delta} T\/(K)}$',xy=(0.5,1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('(a)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
plot_shade_plot(ax,lat,lev_new,T_total,T_1970,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(4,3,4);

ax.annotate('GHGs',xy=(-0.15,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)
ax.annotate('(b)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
plot_shade_plot(ax,lat,lev_new,T_GHG,T_1970,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(4,3,7);

ax.annotate('AAs',xy=(-0.15,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)
ax.annotate('(c)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
plot_shade_plot(ax,lat,lev_new,T_AAs,T_1970,colormap,colorbar_min,colorbar_max)

ax = plt.subplot(4,3,10);

ax.annotate('O3',xy=(-0.15,0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='vertical',fontsize=10)
ax.annotate('(d)',xy=(0.02,0.90), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)					
				
colormesh1 = plot_shade_plot(ax,lat,lev_new,T_O3,T_1970,colormap,colorbar_min,colorbar_max)
cbar_ax = fig.add_axes([0.09, 0.03, 0.28, 0.015])
char = fig.colorbar(colormesh1,orientation='horizontal',extend='both',cax=cbar_ax,ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2))

plt.subplots_adjust(left=0.08, bottom=0.10, right=0.98, top=0.95, wspace=0.1, hspace=0.15); 
plt.savefig('Zonal_mean_VQ.png', format='png', dpi=1000)
"""