# -*- coding: utf-8 -*-
"""
Created on Jan 25 2017
This is to show The temporal-spatial evelution of Aerosol, cloud, dynamic and thermodynamic related fields
This could include the aerosol emissions (BC,OC and SO2) from the emission Scenario (here refers to the RCP 85)
The model simulated :
	aerosol burden, AOD at 550 mm
	The Vertically-integrated droplet concentration CDNUMC
	The Average droplet effective radius
	The net clear/all sky shortwave flux
	
	THe wind, temperature, pressure, humidity, flux...
@author: Alcide.Zhao
"""
import scipy

import netCDF4 as nc4
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy import stats
import site
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
# from scipy.interpolate import interp2d

lib_path = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)
    ), 
    os.path.pardir, 
    # 'lib'
)

site.addsitedir(lib_path)

from lib import *


#######################################
# 0.1 varibale and function definition
#######################################
time_s = 2006
time_e = 2100
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],\
'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}

region = rergion_dic['ASIA'][:]


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name):
	'''
	This function is to process thw diagnostic field into JJA mean from the CESM
	monthly data
	The input is the file which contains both rcp85 and fixA simulations
	The output is should be normaly lons, lats, time_series and range clipped data
	the lev information willl also be returned for 4D field
	'''
	nc_fid = nc4.Dataset(file_name,mode='r')
	time = nc_fid.variables['time'][:]
	lat = nc_fid.variables['lat'][:]
	lon = nc_fid.variables['lon'][:]
	rcp85 = nc_fid.variables['rcp85'][:] 
	RCP85_MV = nc_fid.variables['rcp85'].missing_value
	rcp85[rcp85 == RCP85_MV] = np.nan
	rcp85_fixA = nc_fid.variables['rcp85_fixA'][:] 
	RCP85_fixA_MV = nc_fid.variables['rcp85_fixA'].missing_value
	rcp85_fixA[rcp85_fixA == RCP85_fixA_MV] = np.nan
	units =  nc_fid.variables['rcp85_fixA'].units
	# long_name =  nc_fid.variables['rcp85_fixA'].long_name
	long_name= 'Surface Temperature (radiative)'
	diff = rcp85-rcp85_fixA
	lons,lats,diff_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,diff)
	
	size = np.shape(diff_clipped)
	year= [ iyear/10000 for iyear in map(int,time)]
	year_series = [value for value in np.unique(year) if value <=2100]
	###JJA mean
	if (np.rank(diff_clipped) == 4):
		levs = nc_fid.variables['lev'][:]
		variable_JJAmean = np.empty((95,30,size[2],size[3])); variable_JJAmean[:] =np.nan
		
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = diff_clipped[layer_b:layer_e+1,:,:,:]
			variable_JJAmean[iyear-2006,:,:,:] = stats.nanmean(cache,axis=0)
	else:
		levs=np.nan
		variable_JJAmean = np.empty((95,size[1],size[2])); variable_JJAmean[:] =np.nan
		layer_e = -6  # CESMN MONTHLY DATA BEGINS FROM FEB OF 2006
		# print year_series
		for iyear in year_series:
			layer_b = layer_e + 10
			layer_e = layer_b + 2
			cache = diff_clipped[layer_b:layer_e+1,:,:]
			variable_JJAmean[iyear-2006,:,:] = stats.nanmean(cache,axis=0)
	nc_fid.close()
	return lons, lats, levs, year_series, variable_JJAmean,units,long_name


##############################################
# 850hpa moisture flux convergence :cnvergence and advection
##############################################
def numerical_dif_2D(dependetnt,variable,ax):
	if (ax == 1):
		variable_gradient = np.gradient(variable)[0]
		dependetnt_gradient = np.gradient(dependetnt)[1]
	elif (ax == 2 ):
		variable_gradient = np.gradient(variable)[1]
		dependetnt_gradient = np.gradient(dependetnt)[2]
	deriviative = np.divide(dependetnt_gradient,variable_gradient)
	return deriviative
	
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
variable = 'Q'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, scaler_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)


variable = 'U'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, u_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

variable = 'V'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
_, _, _, _, v_JJAmean,_,_= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

# meshgrid the lat and lon into regular dependent shape
lonm, latm = np.meshgrid(lons, lats)

dudlon = numerical_dif_2D(u_JJAmean,lonm,ax=2); dvdlat = numerical_dif_2D(v_JJAmean,latm,ax=1);
dqdlon = numerical_dif_2D(scaler_JJAmean,lonm,ax=2); dqdlat = numerical_dif_2D(scaler_JJAmean,latm,ax=1);
HMFCc = np.multiply(-1*scaler_JJAmean, (dudlon+dvdlat))*10**5
HMFCa = -1*(np.multiply(u_JJAmean,dqdlon)+np.multiply(v_JJAmean,dqdlat))*10**5
HMFC = HMFCc+ HMFCa


# Now plot
year_b = 2010 #beginning point of the time slice 
year_e = 2100
year_step = 10 # moving bins every tenyears
year_span = 20 #length of the time slice]


year_s = range(year_b,year_e-year_span+year_step,year_step)
figure_no = np.shape(year_s)[0]
spst = [3,int(figure_no)] # subplot_style define the rows and the columns of the s being plotted

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[20,10])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(HMFC[year-2006:year-2006+year_span,:,:],axis=0)
	
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='BrBG' #RdYlGn PRGn PiYG bwr seismic jet
	colorbar_min = -3
	colorbar_max = 3
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.65, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label('10^-5*kg/(s*kg)',fontsize=10)

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[20,10])	
	ax = plt.subplot(spst[0],spst[1],subplot_no+8)
	variable_JJAmean_Slice = stats.nanmean(HMFCc[year-2006:year-2006+year_span,:,:],axis=0)

	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='BrBG' #RdYlGn PRGn PiYG bwr seismic jet RdGy gist_rainbow
	colorbar_min = -2
	colorbar_max =  2
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.35, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label('10^-5*kg/(s*kg)',fontsize=10)

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[20,10])	
	ax = plt.subplot(spst[0],spst[1],subplot_no+2*8)
	variable_JJAmean_Slice = stats.nanmean(HMFCa[year-2006:year-2006+year_span,:,:],axis=0)

	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='BrBG' #RdYlGn PRGn PiYG bwr seismic jet
	colorbar_min = -2
	colorbar_max = 2
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)
	
fig.suptitle('Horizontal Moisture Convergence (convergence + advection) (RCP8.5-FixA)', fontsize="x-large")

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.05, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label('10^-5*kg/(s*kg)',fontsize=10)


###############################################
# CESM scaler FIELDS
###############################################

year_b = 2010 #beginning point of the time slice 
year_e = 2100
year_step = 10 # moving bins every tenyears
year_span = 20 #length of the time slice]

year_s = range(year_b,year_e-year_span+year_step,year_step)
figure_no = np.shape(year_s)[0]
spst = [2,int(figure_no/2)] # subplot_style define the rows and the columns of the s being plotted


variable = 'TS' # CDNUMC ADVIS AERL CLOUD
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, variable_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[15,5])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)
	
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='jet' #RdYlGn PRGn PiYG bwr seismic jet
	colorbar_min = None
	colorbar_max = None
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)
#Vertically-integrated droplet concentration 
#CLOUD FRACTION
#Average droplet effective radius
#Clear-sky net solar flux at surface#
#Net solar flux at top of atmosphere
#Surface latent heat flux
#Surface temperature
fig.suptitle('Surface temperature (RCP8.5-FixA)', fontsize="x-large")

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.1, 0.01, 0.8])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')

char.set_label(units,fontsize=10)


file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/ensumble_mean_PEI_global_2006_2100_rcp85.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
rcp85_meanPrecip = nc_fid.variables['mean_precip'][:] 
file_name = '/exports/csce/datastore/geos/users/s1667168/CESM/extremes_indices/prep/ensumble_mean_PEI_global_2006_2100_fixa.nc'
nc_fid = nc4.Dataset(file_name,mode='r')
fixa_meanPrecip = nc_fid.variables['mean_precip'][:] 
meanPrecip_diff = rcp85_meanPrecip -fixa_meanPrecip
lons,lats,diff_clipped = range_clip(region[0],region[1],region[2],region[3],lon,lat,meanPrecip_diff)
variable_JJAmean_ts = stats.nanmean(variable_JJAmean,axis = 0)
meanPrecip_diff_ts = stats.nanmean(meanPrecip_diff,axis = 0)
fig=plt.figure(figure_n, facecolor='White',figsize=[15,5])
ax = plt.subplot(111)
ax.plot(year_series,meanPrecip_diff_ts)
ax.plot(year_series,variable_JJAmean_ts)
plt.show()

###############################################
# CESM scaler FIELDS with vector fields overlapped
###############################################
# this function is specifically for plotting sclaler fields with vector fields embbed
def spatial_scaler_vector(lons,lats,subplot_title,colormap,colorbar_min,colorbar_max,scaler,vector_u,vector_v,projection='cyl'):
	# calculate the origin of the map
	lon_0 = lons.mean();lat_0 = lats.mean()
	
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = int((lon_e -lon_b)/3); lat_bin = int((lat_e -lat_b)/3)
	
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# plot SLP contours.
	# set desired contour levels.
	clevs = np.arange(np.nanmin(scaler[:]),np.nanmax(scaler[:]),(np.nanmax(scaler[:])-np.nanmin(scaler[:]))/20)
	cmap = discrete_cmap(20,colormap)
	# cmap= reverse_colourmap(cmap)
	cmap.set_bad('w',alpha = 1.0); cmap.set_over('k'); cmap.set_under('w')
	masked_obj = np.ma.masked_where(np.isnan(scaler), scaler)
	# norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max)
	CS1 = map.contourf(xi,yi,masked_obj,clevs,linewidths=0.5,colors='k',animated=True)
	CS2 = map.contourf(xi,yi,masked_obj,clevs,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max)
	# # plot wind vectors on projection grid.
	# # first, shift grid so it goes from -180 to 180 (instead of 0 to 360
	# # in longitude).  Otherwise, interpolation is messed up.
	# ugrid,newlons = shiftgrid(180.,vector_u,lons,start=False)
	# vgrid,newlons = shiftgrid(180.,vector_v,lons,start=False)
	# # transform vectors to projection grid.
	# uproj,vproj,xx,yy = \
	# map.transform_vector(ugrid,vgrid,newlons,lats,31,31,returnxy=True,masked=True)
	# now plot.
	# urot,vrot,x,y = map.rotate_vector(vector_u,vector_v,lons,lats],returnxy=True)
	Q = map.quiver(xi,yi,vector_u,vector_v,scale=None)
	# make quiver key.
	qk = plt.quiverkey(Q, 0.9, 1.1, 2, '2g/kg m/s', labelpos='W') #2 m/s
	
	# Add Grid Lines
	map.drawparallels(np.arange(round(lat_b,1), round(lat_e,1), lat_bin), labels=[1,0,0,0],linewidth=0.0, fontsize=10)
	map.drawmeridians(np.arange(round(lon_b,1), round(lon_e,1), lon_bin), labels=[0,0,0,1],linewidth=0.0, fontsize=10)
	map.drawcoastlines(); map.drawcountries()
	ax.set_title(subplot_title) 	# Add Title
	return CS2
	
	
year_b = 2010 #beginning point of the time slice 
year_e = 2100
year_step = 10 # moving bins every tenyears
year_span = 20 #length of the time slice]

year_s = range(year_b,year_e-year_span+year_step,year_step)
figure_no = np.shape(year_s)[0]
spst = [2,int(figure_no/2)] # subplot_style define the rows and the columns of the s being plotted



input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
# variable = 'T'
# file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
# lons, lats, levs, year_series, scaler,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
# scaler_JJAmean = stats.nanmean(scaler[:,12:18,:,:],axis=1) # 12-200hpa, 18-500hpa, 19-600, 23-850hpa, 29-1000
# scaler_JJAmean = scaler[:]
variable = 'UQ'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, u,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
u_JJAmean = u[:,23,:,:]*1000# 12-200hpa, 18-500hpa, 19-600, 23-850hpa, 29-1000
# u_JJAmean = stats.nanmean(u[:,12:18,:,:],axis=1)# 12-200hpa, 18-500hpa, 19-600, 23-850hpa, 29-1000
variable = 'VQ'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
_, _, _, _, v,_,_= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
v_JJAmean = v[:,23,:,:]*1000# 12-200hpa, 18-500hpa, 19-600, 23-850hpa, 29-1000
# v_JJAmean = stats.nanmean(v[:,12:18,:,:],axis=1)# 12-200hpa, 18-500hpa, 19-600, 23-850hpa, 29-1000
scaler_JJAmean = np.sqrt(np.power(u_JJAmean,2)+np.power(v_JJAmean,2))

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[50,20])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	scaler_JJAmean_Slice = stats.nanmean(scaler_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)
	u_JJAmean_Slice = stats.nanmean(u_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)
	v_JJAmean_Slice = stats.nanmean(v_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap = plt.cm.RdBu_r #RdYlGn PRGn PiYG bwr seismic jet PuOr RdBu_r
	colorbar_min = None
	colorbar_max = None
	colormesh = spatial_scaler_vector(lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max,scaler_JJAmean_Slice,u_JJAmean_Slice,v_JJAmean_Slice)

fig.suptitle('850hpa Water Transport(RCP8.5-FixA)', fontsize="x-large")
# Surface Pressure and 850hpa winds
# 500-200hpa averaged temperature and winds
# 850hpa winds and 1000-600 hpa integrated specific humidity
# 200hpa winds
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.1, 0.01, 0.8])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')

char.set_label(units,fontsize=10)



'''
#########################################
# CESM AEROSOL BURDEN
#########################################
variable = 'BURDENBC'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, variable_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
year_b = 2010 #beginning point of the time slice 
year_e = 2100
year_step = 10 # moving bins every tenyears
year_span = 20 #length of the time slice]

year_s = range(year_b,year_e-year_span+year_step,year_step)
figure_no = np.shape(year_s)[0]
spst = [3,figure_no] # subplot_style define the rows and the columns of the s being plotted
for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[15,5])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)	
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='PRGn'
	colorbar_min = -2.5
	colorbar_max = 0.5
	colormesh = spatial_figure(variable_JJAmean_Slice*10**6,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.65, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label('mg/m2',fontsize=10)

variable = 'BURDENSOA'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, SOA_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

variable = 'BURDENPOM'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, POM_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)
variable_JJAmean_Slice = SOA_JJAmean + POM_JJAmean

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+figure_no+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[30,10])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)	
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='PRGn'
	colorbar_min = -2.5
	colorbar_max = 0.5
	colormesh = spatial_figure(variable_JJAmean_Slice*10**6,lons,lats,subplot_title,colormap,colorbar_min,colorbar_max)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.35, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label('mg/m2',fontsize=10)


variable = 'BURDENSO4'
input_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Aerosol_cloud_dynamic_thermodynamic_hydro/'
file_name = input_path+'ensumble_mean_'+variable+'_200602_210101.nc'
lons, lats, levs, year_series, variable_JJAmean,units,long_name= Aerosol_cloud_dynamic_thermodynamic_hydro_prepro_JJAmean(file_name)

for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	figure_n = 1; subplot_no = year_slice+figure_no*2+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[30,10])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2006:year-2006+year_span,:,:],axis=0)	
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='PRGn'
	colorbar_min = -10
	colorbar_max = 10
	colormesh = spatial_figure(variable_JJAmean_Slice*10**6,lons,lats,subplot_title,colormap,colorbar_min,colorbar_max)

fig.suptitle('Aerosol Burden (RCP8.5-FixA)', fontsize="x-large")
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.05, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label('mg/m2',fontsize=10)


#########################################
# RCP EMISSIONS
#########################################

emission_dic = {'BC':np.zeros((96,141,181)),'OC':np.zeros((96,141,181)),'SO2':np.zeros((96,141,181))}

input_path = '/exports/csce/datastore/geos/users/s1667168/RCP/'
for key in emission_dic.keys():
	file_name ='accmip_interpolated_emissions_RCP85_'+key+'_2005_2100_0.5x0.5.nc'
	file = input_path+file_name
	nc_fid = nc4.Dataset(input_path+file_name,mode='r')
	time_rcp = nc_fid.variables['time'][:]
	lon_rcp = nc_fid.variables['lon'][:]
	lat_rcp = nc_fid.variables['lat'][:]
	ships = nc_fid.variables['ships'][:]
	units = nc_fid.variables['ships'].units
	anthropogenic = nc_fid.variables['anthropogenic'][:]
	biomass_burning = nc_fid.variables['biomass_burning'][:]
	value = ships+anthropogenic+biomass_burning
	lons,lats,value_clipped = range_clip(region[0],region[1],region[2],region[3],lon_rcp,lat_rcp,value)
	
	layer_e = -5   # RCP hae the record for 2005 and begins from Jan
	for iyear in range(2005,2101):
		layer_b = layer_e + 10
		layer_e = layer_b + 2
		cache = value_clipped[layer_b:layer_e+1,:,:]
		emission_dic[key][iyear-2005,:,:] = stats.nanmean(cache,axis=0)

year_b = 2010 #beginning point of the time slice 
year_e = 2100
year_step = 10 # moving bins every tenyears
year_span = 20 #length of the time slice]




year_s = range(year_b,year_e-year_span+year_step,year_step)
figure_no = np.shape(year_s)[0]
spst = [3,figure_no] # subplot_style define the rows and the columns of the s being plotted


variable_JJAmean = emission_dic['BC']*10**12
for year_slice in range(figure_no):
	year = year_s[year_slice]
	print year
	# spatial evelution under the RCP85 scenario
	figure_n = 2; subplot_no = year_slice+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[15,5])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2005:year-2005+year_span,:,:],axis=0)-	variable_JJAmean[0,:,:]
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='PRGn'
	colorbar_min = -20
	colorbar_max = 5
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.65, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label(units+'10^-12',fontsize=10)


variable_JJAmean = emission_dic['OC']*10**12
for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	subplot_no = year_slice+figure_no+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[15,5])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2005:year-2005+year_span,:,:],axis=0)-	variable_JJAmean[0,:,:]
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='PRGn'
	colorbar_min = -50
	colorbar_max = 12
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.35, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label(units+'*10^-12',fontsize=10)

variable_JJAmean = emission_dic['SO2']*10**12
for year_slice in range(figure_no):
	year = year_s[year_slice]
	
	# spatial evelution under the RCP85 scenario
	subplot_no = year_slice+figure_no*2+1
	fig=plt.figure(figure_n, facecolor='White',figsize=[15,5])
	ax = plt.subplot(spst[0],spst[1],subplot_no)
	variable_JJAmean_Slice = stats.nanmean(variable_JJAmean[year-2005:year-2005+year_span,:,:],axis=0)-	variable_JJAmean[0,:,:]
	subplot_title = str(year)+'_'+str(year+year_span)
	colormap ='PRGn'
	colorbar_min = -100
	colorbar_max = 60
	colormesh = spatial_figure(variable_JJAmean_Slice,lons,lats,subplot_title,colormap,colorbar_min ,colorbar_max)

fig.suptitle('Emissions (RCP8.5-2005)', fontsize="x-large")
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.95, 0.05, 0.01, 0.25])
char = fig.colorbar(colormesh, cax=cbar_ax,extend='both')
char.set_label(units+'*10^-12',fontsize=10)
'''
plt.show()
