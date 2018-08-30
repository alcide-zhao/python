
"""
This is to show the variation in  
		1. The surface temperature
		2. the sea level pressure with 850pa winds overlapped
		3. The 200hpa winds 
"""
import numpy as np
from scipy import stats
import math
import netCDF4 as nc4
import os; import site
from scipy.stats import ttest_ind as test
from scipy.interpolate import interp2d  as interp2d
from scipy import stats
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

import scipy.io as sio
oceanmask=sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan
def numerical_dif_2D(dependetnt,lon,lat,ax):
	if (ax == 1):
		lon_interval = 111.320*1000*np.cos(np.pi*lat/180)
		variable_gradient = np.gradient(lon)[1]*lon_interval
		# colormap='gray';ax=plt.imshow(variable_gradient,colormap);plt.show();plt.colorbar(ax)
		# print variable_gradient
		dependetnt_gradient = np.gradient(dependetnt)[1]
	elif (ax == 0):
		variable_gradient = np.gradient(lat)[0]*110.574*1000
		dependetnt_gradient = np.gradient(dependetnt)[0]
	deriviative = np.divide(dependetnt_gradient,variable_gradient)
	return deriviative
# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55.1],'EA':[100,145,20,50],'SA':[65,100,5,30]}
region = rergion_dic['ASIA'][:]

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_TS_PSL_Q_UV(Q)_LHFLX_FSNS(C).mat')
lon = data['lon'][0,:];lat = data['lat'][0,:];#lev = data['lev'][0,:];

# Q_rcp85 = np.nanmean(np.nanmean(data['Q850_rcp85'][0:20,:,:]*10**3;
UQ850 = np.nanmean(data['UQ850_rcp85'][0:20,:,:]*10**3,axis=0);
VQ850 = np.nanmean(data['VQ850_rcp85'][0:20,:,:]*10**3,axis=0);
LHFLX = np.nanmean(data['LHFLX_rcp85'][0:20,:,:]*0.03456,axis=0);
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_UQ_VQ.mat')
UQI =  np.nanmean(data['UQ_rcp85'][0:20,:,:,:]*10**3,axis=0);
VQI =  np.nanmean(data['VQ_rcp85'][0:20,:,:,:]*10**3,axis=0);

del data

def ps_readin():
	input_path ='/exports/csce/datastore/geos/users/s1667168/CESM/'
	var_path = input_path+'ensumble_mean_PS_200602_210101.nc'
	# print var_path
	nc_fid = nc4.Dataset(var_path,mode='r')
	time =nc_fid.variables['time']
	PS = nc_fid.variables['rcp85'][0:240,:,:]/100
	cache = np.empty((20,3,192,288))
	for iyear in range(2006,2026):
		layer_b = [layer for layer in range(len(time)) if time[layer]//100 == iyear*100+6][0]  #June01
		layer_e = [layer for layer in range(len(time)) if time[layer]//100  == iyear*100+9][0]  #August 31
		cache[iyear-2006,:,:,:] =PS[layer_b:layer_e,:,:]
	
	return np.nanmean(np.nanmean(cache,axis=1),axis=0)

def spatial_scaler_vector(axs,lons,lats,colormap,colorbar_min,colorbar_max,scaler,vector_u,vector_v,qk_scale=None,qk_caption=None,qk_is = False,tb_lef=True,tb_bot=True ):
	# calculate the origin of the map
	lon_0 = lons.mean();lat_0 = lats.mean()
	lon_b = np.min(lons); lon_e = np.max(lons)
	lat_b = np.min(lats); lat_e = np.max(lats)	
	lon_bin = 20; lat_bin = 15
	# lon_bins = np.arange(60.0,150.0,20.0); lat_bins = np.arange(-10.0,60.0,20.0)
	map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e,ax=axs)
	lon, lat = np.meshgrid(lons, lats)
	xi, yi = map(lon, lat)
	# plot contours.
	# set desired contour levels.
	# clevs = np.arange(np.nanmin(scaler[:]),np.nanmax(scaler[:]),(np.nanmax(scaler[:])-np.nanmin(scaler[:]))/20)
	cmap = discrete_cmap(20,colormap)
	cmap.set_bad([1,1,1],alpha = 1.0); #cmap.set_over([0,0,0]); cmap.set_under('b')
	masked_obj = np.ma.masked_where(np.isnan(scaler), scaler)
	norm = MidpointNormalize(midpoint=0,vmin = colorbar_min, vmax = colorbar_max)
	# CS1 = map.contourf(xi,yi,masked_obj,linewidths=0.5,colors='k',animated=True)
	# CS2 = map.contourf(xi,yi,masked_obj,cmap=cmap,vmin=colorbar_min,vmax=colorbar_max,norm=norm) #
	colormesh = map.pcolormesh(xi, yi, masked_obj,cmap=cmap,vmin=colorbar_min, vmax=colorbar_max,latlon=True)  #norm=norm,
	Q = map.quiver(xi,yi,vector_u,vector_v,scale=None)
	def plot_rectangle(bmap, lonmin,lonmax,latmin,latmax):
		xs = [lonmin,lonmax,lonmax,lonmin,lonmin]
		ys = [latmin,latmin,latmax,latmax,latmin]
		bmap.plot(xs, ys,latlon = True,linewidth = 1,color='k')
	# plot_rectangle(map, 100,145,20,55); x,y = map(102,48);plt.text(x,y,'East Asia',color='k')
	# plot_rectangle(map, 65,100,5,30); x,y = map(67,7);plt.text(x,y,'South Asia',color='k')
	
	# make quiver key.
	if qk_is:
		qk = plt.quiverkey(Q,1.11,0.05,qk_scale, qk_caption, labelpos='E') #10 m/s ,coordinates='data'   1.12,0.1 1.16,0.05  0.80,-0.20
		qk.text.set_backgroundcolor('w');qk.text.set_fontsize(10);
		qk.text.set_horizontalalignment('left')
	# Add Grid Lines
	map.drawcoastlines(); map.drawcountries()
	return colormesh	
PS =  ps_readin()
lons,lats,PS = range_clip(region[0],region[1],region[2],region[3],lon,lat,PS)
lons,lats,UQ_rcp85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UQI)
lons,lats,VQ_rcp85 = range_clip(region[0],region[1],region[2],region[3],lon,lat,VQI)
lons,lats,UQ850 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UQ850)
lons,lats,VQ850 = range_clip(region[0],region[1],region[2],region[3],lon,lat,VQ850)
# lons,lats,UQI8100 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UQI8100)
lons,lats,LHFLX = range_clip(region[0],region[1],region[2],region[3],lon,lat,LHFLX)



	
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
	print Pnew
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

Plwv = hybrid2p_interp(PS)
print PS
lev = np.array([3.64,7.59,14.36,24.61,39.27,54.60,72.01,87.82,103.31,121.55,142.89,168.21,297.91,232.83,273.91,322.24,379.10,445.99,524.68,609.77,691.39,763.40,820.86,859.53,887.02,912.64,936.20,957.49,976.33,992.56])
lev_gradient = np.gradient(lev,axis=0);  # print lev_gradient;
for x in range(65):
	for y in range(73):
		if PS[x,y]<925:
			dif = abs(lev-850)
			index = [layer for layer in range(len(lev)) if dif[layer] == np.min(dif) ][0]
			# Kprint index
			UQ_rcp85[index:30,x,y]=np.nan
			VQ_rcp85[index:30,x,y]=np.nan
		
p_g =10 # p over g =10 kg /m2
factor_c2mmd = 86.4*1.2      # 24*60*60*10/10000
lonm, latm = np.meshgrid(lons, lats)
QC_0615 = np.empty([30,65,73]);  #QC_8100 = np.empty([30,65,73]); 
for ilev in range(0,30):
	UQ_TMP = UQ_rcp85[ilev,:,:];VQ_TMP= VQ_rcp85[ilev,:,:];
	dvdlat = numerical_dif_2D(VQ_TMP,lonm, latm,ax=0); dudlon = numerical_dif_2D(UQ_TMP,lonm, latm,ax=1); 
	QC_0615[ilev-0,:,:]=(dvdlat+dudlon)*lev_gradient[ilev]*factor_c2mmd*p_g

QC_0615 =np.nansum(QC_0615,axis=0);

lon_interp_AM = np.arange(60,150,2);lat_interp_AM = np.arange(-5,55.1,2);
f = interp2d(lons,lats,QC_0615,kind='linear')
QC_0615 = f(lon_interp_AM, lat_interp_AM)
f = interp2d(lons,lats,UQ850,kind='linear')
UQ850 = f(lon_interp_AM, lat_interp_AM)
f = interp2d(lons,lats,VQ850,kind='linear')
VQ850 = f(lon_interp_AM, lat_interp_AM)



# NOW PLOTTING

fig = plt.figure(facecolor='White',figsize=(6.5, 2.4));plot_setup();	pad= 5 
colorbar_min=-5;  colorbar_max = 5 ;colormap ='RdBu_r'#;colormap= reverse_colourmap(colormap)
p_value=np.zeros((np.shape(LHFLX)))

ax = plt.subplot(1,2,1)

ax.annotate('Evap',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
colormesh1 = spatial_figure_norm(ax,LHFLX,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=False,tb_bot=False )
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.06, 0.07, 0.40, 0.03])
char = fig.colorbar(colormesh1,cax=cbar_ax,orientation ='horizontal',extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2)); # 
char.set_label(' mm/day')

colorbar_min= -10;  colorbar_max =10; 	
ax = plt.subplot(1,2,2)
qk_caption = r'$\mathrm{\mathsf{10\/g\/kg^{-1}\/m\/s^{-1}}}$'; qk_scale = 20
colormesh2 = spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,-QC_0615, UQ850, VQ850, qk_scale,qk_caption,qk_is = False, tb_lef=False,tb_bot=False)
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('Moisture',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
cbar_ax = fig.add_axes([0.56, 0.07, 0.40, 0.03])
char = fig.colorbar(colormesh2,orientation ='horizontal',cax=cbar_ax,extend='both',ticks=np.round(np.arange(-1,1.1,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2)); # 
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}}$')  #r'$\mathrm{\mathsf{ug\/kg^{-1}s^{-1}}}}$'


plt.subplots_adjust(left=0.05, bottom=0.09, right=0.96, top=0.95, wspace=0.2, hspace=0.07);
plt.savefig('Fig14_g.png', format='png', dpi=1200)
plt.show()
