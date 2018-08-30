
"""
This is to show the variation in  
		1. The surface temperature
		2. the sea level pressure with 850pa winds overlapped
		3. The 200hpa winds 
"""
import numpy as np
from scipy import stats
import math
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
oceanmask=sio.loadmat('//home/s1667168/coding/python/external_data/landoceanmask_CESM.mat')['landoceanmask']	
oceanmask[oceanmask==0]=np.nan

# region dictionary storing the resgion boundaey by lon_b.lon_e.lat_b and lat_e.
rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-5,55.1],'EA':[100,145,20,50],'SA':[65,100,5,30]}
region = rergion_dic['ASIA'][:]

# input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
# data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_T_4D.mat')

# T_rcp85 = data['T_rcp85'];T_fixa = data['T_fixa'];T = T_rcp85-T_fixa;

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_P.mat')
PS_rcp85 = data['PS_rcp85'];PS_fixa = data['PS_fixa']; 
del data

data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_TS_PSL_Q_UV(Q)_LHFLX_FSNS(C).mat')
lon = data['lon'][0,:];lat = data['lat'][0,:];#lev = data['lev'][0,:];
TS_rcp85 = data['TS_rcp85'];TS_fixa = data['TS_fixa'];TS = TS_rcp85-TS_fixa;
PSL_rcp85 = data['PSL_rcp85'];PSL_fixa = data['PSL_fixa']; 
U_rcp85 = data['U850_rcp85'];U_fixa = data['U850_fixa'];
V_rcp85 = data['V850_rcp85'];V_fixa = data['V850_fixa'];

Q_rcp85 = data['Q850_rcp85']*10**3;Q_fixa = data['Q850_fixa']*10**3; Q = Q_rcp85-Q_fixa;
UQ_rcp85 = data['UQ850_rcp85']*10**3;UQ_fixa = data['UQ850_fixa']*10**3;
VQ_rcp85 = data['VQ850_rcp85']*10**3;VQ_fixa = data['VQ850_fixa']*10**3;
del data

PSL = np.multiply(PSL_rcp85-PSL_fixa,1);
U = np.multiply(U_rcp85-U_fixa,1);
V = np.multiply(V_rcp85-V_fixa,1);
UQ = np.multiply(UQ_rcp85-UQ_fixa,1);
VQ = np.multiply(VQ_rcp85-VQ_fixa,1);

	
# p_mask is created to mask off the places where have SLP larger than reality. Thi is
#due to the vertical coordinate thet models use are hybrid p/p0
p_mask = stats.nanmean(PSL_rcp85-PS_rcp85,axis=0)/100;  
p_mask[p_mask<=300] =1 ;p_mask[p_mask>300] =0 
_,_,p_mask = range_clip(region[0],region[1],region[2],region[3],lon,lat,p_mask)
# print p_mask
# ax=plt.imshow(p_mask/100,origin='lower');plt.show();plt.colorbar(ax)

def spatial_diff_sig(region,lon,lat,variable_D):
	lons,lats,variable = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable_D)
	lon_interp_AM = np.arange(60,150,2);lat_interp_AM = np.arange(-5,55.1,2);
	f = interp2d(lons,lats,p_mask,kind='linear')
	p_ma = f(lon_interp_AM, lat_interp_AM); p_ma[p_ma==0] =np.nan
	variable_m =stats.nanmean(variable[0:10,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_0615_D = np.multiply(f(lon_interp_AM, lat_interp_AM),p_ma)
	variable_m = stats.nanmean(variable[25:45,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_3150_D = np.multiply(f(lon_interp_AM, lat_interp_AM),p_ma)
	variable_m = stats.nanmean(variable[75:95,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_8100_D = np.multiply(f(lon_interp_AM, lat_interp_AM),p_ma)
	return lon_interp_AM,lat_interp_AM,att_0615_D,att_3150_D,att_8100_D

lons,lats,PSL_0615_D,PSL_3150_D,PSL_8100_D = spatial_diff_sig(region,lon,lat,PSL)

# _,_,Q_0615_D,Q_3150_D,Q_8100_D = spatial_diff_sig(region,lon,lat,QC)
_,_,U_0615_D,U_3150_D,U_8100_D = spatial_diff_sig(region,lon,lat,U)
_,_,V_0615_D,V_3150_D,V_8100_D = spatial_diff_sig(region,lon,lat,V)
_,_,UQ_0615_D,UQ_3150_D,UQ_8100_D = spatial_diff_sig(region,lon,lat,UQ)
_,_,VQ_0615_D,VQ_3150_D,VQ_8100_D = spatial_diff_sig(region,lon,lat,VQ)

def spatial_diff_sig(region,lon,lat,variable_D):
	lons,lats,variable = range_clip(region[0],region[1],region[2],region[3],lon,lat,variable_D)
	lon_interp_AM = np.arange(60,150,2);lat_interp_AM = np.arange(-5,55.1,2);
	variable_m =stats.nanmean(variable[0:10,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_0615_D = np.multiply(f(lon_interp_AM, lat_interp_AM),1)
	variable_m = stats.nanmean(variable[25:45,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_3150_D = np.multiply(f(lon_interp_AM, lat_interp_AM),1)
	variable_m = stats.nanmean(variable[75:95,:,:],axis=0)
	f = interp2d(lons,lats,variable_m,kind='linear')
	att_8100_D = np.multiply(f(lon_interp_AM, lat_interp_AM),1)
	return lon_interp_AM,lat_interp_AM,att_0615_D,att_3150_D,att_8100_D

_,_,TS_0615_D,TS_3150_D,TS_8100_D = spatial_diff_sig(region,lon,lat,TS)


##########################################
# VERTICALLY INTEGRATED WATER CONVERGENCE
##########################################
# deriviative
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
	
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_UQ_VQ.mat')
UQI_rcp85 = data['UQ_rcp85']*10**3;UQI_fixa = data['UQ_fixa']*10**3;UQI = UQI_rcp85-UQI_fixa; 
VQI_rcp85 = data['VQ_rcp85']*10**3;VQI_fixa = data['VQ_fixa']*10**3;VQI = VQI_rcp85-VQI_fixa;
del data


UQI3150 = stats.nanmean(UQI[25:45,:,:,:],axis=0);UQI8100 = stats.nanmean(UQI[75:95,:,:,:],axis=0)
VQI3150 = stats.nanmean(VQI[25:45,:,:,:],axis=0);VQI8100 = stats.nanmean(VQI[75:95,:,:,:],axis=0)

lons,lats,UQI3150 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UQI3150)
lons,lats,VQI3150 = range_clip(region[0],region[1],region[2],region[3],lon,lat,VQI3150)
lons,lats,UQI8100 = range_clip(region[0],region[1],region[2],region[3],lon,lat,UQI8100)
lons,lats,VQI8100 = range_clip(region[0],region[1],region[2],region[3],lon,lat,VQI8100)

lev = np.array([3.64,7.59,14.36,24.61,39.27,54.60,72.01,87.82,103.31,121.55,142.89,168.23,197.91,232.83,273.91,322.24,379.10,445.99,524.68,609.77,691.39,763.40,820.86,859.53,887.02,912.64,936.20,957.49,976.33,992.56])
lev_gradient = np.gradient(lev);  # print lev_gradient;
p_g =10 # p over g =10 kg /m2
factor_c2mmd = 86.4*1.2      # 24*60*60*10/10000
lonm, latm = np.meshgrid(lons, lats)
QC_3150 = np.empty([30,65,73]); QC_8100 = np.empty([30,65,73]); 
for ilev in range(0,30):
	#2031_2050
	UQ_TMP = UQI3150[ilev,:,:];VQ_TMP= VQI3150[ilev,:,:];
	dvdlat = numerical_dif_2D(VQ_TMP,lonm, latm,ax=0); dudlon = numerical_dif_2D(UQ_TMP,lonm, latm,ax=1); 
	QC_3150[ilev-0,:,:]=(dvdlat+dudlon)*lev_gradient[ilev]*factor_c2mmd*p_g
	##2081_2100
	UQ_TMP = UQI8100[ilev,:,:];VQ_TMP = VQI8100[ilev,:,:];
	dvdlat = numerical_dif_2D(VQ_TMP,lonm, latm,ax=0); dudlon = numerical_dif_2D(UQ_TMP,lonm, latm,ax=1); 
	QC_8100[ilev-0,:,:]=(dvdlat+dudlon)*lev_gradient[ilev]*factor_c2mmd*p_g

QC_3150 =np.nansum(QC_3150,axis=0)	;QC_8100 =np.nansum(QC_8100,axis=0)

lon_interp_AM = np.arange(60,150,2);lat_interp_AM = np.arange(-5,55.1,2);
f = interp2d(lons,lats,QC_3150,kind='linear')
QC_3150 = f(lon_interp_AM, lat_interp_AM)
f = interp2d(lons,lats,QC_8100,kind='linear')
QC_8100 = f(lon_interp_AM, lat_interp_AM)



# NOW PLOTTING

fig = plt.figure(facecolor='White',figsize=(6.5, 5.1));plot_setup();	pad= 5 
colorbar_min=-1.5;  colorbar_max =1.5; colormap ='RdYlBu';colormap= reverse_colourmap(colormap);p_value=np.zeros((np.shape(TS_0615_D)))

# ax = plt.subplot(3,2,1)
# spatial_figure_norm(ax,TS_0615_D,lons,lats,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=False )
# ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
# ax.annotate('2006-2015',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')

ax = plt.subplot(3,2,1)
spatial_figure_norm(ax,TS_3150_D,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=False )
ax.annotate('SAT',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(a)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax.annotate('2031-2050',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
ax = plt.subplot(3,2,2)
ax.annotate('2081-2100',xy=(0.5, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center')
colormesh1 = spatial_figure_norm(ax,TS_8100_D,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,p_value, tb_lef=False,tb_bot=False )
ax.annotate('(b)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.86, 0.67, 0.01, 0.26])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.2,0.2)*(colorbar_max-colorbar_min)+colorbar_min,2)); # 
char.set_label('K')

# PSL+winds 850
colorbar_min=-50.1;  colorbar_max =50; 	
# ax = plt.subplot(3,2,4)
# spatial_scaler_vector(ax,lons,lats,colormap,colorbar_min,colorbar_max,PSL_0615_D, U_0615_D, V_0615_D, qk_scale=None,qk_caption=None,qk_is = False, tb_lef=True,tb_bot=False )
# ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)
# ax.annotate('2006-2015',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')

ax = plt.subplot(3,2,3)
spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,PSL_3150_D, U_3150_D, V_3150_D,qk_scale=None,qk_caption=None,qk_is = False, tb_lef=True,tb_bot=False )
ax.annotate('SLP + Winds850',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(c)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
# ax.annotate('2031-2050',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')
ax = plt.subplot(3,2,4)
# ax.annotate('2081-2100',xy=(0.5, 1.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # size=10, ha='center', va='center')
qk_caption = r'$\mathrm{\mathsf{1\/m\/s^{-1}}}$'; qk_scale =1
colormesh1 = spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,PSL_8100_D, U_8100_D, V_8100_D,qk_scale,qk_caption,qk_is = True, tb_lef=False,tb_bot=False )
ax.annotate('(d)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
cbar_ax = fig.add_axes([0.86, 0.39, 0.01, 0.24])
char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(-50,50.1,20))); # 
char.set_label('Pa')

# 850Q+MOISTURE 850
colorbar_min=-1.01;  colorbar_max =1.01; 	
# ax = plt.subplot(3,2,7)
# spatial_scaler_vector(ax,lons,lats,colormap,colorbar_min,colorbar_max,Q_0615_D, UQ_0615_D, VQ_0615_D,qk_scale=None,qk_caption=None,qk_is = False, tb_lef=True,tb_bot=True )
# ax.annotate('(g)',xy=(0.02, 0.01), xytext=(0, pad),
                # xycoords='axes fraction', textcoords='offset points',
                # ha='left', va='baseline',rotation='horizontal',fontsize=10)

ax = plt.subplot(3,2,5)
QC_3150[QC_3150>=0.1] =QC_3150[QC_3150>=0.1]+0.1
spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,-QC_3150, UQ_3150_D, VQ_3150_D,qk_scale=None,qk_caption=None,qk_is = False, tb_lef=True,tb_bot=True )
ax.annotate('Moisture',xy=(-0.2, 0.5), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)
ax.annotate('(e)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)
ax = plt.subplot(3,2,6)
QC_8100[QC_8100>=0] =QC_8100[QC_8100>=0]+0.2
qk_caption = r'$\mathrm{\mathsf{10\/g\/kg^{-1}\/m\/s^{-1}}}$'; qk_scale = 20
colormesh2 = spatial_scaler_vector(ax,lon_interp_AM,lat_interp_AM,colormap,colorbar_min,colorbar_max,-QC_8100, UQ_8100_D, VQ_8100_D, qk_scale,qk_caption,qk_is = True, tb_lef=False,tb_bot=True)
ax.annotate('(f)',xy=(0.02, 0.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='baseline',rotation='horizontal',fontsize=10)

cbar_ax = fig.add_axes([0.86, 0.09, 0.01, 0.24])
char = fig.colorbar(colormesh2,cax=cbar_ax,extend='both',ticks=np.arange(-1,1.1,0.4)); # 
char.set_label(r'$\mathrm{\mathsf{mm\/day^{-1}}}}$')  #r'$\mathrm{\mathsf{ug\/kg^{-1}s^{-1}}}}$'


plt.subplots_adjust(left=0.1, bottom=0.05, right=0.81, top=0.95, wspace=0.000000000001, hspace=0.07);
plt.savefig('Fig14.png', format='png', dpi=1200)
plt.show()
