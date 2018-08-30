import matplotlib.pyplot as plt
import numpy as np

import os; import site
lib_path = os.path.join(
	os.path.realpath(
       os.path.dirname(__file__)
	), 
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

import scipy.io as sio

ocean_mask_CESM = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_360_720.mat')['landoceanmask']
ocean_mask_CESM[ocean_mask_CESM==0]=np.nan;#ocean_mask_CESM[0:27,:]=np.nan
print np.shape(ocean_mask_CESM)
ocean_mask_CESM[0:221,:] =0;ocean_mask_CESM[234:,:] =0
ocean_mask_CESM[:,0:241] =0;ocean_mask_CESM[:,248:] =0
Taiwan =ocean_mask_CESM
data_path ='/home/s1667168/coding/python/shapefiles/'

fig = plt.figure(facecolor='White',figsize=(7, 5));plot_setup();pad= 5 
ax = plt.subplot(1,1,1);colorbar_min=1;colorbar_max=8;cmap ='jet';


# data_map = sio.loadmat(data_path+'Euro_USA_AUS_BRICS_720_360.mat')

# lon = data_map['lon'][0,:];lat = data_map['lat'][0,:];
# Europe = data_map['Europe'][:]*1;USA = data_map['USA'][:]*2
# Australia = data_map['Australia'][:]*3;India = data_map['India'][:]*4;
# Brazil = data_map['Brazil'][:]*5;Russia = data_map['Russia'][:]*6;
# SA = data_map['South Africa'][:]*7;China = data_map['China'][:]*8;

# data_map = sio.loadmat(data_path+'Russia_720_360.mat')
# Russia = data_map['Russia'][:]
# data = Russia
# data=Europe+USA+Australia+India+Brazil+Russia+SA+China

# p_value=np.empty((np.size(data)));p_value[:]=np.nan
# colormesh1=spatial_figure(ax,data,lon,lat,cmap,colorbar_min,colorbar_max,p_value, tb_lef=True,tb_bot=True );
# cbar_ax = fig.add_axes([0.95, 0.1, 0.01, 0.8])
# char = fig.colorbar(colormesh1,cax=cbar_ax,extend='both',ticks=np.round(np.arange(0,1.25,0.125)*(colorbar_max-colorbar_min)+colorbar_min,2)); #
# plt.show()


data_map = sio.loadmat(data_path+'Southern_Afirca_720_360.mat')
Southern_Africa = data_map['SouthernAfrica'][:]

data_map = sio.loadmat(data_path+'Euro_StAf_USA_AUS_BRICS_720_360.mat')
lon = data_map['lon'][0,:];lat = data_map['lat'][0,:];
India = data_map['India'][:] ;Brazil = data_map['Brazil'][:];
SA = data_map['South Africa'][:];China = data_map['China'][:]+Taiwan;
Russia = data_map['Russia'][:];Australia = data_map['Australia'][:]
Europe = data_map['Europe'][:];USA = data_map['USA'][:]

Gloland = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_360_720.mat')['landoceanmask'];

All = sio.loadmat('/home/s1667168/coding/python/external_data/landoceanmask_360_720.mat')['landoceanmask'];All[:]=1
#globe[0:48,:]=0


import scipy.io as sio
data_path = './'
sio.savemat(data_path+'Euro_USA_AUS_BRICS_STA_720_360.mat',{'lon':lon,'lat':lat,\
'All':All,'GloLand':Gloland,'Europe':Europe,'USA':USA,'India':India,'Brazil':Brazil,'China':China,'Russia':Russia,'South Africa':SA,'Australia':Australia,'Southern African':Southern_Africa})   

