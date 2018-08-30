import scipy.io as spio
import matplotlib.pyplot as plt
import numpy as np
from lib import *
# import matplotlib.patches as mpatch
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
import pickle

file_path = '//home/s1667168/scratch/extremes_indices/'


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=10):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


oceanmask=spio.loadmat('//home/s1667168/coding/python/climate_extremes_cesm/external_data/world.oceanmask.1440x720_0_360_-90_90.mat')['landoceanmask']
oceanmask[oceanmask==255]=1
lats = np.arange(-90,90.25,0.25)
lons = np.arange(0,360.25,0.25)
# region_dic = {'SWA':[50.,75.,24.,32.],'SI':[75.,82.,5.,15.],'CI':[72.,90.,15.,22.],\
# 'CEA':[100.,122.,34.,40.],'SNA':[70.,115.,45.,55.],'SASEA':[70.,120.,10.,30.],}

asia_dic={'GLOBE':[0,360,-90,90],'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'SAH':[0,65,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}
# 'EA':[100.,145.,20.,54.],'SA':[65.,100.,0.,30.]
value_dic={'GLOBE':1,'SEA':1,'CA1':3,'CA2':3,'NA':4,'TBT':3,\
'EAF':3,'SAH':4,'WAF':2,'SAF':4,'NEU':2,'MED':1,'AUS':4,\
'SI':1,'CI':2,'WSA':4,'CC':1,'EC':4,'NWC':2,'SC':2,'NC':3,'NEC':2,'NI':1}

# value_dic = {'GLOBE':1.5,'SEA':1,'EA':2,'SA':3,'CA':4,'NA':5,'TBT':6,'EAF':7,'SAH':8,'WAF':9,'SC':10,'NEU':11,'NC':12,'AUS':13}

fig = plt.figure(1, facecolor='White',figsize=[9,7.5])
# ax = fig.add_subplot(111)

# clip the region out and asssign them with different values to discriminate
map_subregion = np.empty((720,1440))
map_subregion[:] = 0.
for region in asia_dic:
	row_s = [row for row in range(len(lats)) if lats[row] == asia_dic[region][2]]
	row_s = row_s[0]
	row_e = [row for row in range(len(lats)) if lats[row] == asia_dic[region][3]]
	row_e = row_e[0]
	col_s = [col for col in range(len(lons)) if lons[col] == asia_dic[region][0]]
	col_s = col_s[0]
	col_e = [col for col in range(len(lons)) if lons[col] == asia_dic[region][1]]
	col_e = col_e[0]
	map_temp = np.empty((720,1440))
	map_temp[:] = 0.
	map_temp[row_s:row_e,col_s:col_e] = value_dic[region]
	map_subregion = map_subregion + map_temp

	# ax.add_patch(mpatch.Rectangle((region_dic[region][0],region_dic[region][2]),region_dic[region][1]-region_dic[region][0],\
	# region_dic[region][3]-region_dic[region][2],edgecolor="k",linewidth=0,fill=True,facecolor="k",alpha = alpha_value))
	# alpha_value = alpha_value +0.05
# print map_subregion
map_subregion = np.multiply(map_subregion,oceanmask)	
lons,lats,map_subregion=range_clip(0,180,-60,90,lons,lats,map_subregion)	

# define the whole region and add the basemap	
lon_b = 0; lon_e = 180; lat_b = -60;lat_e = 90
plt.xlim([lon_b, lon_e])
plt.ylim([lat_b, lat_e])
lon_0 = lons.mean()
lat_0 = lats.mean()
map = Basemap(lat_0=lat_0, lon_0=lon_0,llcrnrlon=lon_b,llcrnrlat=lat_b,urcrnrlon=lon_e,urcrnrlat=lat_e)	
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)
cmap = plt.get_cmap('Greys')
# selext a subset of the colormap, from 0 to 1 for example, and divide the colormap into five slices
# set bad values as color white; set over values as color black ans under as white 
new_cmap = truncate_colormap(cmap, 0.1,1,5)
new_cmap.set_bad('w',alpha = 1)
new_cmap.set_over('k')
new_cmap.set_under('w')
s = map.pcolor(xi, yi, map_subregion,cmap=new_cmap,vmin=1.5, vmax=5,alpha = 0.5)
map.drawparallels(np.arange(lat_b, lat_e, 30.), labels=[1,0,0,0],linewidth=0.5, fontsize=20)
map.drawmeridians(np.arange(lon_b, lon_e, 30.), labels=[0,0,0,1],linewidth=0.5, fontsize=20)
map.drawcoastlines(linewidth=2)
map.drawcountries(linewidth=1)


# add text indicating the subregion name
region_dic={'SEA':[94.,155.,-11.,20.],'CA':[40.,74.,38.,53.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
'SAH':[0,65,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],\
'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}
for region in region_dic:
	# patch= plt.add_patch(patches.Rectangle((value_dic[region][0], value_dic[region][2]), value_dic[region][1]-value_dic[region][0], value_dic[region][3]-value_dic[region][2],
    # alpha=0,facecolor='none',label='Label',linewidth = 0))
	centerx =region_dic[region][0] + (region_dic[region][1]-region_dic[region][0])/2
	centery =region_dic[region][2] + (region_dic[region][3]-region_dic[region][2])/2
	plt.annotate(region, xy=(centerx, centery),ha="center",fontsize=18,color = 'r')

"""
# pickle a figure to disk
fig_name = file_path +'subregion_asia.pickle'	
pickle.dump(fig,file(fig_name,'w'))
fig_object = pickle.load(open(fig_name,'rb'))
fig_object.show()
'''

# rergion_dic={'GLOBE':[0,360,-90,90],'ASIA':[60,150,-15,55],'EA':[100,145,20,50],'SA':[65,100,5,30],\
# 'SEA':[94.,155.,-11.,20.],'CA1':[40.,74.,38.,53.],'CA2':[40.,65.,30.,38.],'NA':[40.,180.,53.,70.],'WSA':[90.,100.,20.,28.],\
# 'ESAH':[40,65,18,30],'WSAH':[0,40,18,30],'EAF':[22,60,-12,18],'WAF':[0.,22.,-12.,18.],'SAF':[0.,52.,-35.,-12.],\
# 'NEU':[0.,40.,48,75],'MED':[0.,40.,30.,48.],'AUS':[110.,155.,-45.,-11.],'SI':[73.,82.,6.,17.],'CI':[65.,90.,17.,28.],'NI':[65.,79.,28.,38.],'SCI':[65.,90.,5.,28.],\
# 'TBT':[79.,105.,28.,38.],'CC':[105.,115.,28.,38.],'EC':[115.,122.,28.,38.],'NWC':[74.,110.,38.,53.],'SC':[100.,124.,20.,28.],'NC':[110.,122.,38.,53.],'NEC':[122.,144.,30.,53.]}



