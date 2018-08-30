# -*- coding: utf-8 -*-
"""
Created by Alcide Zhao in May 10 2017, at UoE
This is to produce mask maps for the countries you are interested
For usig this script, all what you need to do is to add/remove 
the country you would like to/from the country_dic in line 113
The countries which can be now produced by this script include China, India and Australia
for the exacrt long name of the countries, refer https://unstats.un.org/unsd/methodology/m49/
"""
import numpy as np
# import numpy.concentrate as concentrate
from cartopy.io.shapereader import natural_earth, Reader

def get_geometries(country_name):
	"""
	Get an iterable of Shapely geometries corrresponding to given country name.

	"""
	# Using the Natural Earth feature interface provided by cartopy 
	# see www.naturalearthdata.com/features/ for the name of the shapefiles. e.g. admin_0_countries as below 
	# You could use a different source, all you need is the geometries.
	shape_records = Reader(natural_earth(resolution='110m', # note you can choose 50m or 10m instead of 110m here
										 category='Cultural',
										 name='admin_0_countries')).records()
	geoms = []
	for country in shape_records:
		if country.attributes['name_long'] in country_name:
			try:
				geoms += country.geometry.boundary
			except TypeError:
				geoms.append(country.geometry.boundary)
	lat = np.array([]);lon = np.array([])
	for patch in geoms:
		x, y = patch.xy
		lon = np.append(lon,x)
		lat = np.append(lat,y)
	# the administrative boundary of the country in stored as longitude-lattitude coordinates
	boundary_cords = np.transpose(np.array([lon,lat]))
	return boundary_cords


def outline_to_mask(line, x, y):
	"""Create mask from outline contour

	# Parameters
	# ----------
	# line: array-like (N, 2)
	# x, y: 1-D grid coordinates (input for meshgrid)

	# Returns
	# -------
	# mask : 2-D 0 and 1 array (1 for land)
	"""
	import matplotlib.path as mplp
	mpath = mplp.Path(line)
	X, Y = np.meshgrid(x, y)
	points = np.array((X.flatten(), Y.flatten())).T
	mask = (mpath.contains_points(points).reshape(X.shape))*1
	return mask


def countries_mask(country_dic):
	"""
	The mask of each country is stored seperately in the dictionary, 
	while the mask of all the countries is in the  variable mask
	NB because fome counties have disconnected areas, but by default they will be recognized as one patch. As a consequence.
		there would be straight lines connecting the two (or more) patches. To get rid of this, 5deg of distance (This can be 
		tweasted, but I found 5 is the best one ) is used as a threshold to identify those closed patches for one country, 
		and then the boundary coordinates for each closed patch is detected seperately.
	"""
	lon= np.arange(-180,180,0.5); lat = np.arange(-90,90,0.5)
	for country in country_dic.keys():
		print country
		country_mask = np.zeros([360,720])
		boundary_cords= get_geometries(country)#[0:280,:]
		if (np.shape(boundary_cords)[0] <=2):
			country_mask =  outline_to_mask(boundary_cords, lon, lat)
			mask = np.zeros([360,720]); mask[:,0:360] = country_mask[:,360:720];mask[:,360:720] = country_mask[:,0:360]
			country_dic[country] = mask
		else:
			x = boundary_cords[:,0]; y = boundary_cords[:,1];
			# calculate the distance between consecutive boundary coordinate
			x_gradient = np.gradient(x);y_gradient = np.gradient(y); 
			distance = np.sqrt(x_gradient*x_gradient +y_gradient*y_gradient); 
			# identify disconnected and closed patched 
			g_boundary = [i for i in range(len(distance)) if distance[i] >=5];
			g_b_patch =np.empty([len(g_boundary)+2]); # g_b_patch stores the staring and ending coordinate of each patch
			g_b_patch[0]=0;g_b_patch[-1]=len(distance);g_b_patch[1:-1]=g_boundary;
			for k_patch in range(0,len(g_b_patch)-1):
				print k_patch
				patch_start = int(g_b_patch[k_patch]);patch_end = int(g_b_patch[k_patch+1])
				print patch_start,patch_end
				# detect boundary coordinate for each patch
				mask = outline_to_mask(boundary_cords[patch_start:patch_end,:], lon, lat)
				# add all patches up for the country, Russia is the most annoying one, holy
				country_mask=country_mask+mask
			# convert to 0-360 -90-90
			mask = np.zeros([360,720]); mask[:,0:360] = country_mask[:,360:720];mask[:,360:720] = country_mask[:,0:360]
			country_dic[country] = mask
	country_dic.update(all=sum(country_dic.values()))
	lon = np.arange(0,360,0.5);
	return lon,lat,country_dic

#############
#  Main     #
#############
# country_dic = {'China':[],'India':[],'Australia':[]}


#european counties
# country_dic = {'Albania':[],'Andorra':[],'Armenia':[],'Austria':[],'Azerbaijan':[],'Belarus':[],'Belgium':[],'Bosnia and Herzegovina':[],'Bulgaria':[],\
# 'Croatia':[],'Cyprus':[],'Czech Republic':[],'Denmark':[],'Estonia':[],'Finland':[],'Georgia':[],'Germany':[],'Ireland':[],\
# 'Greece':[],'Hungary':[],'Italy':[],'Sweden':[],'Kosovo':[],'Latvia':[],'Liechtenstein':[],\
# 'Lithuania':[],'Luxembourg':[],'The former Yugoslav Republic of Macedonia':[],'Malta':[],'Moldova':[],'Monaco':[],'Montenegro':[],'Netherlands':[],\
# 'Poland':[],'Portugal':[],'Romania':[],'Russia':[],'San Marino':[],'Serbia':[],'Slovakia':[],'Spain':[],'Switzerland':[],\
# 'Turkey':[],'Ukraine':[],'United Kingdom':[],'Slovenia':[],'French Guiana':[],'Russian Federation':[]}

## souther African countries
country_dic = {'Botswana':[],'Lesotho':[],'Zimbabwe':[],'Angola':[],'Malawi':[],'Mozambique':[],'South Africa':[],'Swaziland':[],'Zambia':[],'Namibia':[]}  #



lon,lat,country_dic = countries_mask(country_dic)



import scipy.io as sio
data_path = './'
sio.savemat(data_path+'Southern_Afirca_720_360.mat',{'lon':lon,'lat':lat,\
'SouthernAfrica':country_dic['all']})    #,'China':country_dic['China'],'India':country_dic['India'],'Australia':country_dic['Australia']})

# read out 
data = sio.loadmat(data_path+'Southern_Afirca_720_360.mat')
lon = data['lon'][0,:];lat = data['lat'][0,:];STA = data['SouthernAfrica'][:]
import matplotlib.pyplot as plt
plt.imshow(STA,origin='lower');plt.show()


