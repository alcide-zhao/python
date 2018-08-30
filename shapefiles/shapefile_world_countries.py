#   -- import --
import numpy as np
import shapefile
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

#   -- input --
import cartopy.io.shapereader as shpreader

shpfilename = shpreader.natural_earth(resolution='50m',
                                      category='cultural',
                                      name='admin_0_countries')								  
sf = shapefile.Reader(shpfilename)
recs    = sf.records()
shapes  = sf.shapes()
Nshp    = len(shapes)
print len(recs)
print len(shapes)

#   -- plot --
fig     = plt.figure()
ax      = fig.add_subplot(111)
cm= plt.cm.get_cmap('Dark2')
cccol = cm(1.*np.arange(Nshp)/Nshp)

for nshp in xrange(Nshp):
	ptchs   = []
	pts     = np.array(shapes[nshp].points)
	prt     = shapes[nshp].parts
	par     = list(prt) + [pts.shape[0]]
	for pij in xrange(len(prt)):
		ptchs.append(Polygon(pts[par[pij]:par[pij+1]]))
		ax.add_collection(PatchCollection(ptchs,facecolor=cccol[nshp,:],edgecolor='k', linewidths=.1))
print ptchs.get_closed()
ax.set_aspect(1)
ax.set_xlim(-180,+180)
ax.set_ylim(-90,90)
fig.savefig('test.png',dpi=1000)
