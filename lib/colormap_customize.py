from matplotlib.pyplot import *
import matplotlib
import numpy as np

colors = [(0,178/255.0,238/255.0),(152/255.0,245/255.0,255/255.0),
		 (248/255.0,248/255.0,255/255.0), (255/255.0,246/255.0,143/255.0),
		 (255/255.0,255/255.0,0), (255/255.0,128/255.0,0),
		 (255/255.0,0,0), (139/255.0,37/255.0,0),
		 ]  
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

pcolor(np.random.rand(10,10),cmap=my_cmap)
colorbar()
show()