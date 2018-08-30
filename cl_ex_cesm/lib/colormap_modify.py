# discrete colormap into slices
import matplotlib.pyplot as plt
	
def discrete_cmap(N, base_cmap=None):
	"""Create an N-bin discrete colormap from the specified input map"""
	import numpy as np
	# Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
	base = plt.cm.get_cmap(base_cmap)
	color_list = base(np.linspace(0, 1, N)) #
	cmap_name = base.name + str(N)
	return base.from_list(cmap_name, color_list, N)
	
# reverse colorbar
import matplotlib as mpl
def reverse_colourmap(cmap, name = 'my_cmap_r'):      
	reverse = []
	k = []   
	base = plt.cm.get_cmap(cmap)
	for key in base._segmentdata:    
		k.append(key)
		channel = base._segmentdata[key]
		data = []
		
		for t in channel:                    
			data.append((1-t[0],t[2],t[1]))            
		reverse.append(sorted(data))
		
	LinearL = dict(zip(k,reverse))
	my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
	return my_cmap_r