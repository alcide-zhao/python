from scipy.integrate import quad
from math import pi, exp
import numpy as np
import matplotlib.pyplot as plt
mean = (np.log(3.65)-np.log(0.8))/2+np.log(0.8); 
sd = np.log(1.8); 
# ll = np.log(0)


def get_accumulated_pb(ul):
	loop = np.arange(0,ul,0.0001); probability = np.empty((len(loop))); volume = np.empty((len(loop)));
	for i in range(len(loop)-1):
		ll = np.log(loop[i])
		ul = np.log(loop[i+1])
		probability[i],_= quad(lambda x: 1 / ( sd * ( 2 * pi ) ** 0.5 ) * exp( (x-mean) ** 2 / (-2 * sd ** 2) ), ll, ul)
		volume[i] = (loop[i]+0.005)**3
	plt.plot(loop,probability);plt.show()
	return np.sum(np.multiply(volume,probability))
	

total = get_accumulated_pb(15);print total;
pm = get_accumulated_pb(2.5);print pm;print pm/total*100