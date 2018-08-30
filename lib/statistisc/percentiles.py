from numpy import linspace
from scipy.stats.kde import gaussian_kde
import numpy as np
# import matplotlib.pyplot as plt
def percentile(data,percent):
	if len(data) <=1:
		bins = np.nan; kde = np.nan;percentile=data;
	else:
		# data [data <=0]= np.nan; mask = ~np.isnan(data); data = data[mask]
		bins = linspace(min(data), max(data), 100); 
		gkde=gaussian_kde(data)
		kdepdf = gkde.evaluate(bins)
		kde = 100*kdepdf/np.sum(kdepdf)
		percentile = [bins[i] for i in range(0,99) if (sum(kde[0:i]) == percent or (sum(kde[0:i])<=percent and sum(kde[0:i+1])>=percent))][0]
		# kde_pdf = np.zeros((500,1))
		# print kde
		# print np.shape(kde_pdf)
		# for k in range(len(kde)):
			# kde_pdf[k,0] = sum(kde[0:k])
		# plt.plot(bins,kde_pdf); plt.show()
		# print np.mean(bins)
	return percentile,bins,kde