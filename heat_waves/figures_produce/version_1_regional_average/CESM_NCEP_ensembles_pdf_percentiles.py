"""
This is to plot the pdfs of CESM ensembles and NCEP Tn and Tm 
and also the 90 and 95 percentiles
The CESM ensemble pdfs are savedinto the "TNX_BINS_PDF_1971_2000_JJAS_CESM30.mat"
	0 dimenssion: ensemble_no
	1 dimenssion: Tn bins, Tn pdfs, Tm bins and Tm pdfs from 0 to 3
	2 dimenssion: length of 100 corresponding to the bins and pdfs
The NCEP pdfs are saved into the "TNX_BINS_PDF_1971_2000_JJAS_NCEP.mat"
	0 dimenssion: Tn bins, Tn pdfs, Tm bins and Tm pdfs from 0 to 3
	1 dimenssion: length of 100 corresponding to the bins and pdfs
Region keys are respectively the
	SEC_bins_pdfs:  'SouthEastChina':[105,125,20,30]
	CEC_bins_pdfs:  'CentralEastChina':[105,125,30,40]
	SID_bins_pdfs	'SouthIndia':[70,90,8,22]
	NID_bins_pdfs:  'NorthIndia':[70,90,22,33]
"""
import numpy as np
from scipy import stats
import scipy.io as sio

import os; import site
lib_path = os.path.join(
	os.path.realpath(
        os.path.dirname(__file__)
	), 
	os.path.pardir,
	os.path.pardir,
)
site.addsitedir(lib_path)
from lib import *

###dta readin
data_path = '/exports/csce/datastore/geos/users/s1667168/CESM/Temp/Temp_pp/'
data=sio.loadmat(data_path+'TNX_BINS_PDF_1971_2000_JJAS_CESM30.mat')
data_NCEP = sio.loadmat(data_path+'TNX_BINS_PDF_1971_2000_JJAS_NCEP.mat')

############
##plotting##
############
import matplotlib.pyplot as plt
fig = plt.figure(facecolor='White',figsize=[8,4]);plot_setup();

######################
## Central East China
region = 'CEC_bins_pdfs'
sec_CESM = data[region][:]; sec_NCEP = data_NCEP[region]
bin_90 = np.empty((28));bin_95 = np.empty((28))
kde_90 = np.empty((28));kde_95 = np.empty((28))

ax = plt.subplot(2,2,1);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_ylim([-1,4.5]);ax.set_yticks(np.arange(0,4.6,1.5));ax.set_xlim([0,35]); ax.set_xticklabels([]);
ax.annotate('(a) CentralEast China',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='left', va='center',rotation = 0)	
#############Tn
for en_no in range(28):
	bins = sec_CESM[en_no,0,:]; kde = sec_CESM[en_no,1,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'k',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
print 'CEC TN'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[0,:]; kde = sec_NCEP[1,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'k',linewidth =1,label='Tmin')


#############Tm
for en_no in range(28):
	bins = sec_CESM[en_no,2,:]; kde = sec_CESM[en_no,3,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'r',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
print 'CEC TX'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[2,:]; kde = sec_NCEP[3,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'r',linewidth =1,label='Tmin')

######################
## South East China
region = 'SEC_bins_pdfs'
sec_CESM = data[region][:]; sec_NCEP = data_NCEP[region]
bin_90 = np.empty((28));bin_95 = np.empty((28))
kde_90 = np.empty((28));kde_95 = np.empty((28))

ax = plt.subplot(2,2,3);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_ylim([-1,6]);ax.set_yticks(np.arange(0,7.0,2.0));ax.set_xlim([0,35]); #ax.set_xticklabels([]);
ax.annotate('(b) SouthEast China',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='left', va='center',rotation = 0)
ax.annotate('Probability density function (%)',xy=(-0.10,1.1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 90)	
#############Tn
for en_no in range(28):
	bins = sec_CESM[en_no,0,:]; kde = sec_CESM[en_no,1,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'k',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
print 'SEC TN'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[0,:]; kde = sec_NCEP[1,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'k',linewidth =1,label='Tmin')
	
#############Tm
for en_no in range(28):
	bins = sec_CESM[en_no,2,:]; kde = sec_CESM[en_no,3,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'r',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
print 'SEC TX'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[2,:]; kde = sec_NCEP[3,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'r',linewidth =1,label='Tmin')

######################
## Northern India
region = 'NID_bins_pdfs'
sec_CESM = data[region][:]; sec_NCEP = data_NCEP[region]
bin_90 = np.empty((28));bin_95 = np.empty((28))
kde_90 = np.empty((28));kde_95 = np.empty((28))

ax = plt.subplot(2,2,2);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_ylim([-1,4.5]);ax.set_yticks(np.arange(0,4.6,1.5));ax.set_xlim([10,45]); ax.set_xticklabels([]);
ax.annotate('(c) Northern India',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='left', va='center',rotation = 0)
	
#############Tn
for en_no in range(28):
	bins = sec_CESM[en_no,0,:]; kde = sec_CESM[en_no,1,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'k',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
print 'nid TN'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[0,:]; kde = sec_NCEP[1,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'k',linewidth =1,label='Tmin')

#############Tm
for en_no in range(28):
	bins = sec_CESM[en_no,2,:]; kde = sec_CESM[en_no,3,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'r',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
print 'nid TX'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[2,:]; kde = sec_NCEP[3,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'r',linewidth =1,label='Tmin')

######################
## sOUTH iNDIA
region = 'SID_bins_pdfs'
sec_CESM = data[region][:]; sec_NCEP = data_NCEP[region]
bin_90 = np.empty((28));bin_95 = np.empty((28))
kde_90 = np.empty((28));kde_95 = np.empty((28))

ax = plt.subplot(2,2,4);pad= 5;
ax.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax.set_ylim([-1,4]);ax.set_yticks(np.arange(0,5,1.5));ax.set_xlim([10,45]);  
ax.annotate('(d) Southern India',xy=(0.01, 1.01), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='left', va='center',rotation = 0)
ax.annotate('Temperature (degree)',xy=(-0.1, -0.3), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=10, ha='center', va='center',rotation = 0)					
#############Tn
for en_no in range(28):
	bins = sec_CESM[en_no,0,:]; kde = sec_CESM[en_no,1,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'k',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='k',markersize=5,markeredgecolor='none',alpha=0.2)
print 'sid TN'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[0,:]; kde = sec_NCEP[1,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'k',linewidth =1,label='Tmin')
	
#############Tm
for en_no in range(28):
	bins = sec_CESM[en_no,2,:]; kde = sec_CESM[en_no,3,:]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
	bin_90[en_no]=bins[i];kde_90[en_no]=kde[i]
	i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
	bin_95[en_no]=bins[i];kde_95[en_no]=kde[i]
	ax.plot(bins,kde,'r',linewidth=1,alpha=0.1)
	ax.plot(bin_90[en_no],-0.5,marker ="<",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
	ax.plot(bin_95[en_no],-0.5,marker =">",color='r',markersize=5,markeredgecolor='none',alpha=0.2)
print 'sid TX'
print np.median(bin_90),np.median(bin_95)
bins =sec_NCEP[2,:]; kde = sec_NCEP[3,:]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 90 or (sum(kde[0:i])<=90 and sum(kde[0:i+1])>=90))][0]
ax.plot(bins[i],-0.5,marker ="<",color='k',markersize=5,markeredgecolor='r')
print bins[i]
i = [i for i in range(0,95) if (sum(kde[0:i]) == 95 or (sum(kde[0:i])<=95 and sum(kde[0:i+1])>=95))][0]
print bins[i]
ax.plot(bins[i],-0.5,marker =">",color='k',markersize=5,markeredgecolor='r')
ax.plot(bins,kde,'r',linewidth =1,label='Tmin')



plt.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.92, wspace=0.1, hspace=0.2);
plt.savefig('Fig1_dynamics.png', format='png', dpi=1200)
plt.show()








	