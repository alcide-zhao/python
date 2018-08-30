'''
plotting out the number of so2/so4 molecules in a1 and a2, surface and elevated level
'''


import scipy
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats
import math
import site
import scipy.io as sio


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



################################################
# Loading the so2/so4 molecule data
################################################

# emissions
file_path = '/exports/csce/datastore/geos/users/s1667168/EDGAR/'
file = file_path+'aerosol_particles.mat'
NA=6.002*10**(26)
data = sio.loadmat(file)
SO4_a1_surf = data['SO4_a1_surf'][:]/NA/10**(15);SO4_a1_elev = data['SO4_a1_elev'][:]/NA/10**(15);SO4_a2_surf = data['SO4_a2_surf'][:]/NA/10**(15); 
OC_a1_surf = data['OC_a1_surf'][:]/NA/10**(15);BC_a1_surf = data['BC_a1_surf'][:]/NA/10**(15);

time_series=range(1,13)

fig1 = plt.figure(facecolor='White',figsize=[6,6]);plot_setup();

ax1 = plt.subplot(3,2,1);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.08,0.17]);ax1.set_yticks(np.arange(0.08,0.171,0.03))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(a) SO4_a1_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,SO4_a1_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,SO4_a1_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,SO4_a1_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,SO4_a1_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)


ax1 = plt.subplot(3,2,2);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0,1.2]);ax1.set_yticks(np.arange(0,1.21,0.4))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(b) SO4_a1_elev',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,SO4_a1_elev[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,SO4_a1_elev[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,SO4_a1_elev[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,SO4_a1_elev[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)


ax1 = plt.subplot(3,2,3);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([5,38]);ax1.set_yticks(np.arange(5,38.1,11))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(c) SO4_a2_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,SO4_a2_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,SO4_a2_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,SO4_a2_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,SO4_a2_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

legend = ax1.legend(shadow=False,ncol=1,bbox_to_anchor=(2.0, 1.0))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.1)



ax1 = plt.subplot(3,2,5);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.6,1.8]);ax1.set_yticks(np.arange(0.6,1.81,0.4))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(d) BC_a1_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,BC_a1_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,BC_a1_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,BC_a1_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,BC_a1_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

ax1 = plt.subplot(3,2,6);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([4,10.3]);ax1.set_yticks(np.arange(4,10.4,2.1))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(e) POM_a1_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,OC_a1_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,OC_a1_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,OC_a1_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,OC_a1_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

plt.subplots_adjust(left=0.1, bottom=0.03, right=0.9, top=0.95, wspace=0.2, hspace=0.3);
plt.savefig('Fig3.png', format='png', dpi=1000)
# plt.show()


