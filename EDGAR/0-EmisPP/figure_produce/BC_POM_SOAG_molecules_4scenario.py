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
file = file_path+'BC_OC_soag_molecules.mat'
NA=6.002*10**(23)
data = sio.loadmat(file)
bc_surf = data['bc_surf'][:]/NA;oc_surf = data['oc_surf'][:]/NA;soag = data['soag'][:]/NA; 

time_series=range(1,13)

fig1 = plt.figure(facecolor='White',figsize=[5.5,5]);plot_setup();

ax1 = plt.subplot(2,2,1);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.5,2.0]);ax1.set_yticks(np.arange(0.5,2.1,0.5))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(a) BC_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,bc_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,bc_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,bc_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,bc_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

legend = ax1.legend(shadow=False,ncol=1,bbox_to_anchor=(2.0, 1.0))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.1)


ax1 = plt.subplot(2,2,3);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([2,6.5]);ax1.set_yticks(np.arange(2,6.51,1.5))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(b) POM_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,oc_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,oc_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,oc_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,oc_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)


ax1 = plt.subplot(2,2,4);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([15.5,24.5]);ax1.set_yticks(np.arange(15.5,24.6,3))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(c) SOAG',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,soag[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,soag[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,soag[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,soag[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.3);
plt.savefig('Fig2.png', format='png', dpi=1000)
# plt.show()


