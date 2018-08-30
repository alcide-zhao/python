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
file = file_path+'sulphate_molecules.mat'
NA=6.002*10**(23)
data = sio.loadmat(file)
so2_surf = data['so2_surf'][:]/NA;so2_elev = data['so2_elev'][:]/NA;so4_a1_surf = data['so4_a1_surf'][:]/NA; 
so4_a1_elev = data['so4_a1_elev'][:]/NA;so4_a2_surf = data['so4_a2_surf'][:]/NA;
# print so4_a1_surf[0,:]
# print so4_a1_surf[1,:]
# print so4_a1_surf[2,:]
# print so4_a1_surf[3,:]
time_series=range(1,13)

fig1 = plt.figure(facecolor='White',figsize=[6,6]);plot_setup();

ax1 = plt.subplot(3,2,1);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.5,2.0]);ax1.set_yticks(np.arange(0.5,2.1,0.5))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(a) so2_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,so2_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,so2_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,so2_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,so2_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)


ax1 = plt.subplot(3,2,2);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([1,7]);ax1.set_yticks(np.arange(1,8,2))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(b) so2_elev',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,so2_elev[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,so2_elev[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,so2_elev[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,so2_elev[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)


ax1 = plt.subplot(3,2,3);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.008,0.016]);ax1.set_yticks(np.arange(0.006,0.0181,0.004))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(c) so4_a1_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,so4_a1_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,so4_a1_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,so4_a1_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,so4_a1_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

ax1 = plt.subplot(3,2,4);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.02,0.17]);ax1.set_yticks(np.arange(0.02,0.171,0.05))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(d) so4_a1_elev',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,so4_a1_elev[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,so4_a1_elev[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,so4_a1_elev[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,so4_a1_elev[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)


ax1 = plt.subplot(3,2,5);pad= 5;
ax1.grid(True,color='k', linestyle='-',alpha=0.2,linewidth = 0.5) 
ax1.set_ylim([0.005,0.041]);ax1.set_yticks(np.arange(0.005,0.042,0.012))
ax1.set_xticks(np.arange(0,12.1,3));ax1.set_xticklabels(['Jan','Mar','June','Sep','Dec']);ax1.set_xlim([1,12]);
ax1.annotate('(e) so4_a2_surf',xy=(0.1, 1.05), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='center',rotation='horizontal')
# ax1.set_ylabel('Aerosols and their Precursor emissions'+' ('+r'$\mathrm{\mathsf{10^{3}\times\/kg\/s^{-1}}}$'+')', color='k')

lins1 = ax1.plot(time_series,so4_a2_surf[0,:],'-',color="k",label='REF1970',linewidth=1.5)
lins2 = ax1.plot(time_series,so4_a2_surf[1,:],'-',color="b",label='REF2010',linewidth=1.5)
lins3 = ax1.plot(time_series,so4_a2_surf[2,:],'-',color="r",label='STAG_TECH',linewidth=1.5)
lins4 = ax1.plot(time_series,so4_a2_surf[3,:],'-',color="g",label='STAG_ENE',linewidth=1.5)

legend = ax1.legend(shadow=False,ncol=1,bbox_to_anchor=(2.0, 1.0))	 
legend.get_frame().set_facecolor('white');legend.get_frame().set_edgecolor('None');legend.get_frame().set_alpha(0.1)

plt.subplots_adjust(left=0.1, bottom=0.03, right=0.9, top=0.95, wspace=0.2, hspace=0.3);
plt.savefig('Fig1.png', format='png', dpi=1000)
# plt.show()


