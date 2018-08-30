"""
This is to produce the bar plotof radiative forcing 
following the IPCC AR5 cormat
"""


import numpy as np
import matplotlib.pyplot as plt

RF={'total':[0.54,0.12],'GHG':[0.72,0.11],'ozone':[-0.07,0.14],'BEoA':[-0.11,0.15],'ENE':[-0.27,0.07],'TECH':[0.04,0.07]}

fig = plt.figure(facecolor='White',figsize=[12,5]);pad= 5;

ax = plt.subplot(1,1,1);
ax.tick_params(axis='both', which='major', direction = 'out',left=False,right=False,bottom=True,top=False, pad=5)
# ax.axvline(x=3,color='k',linewidth = 2);
ax.axhline(y=0,color='k',linewidth = 2);ax.axhline(y=7,color='k',linewidth = 2);
ax.set_xlim([-2.5,1.8]);ax.set_xticks(np.arange(-0.6,1.0,0.3));
ax.set_ylim([-6,6.0]);ax.set_yticks(np.arange(-6,6.01,2))
ax.axvline(x=-0.6,color='k',linewidth = 2);ax.axvline(x=0,ymax =6/7.0,color='k',linewidth = 2);
ax.axvline(x=-0.3,ymax =6/7.0,color='k',ls='--',lw = 0.5);ax.axvline(x=0.3,ymax =6/7.0,color='k',ls='--',lw = 0.5);
ax.axvline(x=0.6,ymax =6/7.0,color='k',ls='--',lw = 0.5);
ax.axvline(x=-2.5,color='k',linewidth = 2);ax.axvline(x=0.9,color='k',linewidth = 2);ax.axvline(x=1.8,color='k',linewidth = 2);


ax.set_ylim([0,7.0]);ax.set_yticks([])

ax.axhspan(1, 2, alpha=0.1, color='y',edgecolor='none');
ax.axhspan(3, 4, alpha=0.1, color='y',edgecolor='none')
ax.axhspan(5, 6, alpha=0.1, color='y',edgecolor='none')
ax.axhline(y=1,color='k',linewidth = 2);
ax.axhline(y=6,color='k',linewidth = 2);



# X_pos=np.arange(0,7.5,2.5);dis=0.5
ax.annotate('RF terms (2010 -1970)', xy=(-2.4,6.3),fontsize=15)
ax.annotate(r'  RF values $(W m^{-2})$', xy=(1.0,6.3),fontsize=15)

key = 'GHG'
ax.barh(5.5, RF[key][0], xerr=RF[key][1], align='center',color='g', ecolor='k',height = 0.8,lw=0,capsize=10)
ax.annotate('Greenhouse gases (GHGs)', xy=(-2.4,5.3),fontsize=15)
ax.annotate('       0.72+/-0.12', xy=(1.0,5.3),fontsize=15)

key = 'ozone'
ax.barh(4.5, RF[key][0], xerr=RF[key][1], align='center',color='m', ecolor='k',height = 0.8,lw=0,capsize=10)
ax.annotate('Ozone', xy=(-2.4,4.3),fontsize=15)
ax.annotate('      -0.07+/-0.11', xy=(1.0,4.3),fontsize=15)

key = 'BEoA'
ax.barh(3.5, RF[key][0], xerr=RF[key][1], align='center',color='b', ecolor='k',height = 0.8,lw=0,capsize=10)
ax.annotate('Anthropogenic aersosols (AAs)', xy=(-2.4,3.3),fontsize=15)
ax.annotate('      -0.11+/-0.14', xy=(1.0,3.3),fontsize=15)

key = 'ENE'
ax.barh(2.5, RF[key][0], xerr=RF[key][1], align='center',color='y', ecolor='k',height = 0.8,lw=0,capsize=10)
ax.annotate('Anthropogenic aerosol increases\ndue to energy consumption', xy=(-2.4,2.3),fontsize=15)
ax.annotate('      -0.27+/-0.15', xy=(1.0,2.3),fontsize=15)

key = 'TECH'
ax.barh(1.5, RF[key][0], xerr=RF[key][1], align='center',color='orange', ecolor='k',height = 0.8,lw=0,capsize=10)
ax.annotate('Anthropogenic aerossol reduction\ndue to technology advancements', xy=(-2.4,1.3),fontsize=15)
ax.annotate('       0.04+/-0.07', xy=(1.0,1.3),fontsize=15)


key = 'total'
ax.barh(0.5, RF[key][0], xerr=RF[key][1], align='center',color='r', ecolor='k',height = 0.8,lw=0,capsize=10)
ax.annotate('Total net (GHGs+AAs+Ozone)', xy=(-2.4,0.3),fontsize=15)
ax.annotate('       0.54+/-0.07', xy=(1.0,0.3),fontsize=15)
# plt.show()


plt.subplots_adjust(left=0.03, bottom=0.05, right=0.97, top=0.98, wspace=0.04, hspace=0.15); 
plt.savefig('RF.png', format='png', dpi=1000)