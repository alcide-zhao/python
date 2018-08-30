import numpy as np
import matplotlib.pyplot as plt

AtmAbsp = np.array([[0.77,-2.08,-1.31],[0.54,-0.48,0.06],[0.11,-0.30,-0.19],[-0.09,0.01,-0.08]])
Precip = np.array([[-9.78,30.72,20.94],[-5.62,0.50,-5.12],[-0.81,4.16,3.35],[0.67,-0.86,-0.19]])
SAT = np.array([[0.12,0.9,1.06],[0.02,0.07,0.09],[0.04,0.09,0.13],[-0.03,0.01,-0.02]])
ERF = np.array([[1.50,-1.24,0.26	],[-0.11,0.13,0.02],[0.24,-0.23,0.01],[-0.07,0.07,0.00]])

print SAT[:,0]
fig = plt.figure(facecolor='White',figsize=[8,4]);pad= 5;colormap='BrBG';
####################ERF###########################
ax = plt.subplot(1,2,1);
# ax.plot(ERF[:,0],SAT[:,0],marker='o',ls='none',markersize=20)
ax.plot(Precip[:,0],AtmAbsp[:,0],marker='s',ls='none',markersize=20)
ax2=ax.twinx()
ax2.plot(Precip[:,0],SAT[:,0],marker='D',ls='none',markersize=20,color='r')

ax = plt.subplot(1,2,2);
# ax.plot(ERF[:,0],SAT[:,0],marker='o',ls='none',markersize=20)
ax.plot(Precip[:,1],AtmAbsp[:,1],marker='s',ls='none',markersize=20)
ax2=ax.twinx()
ax2.plot(Precip[:,1],SAT[:,1],marker='D',ls='none',markersize=20,color='r')


plt.show()