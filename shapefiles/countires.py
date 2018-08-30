import scipy.io as sio
import matplotlib.pyplot as plt
data_path = './'
country_mask = 'Euro_USA_AUS_BRICS_720_360.mat'
data_map=sio.loadmat(data_path+'Euro_USA_AUS_BRICS_720_360.mat')
lon = data_map['lon'][0,:];lat = data_map['lat'][0,:];
Europe = data_map['Europe'][:]*1;USA = data_map['USA'][:]*2
Australia = data_map['Australia'][:]*3;India = data_map['India'][:]*4;
Brazil = data_map['Brazil'][:]*5;Russia = data_map['Russia'][:]*6;
SA = data_map['South Africa'][:]*7;China = data_map['China'][:]*8;

data=Europe+USA+Australia+India+Brazil+Russia+SA+China

plt.imshow(data,origin='lower');plt.show()