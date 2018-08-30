import scipy.io as sio

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_TS_PS_850UV(UQVQ)_FSNS(C)_L(S)HFLX.mat')
time = data['time'][0,:];lon = data['lon'][0,:];lat = data['lat'][0,:];
TS_rcp85 = data['TS_rcp85'];TS_fixa = data['TS_fixa']; 
PSL_rcp85 = data['PSL_rcp85'];PSL_fixa = data['PSL_fixa']; 
U850_rcp85 = data['U850_rcp85'];U850_fixa = data['U850_fixa'];
V850_rcp85 = data['V850_rcp85'];V850_fixa = data['V850_fixa'];
Q850_rcp85 = data['Q850_rcp85'];Q850_fixa = data['Q850_fixa'];
FSNS_rcp85 = data['FSNS_rcp85'];FSNS_fixa = data['FSNS_fixa'];
FSNSC_rcp85 = data['FSNSC_rcp85'];FSNSC_fixa = data['FSNSC_fixa'];
LHFLX_rcp85 = data['LHFLX_rcp85'];LHFLX_fixa = data['LHFLX_fixa'];

input_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
data = sio.loadmat(input_path+'CESM_2006_2100_RCP_FIXA_PSL_UV(Q).mat')
UQ850_rcp85 = data['UQ850_rcp85'];UQ850_fixa = data['UQ850_fixa'];
VQ850_rcp85 = data['VQ850_rcp85'];VQ850_fixa = data['VQ850_fixa'];

output_path = '/exports/csce/datastore/geos/users/s1667168/PP/'
sio.savemat(output_path+'CESM_2006_2100_RCP_FIXA_TS_PSL_Q_UV(Q)_LHFLX_FSNS(C).mat',\
{'time':time,'lon':lon,'lat':lat,\
'TS_rcp85':TS_rcp85,'TS_fixa':TS_fixa,'PSL_rcp85':PSL_rcp85,'PSL_fixa':PSL_fixa,\
'U850_rcp85':U850_rcp85,'U850_fixa':U850_fixa,'V850_rcp85':V850_rcp85,'V850_fixa':V850_fixa,\
'UQ850_rcp85':UQ850_rcp85,'UQ850_fixa':UQ850_fixa,'VQ850_rcp85':VQ850_rcp85,'VQ850_fixa':VQ850_fixa,\
'Q850_rcp85':Q850_rcp85,'Q850_fixa':Q850_fixa,\
'FSNS_rcp85':FSNS_rcp85,'FSNS_fixa':FSNS_fixa,'FSNSC_rcp85':FSNSC_rcp85,'FSNSC_fixa':FSNSC_fixa,\
'LHFLX_rcp85':LHFLX_rcp85,'LHFLX_fixa':LHFLX_fixa})