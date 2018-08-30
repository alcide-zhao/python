import numpy as np
import netCDF4 as nc4
import time as clock
import os

input_path_0680 = '/exports/csce/datastore/geos/users/s1667168/CESM/temp/0680'
input_path_8000 = '/exports/csce/datastore/geos/users/s1667168/CESM/temp/8000'
file_path = '/exports/csce/datastore/geos/users/s1667168/CESM/temp/rcp8.5/'

# print os.path.abspath(input_path_0680)
os.system('find ' + os.path.abspath(input_path_0680) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path_0680 + '/file_list.txt')
text_file_0680 = open(input_path_0680 + '/file_list.txt', "r")

text_content_0680 = text_file_0680.readlines()
# print text_content_0680
os.system('find ' + os.path.abspath(input_path_8000) + ' -name "*' + '.nc' + '" -print | sort > ' + input_path_8000 + '/file_list.txt')
text_file_8000 = open(input_path_8000 + '/file_list.txt', "r")
text_content_8000 = text_file_8000.readlines()
# print text_content_8000


for ensumble_member in range(0,30): 
	nc_f_0680 = text_content_0680[ensumble_member][:-1]
	nc_f_8000 = text_content_8000[ensumble_member][:-1]
	fie_output = file_path+'Rcp8.5_2006_2100__Temp_M'+str(['_0'+str(ensumble_member+1) if ensumble_member<9 else '_'+str(ensumble_member+1)])[2:5]+'_EI_JJA.nc'
	nc_fid0680 = nc4.Dataset(nc_f_0680,mode='r')
	nc_fid8000 = nc4.Dataset(nc_f_8000,mode='r')
	att_list = ['lat','lon','time','txx', 'txn', 'tnx', 'tnn', 'dtr', 'fd0', 'su25','id0','tr20', 'tn10p', 'tn90p','tx10p','tx90p'] #, 
	att_dic = {'lat':[],'lon':[],'time':[],'txx':[], 'txn':[], 'tnx':[], 'tnn':[], 'dtr':[], 'fd0':[], 'su25':[],'id0':[], 'tr20':[], 'tn10p':[], 'tn90p':[],'tx10p':[],'tx90p':[]}
	# att_dic = {'total_precip':[],'lat':[],'lon':[],'time':[],'rx1day':[], 'rx5day':[], 'sdii':[], 'r10':[], 'r20':[], 'rnm':[], 'cdd':[],'cwd':[], 'r95p':[], 'r99p':[], 'precptot':[]}
	# att_list = ['total_precip','lat','lon','time','rx1day', 'rx5day', 'sdii', 'r10', 'r20', 'rnm', 'cdd','cwd', 'r95p', 'r99p', 'precptot']
	for att_name in att_list:
		# if (att_name == 'lat' or att_name == 'lon'):
		if (att_name == 'ensumble_no'):
			# 
			continue
		else:
			att_value0680 = nc_fid0680.variables[att_name][:]
			att_value8000 = nc_fid8000.variables[att_name][:]
			if (att_name == 'lat' or att_name == 'lon'):
				att_dic[att_name]=att_value0680
			else:
				att_dic[att_name] = np.concatenate((att_value0680, att_value8000), axis=0)
		

	# print att_dic

	f = nc4.Dataset(fie_output,'w', format='NETCDF4') #'w' stands for write
	f.createDimension('lat', len(att_dic['lat']))
	f.createDimension('lon', len(att_dic['lon']))
	f.createDimension('time', len(att_dic['time']))

	times = f.createVariable('time',np.float64, ('time'))
	lats = f.createVariable('lat',np.float32, ('lat'))
	lons = f.createVariable('lon',np.float32, ('lon'))

	"""
	## prep
	rx1days = f.createVariable('rx1day',np.float32,('time','lat','lon'))
	rx5days = f.createVariable('rx5day',np.float32,('time','lat','lon'))
	sdiis = f.createVariable('sdii',np.float32,('time','lat','lon'))
	r10s = f.createVariable('r10',np.float32,('time','lat','lon'))
	r20s = f.createVariable('r20',np.float32,('time','lat','lon'))
	rnms = f.createVariable('rnm',np.float32,('time','lat','lon'))
	cdds = f.createVariable('cdd',np.float32,('time','lat','lon'))
	cwds = f.createVariable('cwd',np.float32,('time','lat','lon'))
	r95ps = f.createVariable('r95p',np.float32,('time','lat','lon'))
	r99ps = f.createVariable('r99p',np.float32,('time','lat','lon'))
	prcptots = f.createVariable('precptot',np.float32,('time','lat','lon'))
	total_precips = f.createVariable('total_precip',np.float32,('time','lat','lon'))
		
	times[:] = att_dic['time']
	lats[:] = att_dic['lat']
	lons[:] = att_dic['lon']
	rx1days[:] = att_dic['rx1day']
	rx5days[:] = att_dic['rx5day']
	sdiis[:] = att_dic['sdii']
	r10s[:] = att_dic['r10']
	r20s[:] = att_dic['r20']
	rnms[:] = att_dic['rnm']
	cdds[:] = att_dic['cdd']
	cwds[:] = att_dic['cwd']
	r95ps[:] = att_dic['r95p']
	r99ps[:] =  att_dic['r99p']
	prcptots[:] =  att_dic['precptot']
	total_precips[:] = att_dic['total_precip']

	f.description = 'time series of precipitation extrme indeces calculated using the CESM model outputs'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'

	times.long_name = 'Year'
	times.units = '1'
		
	lats.long_name = 'latitude'
	lats.units = 'degree_north'
	lons.long_name = 'longitude'
	lons.units = 'degree_east'
		
	rx1days.standard_name = 'maximum 1-day precip amount'
	rx1days.long_name = 'maximum 1-day precipitation amount'
	rx1days.units = 'mm'
		
	rx5days.standard_name = 'maximum 5-day precip amount'
	rx5days.long_name = 'maximum 5-day precipitation amount'
	rx5days.units = 'mm'
		
	sdiis.standard_name = 'simple daily intensity index'
	sdiis.long_name = 'annual precip amount devided by the number of wet days (> 1mm)'
	sdiis.units = 'm/s/day'
		
	r10s.standard_name = 'number of heavy precip days'
	r10s.long_name = 'nnual count of days woth precipitation >= 10mm'
	r10s.units = 'day'
		
	r20s.standard_name = 'number of very heavy precip days'
	r20s.long_name = 'nnual count of days woth precipitation >= 20mm'
	r20s.units = 'day'	

	rnms.standard_name = 'number of days above nn mm here defines as 30 mm'
	rnms.long_name = 'annual count of days woth precipitation >= 30mm'
	rnms.units = 'day'	
		
	cdds.standard_name = 'consective dry days'
	cdds.long_name = 'maximum consecive days with precipitation<1mm'
	cdds.units = 'day'	
		
	cwds.standard_name = 'consective wet days'
	cwds.long_name = 'maximum consecive days with precipitation>=1mm'
	cwds.units = 'day'	
		
	r95ps.standard_name = 'very wet days'
	r95ps.long_name = 'annual total precipitation >95th percetile'
	r95ps.units = 'day'	
		
	r99ps.standard_name = 'extremely wet days'
	r99ps.long_name = 'annual total precipitation >99th percetile'
	r99ps.units = 'day'	

	prcptots.standard_name = 'annual total wet day prrecip'
	prcptots.long_name = 'annual total precipitation in wet days (rr>1mm)'
	prcptots.units = 'mm'		
		
	total_precips.standard_name = 'annual total prrecip'
	total_precips.long_name = 'annual total precipitation'
	total_precips.units = 'mm'
	f.close()
	"""
	
	# temp
	txxs = f.createVariable('txx',np.float32,('time','lat','lon'))
	txns = f.createVariable('txn',np.float32,('time','lat','lon'))
	tnxs = f.createVariable('tnx',np.float32,('time','lat','lon'))
	tnns = f.createVariable('tnn',np.float32,('time','lat','lon'))
	dtrs = f.createVariable('dtr',np.float32,('time','lat','lon'))
	fd0s = f.createVariable('fd0',np.float32,('time','lat','lon'))
	su25s = f.createVariable('su25',np.float32,('time','lat','lon'))
	id0s = f.createVariable('id0',np.float32,('time','lat','lon'))
	tr20s = f.createVariable('tr20',np.float32,('time','lat','lon'))
	tn10ps = f.createVariable('tn10p',np.float32,('time','lat','lon'))
	tn90ps = f.createVariable('tn90p',np.float32,('time','lat','lon'))
	tx10ps = f.createVariable('tx10p',np.float32,('time','lat','lon'))
	tx90ps = f.createVariable('tx90p',np.float32,('time','lat','lon'))


	times[:] = att_dic['time']
	lats[:] = att_dic['lat']
	lons[:] = att_dic['lon']
	txxs[:] = att_dic['txx']
	txns[:] = att_dic['txn']
	tnxs[:] = att_dic['tnx']
	tnns[:] = att_dic['tnn']
	dtrs[:] = att_dic['dtr']
	fd0s[:] = att_dic['fd0']
	su25s[:] = att_dic['su25']
	id0s[:] = att_dic['id0']
	tr20s[:] = att_dic['tr20']
	tn10ps[:] =  att_dic['tn10p']
	tn90ps[:] =  att_dic['tn90p']
	tx10ps[:] =  att_dic['tx10p']
	tx90ps[:] =  att_dic['tx90p']
	
	f.description = 'Temperature extrme indeces calculated using the CESM model outputs'
	f.history = 'Created at ' + clock.asctime( clock.localtime(clock.time()))
	f.institution = 'Alcide Zhao at the university of Edinburgh'
	
	times.long_name = 'Year'
	times.units = '1'
	
	lats.long_name = 'latitude'
	lats.units = 'degree_north'
	lons.long_name = 'longitude'
	lons.units = 'degree_east'
	
	txxs.standard_name = 'Max Tmax'
	txxs.long_name = 'JJA maximum of daily temperature maximam'
	txxs.units = 'Degrees'
	
	txns.standard_name = 'Min Tmax'
	txns.long_name = 'JJA minimum of daily temperature maxmam'
	txns.units = 'Degrees'
	
	tnxs.standard_name = 'Max Tmin'
	tnxs.long_name = 'JJA maximum of daily temperature minimam'
	tnxs.units = 'Degrees'
	
	tnns.standard_name = 'Min Tmin'
	tnns.long_name = 'JJA minimum of daily temperature minimam'
	tnns.units = 'Degrees'
	
	dtrs.standard_name = 'Diurnal temperature range'
	dtrs.long_name = 'JJA mean between daily TX and TM'
	dtrs.units = 'Degrees'	

	fd0s.standard_name = 'Frost day '
	fd0s.long_name = 'JJA count when daily minimum temperature <0 degree'
	fd0s.units = 'days'	
	
	su25s.standard_name = 'Summer days'
	su25s.long_name = 'JJA count when daily minimum temperature >25 degree'
	su25s.units = 'days'	
	
	id0s.standard_name = 'Ice day'
	id0s.long_name = 'JJA count when daily maximum temperature <0 degree'
	id0s.units = 'days'	
	
	tr20s.standard_name = 'Tropical night'
	tr20s.long_name = 'JJA count when daily minimum temperature >25 degree'
	tr20s.units = 'days'	
	
	tn10ps.standard_name = 'Cool night'
	tn10ps.long_name = 'Percent of days when TN <10th percentile'
	tn10ps.units = 'Percntile'	

	tn90ps.standard_name = 'warm night'
	tn90ps.long_name = 'Percent of days when TN >90th percentile'
	tn90ps.units = 'Percntile'		
	
	tx10ps.standard_name = 'cool day'
	tx10ps.long_name = 'Percent of days when TX <10th percentile'
	tx10ps.units = 'Percntile'
	
	tx90ps.standard_name = 'Warm days'
	tx90ps.long_name = 'Percent of days when TX >90th percentile'
	tx90ps.units = 'Percntile'

