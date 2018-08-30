
def datetime_second_datetime_integer(second_from_reference,reference = '1970-01-01 00:00:00'):
	'''
	Thus module is to conver the time from second to integer datetime Year*1000+Month*100+day*1
	If the the time is given in hours, please convert it to seconds before being inputted into this module
	tHE lINUX refers the 1970-01-01 :00:00:00 as the epoch, you need therefore to confirm the reference  
	time since which the hours/seconds are stored. If the reference time
	note the format of the reference time : %d/%m/%Y %H:%M:%S
	'''
	import time,datetime,calendar
	import numpy as np
	# the Linux epoch refers to 1970-01-01 00:00:00 fuck you, waste me a lot of time 
	second_refer_to_toepoch= calendar.timegm(time.strptime(reference,"%Y-%m-%d %H:%M:%S"))
	second_from_epoch = second_refer_to_toepoch +second_from_reference
	datetime = [time.gmtime(value) for value in second_from_epoch]
	datetime_float= np.array([value[0]*10000+value[1]*100+value[2] for value in datetime])[:]
	return datetime_float